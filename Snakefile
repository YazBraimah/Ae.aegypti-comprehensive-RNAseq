"""
Author: Y. Ahmed-Braimah
--- Snakemake workflow to process public Ae. aegypti RNAseq data
"""

import json
from os.path import join, basename, dirname
from os import getcwd
from subprocess import check_output
import os

##--------------------------------------------------------------------------------------##
## Functions
##--------------------------------------------------------------------------------------##

# To print process messages
def message(x):
  print()

# To remove suffix from a string
def rstrip(text, suffix):
    if not text.endswith(suffix):
        return text
    return text[:len(text)-len(suffix)]

##--------------------------------------------------------------------------------------##
## Global config parameters: 
##--------------------------------------------------------------------------------------##

configfile: 'config.yml'

# Full path to an uncompressed FASTA file with all chromosome sequences.
DNA = config['DNA']
INDEX = config['INDEX']

# Full path to an uncompressed GTF file with all gene annotations.
GTF = config['GTF']

# Full path to a folder where final output files will be deposited.
OUT_DIR = config['OUT_DIR']
WORK_DIR = config['WORK_DIR']
HOME_DIR = config['HOME_DIR']  # the "launch_snakemake.sh" and "config.yml" files should be here

## set the usr and job environments for each job (specific for CBSU qsub jobs)
USER = os.environ.get('USER')
JOB_ID = os.environ.get('JOB_ID')

# Samples and their corresponding filenames.
# single-end
seFILES = json.load(open(config['SE_SAMPLES_JSON'])) 
seSAMPLES = sorted(seFILES.keys())                  
# paired-end:
peFILES = json.load(open(config['PE_SAMPLES_JSON'])) 
peSAMPLES = sorted(peFILES.keys())           

# read both
# FILES = json.load(open(config['SAMPLES_JSON']))
combinedSam = [peSAMPLES, seSAMPLES]
SAMPLES = [y for x in combinedSam for y in x]  

## Create the final output directory if it doesn't already exist
if not os.path.exists(OUT_DIR):
            os.makedirs(OUT_DIR)

## In case fsrv5 needs to be mounted
# os.system("/programs/bin/labutils/mount_server cbsufsrv5 /data1")

##--------------------------------------------------------------------------------------##
## RULES
##--------------------------------------------------------------------------------------##


## Final expected output(s)
rule all: 
    input: 
        join(OUT_DIR, 'MultiQC', 'multiqc_report.html'),
        join(OUT_DIR, 'StringTie', 'gffcmp.annotated.gtf'),
        join(OUT_DIR, 'ballgown', 'gene_counts.csv'), 
        join(OUT_DIR, 'ballgown', 'transcript_counts.csv')
        # join(OUT_DIR, 'cleanUp_complete')
        
##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

## Rule to check raw SE read quality
rule fastqcSE:
    input:
        r1 = lambda wildcards: seFILES[wildcards.sample]['R1']
    output:
        r1 = join(OUT_DIR, 'fastQC', '{sample}' + '.R1_fastqc.html')
    log:
        join(OUT_DIR, 'fastQC', '{sample}' + '.fastQC_se.log')
    benchmark:
        join(OUT_DIR, 'fastQC', '{sample}' + '.fastQC_se.benchmark.tsv')
    message: 
        """--- Checking read quality of SE sample "{wildcards.sample}" with FastQC """
    run:
        if not os.path.exists(join(OUT_DIR, 'fastQC')):
            os.makedirs(join(OUT_DIR, 'fastQC'))
        shell('/programs/bin/labutils/mount_server cbsufsrv5 /data1')
        shell('mkdir -p ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cp {input.r1} ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cd ' + join(WORK_DIR, USER, JOB_ID) + 
                ' && fastqc {wildcards.sample}.R1.fq.gz' 
                ' > {log} 2>&1'
                ' && mv ' + join(WORK_DIR, USER, JOB_ID) + '/*fastqc* ' + join(OUT_DIR, 'fastQC'))
        shell('rm -r ' + join(WORK_DIR, USER, JOB_ID))

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

## Rule to check raw PE read quality
rule fastqcPE:
    input:
        r1 = lambda wildcards: peFILES[wildcards.sample]['R1'],
        r2 = lambda wildcards: peFILES[wildcards.sample]['R2']
    output:
        r1 = join(OUT_DIR, 'fastQC', '{sample}' + '.R1_fastqc.html'),
        r2 = join(OUT_DIR, 'fastQC', '{sample}' + '.R2_fastqc.html')
    log:
        join(OUT_DIR, 'fastQC', '{sample}' + '.fastQC_init_pe.log')
    benchmark:
        join(OUT_DIR, 'fastQC', '{sample}' + '.fastQC_init_pe.benchmark.tsv')
    message: 
        """--- Checking read quality of PE sample "{wildcards.sample}" with FastQC """
    run:
        if not os.path.exists(join(OUT_DIR, 'fastQC')):
            os.makedirs(join(OUT_DIR, 'fastQC'))
        shell('/programs/bin/labutils/mount_server cbsufsrv5 /data1')
        shell('mkdir -p ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cp {input.r1} {input.r2} ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cd ' + join(WORK_DIR, USER, JOB_ID) + 
                ' && fastqc {wildcards.sample}.R1.fq.gz {wildcards.sample}.R2.fq.gz' 
                ' > {log} 2>&1'
                ' && mv ' + join(WORK_DIR, USER, JOB_ID) + '/*fastqc* ' + join(OUT_DIR, 'fastQC'))
        shell('rm -r ' + join(WORK_DIR, USER, JOB_ID))

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

## Rule to map PE reads with HISAT2
rule hisat2_se_mapping:
    input:
        r1 = lambda wildcards: seFILES[wildcards.sample]['R1']
    output:
        bam = join(OUT_DIR, 'HISAT-2', '{sample}', '{sample}' + '.csorted.bowtie2.bam')
    log:
        join(OUT_DIR, 'HISAT-2', '{sample}', 'hs2_map_se.log')
    benchmark:
        join(OUT_DIR, 'HISAT-2', '{sample}', 'hs2_map_se.benchmark.tsv')
    message: 
        """--- Mapping SE reads for sample {wildcards.sample} to genome with HISAT-2 """
    run:
        shell('/programs/bin/labutils/mount_server cbsufsrv5 /data1')
        shell('mkdir -p ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cp ' + join(dirname(DNA), rstrip(os.path.basename(DNA), '.fa') + '*') + ' ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cp ' + join(INDEX, '*') + ' ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cp {input.r1} ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cd ' + join(WORK_DIR, USER, JOB_ID) + 
                ' && (hisat2'
                ' -p 16'
                ' --dta'
                ' -x ' + join(rstrip(os.path.basename(DNA), '.fa') + '_tran') +
                ' -U {wildcards.sample}.R1.fq.gz) 2>{log}'
                ' | samtools sort -@ 8 -o csorted.bowtie2.bam -')
        shell('mv ' + join(WORK_DIR, USER, JOB_ID, 'csorted.bowtie2.bam') + ' ' + join(OUT_DIR, 'HISAT-2', '{wildcards.sample}', '{wildcards.sample}' + '.csorted.bowtie2.bam'))
        shell('rm -r ' + join(WORK_DIR, USER, JOB_ID))

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

## Rule to check raw PE read quality
rule hisat2_pe_mapping:
    input:
        r1 = lambda wildcards: peFILES[wildcards.sample]['R1'],
        r2 = lambda wildcards: peFILES[wildcards.sample]['R2']
    output:
        bam = join(OUT_DIR, 'HISAT-2', '{sample}', '{sample}' + '.csorted.bowtie2.bam')
    log:
        join(OUT_DIR, 'HISAT-2', '{sample}', 'hs2_map_pe.log')
    benchmark:
        join(OUT_DIR, 'HISAT-2', '{sample}', 'hs2_map_pe.benchmark.tsv')
    message: 
        """--- Mapping PE reads for sample {wildcards.sample} to genome with HISAT-2 """
    run:
        shell('/programs/bin/labutils/mount_server cbsufsrv5 /data1')
        shell('mkdir -p ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cp ' + join(dirname(DNA), rstrip(os.path.basename(DNA), '.fa') + '*') + ' ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cp ' + join(INDEX, '*') + ' ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cp {input.r1} {input.r2} ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cd ' + join(WORK_DIR, USER, JOB_ID) + 
                ' && (hisat2'
                ' -p 16'
                ' --dta'
                ' -x ' + join(rstrip(os.path.basename(DNA), '.fa') + '_tran') +
                ' -1 {wildcards.sample}.R1.fq.gz'
                ' -2 {wildcards.sample}.R2.fq.gz)'
                ' 2>{log}'
                ' | samtools sort -@ 8 -o csorted.bowtie2.bam -')
        shell('mv ' + join(WORK_DIR, USER, JOB_ID, 'csorted.bowtie2.bam') + ' ' + join(OUT_DIR, 'HISAT-2', '{wildcards.sample}', '{wildcards.sample}' + '.csorted.bowtie2.bam'))
        shell('rm -r ' + join(WORK_DIR, USER, JOB_ID))

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

## Rule to assemble transcripts with StringTie
rule stringtie_assembly:
    input:
        bam = join(OUT_DIR, 'HISAT-2', '{sample}', '{sample}' + '.csorted.bowtie2.bam')
    output:
        asmbly = join(OUT_DIR, 'StringTie', '{sample}', '{sample}' + '.gtf')
    params:
        gtf = GTF
    log:
        join(OUT_DIR, 'StringTie', '{sample}', 'st_asmbly.log')
    benchmark:
        join(OUT_DIR, 'StringTie', '{sample}', 'st_asmbly.benchmark.tsv')
    message: 
        """--- Assembling transcripts for sample {wildcards.sample} with StringTie """
    run:
        shell('/programs/bin/labutils/mount_server cbsufsrv5 /data1')
        shell('mkdir -p ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cp {params.gtf} ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cp {input.bam} ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cd ' + join(WORK_DIR, USER, JOB_ID) + 
                ' && stringtie'
                ' {wildcards.sample}.csorted.bowtie2.bam'
                ' -p 16'
                ' -G ' + os.path.basename(GTF) +
                ' -o {wildcards.sample}.gtf'
                ' -l {wildcards.sample} > {log}')
        shell('mv ' + join(WORK_DIR, USER, JOB_ID, '{wildcards.sample}.gtf') + ' ' + join(OUT_DIR, 'StringTie', '{wildcards.sample}'))
        shell('rm -r ' + join(WORK_DIR, USER, JOB_ID))

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

## Rule to merge StringTie assemblies
rule merge_assemblies:
    input:
        assemblies = expand(join(OUT_DIR, 'StringTie', '{sample}', '{sample}' + '.gtf'), sample = SAMPLES)
    output:
        asmbly = join(OUT_DIR, 'StringTie', 'stringtie_merged.gtf')
    params:
        gtf = GTF
    log:
        join(OUT_DIR, 'StringTie', 'st_mrg.index.log')
    benchmark:
        join(OUT_DIR, 'StringTie', 'st_mrg.index.benchmark.tsv')
    message: 
        """--- Merging StringTie transcripts """
    run:
        shell('/programs/bin/labutils/mount_server cbsufsrv5 /data1')
        shell('mkdir -p ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cp {params.gtf} ' + join(WORK_DIR, USER, JOB_ID) +
                ' && ls -1 ' + join(OUT_DIR) + '/StringTie/*/*.gtf > ' + join(OUT_DIR, 'StringTie', 'assemblies.txt') +
                ' && cd ' + join(WORK_DIR, USER, JOB_ID) + 
                ' && stringtie'
                ' --merge'
                ' -p 8'
                ' -G ' + os.path.basename(GTF) +
                ' -o stringtie_merged.gtf ' + join(OUT_DIR, 'StringTie', 'assemblies.txt'))
        shell('mv ' + join(WORK_DIR, USER, JOB_ID, 'stringtie_merged.gtf') + ' ' + join(OUT_DIR, 'StringTie'))
        shell('rm -r ' + join(WORK_DIR, USER, JOB_ID))

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

## Rule to merge StringTie assemblies
rule compare_gtf:
    input:
        STasmbly = rules.merge_assemblies.output.asmbly
    output:
        asmbly = join(OUT_DIR, 'StringTie', 'gffcmp.annotated.gtf')
    params:
        gtf = GTF
    message: 
        """--- Comparing StringTie merged assembly with reference GTF """
    run:
        shell('/programs/bin/labutils/mount_server cbsufsrv5 /data1')
        shell('mkdir -p ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cp {params.gtf} {input.STasmbly} ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cd ' + join(WORK_DIR, USER, JOB_ID) + 
                ' && gffcompare'
                ' -G'
                ' -r ' + os.path.basename(GTF) +
                ' stringtie_merged.gtf')
        shell('mv ' + join(WORK_DIR, USER, JOB_ID, 'gffcmp.*') + ' ' + join(OUT_DIR, 'StringTie'))
        shell('rm -r ' + join(WORK_DIR, USER, JOB_ID))

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

## Rule to measure transcript abundances with Stringtie
rule abundances:
    input:
        bam = join(OUT_DIR, 'HISAT-2', '{sample}', '{sample}' + '.csorted.bowtie2.bam'),
        mrgd = rules.merge_assemblies.output.asmbly
    output:
        abundance = join(OUT_DIR, 'ballgown', '{sample}', '{sample}' + '_abundance.gtf')
    log:
        join(OUT_DIR, 'ballgown', '{sample}', 'st_abnd.log')
    benchmark:
        join(OUT_DIR, 'ballgown', '{sample}', 'st_abnd.benchmark.tsv')
    message: 
        """--- Estimating transcript abundances for sample {wildcards.sample} with StringTie"""
    run:
        shell('/programs/bin/labutils/mount_server cbsufsrv5 /data1')
        shell('mkdir -p ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cp {input.bam} {input.mrgd} ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cd ' + join(WORK_DIR, USER, JOB_ID) + 
                ' && stringtie'
                ' -e -B -p 8'
                ' -G stringtie_merged.gtf' 
                ' -o {wildcards.sample}_abundance.gtf'
                ' {wildcards.sample}.csorted.bowtie2.bam'
                ' > {log}')
        shell('mv ' + join(WORK_DIR, USER, JOB_ID, '{wildcards.sample}_abundance.gtf') + ' ' + join(OUT_DIR, 'ballgown', '{wildcards.sample}'))
        shell('rm -r ' + join(WORK_DIR, USER, JOB_ID))

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

## Rule to combine abundance counts for downstream analysis
rule collate_counts:
    input:
        abundances = expand(join(OUT_DIR, 'ballgown', '{sample}', '{sample}' + '_abundance.gtf'), sample = SAMPLES)
    output:
        geneCounts = join(OUT_DIR, 'ballgown', 'gene_counts.csv'),
        transcriptCounts = join(OUT_DIR, 'ballgown', 'transcript_counts.csv')
    message: 
        """--- Outputting count matrices """
    run:
        shell('/programs/bin/labutils/mount_server cbsufsrv5 /data1')
        shell('prepDE.py'
                ' -i ' + join(OUT_DIR, 'ballgown') + 
                ' -g ' + join(OUT_DIR, 'ballgown', 'gene_counts.csv') +
                ' -t ' + join(OUT_DIR, 'ballgown', 'transcript_counts.csv'))

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

## Rule to collate fastQC and HISAT2 outputs with multiQC
rule multiQC:
    input:
        expand(join(OUT_DIR, 'HISAT-2', '{sample}', '{sample}' + '.csorted.bowtie2.bam'), sample = SAMPLES),
        expand(join(OUT_DIR, 'fastQC', '{sample}' + '.R1_fastqc.html'), sample = SAMPLES),
        expand(join(OUT_DIR, 'fastQC', '{sample}' + '.R2_fastqc.html'), sample = peSAMPLES)
    output:
        join(OUT_DIR, 'MultiQC', 'multiqc_report.html')
    log:
        join(OUT_DIR, 'MultiQC', 'multiQC.log')
    benchmark:
        join(OUT_DIR, 'MultiQC', 'multiQC.benchmark.tsv')
    message: 
        """--- Running MultiQC """
    run:
        shell('/programs/bin/labutils/mount_server cbsufsrv5 /data1')
        shell('ls -1 ' + join(OUT_DIR) + '/HISAT-2/*/*log > ' + join(OUT_DIR, 'MultiQC', 'summary_files.txt'))
        shell('ls -1 ' + join(OUT_DIR) + '/fastQC/*fastqc.zip >> ' + join(OUT_DIR, 'MultiQC', 'summary_files.txt'))
        shell('multiqc'
                ' -f'
                ' -o ' + join(OUT_DIR, 'MultiQC') + ' -d -dd 2 -l ' + join(OUT_DIR, 'MultiQC', 'summary_files.txt') +
                ' > {log} 2>&1')

# ## Rule to collate fastQC and HISAT2 outputs with multiQC
# rule clean_up:
#     input:
#         expand(join(OUT_DIR, 'HISAT-2', '{sample}', '{sample}' + '.csorted.bowtie2.bam'), sample = SAMPLES),
#         expand(join(OUT_DIR, 'fastQC', '{sample}' + '.R1_fastqc.html'), sample = SAMPLES),
#         expand(join(OUT_DIR, 'fastQC', '{sample}' + '.R2_fastqc.html'), sample = peSAMPLES),
#         join(OUT_DIR, 'ballgown', 'gene_counts.csv')
#     output:
#         join(OUT_DIR, 'cleanUp_complete')
#     message: 
#         """--- Cleaning Up """
#     run:
#         shell('rm -r ' + join(OUT_DIR, 'HISAT-2', '*', '*.csorted.bowtie2.bam'))
#         shell('touch ' + join(OUT_DIR) + 'cleanUp_complete')
