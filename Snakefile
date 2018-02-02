"""
Author: Y. Ahmed-Braimah
--- Snakemake workflow to process public Ae. aegypti RNAseq data
"""

import json
import pandas as pd
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

# Load SRA runs file
RUN_INFO = pd.read_csv(config['SRA_INFO'])

# define SRA sample list
sraSAMPLES = RUN_INFO["Sample_Name"].tolist()

# split SRA sample list by library type (pe vs. se)
SRAseSAMPLES = list(RUN_INFO[RUN_INFO["LibraryLayout"] == "SINGLE"]["Sample_Name"])
SRApeSAMPLES = list(RUN_INFO[RUN_INFO["LibraryLayout"] == "PAIRED"]["Sample_Name"]) 

# Local samples and their corresponding filenames.
# single-end
local_seFILES = json.load(open(config['SE_SAMPLES_JSON'])) 
LOCseSAMPLES = sorted(local_seFILES.keys())                  
# paired-end:
local_peFILES = json.load(open(config['PE_SAMPLES_JSON'])) 
LOCpeSAMPLES = sorted(local_peFILES.keys())           

# read both
combinedLocSam = [LOCpeSAMPLES, LOCseSAMPLES]
locSAMPLES = [y for x in combinedLocSam for y in x]  

# also combine the two different file sources by library type:
combinedAllSam = [locSAMPLES, sraSAMPLES]
SAMPLES = [y for x in combinedAllSam for y in x]

# combine local and SRA samples:
combinedAllpeSam = [LOCpeSAMPLES, SRApeSAMPLES]
peSAMPLES = [y for x in combinedAllpeSam for y in x]

combinedAllseSam = [LOCseSAMPLES, SRAseSAMPLES]
seSAMPLES = [y for x in combinedAllseSam for y in x]

# define run and ID
RUNS = RUN_INFO["Run"].tolist()
LIBL = RUN_INFO["LibraryLayout"].tolist()

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


## Create the final output directory if it doesn't already exist
if not os.path.exists(OUT_DIR):
            os.makedirs(OUT_DIR)


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

## Rule to fetch raw SRA data
rule get_sra_data:
    output: 
        raw_sra = join(OUT_DIR, 'Reads', 'sra', '{lib}', '{sample}', '{run}', '{run}.sra')
    log:
        join(OUT_DIR, 'Reads', 'sra', '{lib}', '{sample}', '{run}', '{run}.log')
    benchmark:
        join(OUT_DIR, 'Reads', 'sra', '{lib}', '{sample}', '{run}', '{run}.benchmark.tsv')
    message: 
        """--- Downloading raw SRA record for "{wildcards.sample}" """
    params: 
        run_prefix=lambda wildcards: wildcards.run[:6], 
        sra_prefix=lambda wildcards: wildcards.run[:3]
    run:
        if not os.path.exists(join(OUT_DIR, 'fastQC')):
            os.makedirs(join(OUT_DIR, 'fastQC'))
        # shell('/programs/bin/labutils/mount_server cbsufsrv5 /data1')
        shell('mkdir -p ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cd ' + join(WORK_DIR, USER, JOB_ID) + 
                ' && wget -O ' + join(WORK_DIR, USER, JOB_ID, '{run}.sra') + 
                ' ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/{params.sra_prefix}/{params.run_prefix}/{wildcards.run}/{wildcards.run}.sra')
        shell('mv ' + join(WORK_DIR, USER, JOB_ID, '{run}.sra') + ' ' + join(OUT_DIR, 'Reads', 'sra', '{lib}', '{sample}', '{run}'))
        shell('rm -r ' + join(WORK_DIR, USER, JOB_ID))

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

## Rule to convert SRA format raw file to fastq
rule get_fastq_files_from_sra_file:
    input: 
        raw_sra = rules.get_sra_data.output.raw_sra
    output: 
        raw_fastq = join(OUT_DIR, 'Reads', 'fastq', '{lib}', '{sample}', '{run}', '{run}_1.fastq.gz')
    log:
        join(OUT_DIR, 'Reads', 'fastq', '{lib}', '{sample}', '{run}', '{run}.log')
    benchmark:
        join(OUT_DIR, 'Reads', 'fastq', '{lib}', '{sample}', '{run}', '{run}.benchmark.tsv')
    message: 
        """--- Downloading raw SRA record for "{wildcards.sample}" """
    run:
        shell('mkdir -p ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cd ' + join(WORK_DIR, USER, JOB_ID) + 
                ' && fastq-dump --defline-seq \'@[$ac_]$sn/$ri\' --defline-qual \'+\' --split-files --gzip --outdir {wildcards.run}.output {wildcards.run}.sra')
        shell('mv ' + join(WORK_DIR, USER, JOB_ID, '{wildcards.run}.output', '{wildcards.run}') + '*fastq.gz ' + join(OUT_DIR, 'Reads', 'fastq', '{lib}', '{sample}', '{run}'))
        shell('rm -r ' + join(WORK_DIR, USER, JOB_ID))

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

## Rule to merge multiple se fastq files from teh same sample
rule merge_se_fastq_files:
    output:
        join(OUT_DIR, 'Reads', 'fastq', 'SINGLE', '{sample}', 'merged', '{sample}.R1.fastq.gz')
    message: 
        """--- Merging SE fastq files for sample "{wildcards.sample}" """
    run:
        shell('mkdir -p ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cd ' + join(WORK_DIR, USER, JOB_ID) + 
                ' && cp ' + join(OUT_DIR, 'Reads', 'fastq', 'SINGLE', '{wildacrds.sample}', '*', '*_1.fastq.gz') + ' .' +        
                ' && zcat *_1.fastq.gz | gzip - > {wildcards.sample}.R1.fastq.gz')
        shell('mv ' + join(WORK_DIR, USER, JOB_ID, '{wildcards.sample}.R1.fastq.gz') + ' ' + join(OUT_DIR, 'Reads', 'fastq', 'SINGLE', '{wildcards.sample}', 'merged', '{wildcards.sample}.R1.fastq.gz'))
        shell('rm -r ' + join(WORK_DIR, USER, JOB_ID))


##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

## Rule to merge multiple se fastq files from teh same sample
rule merge_left_pe_fastq_files:
    output:
        join(OUT_DIR, 'Reads', 'fastq', 'PAIRED', '{sample}', 'merged', '{sample}.R1.fastq.gz')
    message: 
        """--- Merging left fastq PE reads for sample "{wildcards.sample}" """
    run:
        shell('mkdir -p ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cd ' + join(WORK_DIR, USER, JOB_ID) + 
                ' && cp ' + join(OUT_DIR, 'Reads', 'fastq', 'PAIRED', '{wildacrds.sample}', '*', '*_1.fastq.gz') + ' .' +        
                ' && zcat *_1.fastq.gz | gzip - > {wildcards.sample}.R1.fastq.gz')
        shell('mv ' + join(WORK_DIR, USER, JOB_ID, '{wildcards.sample}.R1.fastq.gz') + ' ' + join(OUT_DIR, 'Reads', 'fastq', 'PAIRED', '{wildcards.sample}', 'merged', '{wildcards.sample}.R1.fastq.gz'))
        shell('rm -r ' + join(WORK_DIR, USER, JOB_ID))


##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

## Rule to merge multiple se fastq files from teh same sample
rule merge_right_pe_fastq_files:
    output:
        join(OUT_DIR, 'Reads', 'fastq', 'PAIRED', '{sample}', 'merged', '{sample}.R2.fastq.gz')
    message: 
        """--- Merging right fastq PE reads for sample "{wildcards.sample}" """
    run:
        shell('mkdir -p ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cd ' + join(WORK_DIR, USER, JOB_ID) + 
                ' && cp ' + join(OUT_DIR, 'Reads', 'fastq', 'PAIRED', '{wildacrds.sample}', '*', '*_2.fastq.gz') + ' .' +        
                ' && zcat *_2.fastq.gz | gzip - > {wildcards.sample}.R2.fastq.gz')
        shell('mv ' + join(WORK_DIR, USER, JOB_ID, '{wildcards.sample}.R2.fastq.gz') + ' ' + join(OUT_DIR, 'Reads', 'fastq', 'PAIRED', '{wildcards.sample}', 'merged', '{wildcards.sample}.R2.fastq.gz'))
        shell('rm -r ' + join(WORK_DIR, USER, JOB_ID))


##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

## Rule to check raw SE read quality
rule fastqcSE:
    input:
        r1 = join(OUT_DIR, 'Reads', 'fastq', 'SINGLE', '{sample}', 'merged', '{sample}.R1.fastq.gz')
    output:
        r1 = join(OUT_DIR, 'fastQC', 'SINGLE', '{sample}' + '.R1_fastqc.html')
    log:
        join(OUT_DIR, 'fastQC', '{sample}' + '.fastQC_se.log')
    benchmark:
        join(OUT_DIR, 'fastQC', '{sample}' + '.fastQC_se.benchmark.tsv')
    message: 
        """--- Checking read quality of SE sample "{wildcards.sample}" with FastQC """
    run:
        if not os.path.exists(join(OUT_DIR, 'fastQC')):
            os.makedirs(join(OUT_DIR, 'fastQC'))
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
        r1 = join(OUT_DIR, 'Reads', 'fastq', 'PAIRED', '{sample}', 'merged', '{sample}.R1.fastq.gz'),
        r2 = join(OUT_DIR, 'Reads', 'fastq', 'PAIRED', '{sample}', 'merged', '{sample}.R2.fastq.gz')
    output:
        r1 = join(OUT_DIR, 'fastQC', 'PAIRED', '{sample}' + '.R1_fastqc.html'),
        r2 = join(OUT_DIR, 'fastQC', 'PAIRED', '{sample}' + '.R2_fastqc.html')
    log:
        join(OUT_DIR, 'fastQC', '{sample}' + '.fastQC_init_pe.log')
    benchmark:
        join(OUT_DIR, 'fastQC', '{sample}' + '.fastQC_init_pe.benchmark.tsv')
    message: 
        """--- Checking read quality of PE sample "{wildcards.sample}" with FastQC """
    run:
        if not os.path.exists(join(OUT_DIR, 'fastQC')):
            os.makedirs(join(OUT_DIR, 'fastQC'))
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
        r1 = join(OUT_DIR, 'Reads', 'fastq', 'SINGLE', '{sample}', 'merged', '{sample}.R1.fastq.gz')
    output:
        bam = join(OUT_DIR, 'HISAT-2', 'SINGLE', '{sample}', '{sample}' + '.csorted.bowtie2.bam')
    log:
        join(OUT_DIR, 'HISAT-2', 'SINGLE', '{sample}', 'hs2_map_se.log')
    benchmark:
        join(OUT_DIR, 'HISAT-2', 'SINGLE', '{sample}', 'hs2_map_se.benchmark.tsv')
    message: 
        """--- Mapping SE reads for sample {wildcards.sample} to genome with HISAT-2 """
    run:
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
        shell('mv ' + join(WORK_DIR, USER, JOB_ID, 'csorted.bowtie2.bam') + ' ' + join(OUT_DIR, 'HISAT-2', 'SINGLE', '{wildcards.sample}', '{wildcards.sample}' + '.csorted.bowtie2.bam'))
        shell('rm -r ' + join(WORK_DIR, USER, JOB_ID))

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

## Rule to check raw PE read quality
rule hisat2_pe_mapping:
    input:
        r1 = join(OUT_DIR, 'Reads', 'fastq', 'PAIRED', '{sample}', 'merged', '{sample}.R1.fastq.gz'),
        r2 = join(OUT_DIR, 'Reads', 'fastq', 'PAIRED', '{sample}', 'merged', '{sample}.R2.fastq.gz')
    output:
        bam = join(OUT_DIR, 'HISAT-2', 'PAIRED', '{sample}', '{sample}' + '.csorted.bowtie2.bam')
    log:
        join(OUT_DIR, 'HISAT-2', 'PAIRED', '{sample}', 'hs2_map_pe.log')
    benchmark:
        join(OUT_DIR, 'HISAT-2', 'PAIRED', '{sample}', 'hs2_map_pe.benchmark.tsv')
    message: 
        """--- Mapping PE reads for sample {wildcards.sample} to genome with HISAT-2 """
    run:
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
        shell('mv ' + join(WORK_DIR, USER, JOB_ID, 'csorted.bowtie2.bam') + ' ' + join(OUT_DIR, 'HISAT-2', 'PAIRED', '{wildcards.sample}', '{wildcards.sample}' + '.csorted.bowtie2.bam'))
        shell('rm -r ' + join(WORK_DIR, USER, JOB_ID))

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

## Rule to assemble transcripts with StringTie
rule stringtie_pe_assembly:
    input:
        bam = join(OUT_DIR, 'HISAT-2', 'PAIRED', '{sample}', '{sample}' + '.csorted.bowtie2.bam')
    output:
        asmbly = join(OUT_DIR, 'StringTie', 'PAIRED', '{sample}', '{sample}' + '.gtf')
    params:
        gtf = GTF
    log:
        join(OUT_DIR, 'StringTie', 'PAIRED', '{sample}', 'st_asmbly.log')
    benchmark:
        join(OUT_DIR, 'StringTie', 'PAIRED', '{sample}', 'st_asmbly.benchmark.tsv')
    message: 
        """--- Assembling transcripts for sample {wildcards.sample} with StringTie """
    run:
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
        shell('mv ' + join(WORK_DIR, USER, JOB_ID, '{wildcards.sample}.gtf') + ' ' + join(OUT_DIR, 'StringTie', 'PAIRED', '{wildcards.sample}'))
        shell('rm -r ' + join(WORK_DIR, USER, JOB_ID))

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

## Rule to assemble transcripts with StringTie
rule stringtie_se_assembly:
    input:
        bam = join(OUT_DIR, 'HISAT-2', 'SINGLE', '{sample}', '{sample}' + '.csorted.bowtie2.bam')
    output:
        asmbly = join(OUT_DIR, 'StringTie', 'SINGLE', '{sample}', '{sample}' + '.gtf')
    params:
        gtf = GTF
    log:
        join(OUT_DIR, 'StringTie', 'SINGLE', '{sample}', 'st_asmbly.log')
    benchmark:
        join(OUT_DIR, 'StringTie', 'SINGLE', '{sample}', 'st_asmbly.benchmark.tsv')
    message: 
        """--- Assembling transcripts for sample {wildcards.sample} with StringTie """
    run:
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
        shell('mv ' + join(WORK_DIR, USER, JOB_ID, '{wildcards.sample}.gtf') + ' ' + join(OUT_DIR, 'StringTie', 'SINGLE', '{wildcards.sample}'))
        shell('rm -r ' + join(WORK_DIR, USER, JOB_ID))

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

## Rule to merge StringTie assemblies
rule merge_assemblies:
    input:
        assemblies = expand(join(OUT_DIR, 'StringTie', '{lib}', '{sample}', '{sample}' + '.gtf'), sample = SAMPLES, lib=LIBL)
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
        shell('mkdir -p ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cp {params.gtf} ' + join(WORK_DIR, USER, JOB_ID) +
                ' && ls -1 ' + join(OUT_DIR) + '/StringTie/*/*/*.gtf > ' + join(OUT_DIR, 'StringTie', 'assemblies.txt') +
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
rule abundances_pe:
    input:
        bam = join(OUT_DIR, 'HISAT-2', 'PAIRED', '{sample}', '{sample}' + '.csorted.bowtie2.bam'),
        mrgd = rules.merge_assemblies.output.asmbly
    output:
        abundance = join(OUT_DIR, 'ballgown', 'PAIRED', '{sample}', '{sample}' + '_abundance.gtf')
    log:
        join(OUT_DIR, 'ballgown', 'PAIRED', '{sample}', 'st_abnd.log')
    benchmark:
        join(OUT_DIR, 'ballgown', 'PAIRED', '{sample}', 'st_abnd.benchmark.tsv')
    message: 
        """--- Estimating transcript abundances for sample {wildcards.sample} with StringTie"""
    run:
        shell('mkdir -p ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cp {input.bam} {input.mrgd} ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cd ' + join(WORK_DIR, USER, JOB_ID) + 
                ' && stringtie'
                ' -e -B -p 8'
                ' -G stringtie_merged.gtf' 
                ' -o {wildcards.sample}_abundance.gtf'
                ' {wildcards.sample}.csorted.bowtie2.bam'
                ' > {log}')
        shell('mv ' + join(WORK_DIR, USER, JOB_ID, '{wildcards.sample}_abundance.gtf') + ' ' + join(OUT_DIR, 'ballgown', 'PAIRED', '{wildcards.sample}'))
        shell('rm -r ' + join(WORK_DIR, USER, JOB_ID))

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

## Rule to measure transcript abundances with Stringtie
rule abundances_se:
    input:
        bam = join(OUT_DIR, 'HISAT-2', 'SINGLE', '{sample}', '{sample}' + '.csorted.bowtie2.bam'),
        mrgd = rules.merge_assemblies.output.asmbly
    output:
        abundance = join(OUT_DIR, 'ballgown', 'SINGLE', '{sample}', '{sample}' + '_abundance.gtf')
    log:
        join(OUT_DIR, 'ballgown', 'SINGLE', '{sample}', 'st_abnd.log')
    benchmark:
        join(OUT_DIR, 'ballgown', 'SINGLE', '{sample}', 'st_abnd.benchmark.tsv')
    message: 
        """--- Estimating transcript abundances for sample {wildcards.sample} with StringTie"""
    run:
        shell('mkdir -p ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cp {input.bam} {input.mrgd} ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cd ' + join(WORK_DIR, USER, JOB_ID) + 
                ' && stringtie'
                ' -e -B -p 8'
                ' -G stringtie_merged.gtf' 
                ' -o {wildcards.sample}_abundance.gtf'
                ' {wildcards.sample}.csorted.bowtie2.bam'
                ' > {log}')
        shell('mv ' + join(WORK_DIR, USER, JOB_ID, '{wildcards.sample}_abundance.gtf') + ' ' + join(OUT_DIR, 'ballgown', 'SINGLE', '{wildcards.sample}'))
        shell('rm -r ' + join(WORK_DIR, USER, JOB_ID))

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

## Rule to combine abundance counts for downstream analysis
rule collate_counts:
    input:
        abundances = expand(join(OUT_DIR, 'ballgown', '{lib}', '{sample}', '{sample}' + '_abundance.gtf'), sample = SAMPLES, lib=LIBL)
    output:
        geneCounts = join(OUT_DIR, 'ballgown', 'gene_counts.csv'),
        transcriptCounts = join(OUT_DIR, 'ballgown', 'transcript_counts.csv')
    message: 
        """--- Outputting count matrices """
    run:
        shell('prepDE.py'
                ' -i ' + join(OUT_DIR, 'ballgown') + 
                ' -g ' + join(OUT_DIR, 'ballgown', 'gene_counts.csv') +
                ' -t ' + join(OUT_DIR, 'ballgown', 'transcript_counts.csv'))

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

## Rule to collate fastQC and HISAT2 outputs with multiQC
rule multiQC:
    input:
        expand(join(OUT_DIR, 'HISAT-2', '{lib}', '{sample}', '{sample}' + '.csorted.bowtie2.bam'), sample = SAMPLES, lib = LIBL),
        expand(join(OUT_DIR, 'fastQC', 'SINGLE', '{sample}' + '.R1_fastqc.html'), sample = SRAseSAMPLES),
        expand(join(OUT_DIR, 'fastQC', 'PAIRED', '{sample}' + '.R1_fastqc.html'), sample = SRApeSAMPLES),
        expand(join(OUT_DIR, 'fastQC', 'PAIRED', '{sample}' + '.R2_fastqc.html'), sample = SRApeSAMPLES),
    output:
        join(OUT_DIR, 'MultiQC', 'multiqc_report.html')
    log:
        join(OUT_DIR, 'MultiQC', 'multiQC.log')
    benchmark:
        join(OUT_DIR, 'MultiQC', 'multiQC.benchmark.tsv')
    message: 
        """--- Running MultiQC """
    run:
        shell('ls -1 ' + join(OUT_DIR) + '/HISAT-2/*/*/*log > ' + join(OUT_DIR, 'MultiQC', 'summary_files.txt'))
        shell('ls -1 ' + join(OUT_DIR) + '/fastQC/*/*fastqc.zip >> ' + join(OUT_DIR, 'MultiQC', 'summary_files.txt'))
        shell('multiqc'
                ' -f'
                ' -o ' + join(OUT_DIR, 'MultiQC') + ' -d -dd 2 -l ' + join(OUT_DIR, 'MultiQC', 'summary_files.txt') +
                ' > {log} 2>&1')
