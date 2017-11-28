# *Ae.aegypti* RNAseq data collection

This repository contains an analysis workflow to process Ae. aegypti RNAseq data that is available through the Sequencing Read Archive (SRA). The current workflow uses RNAseq data from the following studies (more will be added in the future):

* [Akbari *et al*. (2013)](http://www.g3journal.org/content/3/9/1493) "The Developmental Transcriptome of the Mosquito *Aedes aegypti*, an Invasive Species and Major Arbovirus Vector"
* [Alfonso-Parra *et al*. (2016)](http://journals.plos.org/plosntds/article?id=10.1371/journal.pntd.0004451) "Mating-Induced Transcriptome Changes in the Reproductive Tract of Female *Aedes aegypti*"
* [Matthews *et al*. (2016)](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-015-2239-0) "The neurotranscriptome of the *Aedes aegypti* mosquito"

The expression data for a given gene can be queried using a Shiny App. This is just a beta version and only works locally (will be ona server in the near future), so you would need to install the Shiny App package and intitate the app using R. For mroe details on Shiny Apps, visit the [Shiny website](https://shiny.rstudio.com). But here're some quick and dirty instructions for running this particular Shiny App on your computer (it is best to use RStudio, but should run on any R console):

1. Install the Shiny package on your computer:

`install.packages("shiny")`

2. Install other requried packages (ggplot2, cowplot, data.table), if not already installed:

3. Download the Shiny App file (`app.R`), the functions file (`functions.R`), and the data files and keep them in the same folder (~370MB). Those can be found [here](https://www.dropbox.com/sh/wiyv3vbc9q069ri/AACpuv4UHhyPIE9-oHcNUnz8a?dl=0)

4. Now simply open the Shiny app file, `app.R`, and run the entire script. In 10-30 seconds a window will pop up with a slot to enter gene names.

For now, the gene names you have to enter are those from the newest genome release (AaegL5).

Again, this is a work in progress. Please email suggestions, comments, issues (there should be lots of those) to yazahmed@gmail.com.
