# scRNA-utils
modules and utility scripts for processing scRNA data

## Installation

Installation requires the use of Nextflow, a workflow description language (WDL) that enables reproducible parallelization of common bioinformatics tasks. Nextflow provides an executable that requires both [Groovy](https://groovy-lang.org/) and [Java/JDK](https://www.oracle.com/java/technologies/downloads/). Installation of the portable executable is as follows:

```
wget -qO- https://get.nextflow.io | bash
```
or 
```
curl -s https://get.nextflow.io | bash
```
More specific installation instructions for Nextflow can be found [here](https://www.nextflow.io/docs/latest/getstarted.html).

## Pipeline installation & Updating

```
git clone https://github.com/GaitiLab/scRNA-utils.git
git checkout main
git pull
```

## Modules

Modules represent individual processes for dedicated single cell tasks, such as executing cellranger to running scrublet doublet detection on a series of matrices. They are designed to be run individually, or as part of a larger workflow/pipeline. 

Modules can be run using the following generic command: 

```nextflow run scRNA-utils/modules/{module_selection}/``` where the module_selection is the name of the specific module to be run. Currently the available modules can be found in the `modules` directory: 

* ***cellranger***: runs `cellranger count` on either a directory (recursive or not) of FASTQ files, or a sample sheet with samples and their file outputs specified. [See below](#running-the-cellranger-count-pipeline-for-10x-genomics-scRNA-data) for more information. 
* ***epiAneuFinder***: NOTE: **experimental**: currently not maintained. Runs an experimental scATAC-seq CNV caller on count matrices using the epiAneuFinder R package. 
* ***fastqc_multiqc***: Given a directory (recursive or not) of FASTQ files, run FastQC and multiQC (optional) on the files. 
* ***kb-python***: NOTE: **experimental**: currently not maintained. Runs kb-python (kallisto bustools) on a set of FASTQ files. 
* ***scrublet***: NOTE: **experimental**: currently not maintained. Runs scrublet for doublet detection on a count matrix output of either split-pipe or cellranger count. 
* ***split-pipe***: Process ParseBio FASTQ files into count matrices and alignment files using `split-pipe --mode all` and/or `split-pipe --mode comb`. [See below](#running-the-split-pipe-pipeline-for-parsebio-data) for more information. 


### Running individual scRNA processing pipelines from modules

Within the **modules** directory are two basic pipelines for processing scRNA-seq data from raw FASTQ files: 

* ParseBio data (to be processed using split-pipe)
* 10X Genomics data (to be processed using cellranger)

Below are the links to the specific user documentation for each type of scRNA-seq data. 

#### Running the split pipe pipeline for ParseBio data

[parseBio split-pipe analysis pipeline](https://support.parsebiosciences.com/hc/en-us/categories/360004765711-Computational-Support) \
[split-pipe instructions using Nextflow](modules/split-pipe/README.md)

#### Running the cellranger count pipeline for 10x Genomics scRNA data

[cellranger count pipeline for 10X scRNA](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger) \
[cellranger instructions using Nextflow](modules/cellranger/README.md)

## Workflows

Workflows represent more complex and linked series of processes. Currently there is one workflow in development for toggling between both ParseBio and 10X scRNA data. The workflow can be found in `workflows` and enables the following behaviour: 

* Specifying the type of input scRNA data with `--method` as either `split-pipe` or `cellranger`.
* For either mode, the pipeline will generate count matrices from FASTQ files, run FastQC (and MultiQC optionally), and run scrublet on the filtered output count matrices. 

### Running combined/toggled scRNA processing pipelines from workflows

**Currently under development**. 

