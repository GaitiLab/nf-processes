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

## Repository structure

### Modules

Modules represent individual processes for dedicated single cell tasks, such as executing cellranger to running scrublet doublet detection on a series of matrices. They are designed to be run individually, or as part of a larger workflow/pipeline. 

Modules can be run using the following generic command: 

```nextflow run scRNA-utils/modules/{module_selection}/``` where the module_selection is the name of the specific module to be run. 

### Workflows

Workflows represent more complex and linked series of processes. Currently there is one workflow in development for split-pipe ParseBio, which in addition to the core split-pipe commands, also runs FastQC, MultiQC, and scrublet. 

## Running the split-pipe pipeline for ParseBio data

[parseBio split-pipe analysis pipeline](https://support.parsebiosciences.com/hc/en-us/categories/360004765711-Computational-Support) \

[split-pipe instructions using Nextflow](modules/split-pipe/README.md)

## Running the cellranger count pipeline for 10x Genomics scRNA data

[cellranger count pipeline for 10X scRNA](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger) \
[cellranger instructions using Nextflow](modules/cellranger/README.md)
