# Processing scRNA samples with ParseBio split-pipe in Nextflow

This module allows for the parallelization of the split-pipe scRNA processing pipeline from [ParseBio](https://support.parsebiosciences.com/hc/en-us/categories/360004765711-Computational-Support)

## Installation

Installation requires the use of Nextflow, a workflow description language (WDL) that enables reproducible parallelization of bioinformatics tasks. 

```
wget -qO- https://get.nextflow.io | bash
```
or 
```
curl -s https://get.nextflow.io | bash
```
More specific installation instructions for Nextflow can be found [here](https://www.nextflow.io/docs/latest/getstarted.html). 


The pipeline requires installation of split-pipe, provided through Parse Biosciences. On the H4H cluster maintained by UHN, the code is available as a loadable module with the following: 

```
module load ParseBiosciences/0.9.6p
```


## Pipeline installation & Updating

```
git clone https://github.com/GaitiLab/scRNA-utils.git
git checkout main
git pull
```

## Inputs

split-pipe operates on the basis of sub-libraries, which are replicates of samples processed across multiple wells. The pipeline requires as input a .txt file of sample names contained in these sub-libraries and the corresponding plate wells for each sample. A suitable .txt file will look as the following with no sample header: 

```
sample_1 A1-A3
sample_2 A4-A6
```

Ensure that there is only one space between the sample name and its corresponding well(s). 

The pipeline allows requires a .txt list of the sub-library names. These names should correspond to the fastq file names generated from the sequencing run. A suitable input would be as follows: 

```
Sublibrary-1
Sublibrary-2
```

Ensure that each sub library name is on its own line in the .txt file input. 

**IMPORTANT**: it is absolutely critical that sub-libraries **ARE NOT** concatenated together. Individual sub-libraries as they were prepared must be kept as distinct fastq files to ensure that the combinatorial barcoding is maintained. However, it **IS APPROPRIATE** to concatenate fastq files from the same sub-library across multiple lanes. 





