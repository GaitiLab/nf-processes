# Processing 10x cellranger samples in Nextflow

This module allows for the parallelization of the cellranger scRNA processing pipeline from [10x Genomics](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger)

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


The pipeline requires installation of cellranger, provided through 10X Genomics. On the H4H cluster maintained by UHN, the code is available as a loadable module with the following: 

```
module load cellranger
```

## Pipeline installation & Updating

```
git clone https://github.com/GaitiLab/scRNA-utils.git
git checkout main
git pull
```

## Inputs

### Sample Input Option #1: Directory

zthere are two possibilities for specifying the input FASTQ files required to run cellranger count. The first allows for specifying a directory that contains 4 FASTQS for every sample, with multiple samples permitted:

* Index 1, as {sample}_I1_001.fastq.gz
* Index 2, as {sample}_I2_001.fastq.gz
* Read 1, {sample}_R1_001.fastq.gz
* Read 2, {sample}_R1_001.fastq.gz

From thes 4 FASTQ files, cellranger is able to identify the sample name (here given as {sample}).

When running cellranger in Nextflow, the sample name will be recognized automatically through a regex pattern. For example, when the following FASTQ files: 

```
144556_S52_L001_I1_001.fastq.gz
144556_S52_L001_I2_001.fastq.gz
144556_S52_L001_R1_001.fastq.gz
144556_S52_L001_R2_001.fastq.gz
144556_S52_L002_I1_001.fastq.gz
144556_S52_L002_I2_001.fastq.gz
144556_S52_L002_R1_001.fastq.gz
144556_S52_L002_R2_001.fastq.gz
```

The sample name will be derived as 144556 (the filenames are split at _S). \
**Note**: FASTQ files that correspond to the same sample, but across multiple lanes, will be collapsed together. In the example above, 144556 is apread out across 2 lanes, and the resulting analysis will combine the FASTQ files for these 2 lanes into one output directory automatically by cellranger, as long as the portion of the FASTQ names before the lane assignment (such as 144556_S52_) is the same.

### Sample Input Option #2: Sample sheet

The second option for specifying input FASTQs is to use a sample sheet with two columns named **sample_name** and **sample_path**, as seen below: 

```
sample_name,sample_path
610-7,/cluster/projects/gaitigroup/ParseBio_data/220215_A00827_0509_BHVT53DSX2_NU210019/fastq/
331R-11,/cluster/projects/gaitigroup/ParseBio_data/220215_A00827_0509_BHVT53DSX2_NU210019/fastq/
```

This permits the running of samples in different directories. Note that these directories must contain the 4 required FASTQ files as described above with the same naming conventions. 

## Configuration

The following input variables can be used with the splitpipe Nextflow module. Default values are set that can be changed at runtime: 

```
params {

input_dir = '.input/'
output_dir = './results'
resursive_search = false

cellranger.fastq_pattern = '*_R{1,2}*.fastq.gz'
cellranger.sample_name_split = '_S'
cellranger.ref = ''
cellranger.sample_sheet = ''
cellranger.include_introns = true
cellranger.expected_cells = 3000
}

}
```

***input_dir***: The absolute path of the input directory where the raw FASTQ files are held. All FASTQ files should be contained in this directory unless **recursive_search** is set to true. \
***recursive_search***: Whether or not the FASTQ file search should be recursive for the input directory. If set to yes, the runner will search recursively for any paired FASTQ files in the main durectory and all sub-directories specified by ***input_dir**. Default is false. \
***output_dir***: The absolute path of the output directory where the results are to be written. The module will create the output directory if it does not exist. \
***cellranger.fastq_pattern***: A glob pattern that is combined with the input directory to detect pairs of FASTQ files for input. This pattern is used to get the final sample names. \
***cellranger.sample_name_split***: The FASTQ name pattern at which the split should occur to derive the unique sample name. The default is set to **_S**, which for Illumina sequencing corresponds to the beginning of the index name. Splitting at the index name makes the unique name compatible with cellranger. \
***cellranger.ref***: A compatible RNA-based reference for running cellranger. The HPC Bioinformatics core is able to provide the path of various reference genomes at `/cluster/tools/data/commondata/cellranger/`. \
***cellranger.sample_sheet***: An alternative way of specifying the input files and sample names using a sheet. Overrides the input_dir variable if both are specified. \
***cellranger.include_introns***: Whether or not cellranger should evaluate reads that align to introns or not. Unless specifically stated, this parameter should be set to true. \
***cellranger.expected_cells***: The expected number of cells to retain from the experiment. Typical values are set between 2000-3000.

## Profiles

Nextflow supports the use of profiles for resource allocation configuration. Currently, splitpipe with Nextflow supports two profiles: 

***local***: Recommended only for local analysis (not on a HPC cluster). This specifies 128G of memory across 8 cpus, which are the recommended minimum requirements for splitpipe. \
***slurm_h4h***: The recommended profile for use on the UHN H4H cluster. This profile directly allows Nextflow to parallelize each process as a Slurm batch job using a veryhighmem partition with 8 cpus and 128GB of memory (time limit of 24H). 

Profiles can be specified using ```-profile profile_of_choice``` at execution. Note that multiple configuration profiles can be used at once.

## Running the pipeline

The configuration variables above can either be modified directly in the nextflow.config file (this method is **NOT RECOMMENDED** unless variables are to be used repeatedly), or as command line options at runtime. An example of execution would be as follows for the H4H cluster: 


```

./nextflow run scRNA-utils/modules/cellranger/ \
--input_dir  /cluster/projects/gaitigroup/ParseBio_data/220215_A00827_0509_BHVT53DSX2_NU210019/fastq/ \
--output_dir /cluster/projects/gaitigroup/ParseBio_data/test_cellranger/cellranger/ \
--cellranger.ref /cluster/tools/data/commondata/cellranger/refdata-gex-GRCh38-2020-A/ \
-profile slurm_h4h
```

Each input parameter can be specified using ```--parameter```. Note that profiles must be only one ```-```. 

The command above can be submitted through a slurm job with the following code preceeding it: 

```
#!/bin/bash

#SBATCH -t 23:00:00
#SBATCH --mem=4G
#SBATCH -J cellranger
#SBATCH -p all
#SBATCH -c 2
#SBATCH -N 1     
#SBATCH -o %x-%j.out
```

and the file can be submitted using `sbatch`.



## Outputs

The basic structure of the output directory will be as follows, with a directory named cellranger containing the 3 sub-directories listed below: 

```
cellranger_count/
├── 331R-11
│   ├── 331R-11.mri.tgz
│   ├── _cmdline
│   ├── _filelist
│   ├── _finalstate
│   ├── _invocation
│   ├── _jobmode
│   ├── _log
│   ├── _mrosource
│   ├── outs
│   ├── _perf
│   ├── SC_RNA_COUNTER_CS
│   ├── _sitecheck
│   ├── _tags
│   ├── _timestamp
│   ├── _uuid
│   ├── _vdrkill
│   └── _versions
└── 610-7
    ├── 610-7.mri.tgz
    ├── _cmdline
    ├── _filelist
    ├── _finalstate
    ├── _invocation
    ├── _jobmode
    ├── _log
    ├── _mrosource
    ├── outs
    ├── _perf
    ├── SC_RNA_COUNTER_CS
    ├── _sitecheck
    ├── _tags
    ├── _timestamp
    ├── _uuid
    ├── _vdrkill
    └── _versions
```

In this example there are two samples in the cellranger_count sub-directory: 331R-11 and 610-7. A sub-directory for each sample, by name, will be created in the output directory.