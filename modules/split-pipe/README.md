# Processing scRNA samples with ParseBio split-pipe in Nextflow

This module allows for the parallelization of the split-pipe scRNA processing pipeline from [ParseBio](https://support.parsebiosciences.com/hc/en-us/categories/360004765711-Computational-Support)

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


### Specifying sub-library inputs

The pipeline requires the detection of the sub-libraries corresponding to the FASTQ files for the raw data. There are two methods for providing sub-library information to the pipeline: 

#### Option 1: Sub-library list and FASTQ merging

If the user wishes to concatenate multiple FASTQ lanes OF THE SAME SUB-LIBRARY prior to running splitpipe, then the following input formats will apply: \

* A list of sub-library names contained within an input directory. These names should correspond to the FASTQ file names generated from the sequencing run. A suitable input would be as follows in the form of a .txt file: 

```
Sublibrary-1
Sublibrary-2
```
Alternatively, the user may modify the ```nextflow.config``` file and pass the sub-library names as a Groovy-compatible list of strings: 

```sublibrary = ['Sublibrary-1', 'Sublibrary-2']```

However this method is not recommended. \

Ensure that each sub library name is on its own line in the .txt file input. 

**IMPORTANT**: it is absolutely critical that sub-libraries **ARE NOT** concatenated together. Individual sub-libraries as they were prepared must be kept as distinct fastq files to ensure that the combinatorial barcoding is maintained. However, it **IS APPROPRIATE** to concatenate fastq files from the same sub-library across multiple lanes. For example, with the following fastq input files using paired end sequencing: 

```
Sublibrary-1_S1_L001_R2_001.fastq.gz  Sublibrary-1_S1_L002_R2_001.fastq.gz  Sublibrary-2_S2_L001_R2_001.fastq.gz  Sublibrary-2_S2_L002_R2_001.fastq.gz
Sublibrary-1_S1_L001_R1_001.fastq.gz  Sublibrary-1_S1_L002_R1_001.fastq.gz  Sublibrary-2_S2_L001_R1_001.fastq.gz  Sublibrary-2_S2_L002_R1_001.fastq.gz
```

An appropriate concatenation would produce the following paired fastq files: 


```
Sublibrary-1_R1_001.fastq.gz       Sublibrary-1_R2_001.fastq.gz
Sublibrary-2_R1_001.fastq.gz       Sublibrary-2_R2_001.fastq.gz
```
Note that the module will also concatenate the fastq files for you if you provide a list of the sub libraries as shown above. 

#### Option 2: No FASTQ merging

Alternatively, if the user does not wish to concatenate the FASTQ files prior to running splitpipe, the pipeline will detect the FASTQ pairs from the input directory and a FASTQ pattern. Please look at ```params.fastq_pattern``` in the ```nextflow.config``` file for an example of a Nextflow-compatible glob pattern that will identify pairs of FASTQ files for input. 

## Configuration

The following input variables can be used with the splitpipe Nextflow module: 

```
params {

input_dir = './fastq/'
output_dir = './results/'
merge_fastqs = true
ref = ''
sample_list = ''
kit = 'WT'
mode = 'all'
sublibrary = ''
combine = true
fastq_pattern = '*_R{1,1}*.fastq.gz'
}
```

***input_dir***: The absolute path of the input directory where the raw FASTQ files are held. \
***output_dir***: The absolute path of the output directory where the results are to be written. \
***merge_fastqs***: Setting to false will not concatenate the fastq files, and each pair of fastq files in the input directory will be treated as a sub-library. Default is TRUE. \
***ref***: The absolute path to a suitable reference genome. \
***sample_list***: The absolute path to a list of samples as shown above. \
***kit***: The type of kit used by ParseBio. This will correspond to the number of samples processed. The options currently available are WT_mini, WT, or WT_mega. \
***mode***: Unless custom analysis is required, the mode should always be set to 'all'. \ 
***sublibrary***: The absolute path to a list of sub library names as shown above. ONLY requires if merging fastqs from a specific input directory. \
***combine***: Whether or not the results should be combined by sub library. Default is TRUE. Setting to FALSE will keep all samples across different sub-libraries are separate outputs (NOT RECOMMENDED). \
***fastq_pattern***: A glob pattern that is combined with the input directory to detect pairs of FASTQ files for input. Only used if merge_fastqs is set to False.  

## Profiles 

Nextflow supports the use of profiles for resource allocation configuration. Currently, splitpipe with Nextflow supports two profiles: 

***local***: Recommended only for local analysis (not on a HPC cluster). This specifies 128G of memory across 8 cpus (recommended minimum requirements for splitpipe). \
***slurm_h4h***: The recommended profile for use on the UHN H4H cluster. This profile directly allows Nextflow to parallelize each process as a Slurm batch job using a veryhighmem partition with 8 cpus and 128GB of memory (time limit of 24H). 

Profiles can be specified using ```-profile profile_of_choice``` at execution. 

## Running the pipeline

The configuration variables above can either be modified directly in the nextflow.config file (NOT RECOMMENDED unless variables are to be used repeatedly), or as command line options at runtime. An example of execution would be as follows: 


```
module load ParseBiosciences/0.9.6p

./nextflow run scRNA-utils/modules/split-pipe/ --input_dir /cluster/projects/gaitigroup/ParseBio_data/fastq/ \
--output_dir /cluster/projects/gaitigroup/ParseBio_data/human_mouse \
--sample_list /cluster/projects/gaitigroup/ParseBio_data/test_HH/nextflow_test_pb_human/sample_list_human_mouse.txt \
--ref /cluster/tools/data/commondata/parsebio/genomes/hg38_mm10/ \
--kit WT_mini \
--sublibrary /cluster/projects/gaitigroup/ParseBio_data/test_HH/nextflow_test_pb_human/sublibrary.txt \
-profile slurm_h4h
```

Each input parameter can be specified using ```--parameter```. Note that profiles must be only one ```-```. 

## Outputs

The basic structure of the output directory will be as follows, with a directory named splitpipe containing the 3 sub-directories listed below: 

```
└── split_pipe
    ├── all
    ├── combined
    └── merged_fastqs
```

***merged_fastqs*** with contain symlinks to the merged FASTQ files concatenated by sub-library, if the user has enabled merging of FASTQ files. \
***all*** will contain the splitpipe output for each individual sub-library. These outputs are the outputs corresponding to running ```splitpipe --mode all``` for a set of sub-libraries. \
***combined*** will contain the spitpipe outputs for all samples combined by sub-library. These outputs correspond to the outputs from running ```splitpipe --mode comb```.






