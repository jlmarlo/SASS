# SASS

**S**treamlined **A**nimal dna-**S**equence **S**imulation

SASS is a program built to provide reproducible methods for creating artificial DNA sequencing that represents simulated genomes. This program is possible in large part thanks to the already existing **Ne**xt-generation sequencing **A**nalysis **T**oolkit which is responsible for all of the model creation and artificial sequencing read creation [NEAT](https://github.com/zstephens/neat-genreads). SASS combines the steps of model creation and read generation performed by NEAT with custom code used to select genetic variants for insertion into simulated genomes. This provides a simplified method for creating simulated genomes with preselected genetic variants. 

# Dependencies

**R libraries required to run this code:**
- ggplot2
- tidyverse

(libraries will be automatically installed if not already available)

- BCFTools
- GATK
- SAMTools

# Initial Set Up

1. Collect Necessary Files

	- Reference Genome (FASTA) with index (.fai)
	- VCF with known SNPs
	- VCF with known Indels
	- 1 or more pairs of FASTQ files
	- 1 or more BAMS paired with FASTQ files

2. Create list of empirical samples

	A text file containing a list of unique identifiers that will be used to find fastq and bam files

```
SampleA_R1.fastq.gz
SampleA_R2.fastq.gz
SampleB_R1.fastq.gz
SampleB_R2.fastq.gz
```

File should contain

```
SampleA
SampleB
...
```

3. Download SASS repository

```
git clone https://github.com/jlmarlo/SASS.git
```

4. Download NEAT repository

Original work done with SASS used version 2.0 of NEAT which is since become outdated. Future work will update SASS to work with most recent versions of NEAT

```
git clone https://github.com/zstephens/neat-genreads.git
cd neat-genreads
git checkout 8fd83e0
cd ..
```

5. Edit config.yml file

	There are several fields in the config.yml file that must be modified for the program to be fully functional

	- snpsvcf : pathway to VCF file with known SNPs
	- indelsvcf : pathway to VCF file with known indels
	- simGen : number of genomes to simulate
	- snps : total number of SNPs to be inserted into each genome (average number of SNPs that appear in species of interest)
	- indels : total number of indels to be inserted into each genome (average number of indels that appear in species of interest)
	- seed : random number from 1 - 1,000,000 to set seed for randomization
	- library : pathway to R package library where packages will be retrieved/downloaded
	- output : prefix used for output files

6. Create Environments

	The old version of neat requires python 2 to operate so a separate environment will be created for access to Python2

	If you do not have `conda` already installed it’s recommended to install via [Mambaforge](https://github.com/conda-forge/miniforge#mambaforge).

```
mamba create -c conda-forge -c bioconda -n SASSEnv bcftools gatk4 samtools

mamba create -c conda-forge -c bioconda -n Python2 python=2
```
# Running Instructions
To run SASS on a local login node 
```
conda activate SASSEnv
snakemake -s HowMakeModels.smk --configfile config.yml -c4
```

Due to the high computational demand of the read generating program it is necessary to have access to an HPC service. SASS uses snakemake to parallelize the work which does support job scheduling. Examples of snakemake profiles can be seen here [snakemake profiles(https://github.com/Snakemake-Profiles)

# Current Considerations/Future updates:
- Requires chromosomes to be formated in 'chr1, chr2' format including unnamed contigs
- Y chromosome is named MSY to work with horse reference 
- FASTQ files and bams must begin with identical sample identifiers to facilitate the programs ability to find files. 
- Incorporating slurm profile for HPCs that use the slurm job scheduler
