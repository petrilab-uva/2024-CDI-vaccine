# Development of a Next-Generation Vaccine Against Clostridioides difficile Infection
Farha Naz, Nicholas Hagspiel, Feifan Xu, Brandon Thompson, G. Brett Moreau, Mary Young, Joel Herbein, Christopher B. Fox, William A. Petri, and Mayuresh M. Abhyankar


## Project Summary
The following repository contains the data and code used for analysis of bulk RNA sequencing data from the manuscript *Enhanced immunogenicity of a Clostridioides difficile TcdB vaccine adjuvanted with a synthetic dual-TLR ligand adjuvant*, which has been published in NPJ Vaccines. It is accessable under DOI 10.1038/s41541-025-01075-3 and can be found [here](https://pmc.ncbi.nlm.nih.gov/articles/PMC11836405/). 


## Experimental Setup
This study investigated the use of different adjuvents (alum vs a liposome-based adjuvent GLA 3M-052 LS) in conjunction with a TcdB-based *Clostridioides difficile* vaccine. Mice were vaccinated and then challenged with *C. difficile*. Colonic samples were harvested 7 days post-challenge and RNA isolated to investigate differences in transcriptional profiles across treatment groups.

Three groups were analyzed in the bulk RNA sequencing experiment.

1. **Control (n=4):** Unvaccinated mice as a control group.
2. **TcdB-Alum (n=8):** Mice vaccinated with inactivated toxin B (TcdB) with Alum as an adjuvent.
3. **TcdB-LS (n=9):** Mice vaccinated with inactivated toxin B (TcdB) with liposome adjuvant (GLA 3M-052 LS). 


## Software Implementation.
All source code required for each analysis is included within the `src` directory, while required metadata is included within the `data` directory. 

All data and source code included with this repository can be downloaded by cloning the git repository:
```
git clone https://github.com/petrilab-uva/2024-CDI-vaccine.git

```
Alternatively, you can download a [zip archive of this repository](https://github.com/petrilab-uva/2024-CDI-vaccine/archive/refs/heads/main.zip). 


## Dependencies
Source code for RNA sequencing analysis includes .R scripts as well as a .sh script, which must be run on the command line. These .sh script require a Python environment to run. This environment was set up using [Anaconda](https://www.anaconda.com/download/), which provides the `conda` package manager. Conda virtual environments were used in these analyses to install required packages (FastQC, MultiQC, BBMap, and Kallisto) in isolation. To do this, create a new conda environment using the following code, then install the required packages within this environment.
```
conda create --name ENVIRONMENT_NAME_HERE
```

Before running the .sh script within the terminal, activate the conda environment using the following code.
```
conda activate ENVIRONMENT_NAME_HERE
```


## Setup
There are several files that are required for this analysis but not included. Instructions for downloading these files are outlined below:

1. Sequencing files have been uploaded to the NCBI Sequence Read Archive (SRA) and are available under Bioproject Accession ID PRJNA1139815 (found [here](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1139815)). RNAseq FASTQ files should be placed within a `raw_reads` folder in the `data` directory.

2. Two murine genomic data sets are required for these analyses. These files can be accessed from [Ensembl.org as an FTP download](https://useast.ensembl.org/info/data/ftp/index.html). I used release version 109, but more up-to-date releases may be currently available. These versions can be used, **but ensure that the same release version is used for both files!** Instructions for obtaining these files are outlined below.
     * A *Mus musculus* cDNA file containing all transcript sequences is required to generate an index file used for read mapping. I used release 109, which is available [here](https://ftp.ensembl.org/pub/release-109/fasta/mus_musculus/cdna/). This file should be placed within the `data` directory. Code for generating the .index file from this file is commented out within the .sh script.
     * A *Mus musculus* GTF gene set is required to match transcript and gene identifiers to mapped reads. I used release 109, which is available [here](https://ftp.ensembl.org/pub/release-109/gtf/mus_musculus/). This file should be placed within the `data` directory.

  
## Contact Information
G. Brett Moreau - gbm5pn@virginia.edu

