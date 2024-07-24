# QUALITY CONTROL AND PSEUDOALIGNMENT OF RNA-SEQ FASTQs USING KALLISTO
# NAME: G. Brett Moreau
# DATE: July 8, 2024


# INTRODUCTION -----------------------------------------------------------------

# This script will perform QC and map reads from bulk RNA sequencing FASTQs. QC 
# will be evaluated using FastQC, while reads will be pseudoaligned to the mouse 
# genome using Kallisto. All QC information will be summarized into a QC report 
# using MultiQC.

# Before running this script, raw FASTQ files need to be downloaded and placed 
# into a "raw_reads" subdirectory within the "data" directory. This script 
# should be run using the command line in a conda enviroment containing the 
# packages listed below, which will be used for this script.

# Permissions need to be granted to allow this script to be loaded using the 
# "chmod +x ./01-fastq_processing_and_readMapping.sh" command (no quotes).


# PACKAGE VERSIONS USED --------------------------------------------------------
# The following package versions were use for this analysis.

### python: Version 3.8.16
### conda: Version 23.5.0
### fastqc: Version 0.12.1
### multiqc: Version 1.14
### bbmap: Version 38.18
### kallisto: Version 0.44.0


# SETTING UP THE ENVIRONMENT ---------------------------------------------------

# Change to the project directory (change as necessary)
cd ~/Documents/Github/2024-CDI-Mouse-Vaccine

# Make output directory for fastqc and kallisto report files
mkdir ./results/fastqc
mkdir ./results/fastqc/raw_reads
mkdir ./data/trimmed_reads
mkdir ./results/fastqc/trimmed_reads
mkdir ./results/kallisto
mkdir ./results/kallisto/logs


# QUALITY CONTROL OF PROCESSED FASTQS USING FASTQC -----------------------------

# Run FastQC on processed FASTQ files.
for i in ./data/raw_reads/*fq.gz
do
	fastqc -o ./results/fastqc/raw_reads $i -t 16
done

# Summarize FastQC results using MultiQC
multiqc ./results \
-d -n report_raw_reads.html \
-o ./results


# ADAPTER SEQUENCE TRIMMING USING BBDUK ----------------------------------------

# Trim adapter sequences using BBDuk
for i in $(basename ./data/raw_reads/*.fq.gz | cut -f1,1 -d'_' | uniq)
do 
	bbduk.sh \
	in1=./data/raw_reads/"$i"_1.fq.gz \
	in2=./data/raw_reads/"$i"_2.fq.gz \
	out1=./data/trimmed_reads/"$i"_1_trimmed.fq.gz \
	out2=./data/trimmed_reads/"$i"_2_trimmed.fq.gz \
	ref=./data/adapters.fa \
	ktrim=r \
	mink=11 \
	hdist=1 \
	tpe tbo \

done 2> ./results/bbduk_output.log


# Run FastQC on trimmed FASTQ files
for i in ./data/trimmed_reads/*fq.gz
do
	fastqc -o ./results/fastqc/trimmed_reads $i -t 16
done


# READ MAPPING WITH KALLISTO ---------------------------------------------------

# Build cDNA index from reference FASTA file for aligning reads. I have 
# previously generated this file, but if you haven't run the following lines of
# code.

#kallisto index -i ./data/Mus_musculus.GRCm39.cdna.all.index \
#./data/Mus_musculus.GRCm39.cdna.all.fa

# NOTE: I'm using release 109 of the Mus musculus cDNA FASTA file, which can be
# accessed from the ENSEMBL depository: 
# (https://useast.ensembl.org/info/data/ftp/index.html).


for i in $(basename ./data/trimmed_reads/*.fq.gz | cut -f1,1 -d'_' | uniq)
do
	kallisto quant -i ./data/Mus_musculus.GRCm39.cdna.all.index \
	-o ./results/kallisto/"$i" -t 16 \
	./data/trimmed_reads/"$i"_1_trimmed.fq.gz \
	./data/trimmed_reads/"$i"_2_trimmed.fq.gz \
	&> ./results/kallisto/logs/"$i".log
	echo $i "finished"
done	


# SUMMARIZE REPORTS USING MULTIQC ----------------------------------------------

# Summarize using MultiQC
multiqc ./results \
--ignore ./results/fastqc/raw_reads \
-d -n report_trimmed_reads_and_mapping.html \
-o ./results

echo "Finished"
