## Obtaining Raw Reads

# Creating new virtual environment
$ conda create --name woof-snp-proj
$ conda activate woof-snp-proj

# Setting env to bash scripting: 
$ #!/usr/bin/env bash

# Checking current directory: 
$ pwd
$ mkdir /home/ethan/woof-snp-proj

# Installing sra-toolkit to ubuntu: https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit 

# Installing sra-toolkit 
$ conda install -c bioconda sra-tools

# Using prefetch command to download corgi reads from NCBI SRA
$ prefetch SRR13233412 -O /home/ethan/woof-snp-proj

# Listing files to check if SRA file was downloaded
$ ls ~/woof-snp-proj

# Change directory to current SRA accession
$ cd /home/ethan/woof-snp-proj/SRR13233412

# Extracting FASTQ data from SRA accession using fasterq-dump command
$ fasterq-dump --progress SRR13233412
$ ls 

# Compressing fastq files
$ gzip SRR13233412_1.fastq  SRR13233412_2.fastq
$ ls
