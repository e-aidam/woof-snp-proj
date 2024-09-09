## Data Preparation

# Installing fastqc and Trimmomatic
$ conda install -c bioconda fastqc
$ conda install -c bioconda trimmomatic

# Running fastqc on SRA files
$ fastqc *.fastq.gz
$ ls

# Running Trimmomatic on SRA files
$ trimmomatic PE -phred33 SRR13233412_1.fastq.gz SRR13233412_2.fastq.gz SRR13233412_1.paired.fastq.gz SRR13233412_1.unpaired.fastq.gz SRR13233412_2.paired.fastq.gz SRR13233412_2.unpaired.fastq.gz ILLUMINACLIP:NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:30
