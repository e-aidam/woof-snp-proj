## Alignment and SNP search

# Installing and configuring samtools package
$ conda install -c bioconda bowtie2 samtools
$ conda config --add channels bioconda
$ conda config --add channels conda-forge
$ conda install samtools==1.11

# Aligning pair-end reads to refrence genome and converting SAM file to BAM format.
$ bowtie2 -x corgi -1 SRR13233412_2.paired.fastq.gz -2 SRR13233412_1.paired.fastq.gz > aligned_data.sam
$ samtools view -bS aligned_data.sam > aligned_data.bam
$ rm aligned_data.sam

# Checking alginment statistics
$ samtools flagstat aligned_data.bam

# Indexing reference genome and alignment data for them to be viewed on IGV web application (https://igv.org/)
$ ls
$ samtools faidx Canis_lupus_familiaris.ROS_Cfam_1.0.dna.primary_assembly.31.fa
$ samtools sort aligned_data.bam -o aligned_sorted_data.bam
$ samtools index aligned_sorted_data.bam
$ ls

# Upload the fasta and index reference files (Canis_lupus_familiaris.ROS_Cfam_1.0.dna.primary_assembly.31.fa,Canis_lupus_familiaris.ROS_Cfam_1.0.dna.primary_assembly.31.fa.fai) to the Genome section of IGV 
# Upload sample BAM and index files (aligned_sorted_data.bam, aligned_sorted_data.bam.bai) in the Track section.
