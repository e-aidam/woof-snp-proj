## Refence Genome

# Bowtie2 github: https://github.com/BenLangmead/bowtie2
# Downloading the Canis lupus familiaris reference genome, chromosome 31, from Ensembl database
$ wget -c https://ftp.ensembl.org/pub/release-112/fasta/canis_lupus_familiaris/dna/Canis_lupus_familiaris.ROS_Cfam_1.0.dna.primary_assembly.31.fa.gz
$ ls

# Installing bowtie to align to reference genome
$ conda install python=2.7
$ conda install -c bioconda bowtie2

# Indexing reference genome and naming reference dataset "corgi" 
$ gzip -dk Canis_lupus_familiaris.ROS_Cfam_1.0.dna.primary_assembly.31.fa.gz
$ bowtie2-build Canis_lupus_familiaris.ROS_Cfam_1.0.dna.primary_assembly.31.fa.gz corgi
