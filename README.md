# woof-snp-proj
Bionformatics pipeline to identify single-nucleotide polymorphisms (SNPs) that are associated with the development of Canine Degenerative Myelopathy using Next-Gen Sequencing Data (NGS).

For the project, the Sequence Read Archive (SRA) item with ID SRX9665373 is analyzed.(https://www.ncbi.nlm.nih.gov/sra/SRX9665373[accn]) It includes low-coverage genomic sequencing data of the Canis lupus familiaris sample (dog sample) with a size of 1.8 G bases. The data was sequenced on Illumina NovaSeq 6000 platform, using the whole genome sequencing strategy. The sequencing was Paired-End. The sample has BioSample ID SAMN16778868 , so in the BioSample database we can find more information. For example, the breed of the dog is Pembroke Welsh corgi. Run accessions, whose names start with "SRR", are used to download SRA data. (Current SRR: SRR13233412)

Pipeline Sequence:
1. Obtaining Raw Reads
2. Data Preparation
3. Refrence Genome
4. Alignment and SNP Search
5. Find the Gene
