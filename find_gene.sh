# Dowloading gff file from ensemble to discover gene of interest
$ wget -c https://ftp.ensembl.org/pub/release-112/gff3/canis_lupus_familiaris/Canis_lupus_familiaris.ROS_Cfam_1.0.112.primary_assembly.31.gff3.gz
$ gzip -dk Canis_lupus_familiaris.ROS_Cfam_1.0.112.primary_assembly.31.gff3.gz
$ ls

# Upload gff3 file to Track section
# Gene of interest is SOD1-202, gene_id:ENSCAFG00845022725, description:   superoxide dismutase 1 [Source:NCBI gene (formerly Entrezgene);Acc:403559], Location   31:26,654,409-26,782,319
# Literature analysis: https://www.pnas.org/doi/10.1073/pnas.0812297106 
