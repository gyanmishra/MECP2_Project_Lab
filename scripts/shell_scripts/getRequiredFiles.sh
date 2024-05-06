#!/usr/bin/bash

module load bedops/2.4.14
# First create directory 
mkdir -p ../../data/mm10

## Link to Download files

# 1. Download Allen Brian Map single cell 10X metadata
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M10/gencode.vM10.annotation.gff3.gz \
-P ../../data/mm10/
gunzip ../../data/mm10/gencode.vM10.annotation.gff3.gz 

# 2. Download mm9ToMm10.over.chain.gz from UCSC
wget https://hgdownload.cse.ucsc.edu/goldenpath/mm9/liftOver/mm9ToMm10.over.chain.gz \
-P ../../data/mm10
gunzip ../../data/mm10/mm9ToMm10.over.chain.gz

# 3. Download HiC-loop from 
mkdir -p ../../data/mm10_hic_northwerstern_YueLab
wget http://3dgenome.fsm.northwestern.edu/downloads/loops-mm10.zip \
-P ../../data/mm10_hic_northwerstern_YueLab
unzip ../../data/mm10_hic_northwerstern_YueLab/loops-mm10.zip -d ../../data/mm10_hic_northwerstern_YueLab/

# 4. Download mm10 blacklisted regions
wget http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/mm10-mouse/mm10.blacklist.bed.gz \
-P ../../data/mm10
gunzip ../../data/mm10/mm10.blacklist.bed.gz

# generate geneBody bed file for gencode annotation
awk '{if($3=="gene") print $0}'  ../../data/mm10/gencode.vM10.annotation.gff3 |\
grep 'gene_type=protein_coding' | \
cut -f 1,4,5,7 | awk '{if($4 == "-") {print $1"\t"$2"\t"$3-3000} else {print $1"\t"$2+3000"\t"$3}}' | \
awk '{if($3-$2 >4500) print $0 }' | \
sort-bed - | bedops -n 1 - ../../data/mm10/mm10.blacklist.bed \
>../../data/mm10/gencode.vM10.annotation_proteinCoding_gene_greater_than_4.5kb.bed


#conda create deeptools3.1.2 
#conda activate deeptools3.1.2
#conda install deeptools==3.1.2





