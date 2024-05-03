#!/usr/bin/bash


# Process 
# Boxer LD, Renthal W, Greben AW, Whitwam T et al. 
# MeCP2 Represses the Rate of Transcriptional Initiation of Highly Methylated Long Genes. 
# Mol Cell 2020 Jan 16;77(2):294-309.e9. PMID: 31784358
# GSE139509
mkdir -p ../../data/GSE139509

# Download metadata 
esearch -db sra -query PRJNA580008 | efetch -format runinfo >../../data/GSE67293/GSE139509_SraRunTable.txt

# download raw data
awk -F',' '{ print $1 }' ../../data/GSE139509/GSE139509_SraRunTable.txt  |grep -v 'Run'  | \
parallel --verbose 'prefetch {} -O ../../data/GSE67293/'

# dump fastq files
awk -F',' '{ print $1 }' ../../data/GSE139509/GSE139509_SraRunTable.txt  |grep -v 'Run'  | \
parallel --verbose 'fastq-dump --split-3 {} ../../data/GSE139509/{}/{}.sra -O ../../data/GSE139509/{}'

# Analyze raw data to generate bigWig files
perl Process_MECP2_ChIPseq.pl ../../data/GSE139509/GSE139509_SraRunTable.txt GSE139509 
wait
# wait for this program to complete