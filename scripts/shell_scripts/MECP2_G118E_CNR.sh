#!/bin/bash

#SBATCH -p 512GB                      # partition name
#SBATCH -c 8                         # Request N cores 
#SBATCH -t 0-16:00:00                 # hours:minutes time after which job will be killed 
#SBATCH --mem 501760                  # total amount of memory requested [>100G recommended]
#SBATCH --job-name SRARun             # Job name
#SBATCH -o ../../results/log/G118E.out    # File to which standard out will be written
#SBATCH -e ../../results/log/G118E.err    # File to which standard err will be written


module load UCSC_userApps/v317
module load parallel/20150122
module load bedops/2.4.14
module load deeptools/3.5.0
module load wiggletools
module load bedtools/2.29.2

#ls ../../data/GSE243009/*.bw | parallel --verbose --jobs 5 'bigWigToWig {} {.}.wig'
#ls ../../data/GSE243009/*.wig | parallel --verbose --jobs 5 'wig2bed --zero-indexed < {} > {.}.bed'
#ls ../../data/GSE243009/*.wig | parallel --verbose 'rm {}'

# Download processed bigWig file from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE243009
# using R code given in markdown file 

# GSM7777160	WT IgG rep 1
# GSM7777161	WT MeCP2 rep 1
bigwigCompare \
--bigwig1 ../../data/GSE243009/GSM7777161_2-MeCP2-2.SpikeInNormalized.bw \
--bigwig2 ../../data/GSE243009/GSM7777160_2-IgG-2.SpikeInNormalized.bw \
--outFileName  ../../data/GSE243009/GSM7777161_2-MeCP2-2_over_GSM7777160_2-IgG-2.SpikeInNormalized.bw \
-p "max/2"

# GSM7777162	G118E IgG rep 1
# GSM7777163	G118E MeCP2 rep 1
bigwigCompare \
--bigwig1 ../../data/GSE243009/GSM7777163_3-MeCP2-2.SpikeInNormalized.bw \
--bigwig2 ../../data/GSE243009/GSM7777162_3-IgG-2.SpikeInNormalized.bw \
--outFileName  ../../data/GSE243009/GSM7777163_3-MeCP2-2_over_GSM7777162_3-IgG-2.SpikeInNormalized.bw \
-p "max/2"

# G118E vs WT (Rep1)
bigwigCompare \
--bigwig1 ../../data/GSE243009/GSM7777163_3-MeCP2-2_over_GSM7777162_3-IgG-2.SpikeInNormalized.bw \
--bigwig2 ../../data/GSE243009/GSM7777161_2-MeCP2-2_over_GSM7777160_2-IgG-2.SpikeInNormalized.bw \
--outFileName ../../data/GSE243009/G118E_vs_IgG_Rep1.bw \
-p "max/2" \
--operation subtract

# GSM7777164	WT IgG rep 2
# GSM7777165	WT MeCP2 rep 2
bigwigCompare \
--bigwig1 ../../data/GSE243009/GSM7777165_4-MeCP2-2.SpikeInNormalized.bw \
--bigwig2 ../../data/GSE243009/GSM7777164_4-IgG-2.SpikeInNormalized.bw \
--outFileName  ../../data/GSE243009/GSM7777165_4-MeCP2-2_over_GSM7777164_4-IgG-2.SpikeInNormalized.bw \
-p "max/2"

# GSM7777168	G118E IgG rep 2
# GSM7777169	G118E MeCP2 rep 2
bigwigCompare \
--bigwig1 ../../data/GSE243009/GSM7777169_6-MeCP2-2.SpikeInNormalized.bw \
--bigwig2 ../../data/GSE243009/GSM7777168_6-IgG-2.SpikeInNormalized.bw \
--outFileName  ../../data/GSE243009/GSM7777169_6-MeCP2-2_over_GSM7777168_6-IgG-2.SpikeInNormalized.bw \
-p "max/2"

# G118E vs WT (Rep2)
bigwigCompare \
--bigwig1 ../../data/GSE243009/GSM7777169_6-MeCP2-2_over_GSM7777168_6-IgG-2.SpikeInNormalized.bw \
--bigwig2 ../../data/GSE243009/GSM7777165_4-MeCP2-2_over_GSM7777164_4-IgG-2.SpikeInNormalized.bw \
--outFileName ../../data/GSE243009/G118E_vs_IgG_Rep2.bw \
-p "max/2"  \
--operation subtract


# GSM7777166	WT IgG rep 3
# GSM7777167	WT MeCP2 rep 3
bigwigCompare \
--bigwig1 ../../data/GSE243009/GSM7777167_5-MeCP2-2.SpikeInNormalized.bw \
--bigwig2 ../../data/GSE243009/GSM7777166_5-IgG-2.SpikeInNormalized.bw \
--outFileName  ../../data/GSE243009/GSM7777167_5-MeCP2-2_over_GSM7777166_5-IgG-2.SpikeInNormalized.bw \
-p "max/2"

# GSM7777170	G118E IgG rep 3
# GSM7777171	G118E MeCP2 rep 3
bigwigCompare \
--bigwig1 ../../data/GSE243009/GSM7777171_7-MeCP2-2.SpikeInNormalized.bw \
--bigwig2 ../../data/GSE243009/GSM7777170_7-IgG-2.SpikeInNormalized.bw \
 --outFileName  ../../data/GSE243009/GSM7777171_7-MeCP2-2_over_GSM7777170_7-IgG-2.SpikeInNormalized.bw \
-p "max/2"


# G118E vs WT (Rep3)
bigwigCompare \
--bigwig1 ../../data/GSE243009/GSM7777171_7-MeCP2-2_over_GSM7777170_7-IgG-2.SpikeInNormalized.bw \
--bigwig2 ../../data/GSE243009/GSM7777167_5-MeCP2-2_over_GSM7777166_5-IgG-2.SpikeInNormalized.bw  \
--outFileName ../../data/GSE243009/G118E_vs_IgG_Rep3.bw \
-p "max/2" \
--operation subtract

# convert BigWig to bedGraph
ls ../../data/GSE243009/G118E_vs_IgG_Rep*.bw | parallel --verbose  'bigWigToBedGraph {} {.}.bedGraph'
ls ../../data/GSE243009/*_over_*.bw | parallel --verbose  'bigWigToBedGraph {} {.}.bedGraph'


# generate geneBody bed file for gencode annotation
awk '{if($3=="gene") print $0}'  ../../../CUT_and_RUN/data/mm10/gencode.vM10.annotation.gff3 |\
grep 'gene_type=protein_coding' | \
cut -f 1,4,5,7 | awk '{if($4 == "-") {print $1"\t"$2"\t"$3-3000} else {print $1"\t"$2+3000"\t"$3}}' | \
awk '{if($3-$2 >4500) print $0 }' | \
sort-bed - | bedops -n 1 - ../../data/mm10/mm10.blacklist.bed \
>../../data/mm10/gencode.vM10.annotation_proteinCoding_gene_greater_than_4.5kb.bed


# generate GeneBody Coverage from bedGraph files. 
ls ../../data/GSE243009/*_over_*.bedGraph |\
parallel --verbose 'bedtools map -a \
../../data/mm10/gencode.vM10.annotation_proteinCoding_gene_greater_than_4.5kb.bed \
-b {} -c 4 -o sum \
>../../results/CNR_Coverage/{/.}.bed'

ls ../../data/GSE243009/G118E_vs_IgG_Rep*.bedGraph | parallel --verbose 'grep "chr" {} >{.}.chr.bedGraph'
# generate GeneBody Coverage from bedGraph files (G118E vs WT). 
ls ../../data/GSE243009/G118E_vs_IgG_Rep*.chr.bedGraph |\
parallel --verbose 'bedtools map -a \
../../data/mm10/gencode.vM10.annotation_proteinCoding_gene_greater_than_4.5kb.bed \
-b {} -c 4 -o sum \
>../../results/CNR_Coverage/{/.}.bed'

##########################################################################################################
# Plot Average profile of MECP2 CUT&RUN (WT and G118E)

computeMatrix reference-point \
-S \
../../data/GSE243009/GSM7777160_2-IgG-2.SpikeInNormalized.bw \
../../data/GSE243009/GSM7777161_2-MeCP2-2.SpikeInNormalized.bw \
../../data/GSE243009/GSM7777162_3-IgG-2.SpikeInNormalized.bw \
../../data/GSE243009/GSM7777163_3-MeCP2-2.SpikeInNormalized.bw \
-R  ../../results/MECP2/MECP2WT.vs.KO_overlap_MECP2WT.vs.IgG.noBlacklist_filtered.bed \
--referencePoint center \
-p "max/2" \
-b 2000 -a 2000 \
--skipZeros \
--samplesLabel "IgG WT Rep1" "MECP2 WT Rep1" "IgG G118E Rep1" "MECP2 G118E Rep1" \
-o ../../results/MECP2_CNR_ChIPseq_comparison/GSE243009_MECP2_Rep1_CNR.gz

plotProfile -m ../../results/MECP2_CNR_ChIPseq_comparison/GSE243009_MECP2_Rep1_CNR.gz \
-out ../../results/MECP2_CNR_ChIPseq_comparison/GSE243009_MECP2_Rep1_CNR.pdf \
--perGroup \
--yAxisLabel "Read Density" \
--plotTitle "" \
--plotWidth 8 \
--yMin 0 --yMax 15 

computeMatrix reference-point \
-S \
../../data/GSE243009/GSM7777164_4-IgG-2.SpikeInNormalized.bw \
../../data/GSE243009/GSM7777165_4-MeCP2-2.SpikeInNormalized.bw \
../../data/GSE243009/GSM7777168_6-IgG-2.SpikeInNormalized.bw \
../../data/GSE243009/GSM7777169_6-MeCP2-2.SpikeInNormalized.bw \
-R  ../../results/MECP2/MECP2WT.vs.KO_overlap_MECP2WT.vs.IgG.noBlacklist_filtered.bed \
--referencePoint center \
-p "max/2" \
-b 2000 -a 2000 \
--skipZeros \
--samplesLabel "IgG WT Rep2" "MECP2 WT Rep2" "IgG G118E Rep2" "MECP2 G118E Rep2" \
-o ../../results/MECP2_CNR_ChIPseq_comparison/GSE243009_MECP2_Rep2_CNR.gz

plotProfile -m ../../results/MECP2_CNR_ChIPseq_comparison/GSE243009_MECP2_Rep2_CNR.gz \
-out ../../results/MECP2_CNR_ChIPseq_comparison/GSE243009_MECP2_Rep2_CNR.pdf \
--perGroup \
--yAxisLabel "Read Density" \
--plotTitle "" \
--plotWidth 8


computeMatrix reference-point \
-S \
../../data/GSE243009/GSM7777166_5-IgG-2.SpikeInNormalized.bw \
../../data/GSE243009/GSM7777167_5-MeCP2-2.SpikeInNormalized.bw \
../../data/GSE243009/GSM7777170_7-IgG-2.SpikeInNormalized.bw \
../../data/GSE243009/GSM7777171_7-MeCP2-2.SpikeInNormalized.bw \
-R  ../../results/MECP2/MECP2WT.vs.KO_overlap_MECP2WT.vs.IgG.noBlacklist_filtered.bed \
--referencePoint center \
-p "max/2" \
-b 2000 -a 2000 \
--skipZeros \
--samplesLabel "IgG WT Rep3" "MECP2 WT Rep3" "IgG G118E Rep3" "MECP2 G118E Rep3" \
-o ../../results/MECP2_CNR_ChIPseq_comparison/GSE243009_MECP2_Rep3_CNR.gz

plotProfile -m ../../results/MECP2_CNR_ChIPseq_comparison/GSE243009_MECP2_Rep3_CNR.gz \
-out ../../results/MECP2_CNR_ChIPseq_comparison/GSE243009_MECP2_Rep3_CNR.pdf \
--perGroup \
--yAxisLabel "Read Density" \
--plotTitle "" \
--plotWidth 8




