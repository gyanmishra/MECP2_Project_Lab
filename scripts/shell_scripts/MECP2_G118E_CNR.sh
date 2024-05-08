#!/usr/bin/bash

module load UCSC_userApps/v317
module load parallel/20150122
module load bedops/2.4.14
module load deeptools/3.5.0
module load wiggletools
module load bedtools/2.29.2


bigwig=$1 

##########################################################################################################
# Plot Average profile of MECP2 CUT&RUN (WT and G118E)

computeMatrix reference-point \
-S \
$bigwig/GSM7777160_2-IgG-2.SpikeInNormalized.bw \
$bigwig/GSM7777161_2-MeCP2-2.SpikeInNormalized.bw \
$bigwig/GSM7777162_3-IgG-2.SpikeInNormalized.bw \
$bigwig/GSM7777163_3-MeCP2-2.SpikeInNormalized.bw \
-R  ../../results/MECP2/MECP2WT.vs.KO_overlap_MECP2WT.vs.IgG.noBlacklist_filtered.bed \
--referencePoint center \
-p "max/2" \
-b 2000 -a 2000 \
--skipZeros \
--samplesLabel "IgG WT Rep1" "MECP2 WT Rep1" "IgG G118E Rep1" "MECP2 G118E Rep1" \
-o ../../results/MECP2_CNR_ChIPseq_comparison/GSE243009_MECP2_Rep1_CNR.gz

plotProfile -m ../../results/MECP2_CNR_ChIPseq_comparison/GSE243009_MECP2_Rep1_CNR.gz \
-out ../../Figures/GSE243009_MECP2_Rep1_CNR.pdf \
--perGroup \
--yAxisLabel "Read Density" \
--plotTitle "" \
--plotWidth 8 \
--yMin 0 --yMax 15 

computeMatrix reference-point \
-S \
$bigwig/GSM7777164_4-IgG-2.SpikeInNormalized.bw \
$bigwig/GSM7777165_4-MeCP2-2.SpikeInNormalized.bw \
$bigwig/GSM7777168_6-IgG-2.SpikeInNormalized.bw \
$bigwig/GSM7777169_6-MeCP2-2.SpikeInNormalized.bw \
-R  ../../results/MECP2/MECP2WT.vs.KO_overlap_MECP2WT.vs.IgG.noBlacklist_filtered.bed \
--referencePoint center \
-p "max/2" \
-b 2000 -a 2000 \
--skipZeros \
--samplesLabel "IgG WT Rep2" "MECP2 WT Rep2" "IgG G118E Rep2" "MECP2 G118E Rep2" \
-o ../../results/MECP2_CNR_ChIPseq_comparison/GSE243009_MECP2_Rep2_CNR.gz

plotProfile -m ../../results/MECP2_CNR_ChIPseq_comparison/GSE243009_MECP2_Rep2_CNR.gz \
-out ../../Figures/GSE243009_MECP2_Rep2_CNR.pdf \
--perGroup \
--yAxisLabel "Read Density" \
--plotTitle "" \
--plotWidth 8


computeMatrix reference-point \
-S \
$bigwig/GSM7777166_5-IgG-2.SpikeInNormalized.bw \
$bigwig/GSM7777167_5-MeCP2-2.SpikeInNormalized.bw \
$bigwig/GSM7777170_7-IgG-2.SpikeInNormalized.bw \
$bigwig/GSM7777171_7-MeCP2-2.SpikeInNormalized.bw \
-R  ../../results/MECP2/MECP2WT.vs.KO_overlap_MECP2WT.vs.IgG.noBlacklist_filtered.bed \
--referencePoint center \
-p "max/2" \
-b 2000 -a 2000 \
--skipZeros \
--samplesLabel "IgG WT Rep3" "MECP2 WT Rep3" "IgG G118E Rep3" "MECP2 G118E Rep3" \
-o ../../results/MECP2_CNR_ChIPseq_comparison/GSE243009_MECP2_Rep3_CNR.gz

plotProfile -m ../../results/MECP2_CNR_ChIPseq_comparison/GSE243009_MECP2_Rep3_CNR.gz \
-out ../../Figures/GSE243009_MECP2_Rep3_CNR.pdf \
--perGroup \
--yAxisLabel "Read Density" \
--plotTitle "" \
--plotWidth 8




