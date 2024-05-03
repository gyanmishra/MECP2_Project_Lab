#!/usr/bin/bash

# Reserves  a node before running below code.
# srun --partition=super --nodes=8 --pty --time=5:00:00 /bin/bash

# Load below modules to execute the commands. 
module load parallel/20150122
module load bedtools/2.29.2
module load deeptools/3.5.0
module load bedops/2.4.14
module load UCSC_userApps/v317



######################################################################################################################
# Figure 3B
# Overlap MECP2 peaks with H3K27ac peaks
bedtools intersect -wa -b \
../../results/MECP2/MECP2WT.vs.KO_overlap_MECP2WT.vs.IgG.noBlacklist_filtered.bed \
-a ../../results/H3K27ac/H3K27ac_MECP2WT_filtered_peaks.bed |bedops -m - \
>../../results/MECP2_H3K27ac_overlap/MECP2_H3K27ac_overlap.bed

bedtools intersect -v -wa -b \
../../results/MECP2/MECP2WT.vs.KO_overlap_MECP2WT.vs.IgG.noBlacklist_filtered.bed \
-a ../../results/H3K27ac/H3K27ac_MECP2WT_filtered_peaks.bed |bedops -m - \
>../../results/MECP2_H3K27ac_overlap/MECP2_H3K27ac_H3K27ac_noOverlap.bed

bedtools intersect -v -wa -a \
../../results/MECP2/MECP2WT.vs.KO_overlap_MECP2WT.vs.IgG.noBlacklist_filtered.bed \
-b ../../results/H3K27ac/H3K27ac_MECP2WT_filtered_peaks.bed |bedops -m - \
>../../results/MECP2_H3K27ac_overlap/MECP2_H3K27ac_MECP2_noOverlap.bed

######################################################################################################################
# Figure 3C, Supplementary Figure 3C
bedtools intersect -wa -a \
../../results/MECP2/MECP2WT.vs.KO_overlap_MECP2WT.vs.IgG.noBlacklist_filtered.bed \
-b ../../results/H3K27ac/H3K27ac_MECP2_up_MECP2WT.vs.KO.bed |bedops -m - \
>../../results/MECP2_H3K27ac_overlap/MECP2_H3K27ac_up.bed

bedtools intersect -wa -a \
../../results/MECP2/MECP2WT.vs.KO_overlap_MECP2WT.vs.IgG.noBlacklist_filtered.bed \
-b ../../results/H3K27ac/H3K27ac_MECP2_unchanged_MECP2WT.vs.KO.bed |bedops -m - \
>../../results/MECP2_H3K27ac_overlap/MECP2_H3K27ac_unchanged.bed

bedtools intersect -wa -a \
../../results/MECP2/MECP2WT.vs.KO_overlap_MECP2WT.vs.IgG.noBlacklist_filtered.bed \
-b ../../results/H3K27ac/H3K27ac_MECP2_down_MECP2WT.vs.KO.bed |bedops -m - \
>../../results/MECP2_H3K27ac_overlap/MECP2_H3K27ac_down.bed

ls ../../results/MECP2_H3K27ac_overlap/MECP2_H3K27ac_*.bed | \
parallel --verbose 'computeMatrix reference-point \
-S ../../../CUT_and_RUN/results/BigWig/7wk_CTX_MECP2WT_Mecp2_CNR_HS031219_N031419.bw \
../../../CUT_and_RUN/results/BigWig/7wk_CTX_MECP2KO_Mecp2_CNR_HS031219_N031419.bw \
../../../CUT_and_RUN/results/BigWig/7wk_CTX_MECP2WT_H3K27ac_CNR_HS031219_N031419.bw
../../../CUT_and_RUN/results/BigWig/9wk_CTX_MECP2WT_H3K27ac_CNR_HS032519_N032619.bw \
../../../CUT_and_RUN/results/BigWig/7wk_CTX_MECP2KO_H3K27ac_CNR_HS031219_N031419.bw \
../../../CUT_and_RUN/results/BigWig/9wk_CTX_MECP2KO_H3K27ac_CNR_HS032519_N032619.bw \
-R {} \
--referencePoint center \
       -p 20 \
       -b 2000 -a 2000 \
       --skipZeros \
       --missingDataAsZero \
       -o {.}.gz'

ls ../../results/MECP2_H3K27ac_overlap/MECP2_H3K27ac_*.bed |\
parallel --verbose 'plotHeatmap -m {.}.gz \
--missingDataColor 'white' \
-out {.}.pdf \
--colorMap RdYlBu \
--samplesLabel "7wk MECP2 WT" "7wk MECP2 KO" \
               "9wk MECP2 WT (H3K27ac)" "9wk MECP2 KO (H3K27ac)" \
--heatmapHeight 3'


######################################################################################################################
# Supplementary Figure 3B
cat ../../results/MECP2/MECP2WT.vs.KO_overlap_MECP2WT.vs.IgG.noBlacklist_filtered.bed | \
grep -v 'random' |grep -v 'Un' \
>../../results/MECP2/MECP2WT.vs.KO_overlap_MECP2WT.vs.IgG.noBlacklist_filtered_chr1to19XY.bed

ls ../../data/ENCODE/*.gz | parallel --verbose 'gunzip {}'


