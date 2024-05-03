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
# Figure 2C, Supplementary figure 2E and 2G
# remove blacklisted regions and select only chr1-19,X,Y
ls ../../results/MECP2/MECP2_u*.bed |parallel --verbose 'bedops -n 1 {} ../../data/mm10/mm10.blacklist.bed |\
bedops -m - | grep -v 'chrUn' - |grep -v 'random' >{.}.noBlacklist.bed'

ls ../../results/MECP2/MECP2WT.vs.KO_overlap_MECP2WT.vs.IgG.noBlacklist_filtered_DNMT3AKO_vs_WT.*sort.bed| \
parallel --verbose 'computeMatrix reference-point \
-S  ../../../CUT_and_RUN/results/BigWig/9wk_CTX_DNMT3AWT_Mecp2_CNR_HS031219_N031419.bw \
    ../../../CUT_and_RUN/results/BigWig/Adult_CTX_WT_Mecp2_CNR_HS020619_N020719.bw \
    ../../../CUT_and_RUN/results/BigWig/9wk_CTX_DNMT3AKO_Mecp2_CNR_HS031219_N031419.bw \
    ../../../CUT_and_RUN/results/BigWig/Adult_CTX_DNMT3AKO_Mecp2_CNR_HS020619_N020719.bw \
       -R {} \
       --referencePoint center\
       -p 20 \
       -b 2000 -a 2000 \
       --missingDataAsZero \
       --skipZeros \
       -o {.}.gz'

# pdf output
ls ../../results/MECP2/MECP2WT.vs.KO_overlap_MECP2WT.vs.IgG.noBlacklist_filtered_DNMT3AKO_vs_WT.*sort.gz |\
parallel --verbose 'plotHeatmap -m {} --missingDataColor 'white' -out {.}.pdf --colorMap RdYlBu \
--samplesLabel "9wk DNMT3A WT" "9wk DNMT3A WT" \
               "Adult DNMT3A KO" "Adult DNMT3A KO" \
--heatmapHeight 3' 


# calculate MECP2 coverage 
ls ../../results/MECP2/MECP2WT.vs.KO_overlap_MECP2WT.vs.IgG.noBlacklist_filtered_DNMT3AKO_vs_WT.*sort.bed |\
parallel --verbose 'multiBigwigSummary BED-file  \
-b ../../../CUT_and_RUN/results/BigWig/9wk_CTX_DNMT3AWT_Mecp2_CNR_HS031219_N031419.bw \
    ../../../CUT_and_RUN/results/BigWig/Adult_CTX_WT_Mecp2_CNR_HS020619_N020719.bw \
    ../../../CUT_and_RUN/results/BigWig/9wk_CTX_DNMT3AKO_Mecp2_CNR_HS031219_N031419.bw \
    ../../../CUT_and_RUN/results/BigWig/Adult_CTX_DNMT3AKO_Mecp2_CNR_HS020619_N020719.bw \
-out {.}_multiBigWigOut.npz \
--BED {} \
--labels Dnmt3a_WT_rep1 Dnmt3a_WT_rep2 Dnmt3a_KO_rep1 Dnmt3a_KO_rep2 \
--blackListFileName ../../data/mm10/mm10.blacklist.bed \
--outRawCounts {.}_multiBigWigOut.tab \
-p "max/2"'


######################################################################################################################
# Supplementary figure 2K

computeMatrix reference-point \
-S  ../../../CUT_and_RUN/results/BigWig/9wk_CTX_DNMT3AWT_Mecp2_CNR_HS031219_N031419.bw \
    ../../../CUT_and_RUN/results/BigWig/Adult_CTX_WT_Mecp2_CNR_HS020619_N020719.bw \
    ../../../CUT_and_RUN/results/BigWig/9wk_CTX_DNMT3AKO_Mecp2_CNR_HS031219_N031419.bw \
    ../../../CUT_and_RUN/results/BigWig/Adult_CTX_DNMT3AKO_Mecp2_CNR_HS020619_N020719.bw \
       -R ../../results/MECP2/MECP2_peaks_with_diff_mCG_D3AKO_vs_WT.bed \
       --referencePoint center\
       -p 20 \
       -b 2000 -a 2000 \
       --missingDataAsZero \
       --skipZeros \
       -o ../../results/MECP2/MECP2_peaks_with_diff_mCG_D3AKO_vs_WT.bed.gz


plotProfile -m ../../results/MECP2/MECP2_peaks_with_diff_mCG_D3AKO_vs_WT.bed.gz \
-out ../../results/MECP2/MECP2_peaks_with_diff_mCG_D3AKO_vs_WT.bed.pdf \
--perGroup \
--yAxisLabel "Read Density" \
--samplesLabel "9wk DNMT3A WT" "9wk DNMT3A WT" \
               "Adult DNMT3A KO" "Adult DNMT3A KO" \
--plotWidth 8

