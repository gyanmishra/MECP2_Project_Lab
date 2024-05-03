#!/usr/bin/bash

# Load required tools
module load bedtools/2.29.0
module load bedops/2.4.14
module load parallel/20150122 

mkdir -p ../../results/MECP2

# sort Peak files
ls ../../../CUT_and_RUN/results/macs2_Callpeak/*/*_peaks.narrowPeak | \
parallel --verbose 'sort-bed {} | \bedops -m - >{}.sort.bed'

# process MECP2 peak (MECP2 WT vs MECP2 KO) files
awk '{ if($5 >4) print $1,$2-249,$2+249}' \
../../../CUT_and_RUN/results/macs2_Callpeak/Adult_MECP2_MECP2WT_vs_KO_Macs2/Adult_MECP2_MECP2WT_vs_KO_summits.bed |\
sort-bed - \
>../../../CUT_and_RUN/results/macs2_Callpeak/Adult_MECP2_MECP2WT_vs_KO_Macs2/Adult_MECP2_MECP2WT_vs_KO_summits.pm250.sort.bed 

# process MECP2 peak (MECP2 WT vs MECP2 IgG) files
awk '{ if($5 >4) print $1,$2-249,$2+249}' \
../../../CUT_and_RUN/results/macs2_Callpeak/Adult_MECP2_MECP2WT_vs_IgG_macs2/Adult_MECP2_MECP2WT_vs_IgG_summits.bed |\
sort-bed - \
>../../../CUT_and_RUN/results/macs2_Callpeak/Adult_MECP2_MECP2WT_vs_IgG_macs2/Adult_MECP2_MECP2WT_vs_IgG_summits.pm250.sort.bed


# Filter MECP2 peaks from MECP2WT_vs_KO also called in MECP2WT_vs_IgG comparison 
bedtools intersect -wa \
-a ../../../CUT_and_RUN/results/macs2_Callpeak/Adult_MECP2_MECP2WT_vs_KO_Macs2/Adult_MECP2_MECP2WT_vs_KO_summits.pm250.sort.bed \
-b ../../../CUT_and_RUN/results/macs2_Callpeak/Adult_MECP2_MECP2WT_vs_IgG_macs2/Adult_MECP2_MECP2WT_vs_IgG_summits.pm250.sort.bed \
>../../results/MECP2/MECP2WT.vs.KO_overlap_MECP2WT.vs.IgG.bed

# Filter blacklisted regions
bedops -n ../../results/MECP2/MECP2WT.vs.KO_overlap_MECP2WT.vs.IgG.bed \
          ../../../CUT_and_RUN/data/mm10/mm10.blacklist.bed  |\
bedops -m - >../../results/MECP2/MECP2WT.vs.KO_overlap_MECP2WT.vs.IgG.noBlacklist.bed

# Add peak ID for each peaks
awk '{ print $1,$2,$3,"Peak"NR}' OFS='\t' ../../results/MECP2/MECP2WT.vs.KO_overlap_MECP2WT.vs.IgG.noBlacklist.bed \
>../../results/MECP2/MECP2WT.vs.KO_overlap_MECP2WT.vs.IgG.noBlacklist_peakID.bed

# Extract tag Count for each peaks using HOMER
annotatePeaks.pl ../results/MECP2WT.vs.KO_overlap_MECP2WT.vs.IgG.noBlacklist_peakID.bed mm10 
       -d ../../../CUT_and_RUN/results/7wk_CTX_MECP2WT_Mecp2_CNR_HS031219_N031419.tagDir 
          ../../../CUT_and_RUN/results/9wk_CTX_DNMT3AWT_Mecp2_CNR_HS031219_N031419.tagDir 
          ../../../CUT_and_RUN/results/9wk_CTX_MECP2WT_Mecp2_CNR_HS032519_N032619.tagDir 
          ../../../CUT_and_RUN/results/Adult_CTX_WT_Mecp2_CNR_HS020619_N020719.tagDir 
          ../../../CUT_and_RUN/results/7wk_CTX_MECP2KO_Mecp2_CNR_HS031219_N031419.tagDir 
          ../../../CUT_and_RUN/results/9wk_CTX_DNMT3AKO_Mecp2_CNR_HS031219_N031419.tagDir 
          ../../../CUT_and_RUN/results/9wk_CTX_MECP2KO_Mecp2_CNR_HS032519_N032619.tagDir 
          ../../../CUT_and_RUN/results/Adult_CTX_DNMT3AKO_Mecp2_CNR_HS020619_N020719.tagDir 
       >../../results/MECP2/MECP2WT.vs.KO_overlap_MECP2WT.vs.IgG.noBlacklist_tagCount.tsv
