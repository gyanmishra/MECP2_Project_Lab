#!/usr/bin/bash

# Reserves  a node before running below code.
# srun --partition=super --nodes=8 --pty --time=5:00:00 /bin/bash

# Load below modules to execute the commands. 
module load parallel/20150122
module load bedtools/2.29.2
module load deeptools/3.5.0
module load bedops/2.4.14
module load UCSC_userApps/v317
module load homer/4.10.4 

######################################################################################################################
# Figure 1A 
# calculate mean coverage over geneBody of MECP2 CNR samples
bedtools map -b ../../../CUT_and_RUN/results/BigWig/7wk_CTX_MECP2WT_Mecp2_CNR_HS031219_N031419.bed \
-a ../../data/mm10/gencode.vM10.annotation_proteinCoding_gene_greater_than_4.5kb.bed \
-c 5 -o sum >../../results/MECP2/7wk_CTX_MECP2_gene_Cov.bed

bedtools map -b ../../../CUT_and_RUN/results/BigWig/7wk_CTX_MECP2WT_IgG_CNR_HS031219_N031419.bed \
-a ../../data/mm10/gencode.vM10.annotation_proteinCoding_gene_greater_than_4.5kb.bed \
-c 5 -o sum >../../results/MECP2/7wk_CTX_IgG_gene_Cov.bed

######################################################################################################################
# Figure 1C
computeMatrix reference-point \
-S ../../../CUT_and_RUN/results/BigWig/7wk_CTX_MECP2WT_Mecp2_CNR_HS031219_N031419.bw \
../../../CUT_and_RUN/results/BigWig/7wk_CTX_MECP2KO_Mecp2_CNR_HS031219_N031419.bw \
../../../CUT_and_RUN/results/BigWig/9wk_CTX_MECP2WT_Mecp2_CNR_HS032519_N032619.bw \
../../../CUT_and_RUN/results/BigWig/9wk_CTX_MECP2KO_Mecp2_CNR_HS032519_N032619.bw \
-R  ../../results/MECP2/MECP2WT.vs.KO_overlap_MECP2WT.vs.IgG.noBlacklist_filtered.bed \
--referencePoint center \
       -p 20 \
       -b 2000 -a 2000 \
       --skipZeros \
       --missingDataAsZero \
       -o ../../results/MECP2/MECP2WT.vs.KO_overlap_MECP2WT.vs.IgG.noBlacklist_filtered.gz

plotHeatmap -m ../../results/MECP2/MECP2WT.vs.KO_overlap_MECP2WT.vs.IgG.noBlacklist_filtered.gz \
--missingDataColor 'white' \
-out ../../results/MECP2/MECP2WT.vs.KO_overlap_MECP2WT.vs.IgG.noBlacklist_filtered.pdf \
--colorMap RdYlBu \
--samplesLabel "7wk MECP2 WT" "7wk MECP2 KO" \
               "9wk MECP2 WT" "9wk MECP2 KO" \
--heatmapHeight 3

######################################################################################################################
# Figure 1D
# Calculate coverage of MECP2 within MBHs and H3K9me3 and high mCA regions.
# High mCA group ######################################################################
bedtools map \
-a ../../results/BSMap/8wk_CTX_Stroud2017_mm10_BSmap_top1000_CA_5kbbin.bed \
-b ../../../CUT_and_RUN/results/BigWig/7wk_CTX_MECP2WT_Mecp2_CNR_HS031219_N031419.bed \
-c 5 -o sum >../../results/MECP2/7wk_CTX_MECP2WT_Mecp2_top1000_CA.bed

bedtools map \
-a ../../results/BSMap/8wk_CTX_Stroud2017_mm10_BSmap_top1000_CA_5kbbin.bed \
-b ../../../CUT_and_RUN/results/BigWig/9wk_CTX_MECP2WT_Mecp2_CNR_HS032519_N032619.bed \
-c 5 -o sum >.../../results/MECP2/9wk_CTX_MECP2WT_Mecp2_top1000_CA.bed

# 7wk_CTX_IgG
bedtools map \
-a ../../results/BSMap/8wk_CTX_Stroud2017_mm10_BSmap_top1000_CA_5kbbin.bed \
-b ../../../CUT_and_RUN/results/BigWig/7wk_CTX_MECP2WT_IgG_CNR_HS031219_N031419.bed \
-c 5 -o sum >../../results/MECP2/7wk_CTX_MECP2WT_IgG_top1000_CA.bed

# 9wk_CTX_IgG
bedtools map \
-a ../../results/BSMap/8wk_CTX_Stroud2017_mm10_BSmap_top1000_CA_5kbbin.bed \
-b ../../../CUT_and_RUN/results/BigWig/9wk_CTX_MECP2WT_IgG_CNR_HS032519_N032619.bed \
-c 5 -o sum >../../results/MECP2/9wk_CTX_MECP2WT_IgG_top1000_CA.bed

# MBH group ############################################################################
# 7wk_CTX_MECP2WT
bedtools map \
-a ../../results/MECP2/MECP2WT.vs.KO_overlap_MECP2WT.vs.IgG.noBlacklist_filtered.bed \
-b ../../../CUT_and_RUN/results/BigWig/7wk_CTX_MECP2WT_Mecp2_CNR_HS031219_N031419.bed \
-c 5 -o sum | cut -f 1,2,3,5 >../../results/MECP2/7wk_CTX_MECP2WT_Mecp2_peak_Coverage.bed

# 9wk_CTX_MECP2WT
bedtools map \
-a ../../results/MECP2/MECP2WT.vs.KO_overlap_MECP2WT.vs.IgG.noBlacklist_filtered.bed \
-b ../../../CUT_and_RUN/results/BigWig/9wk_CTX_MECP2WT_Mecp2_CNR_HS032519_N032619.bed \
-c 5 -o sum | cut -f 1,2,3,5 >../../results/MECP2/9wk_CTX_MECP2WT_Mecp2_peak_Coverage.bed

# 7wk_CTX_IgG
bedtools map \
-a ../../results/MECP2/MECP2WT.vs.KO_overlap_MECP2WT.vs.IgG.noBlacklist_filtered.bed \
-b ../../../CUT_and_RUN/results/BigWig/7wk_CTX_MECP2WT_IgG_CNR_HS031219_N031419.bed \
-c 5 -o sum | cut -f 1,2,3,5 >../../results/MECP2/7wk_CTX_MECP2WT_IgG_peak_Coverage.bed

# 9wk_CTX_IgG
bedtools map \
-a ../../results/MECP2/MECP2WT.vs.KO_overlap_MECP2WT.vs.IgG.noBlacklist_filtered.bed \
-b ../../../CUT_and_RUN/results/BigWig/9wk_CTX_MECP2WT_IgG_CNR_HS032519_N032619.bed \
-c 5 -o sum | cut -f 1,2,3,5 >../../results/MECP2/9wk_CTX_MECP2WT_IgG_peak_Coverage.bed


# H3K9me3 ############################################################################

liftOver ../../data/8wk_WT_CTX_H3K9me3_vs_H3_regions.bed \
../../data/mm10/mm9ToMm10.over.chain \
../../data/8wk_WT_CTX_H3K9me3_vs_H3_regions.mm10.bed \
../../data/8wk_WT_CTX_H3K9me3_vs_H3_regions.unMapped

sortBed -i ../../data/8wk_WT_CTX_H3K9me3_vs_H3_regions.mm10.bed \
>../../data/8wk_WT_CTX_H3K9me3_vs_H3_regions.mm10.sort.bed

# 7wk_CTX_MECP2WT
bedtools map \
-a ../../data/8wk_WT_CTX_H3K9me3_vs_H3_regions.mm10.sort.bed \
-b ../../../CUT_and_RUN/results/BigWig/7wk_CTX_MECP2WT_Mecp2_CNR_HS031219_N031419.bed \
-c 5 -o sum >../../results/MECP2/7wk_CTX_MECP2WT_Mecp2_H3K9me3.bed

bedops -n 1 ../../results/MECP2/7wk_CTX_MECP2WT_Mecp2_H3K9me3.bed \
../../data/mm10/mm10.blacklist.bed \
>../../results/MECP2/7wk_CTX_MECP2WT_Mecp2_H3K9me3_noBlacklist.bed

# 9wk_CTX_MECP2WT
bedtools map \
-a ../../data/8wk_WT_CTX_H3K9me3_vs_H3_regions.mm10.sort.bed \
-b ../../../CUT_and_RUN/results/BigWig/9wk_CTX_MECP2WT_Mecp2_CNR_HS032519_N032619.bed \
-c 5 -o sum >../../results/MECP2/9wk_CTX_MECP2WT_Mecp2_H3K9me3.bed

bedops -n 1 ../../results/MECP2/9wk_CTX_MECP2WT_Mecp2_H3K9me3.bed \
../../data/mm10/mm10.blacklist.bed \
>../../results/MECP2/9wk_CTX_MECP2WT_Mecp2_H3K9me3_noBlacklist.bed

# IgG on H3K9me3 sites
# 7wk_CTX_IgG
bedtools map \
-a ../../data/8wk_WT_CTX_H3K9me3_vs_H3_regions.mm10.sort.bed \
-b ../../../CUT_and_RUN/results/BigWig/7wk_CTX_MECP2WT_IgG_CNR_HS031219_N031419.bed \
-c 5 -o sum >../../results/MECP2/7wk_CTX_MECP2WT_IgG_H3K9me3.bed

bedops -n 1 ../../results/MECP2/7wk_CTX_MECP2WT_IgG_H3K9me3.bed \
../../data/mm10/mm10.blacklist.bed \
>../../results/MECP2/7wk_CTX_MECP2WT_IgG_H3K9me3_noBlacklist.bed

# 9wk_CTX_IgG
bedtools map \
-a ../../data/8wk_WT_CTX_H3K9me3_vs_H3_regions.mm10.sort.bed \
-b ../../../CUT_and_RUN/results/BigWig/9wk_CTX_MECP2WT_IgG_CNR_HS032519_N032619.bed \
-c 5 -o sum >../../results/MECP2/9wk_CTX_MECP2WT_IgG_H3K9me3.bed

bedops -n 1 ../../results/MECP2/9wk_CTX_MECP2WT_Mecp2_H3K9me3.bed \
../../data/mm10/mm10.blacklist.bed \
>../../results/MECP2/9wk_CTX_MECP2WT_IgG_H3K9me3_noBlacklist.bed

######################################################################################################################
# Figure 1E
# Extract CA repeat region using homer
mkdir -p ../../results/CAsequence
analyzeRepeats.pl repeats mm10 |grep '(CA)n' >../../results/CAsequence/CA_repeatsHomer.txt
cut -f 2,3,4 ../../results/CAsequence/CA_repeatsHomer.txt >../../results/CAsequence/CA_repeatsHomer.bed

# Overlap repeat regions with MECP2 peaks to check how many MBHs sites overlap with CA repeat regions
bedtools intersect -wa \
-a ../../results/MECP2/MECP2WT.vs.KO_overlap_MECP2WT.vs.IgG.noBlacklist_filtered.bed \
-a ../../results/CAsequence/CA_repeatsHomer.bed \
>../../results/CAsequence/MECP2_CArepeatOverlap.bed

######################################################################################################################
# Figure 1G
# GGTGT, GTGT, GGGTTT and TTTGGG
seq2profile.pl GGTGT 0 GGTGT >../../results/CAsequence/GGTGT.motif
seq2profile.pl GTGT 0 GTGT >../../results/CAsequence/GTGT.motif
seq2profile.pl GGGTTT 0 GGGTTT >../../results/CAsequence/GGGTTT.motif
seq2profile.pl TTTGGG 0 TTTGGG >../../results/CAsequence/TTTGGG.motif
seq2profile.pl CCACA 0 CCACA >../../results/CAsequence/CCACA.motif

annotatePeaks.pl ../../results/MECP2/MECP2WT.vs.KO_overlap_MECP2WT.vs.IgG.noBlacklist_filtered.bed mm10 \
-m ../../results/CAsequence/GGGTTT.motif \
../../results/CAsequence/GGTGT.motif \
../../results/CAsequence/GTGT.motif \
../../results/CAsequence/CCAC.motif \
../../results/CAsequence/TTTGGG.motif \
../../results/CAsequence/CCACA.motif \
-size 4000 -hist 10 \
>../../results/CAsequence/GT_motif_freq_on_MECP2.txt

######################################################################################################################
# Supplementary Figure 1A

# compare genome-wide coverage of MECP2 CUT&RUN and MECP2 ChIPSeq(Boxer et al. 2020) samples.
mkdir -p ../../resutls/MECP2_CNR_ChIPseq_comparison
multiBigwigSummary bins \
-b ../../../CUT_and_RUN/results/BigWig/7wk_CTX_MECP2_over_IgG_CNR_HS031219_N031419.bw \
../../../CUT_and_RUN/results/BigWig/9wk_CTX_MECP2_over_IgG_CNR_HS031219_N031419.bw \
../../data/GSE139509/bigWig/SRR10356997_98_over_SRR10357005_06.bw \
../../data/GSE139509/bigWig/SRR10357001_02_over_SRR10357005_06.bw \
--labels "MECP2 CNR R1" "MECP2 CNR R2" "MECP2 ChIPseq Ab1" "MECP2 ChIPseq Ab2" \
--blackListFileName ../../data/mm10/mm10.blacklist.bed \
-out ../../results/MECP2_CNR_ChIPseq_comparison/MECP2_CNR_ChIPseq_multiBigWigOut.npz \
--outRawCounts ../../results/MECP2_CNR_ChIPseq_comparison/MECP2_CNR_ChIPseq_multiBigWigOut.tab \
-p "max/2"


# xRange and yRange option only works with deeptool v3.1.2
plotCorrelation \
-in ../../results/MECP2_CNR_ChIPseq_comparison/MECP2_CNR_ChIPseq_multiBigWigOut.npz \
--corMethod spearman --skipZeros \
--plotTitle "Pearson Correlation of Average Scores Per bin" \
--whatToPlot scatterplot \
-o ../../results/MECP2_CNR_ChIPseq_comparison/MECP2_CNR_ChIPseq_multiBigWigOut.pdf   \
--outFileCorMatrix ../../results/MECP2_CNR_ChIPseq_comparison/MECP2_CNR_ChIPseq_multiBigWigOut_PearsonCorr_bigwigScores.tab \
--xRange -1 2 \
--yRange -1 2

######################################################################################################################
# Supplementary Figure 1B
# generate computeMatrix for all the MECP2 cut&run and ChIP sample
computeMatrix reference-point \
-S \
        ../../../CUT_and_RUN/results/BigWig/7wk_CTX_MECP2_over_IgG_CNR_HS031219_N031419.bw \
        ../../..//CUT_and_RUN/results/BigWig/9wk_CTX_MECP2_over_IgG_CNR_HS031219_N031419.bw \
        ../../data/GSE139509/bigWig/SRR10356997_98_over_SRR10357005_06.bw \
        ../../data/GSE139509/bigWig/SRR10357001_02_over_SRR10357005_06.bw \
        ../../data/GSE67293/GSM1643934-GSM1643939/GSM1643934_MECP2_over_GSM1643935_INPUT.bw \
-R  ../../results/MECP2/MECP2WT.vs.KO_overlap_MECP2WT.vs.IgG.noBlacklist_filtered.bed \
--referencePoint center \
-p "max/2" \
-b 2000 -a 2000 \
--skipZeros \
-o ../../results/MECP2_CNR_ChIPseq_comparison/GSE139509_MECP2_CNR.gz

plotProfile -m ../../results/MECP2_CNR_ChIPseq_comparison/GSE139509_MECP2_CNR.gz \
-out ../../results/MECP2_CNR_ChIPseq_comparison/GSE139509_MECP2_CNR.AvgProfile.png \
--samplesLabel  "7wk MECP2 CUT&RUN" "9wk MECP2 CUT&RUN" \
"MECP2 ChIPseq Ab1 (Forebrain)" "MECP2 ChIPseq Ab2 (Forebrain)" "MECP2 ChIPseq (FC1)" \
--perGroup \
--yAxisLabel "log2 (CUT&RUN/ChIPseq)/(IgG/Input)" \
--plotTitle ""

######################################################################################################################
# Supplementary Figure 1C
computeMatrix reference-point \
-S \
../../data/GSE213752/GSM6593478_57.SpikeInNormalized.bw \
../../data/GSE213752/GSM6593490_69.SpikeInNormalized.bw \
../../data/GSE213752/GSM6593481_25.SpikeInNormalized.bw \
../../data/GSE213752/GSM6593493_36.SpikeInNormalized.bw \
../../data/GSE213752/GSM6593484_60.SpikeInNormalized.bw \
../../data/GSE213752/GSM6593496_72.SpikeInNormalized.bw \
-R  ../../results/MECP2/MECP2WT.vs.KO_overlap_MECP2WT.vs.IgG.noBlacklist_filtered.bed \
--referencePoint center \
-p "max/2" \
-b 2000 -a 2000 \
--skipZeros \
--samplesLabel "IgG WT" "MECP2 WT" "IgG TG1" "MECP2 TG1" "IgG MECP2 KO" "MECP2 MECP2 KO" \
-o ../../results/MECP2_CNR_ChIPseq_comparison/GSE213752_MECP2_CNR.gz

plotProfile -m ../../results/MECP2_CNR_ChIPseq_comparison/GSE213752_MECP2_CNR.gz \
-out ../../results/MECP2_CNR_ChIPseq_comparison/GSE213752_MECP2_CNR.AvgProfile.pdf \
--perGroup \
--yAxisLabel "Normalized Coverage" \
--plotTitle "" \
--plotWidth 8


######################################################################################################################
# Supplementary Figure 1D
# compare genome-wide coverage of MECP2 CUT&RUN and MECP2 CUT&TAG samples.

multiBigwigSummary bins --binSize 10000 \
-b ../../results/bamCoverage_out/7wk_CTX_MECP2WT_Mecp2_CNR_HS031219_N031419.final.sort.final.RPKM.bw \
../../results/bamCoverage_out/9wk_CTX_MECP2WT_Mecp2_CNR_HS032519_N032619.final.sort.final.RPKM.bw \
../../results/bamCoverage_out/Adult_CTX_4399_0.5_A_CNT_HS010421_N012121.final.sort.RPKM.bw \
../../results/bamCoverage_out/Adult_CTX_MECP2mono_0.5_A_CNT_HS010421_N012121.final.sort.RPKM.bw \
--labels "MECP2 CNR R1" "MECP2 CNR R2" "MECP2 CNT Ab1" "MECP2 CNT Ab2" \
--blackListFileName ../../data/mm10/mm10.blacklist.bed \
-out ../../results/MECP2_CNR_ChIPseq_comparison/MECP2_CNR_CNT_multiBigWigOut.npz \
--outRawCounts ../../results/MECP2_CNR_ChIPseq_comparison/MECP2_CNR_CNT_multiBigWigOut.tab \
-p "max/2"

# xRange and yRange option only works with deeptool v3.1.2)
plotCorrelation \
-in ../../results/MECP2_CNR_ChIPseq_comparison/MECP2_CNR_CNT_multiBigWigOut.npz \
--corMethod spearman --skipZeros \
--plotTitle "Pearson Correlation of Average Scores Per bin" \
--whatToPlot scatterplot \
-o ../../results/MECP2_CNR_ChIPseq_comparison/MECP2_CNR_CNT_multiBigWigOut.pdf   \
--outFileCorMatrix ../../results/MECP2_CNR_ChIPseq_comparison/MECP2_CNR_CNT_multiBigWigOut_PearsonCorr_bigwigScores.tab \
--log1p \
--xRange 0 3 \
--yRange 0 3

######################################################################################################################
# Supplementary Figure 1E

computeMatrix reference-point -S  \
../../../CUT_and_RUN/results/BigWig/Adult_CTX_H3_0.5_Native_CNR_YT042919_N050119.bw \
../../../CUT_and_RUN/results/BigWig/Adult_CTX_IgG_0.5_Native2hr_CNR_YT042919_N050119.bw \
../../../CUT_and_RUN/results/BigWig/Adult_CTX_IgG_0.5_Native_CNR_YT042919_N050119.bw \
../../../CUT_and_RUN/results/BigWig/Adult_CTX_MECP2_0.5_Native2hr_CNR_YT042919_N050119.bw \
../../../CUT_and_RUN/results/BigWig/Adult_CTX_MECP2_0.5_Native_CNR_YT042919_N050119.bw \
../../../CUT_and_RUN/results/BigWig/Adult_CTX_MECP2_10_Native_CNR_YT042919_N050119.bw \
-R ../../results/MECP2/MECP2WT.vs.KO_overlap_MECP2WT.vs.IgG.noBlacklist_filtered.bed \
--referencePoint center \
-p 20 \
-b 2000 -a 2000 \
--missingDataAsZero \
--skipZeros \
 -o ../../results/MECP2/H3_Native_CNR_on_MECP2.gz


plotProfile -m ../../results/MECP2/H3_Native_CNR_on_MECP2.gz \
-out ../../results/MECP2/H3_Native_CNR_on_MECP2.pdf \
--samplesLabel  "H3 (0.5 Native)" "IgG (0.5 Native 2h)" "IgG (0.5 Native)" "IgG (0.5 Native lyse)" \
                "MECP2 (0.5 Native 2h)" "MECP2 (0.5 Native)" "MECP2 (0.5 Native lyse)" "MECP2 (10 Native)" \
--perGroup \
--yAxisLabel "Normalized Coverage" \
--plotTitle "" \
--plotWidth 8

computeMatrix reference-point        \
-S  ../../../CUT_and_RUN/results/BigWig/Adult_CTX_IgG_0.5_Nativelyse_CNR_N050119.bw \
../../../CUT_and_RUN/results/BigWig/Adult_CTX_MECP2_0.5_Nativelyse_CNR_YT042919_N050119.bw        \
-R ../../results/MECP2/MECP2WT.vs.KO_overlap_MECP2WT.vs.IgG.noBlacklist_filtered.bed        \
--referencePoint center \
-p 20 \
-b 2000 -a 2000 \
--skipZeros   \
 -o ../../results/MECP2/MECP2_Nativelyse_CNR_on_MECP2.gz

plotProfile -m ../../results/MECP2/MECP2_Nativelyse_CNR_on_MECP2.gz \
-out ../../results/MECP2/MECP2_Nativelyse_CNR_on_MECP2.pdf \
--samplesLabel  "IgG (0.5 Native lyse)" "MECP2 (0.5 Native lyse)" \
--perGroup \
--yAxisLabel "Normalized Coverage" \
--plotTitle "" \
--plotWidth 8


######################################################################################################################
# Supplementary Figure 1F
computeMatrix reference-point \
-S ../../../CUT_and_RUN/results/BigWig/Adult_CTX_4399_0.5_A_CNT_HS010421_N012121.bw \
../../../CUT_and_RUN/results/BigWig/Adult_CTX_MECP2mono_0.5_A_CNT_HS010421_N012121.bw \
-R  ../../results/MECP2/MECP2WT.vs.KO_overlap_MECP2WT.vs.IgG.noBlacklist_filtered.bed \
--referencePoint center \
-p "max/2" \
-b 2000 -a 2000 \
--skipZeros \
-o ../../results/MECP2_CNR_ChIPseq_comparison/MECP2mono_CUTandTAG.gz

plotProfile -m ../../results/MECP2_CNR_ChIPseq_comparison/MECP2mono_CUTandTAG.gz \
-out ../../results/MECP2_CNR_ChIPseq_comparison/MECP2mono_CUTandTAG.pdf \
--samplesLabel  "MECP2 CUT&TAG (Ab1) Rep1" "MECP2 CUT&TAG (Ab2) Rep1" \
--perGroup \
--yAxisLabel "Noramlized Coverage" \
--plotTitle "" \
--plotWidth 8

######################################################################################################################
# Supplementary Figure 1G,H
bash MECP2_G118E_CNR.sh

# Supplementary Figure 1J
# extract Motif2 and Motif4 from results/MECP2_motifs/MECP2_all_peaks_STREME/streme.txt
annotatePeaks.pl ../../results/MECP2_motifs/MECP2_pm50_CCAC.bed mm10 \
-m ../../results/MECP2_motifs/MECP2_all_peaks_STREME/Motif24.motif -size 2000 -hist 10




