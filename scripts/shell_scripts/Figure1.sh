#!/usr/bin/bash

# Reserve a node before running below code.
# srun --partition=super --nodes=8 --pty --time=5:00:00 /bin/bash

# Load below modules to execute the commands. 
module load parallel/20150122
module load bedtools/2.29.2
#module load deeptools/3.5.0
module load bedops/2.4.14
module load UCSC_userApps/v317
module load homer/4.10.4


bigwig=$1

mkdir -p ../../results/MECP2/
mkdir -p ../../Figures

# Figure 1C
computeMatrix reference-point \
-S $bigwig/7wk_CTX_MECP2WT_Mecp2_CNR_HS031219_N031419.bw \
$bigwig/7wk_CTX_MECP2KO_Mecp2_CNR_HS031219_N031419.bw \
$bigwig/9wk_CTX_MECP2WT_Mecp2_CNR_HS032519_N032619.bw \
$bigwig/9wk_CTX_MECP2KO_Mecp2_CNR_HS032519_N032619.bw \
-R  ../../results/MECP2/MECP2WT.vs.KO_overlap_MECP2WT.vs.IgG.noBlacklist_filtered.bed \
--referencePoint center \
       -p 20 \
       -b 2000 -a 2000 \
       --skipZeros \
       --missingDataAsZero \
       -o ../../results/MECP2/MECP2WT.vs.KO_overlap_MECP2WT.vs.IgG.noBlacklist_filtered.gz

plotHeatmap -m ../../results/MECP2/MECP2WT.vs.KO_overlap_MECP2WT.vs.IgG.noBlacklist_filtered.gz \
--missingDataColor 'white' \
-out ../../Figures/MECP2WT.vs.KO_overlap_MECP2WT.vs.IgG.noBlacklist_filtered.pdf \
--colorMap RdYlBu \
--samplesLabel "7wk MECP2 WT" "7wk MECP2 KO" \
               "9wk MECP2 WT" "9wk MECP2 KO" \
--heatmapHeight 3

######################################################################################################################
# Supplementary Figure 1A

# compare genome-wide coverage of MECP2 CUT&RUN and MECP2 ChIPSeq (Boxer et al. 2020) samples.
mkdir -p ../../results/MECP2_CNR_ChIPseq_comparison

multiBigwigSummary bins \
-b $bigwig/7wk_CTX_MECP2_over_IgG_CNR_HS031219_N031419.bw \
$bigwig/9wk_CTX_MECP2_over_IgG_CNR_HS031219_N031419.bw \
$bigWig/SRR10356997_98_over_SRR10357005_06.bw \
$bigWig/SRR10357001_02_over_SRR10357005_06.bw \
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
-o ../../Figures/MECP2_CNR_ChIPseq_multiBigWigOut.pdf   \
--outFileCorMatrix ../../results/MECP2_CNR_ChIPseq_comparison/MECP2_CNR_ChIPseq_multiBigWigOut_PearsonCorr_bigwigScores.tab \
--xRange -1 2 \
--yRange -1 2

######################################################################################################################
# Supplementary Figure 1B
# generate computeMatrix for all the MECP2 cut&run and ChIP sample
computeMatrix reference-point \
-S \
        $bigwig/7wk_CTX_MECP2_over_IgG_CNR_HS031219_N031419.bw \
        $bigwig/9wk_CTX_MECP2_over_IgG_CNR_HS031219_N031419.bw \
        $bigWig/SRR10356997_98_over_SRR10357005_06.bw \
        $bigWig/SRR10357001_02_over_SRR10357005_06.bw \
        $bigWig/GSM1643934_MECP2_over_GSM1643935_INPUT.bw \
-R  ../../results/MECP2/MECP2WT.vs.KO_overlap_MECP2WT.vs.IgG.noBlacklist_filtered.bed \
--referencePoint center \
-p "max/2" \
-b 2000 -a 2000 \
--skipZeros \
-o ../../results/MECP2_CNR_ChIPseq_comparison/GSE139509_MECP2_CNR.gz

plotProfile -m ../../results/MECP2_CNR_ChIPseq_comparison/GSE139509_MECP2_CNR.gz \
-out ../../Figures/GSE139509_MECP2_CNR.AvgProfile.pdf \
--samplesLabel  "7wk MECP2 CUT&RUN" "9wk MECP2 CUT&RUN" \
"MECP2 ChIPseq Ab1 (Forebrain)" "MECP2 ChIPseq Ab2 (Forebrain)" "MECP2 ChIPseq (FC1)" \
--perGroup \
--yAxisLabel "log2 (CUT&RUN/ChIPseq)/(IgG/Input)" \
--plotTitle ""

######################################################################################################################
# Supplementary Figure 1C
computeMatrix reference-point \
-S \
$bigwig/GSM6593478_57.SpikeInNormalized.bw \
$bigwig/GSM6593490_69.SpikeInNormalized.bw \
$bigwig/GSM6593481_25.SpikeInNormalized.bw \
$bigwig/GSM6593493_36.SpikeInNormalized.bw \
$bigwig/GSM6593484_60.SpikeInNormalized.bw \
$bigwig/GSM6593496_72.SpikeInNormalized.bw \
-R  ../../results/MECP2/MECP2WT.vs.KO_overlap_MECP2WT.vs.IgG.noBlacklist_filtered.bed \
--referencePoint center \
-p "max/2" \
-b 2000 -a 2000 \
--skipZeros \
--samplesLabel "IgG WT" "MECP2 WT" "IgG TG1" "MECP2 TG1" "IgG MECP2 KO" "MECP2 MECP2 KO" \
-o ../../results/MECP2_CNR_ChIPseq_comparison/GSE213752_MECP2_CNR.gz

plotProfile -m ../../results/MECP2_CNR_ChIPseq_comparison/GSE213752_MECP2_CNR.gz \
-out ../../Figures/GSE213752_MECP2_CNR.AvgProfile.pdf \
--perGroup \
--yAxisLabel "Normalized Coverage" \
--plotTitle "" \
--plotWidth 8


######################################################################################################################
# Supplementary Figure 1D
# compare genome-wide coverage of MECP2 CUT&RUN and MECP2 CUT&TAG samples.

multiBigwigSummary bins --binSize 10000 \
-b $bigwig/7wk_CTX_MECP2WT_Mecp2_CNR_HS031219_N031419.final.sort.final.RPKM.bw \
$bigwig/9wk_CTX_MECP2WT_Mecp2_CNR_HS032519_N032619.final.sort.final.RPKM.bw \
$bigwig/Adult_CTX_4399_0.5_A_CNT_HS010421_N012121.final.sort.RPKM.bw \
$bigwig/Adult_CTX_MECP2mono_0.5_A_CNT_HS010421_N012121.final.sort.RPKM.bw \
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
-o ../../Figures/MECP2_CNR_CNT_multiBigWigOut.pdf   \
--outFileCorMatrix ../../results/MECP2_CNR_ChIPseq_comparison/MECP2_CNR_CNT_multiBigWigOut_PearsonCorr_bigwigScores.tab \
--log1p \
--xRange 0 3 \
--yRange 0 3

######################################################################################################################
# Supplementary Figure 1E

computeMatrix reference-point -S  \
$bigwig/Adult_CTX_H3_0.5_Native_CNR_YT042919_N050119.bw \
$bigwig/Adult_CTX_IgG_0.5_Native2hr_CNR_YT042919_N050119.bw \
$bigwig/Adult_CTX_IgG_0.5_Native_CNR_YT042919_N050119.bw \
$bigwig/Adult_CTX_MECP2_0.5_Native2hr_CNR_YT042919_N050119.bw \
$bigwig/Adult_CTX_MECP2_0.5_Native_CNR_YT042919_N050119.bw \
$bigwig/Adult_CTX_MECP2_10_Native_CNR_YT042919_N050119.bw \
-R ../../results/MECP2/MECP2WT.vs.KO_overlap_MECP2WT.vs.IgG.noBlacklist_filtered.bed \
--referencePoint center \
-p 20 \
-b 2000 -a 2000 \
--missingDataAsZero \
--skipZeros \
 -o ../../results/MECP2/H3_Native_CNR_on_MECP2.gz


plotProfile -m ../../results/MECP2/H3_Native_CNR_on_MECP2.gz \
-out ../../Figures/MECP2/H3_Native_CNR_on_MECP2.pdf \
--samplesLabel  "H3 (0.5 Native)" "IgG (0.5 Native 2h)" "IgG (0.5 Native)" "IgG (0.5 Native lyse)" \
                "MECP2 (0.5 Native 2h)" "MECP2 (0.5 Native)" "MECP2 (0.5 Native lyse)" "MECP2 (10 Native)" \
--perGroup \
--yAxisLabel "Normalized Coverage" \
--plotTitle "" \
--plotWidth 8

computeMatrix reference-point        \
-S  $bigwig/Adult_CTX_IgG_0.5_Nativelyse_CNR_N050119.bw \
$bigwig/Adult_CTX_MECP2_0.5_Nativelyse_CNR_YT042919_N050119.bw        \
-R ../../results/MECP2/MECP2WT.vs.KO_overlap_MECP2WT.vs.IgG.noBlacklist_filtered.bed        \
--referencePoint center \
-p 20 \
-b 2000 -a 2000 \
--skipZeros   \
 -o ../../results/MECP2/MECP2_Nativelyse_CNR_on_MECP2.gz

plotProfile -m ../../results/MECP2/MECP2_Nativelyse_CNR_on_MECP2.gz \
-out ../../Figures/MECP2_Nativelyse_CNR_on_MECP2.pdf \
--samplesLabel  "IgG (0.5 Native lyse)" "MECP2 (0.5 Native lyse)" \
--perGroup \
--yAxisLabel "Normalized Coverage" \
--plotTitle "" \
--plotWidth 8


######################################################################################################################
# Supplementary Figure 1F
computeMatrix reference-point \
-S $bigwig/Adult_CTX_4399_0.5_A_CNT_HS010421_N012121.bw \
$bigwig/Adult_CTX_MECP2mono_0.5_A_CNT_HS010421_N012121.bw \
-R  ../../results/MECP2/MECP2WT.vs.KO_overlap_MECP2WT.vs.IgG.noBlacklist_filtered.bed \
--referencePoint center \
-p "max/2" \
-b 2000 -a 2000 \
--skipZeros \
-o ../../results/MECP2_CNR_ChIPseq_comparison/MECP2mono_CUTandTAG.gz

plotProfile -m ../../results/MECP2_CNR_ChIPseq_comparison/MECP2mono_CUTandTAG.gz \
-out ../../Figures/MECP2mono_CUTandTAG.pdf \
--samplesLabel  "MECP2 CUT&TAG (Ab1) Rep1" "MECP2 CUT&TAG (Ab2) Rep1" \
--perGroup \
--yAxisLabel "Noramlized Coverage" \
--plotTitle "" \
--plotWidth 8

######################################################################################################################
# Supplementary Figure 1G, H
bash MECP2_G118E_CNR.sh





