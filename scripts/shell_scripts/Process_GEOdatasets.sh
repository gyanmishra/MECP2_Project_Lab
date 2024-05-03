#!/usr/bin/bash

module load parallel/20150122
module load deeptools/3.5.0
module load bedops/2.4.14
module load UCSC_userApps/v317

# download supplementary file from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE67293
tar -xvf ../../data/GSE67293/GSE67293_RAW.tar

# Convert bigwig to bedGraph
ls ../../data/GSE67293/GSM1643934-GSM1643939/GSM164393*.bw | parallel --verbose  'bigWigToBedGraph {} {.}.mm9.bedGraph'

# Liftover from mm9 to mm10
ls ../../data/GSE67293/GSM1643934-GSM1643939/GSM164393*.mm9.bedGraph | parallel --verbose  \
'liftOver {} ../../data/mm10/mm9ToMm10.over.chain {= s/.mm9.bedGraph// =}.mm10.bedGraph {.}.unMapped'

# Sort the bed file
ls ../../data/GSE67293/GSM1643934-GSM1643939/GSM164393*.mm10.bedGraph | \
parallel --verbose 'sort -k1,1 -k2,2n {} > {.}.sort.bedGraph'

# Fix the overlapping intervals
# The bed file from liftover might have overlapping intervals. 
# You will hit error if directly using bedGraphToBigWig for file conversion. 
# In this case, you need to split the overlapping intervals and assign mean signal to a new bed file.
ls ../../data/GSE67293/GSM1643934-GSM1643939/GSM164393*.mm10.sort.bedGraph | \
parallel --verbose awk -vOFS=\"\\t\" \'{ print \$1, \$2, \$3, \".\", \$4 }\'  {} \> {.}.bed
ls ../../data/GSE67293/GSM1643934-GSM1643939/GSM164393*.mm10.sort.bed | \
parallel --verbose 'bedops --partition {} | bedmap --echo --mean --delim "\t" - {} >  {.}.split.bedGraph'

# Convert bedGraph to bigwig 
ls ../../data/GSE67293/GSM1643934-GSM1643939/GSM164393*.mm10.sort.split.bedGraph | \
parallel --verbose 'bedGraphToBigWig {} ../../data/mm10/mm10.chrom.sizes {.}.bw'

# generate bamComapre bigwig files for MECP2 CUT&RUN data
bigwigCompare --bigwig1 ../../../CUT_and_RUN/results/BigWig/7wk_CTX_MECP2WT_Mecp2_CNR_HS031219_N031419.bw \
              --bigwig2 ../../../CUT_and_RUN/results/BigWig/7wk_CTX_MECP2WT_IgG_CNR_HS031219_N031419.bw \
              --outFileName ../../../CUT_and_RUN/results/BigWig/7wk_CTX_MECP2_over_IgG_CNR_HS031219_N031419.bw \
              -p "max/2"


bigwigCompare --bigwig1 ../../../CUT_and_RUN/results/BigWig/9wk_CTX_MECP2WT_Mecp2_CNR_HS032519_N032619.bw \
              --bigwig2 ../../../CUT_and_RUN/results/BigWig/9wk_CTX_MECP2WT_IgG_CNR_HS032519_N032619.bw \
              --outFileName ../../../CUT_and_RUN/results/BigWig/9wk_CTX_MECP2_over_IgG_CNR_HS031219_N031419.bw  \
              -p "max/2"

# generate bamComapre bigwig files for MECP2 ChIPseq data

bigwigCompare   --bigwig1 ../../data/GSE67293/GSM1643934-GSM1643939/GSM1643934_MeCP2IP_MeCP2_ChIP_FC1.mm10.sort.split.bw \
                --bigwig2 ../../data/GSE67293/GSM1643934-GSM1643939/GSM1643935_INPUT_MeCP2_ChIP_FC1.mm10.sort.split.bw \
                --outFileName GSM1643934_MECP2_over_GSM1643935_INPUT.bw \
                -p "max/2"

bigwigCompare   --bigwig1 ../../data/GSE67293/GSM1643934-GSM1643939/GSM1643936_MeCP2IP_MeCP2_ChIP_FC2.mm10.sort.split.bw \
                --bigwig2 ../../data/GSE67293/GSM1643934-GSM1643939/GSM1643937_INPUT_MeCP2_ChIP_FC2.mm10.sort.split.bw \
                --outFileName GSM1643936_MECP2_over_GSM1643937_INPUT.bw \
                -p "max/2"

bigwigCompare   --bigwig1 ../../data/GSE67293/GSM1643934-GSM1643939/GSM1643938_MeCP2IP_MeCP2_ChIP_CTX1.mm10.sort.split.bw \
                --bigwig2 ../../data/GSE67293/GSM1643934-GSM1643939/GSM1643939_INPUT_MeCP2_ChIP_CTX1.mm10.sort.split.bw \
                --outFileName GSM1643938_MECP2_over_GSM1643939_INPUT.bw \
                -p "max/2"


# Merge repilcates of MECP2 ChIP WT samples from GSE139509 datasets
wiggletools write_bg ../../data/GSE139509/bigWig/SRR10356997_98.Average.bedGraph mean \
../../data/GSE139509/bigWig/SRR10356997.bw ../../data/GSE139509/bigWig/SRR10356998.bw 
bedGraphToBigWig ../../data/GSE139509/bigWig/SRR10356997_98.Average.bedGraph ../../data/mm10/mm10.chrom.sizes \
../../data/GSE139509/bigWig/SRR10356997_98.Average.bw

wiggletools write_bg ../../data/GSE139509/bigWig/SRR10357001_02.Average.bedGraph mean \
../../data/GSE139509/bigWig/SRR10357001.bw ../../data/GSE139509/bigWig/SRR10357002.bw
bedGraphToBigWig ../../data/GSE139509/bigWig/SRR10357001_02.Average.bedGraph ../../data/mm10/mm10.chrom.sizes \
../../data/GSE139509/bigWig/SRR10357001_02.Average.bw

wiggletools write_bg ../../data/GSE139509/bigWig/SRR10357005_06.Average.bedGraph mean \
../../data/GSE139509/bigWig/SRR10357005.bw ../../data/GSE139509/bigWig/SRR10357006.bw
bedGraphToBigWig ../../data/GSE139509/bigWig/SRR10357005_06.Average.bedGraph ../../data/mm10/mm10.chrom.sizes \
../../data/GSE139509/bigWig/SRR10357005_06.Average.bw


# Compare MECP2 ChIPseq with IgG sample (Boxer et al., 2020)
bigwigCompare   --bigwig1 ../../data/GSE139509/bigWig/SRR10356997_98.Average.bw \
                --bigwig2 ../../data/GSE139509/bigWig/SRR10357005_06.Average.bw \
                --outFileName ../../data/GSE139509/bigWig/SRR10356997_98_over_SRR10357005_06.bw \
                -p "max/2"

bigwigCompare   --bigwig1 ../../data/GSE139509/bigWig/SRR10357001_02.Average.bw \
                --bigwig2 ../../data/GSE139509/bigWig/SRR10357005_06.Average.bw \
                --outFileName ../../data/GSE139509/bigWig/SRR10357001_02_over_SRR10357005_06.bw \
                -p "max/2"

