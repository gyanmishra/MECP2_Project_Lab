#!/bin/bash

#SBATCH --job-name macs2
#SBATCH -N 1
#SBATCH --partition=super
#SBATCH -t 0-12:0:0
#SBATCH -o ../../../CUT_and_RUN/results/logs/MECP2_MECP2WT_vs_IgG_and_MECP2KO.out
#SBATCH -o ../../../CUT_and_RUN/results/logs/MECP2_MECP2WT_vs_IgG_and_MECP2KO.err
#SBATCH --mail-type ALL
#SBATCH --mail-user GyanPrakash.Mishra.edu


module load parallel/20150122 
module load macs/2.1.0-20151222 

# MECP2 (MECP2WT_vs_IgG)
macs2 callpeak \
-t ../../../CUT_and_RUN/results/bam/9wk_CTX_MECP2WT_Mecp2_CNR_HS032519_N032619.final.bam \
   ../../../CUT_and_RUN/results/bam/9wk_CTX_DNMT3AWT_Mecp2_CNR_HS031219_N031419.final.bam \
   ../../../CUT_and_RUN/results/bam/7wk_CTX_MECP2WT_Mecp2_CNR_HS031219_N031419.final.bam \
   ../../../CUT_and_RUN/results/bam/Adult_CTX_WT_Mecp2_CNR_HS020619_N020719.final.bam \
-c ../../../CUT_and_RUN/results/bam/7wk_CTX_MECP2WT_IgG_CNR_HS031219_N031419.final.bam \
   ../../../CUT_and_RUN/results/bam/9wk_CTX_DNMT3AWT_IgG_CNR_HS031219_N031419.final.bam \
   ../../../CUT_and_RUN/results/bam/9wk_CTX_MECP2WT_IgG_CNR_HS032519_N032619.final.bam \
   ../../../CUT_and_RUN/results/bam/Adult_CTX_WT_IgG_CNR_HS020619_N020719.final.bam \
-f BAM \
-g mm \
--outdir ../../../CUT_and_RUN/results/macs2_Callpeak/Adult_MECP2_MECP2WT_vs_IgG_macs2 \
--call-summits \
-n Adult_MECP2_MECP2WT_vs_IgG

# MECP2 (MECP2 WT vs MECP2 KO)
macs2 callpeak \
-t ../../../CUT_and_RUN/results/bam/9wk_CTX_MECP2WT_Mecp2_CNR_HS032519_N032619.final.bam \
   ../../../CUT_and_RUN/results/bam/9wk_CTX_DNMT3AWT_Mecp2_CNR_HS031219_N031419.final.bam \
   ../../../CUT_and_RUN/results/bam/7wk_CTX_MECP2WT_Mecp2_CNR_HS031219_N031419.final.bam \
   ../../../CUT_and_RUN/results/bam/Adult_CTX_WT_Mecp2_CNR_HS020619_N020719.final.bam \
-c ../../../CUT_and_RUN/results/bam/7wk_CTX_MECP2KO_Mecp2_CNR_HS031219_N031419.final.bam \
   ../../../CUT_and_RUN/results/bam/9wk_CTX_MECP2KO_Mecp2_CNR_HS032519_N032619.final.bam \
-f BAM \
-g mm \
--outdir ../../../CUT_and_RUN/results/macs2_Callpeak/Adult_MECP2_MECP2WT_vs_KO_Macs2 \
--call-summits \
-n Adult_MECP2_MECP2WT_vs_KO



# MECP2 (MECP2 KO vs MECP2 WT)
macs2 callpeak \
-t ../../../CUT_and_RUN/results/bam/7wk_CTX_MECP2KO_Mecp2_CNR_HS031219_N031419.final.bam \
   ../../../CUT_and_RUN/results/bam/9wk_CTX_MECP2KO_Mecp2_CNR_HS032519_N032619.final.bam \
-c ../../../CUT_and_RUN/results/bam/9wk_CTX_MECP2WT_Mecp2_CNR_HS032519_N032619.final.bam \
   ../../../CUT_and_RUN/results/bam/9wk_CTX_DNMT3AWT_Mecp2_CNR_HS031219_N031419.final.bam \
   ../../../CUT_and_RUN/results/bam/7wk_CTX_MECP2WT_Mecp2_CNR_HS031219_N031419.final.bam \
   ../../../CUT_and_RUN/results/bam/Adult_CTX_WT_Mecp2_CNR_HS020619_N020719.final.bam \
-f BAM \
-g mm \
--outdir ../../../CUT_and_RUN/results/macs2_Callpeak/Adult_MECP2_MECP2KO_vs_WT_Macs2 \
--call-summits \
-n Adult_MECP2_MECP2KO_vs_WT


# MECP2 (MECP2 KO vs IgG)
macs2 callpeak \
-t ../../../CUT_and_RUN/results/bam/7wk_CTX_MECP2KO_Mecp2_CNR_HS031219_N031419.final.bam \
   ../../../CUT_and_RUN/results/bam/9wk_CTX_MECP2KO_Mecp2_CNR_HS032519_N032619.final.bam \
-c ../../../CUT_and_RUN/results/bam/7wk_CTX_MECP2WT_IgG_CNR_HS031219_N031419.final.bam \
   ../../../CUT_and_RUN/results/bam/9wk_CTX_DNMT3AWT_IgG_CNR_HS031219_N031419.final.bam \
   ../../../CUT_and_RUN/results/bam/9wk_CTX_MECP2WT_IgG_CNR_HS032519_N032619.final.bam \
   ../../../CUT_and_RUN/results/bam/Adult_CTX_WT_IgG_CNR_HS020619_N020719.final.bam \
-f BAM \
-g mm \
--outdir ../../../CUT_and_RUN/results/macs2_Callpeak/Adult_MECP2_MECP2KO_vs_IgG_Macs2 \
--call-summits \
-n Adult_MECP2_MECP2KO_vs_IgG





