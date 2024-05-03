#!/bin/bash

#SBATCH --job-name bamCov

# Name of the SLURM partition that this job should run on.
#SBATCH -p super       # partition (queue)
# Number of nodes required to run this job
#SBATCH -N 1

# Memory (RAM) requirement/limit in MB.
#SBATCH --mem 254976      # Memory Requirement (MB)

# Time limit for the job in the format Days-H:M:S
# A job that reaches its time limit will be cancelled.
# Specify an accurate time limit for efficient scheduling so your job runs promptly.
#SBATCH -t 0-2:0:0

# The standard output and errors from commands will be written to these files.
# %j in the filename will be replace with the job number when it is submitted.
#SBATCH -o ../../results/log/bamCoverage.out
#SBATCH -e ../../results/log/bamCoverage.err

# Send an email when the job status changes, to the specfied address.
#SBATCH --mail-type ALL
#SBATCH --mail-user GyanPrakash.Mishra@UTSouthwestern.edu

module load parallel/20150122
module load bedtools/2.29.2
module load deeptools/3.5.0
module load bedops/2.4.14

# Generate bam coverage for CUT&TAG data
ls ../../../CUT_and_RUN/results/bam/*.final.bam | grep 'CNT' | \
parallel --verbose --jobs 1 'samtools sort -@ 20 {} -o {.}.sort.bam'
ls ../../../CUT_and_RUN/results/bam/*.final.sort.bam | grep 'CNT' | \
parallel --verbose --jobs 1 'samtools index {}'

mkdir -p ../../results/bamCoverage_out
# CPM
parallel --verbose --jobs 1 \
'bamCoverage --bam {} -o ../../results/bamCoverage_out/{/.}.CPM.bw \
--binSize 10 \
--effectiveGenomeSize 2652783500 \
--normalizeUsing CPM \
-p "max/2"'

# RPKM 
ls ../../../CUT_and_RUN/results/bam/*.final.sort.bam | grep 'CNT' | \
parallel --verbose \
'bamCoverage --bam {} -o ../../results/bamCoverage_out/{/.}.RPKM.bw \
--binSize 1000 \
--effectiveGenomeSize 2652783500 \
--normalizeUsing RPKM \
-p "max/2"'

# generate bed file coverage from makeUCSC output
ls ../../../CUT_and_RUN/results/BigWig/*.bw | parallel --verbose --jobs 5 'bigWigToWig {} {.}.wig'
ls ../../../CUT_and_RUN/results/BigWig/*.wig | parallel --verbose --jobs 5 'wig2bed --zero--indexed < {} > {.}.bed'
ls ../../../CUT_and_RUN/results/BigWig/*.wig | parallel --verbose 'rm {}'


# generate GeneBody Coverage from bed Coverage files. 
ls ../../results/bamCoverage_out/*.bed |grep 'CNT' |\
parallel --verbose 'bedtools map -a \
../../../CUT_and_RUN/data/mm10/gencode.vM10.annotation_proteinCoding_gene_greater_than_4.5kb.bed \
-b {} -c 5 -o sum \
>../results/CNR_Coverage/{/.}_RPKM.bed'





