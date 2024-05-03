#!/usr/bin/bash


module load parallel/20150122
module load bedtools/2.29.2


mkdir -p ../../results/HiC/

# Jiang.2017
awk '{print $1"\t"$2"\t"$3}' \
    ../../data/mm10_hic_northwerstern_YueLab/mm10/Jiang_2017.Neuron.mm10.peakachu-merged.loops  \
    >../../data/mm10_hic_northwerstern_YueLab/mm10/Jiang_2017.Neuron.mm10.peakachu-merged.loops_start.bed

awk '{print $4"\t"$5"\t"$5}' \
    ../../data/mm10_hic_northwerstern_YueLab/mm10/Jiang_2017.Neuron.mm10.peakachu-merged.loops \
    >>../../data/mm10_hic_northwerstern_YueLab/mm10/Jiang_2017.Neuron.mm10.peakachu-merged.loops_end.bed

cat ../../data/mm10_hic_northwerstern_YueLab/mm10/Jiang_2017.Neuron.mm10.peakachu-merged.loops_start.bed \
../../data/mm10_hic_northwerstern_YueLab/mm10/Jiang_2017.Neuron.mm10.peakachu-merged.loops_end.bed | \
sort -k1,1n -k2,2n \
>../../data/mm10_hic_northwerstern_YueLab/mm10/Jiang_2017.Neuron.mm10.peakachu-merged.loops.bed


# overlap chromatin loops with MBHs
bedtools intersect -wao -a ../../results/MECP2/MECP2WT.vs.KO_overlap_MECP2WT.vs.IgG.noBlacklist_filtered.bed \
-b ../../data/mm10_hic_northwerstern_YueLab/mm10/Jiang_2017.Neuron.mm10.peakachu-merged.loops.bed  \
>../../results/HiC/Jiang_2017.Neuron.mm10.peakachu-merged.loops_closest_MECP2.bed 