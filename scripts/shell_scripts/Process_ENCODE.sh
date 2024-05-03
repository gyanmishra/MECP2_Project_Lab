

module load UCSC_userApps/v317
module load parallel/20150122
module load bedops/2.4.14
module load deeptools/3.5.0
module load wiggletools

# download 8wk cortex encode datasets from UCSC 
# bigWig (Histone)
wget http://hgdownload.soe.ucsc.edu/goldenPath/mm9/encodeDCC/wgEncodeLicrHistone/wgEncodeLicrHistoneCortexH3k27acMAdult8wksC57bl6StdSig.bigWig \
-P ../../data/ENCODE/
wget http://hgdownload.soe.ucsc.edu/goldenPath/mm9/encodeDCC/wgEncodeLicrHistone/wgEncodeLicrHistoneCortexH3k4me3MAdult8wksC57bl6StdSig.bigWig \
-P ../../data/ENCODE/
wget http://hgdownload.soe.ucsc.edu/goldenPath/mm9/encodeDCC/wgEncodeLicrHistone/wgEncodeLicrHistoneCortexH3k4me1MAdult8wksC57bl6StdSig.bigWig \
-P ../../data/ENCODE/
wget http://hgdownload.soe.ucsc.edu/goldenPath/mm9/encodeDCC/wgEncodeLicrHistone/wgEncodeLicrHistoneCortexInputMAdult8wksC57bl6StdSig.bigWig \
-P ../../data/ENCODE/

# broadPeak (Histone)
wget http://hgdownload.soe.ucsc.edu/goldenPath/mm9/encodeDCC/wgEncodeLicrHistone/wgEncodeLicrHistoneCortexH3k27acMAdult8wksC57bl6StdPk.broadPeak.gz \
-P ../../data/ENCODE/
wget http://hgdownload.soe.ucsc.edu/goldenPath/mm9/encodeDCC/wgEncodeLicrHistone/wgEncodeLicrHistoneCortexH3k4me3MAdult8wksC57bl6StdPk.broadPeak.gz \
-P ../../data/ENCODE/
wget http://hgdownload.soe.ucsc.edu/goldenPath/mm9/encodeDCC/wgEncodeLicrHistone/wgEncodeLicrHistoneCortexH3k4me1MAdult8wksC57bl6StdPk.broadPeak.gz \
-P ../../data/ENCODE/

# Ctcf (TFs)
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm9/encodeDCC/wgEncodeLicrTfbs/wgEncodeLicrTfbsCortexCtcfMAdult8wksC57bl6StdSig.bigWig \
-P ../../data/ENCODE/
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm9/encodeDCC/wgEncodeLicrTfbs/wgEncodeLicrTfbsCortexCtcfMAdult8wksC57bl6StdSig.bigWig \
-P ../../data/ENCODE/

# ENCODE 8week wholebrain DHS site
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm9/encodeDCC/wgEncodeUwDnase/wgEncodeUwDnaseWbrainC57bl6MAdult8wksPkRep1.narrowPeak.gz \
-P ../../data/ENCODE/
cut -f 1,2,3 ../../data/ENCODE/wgEncodeUwDnaseWbrainC57bl6MAdult8wksPkRep1.narrowPeak \
>../../data/ENCODE/wgEncodeUwDnaseWbrainC57bl6MAdult8wksPkRep1.narrowPeak.bed
liftOver ../../data/ENCODE/wgEncodeUwDnaseWbrainC57bl6MAdult8wksPkRep1.narrowPeak.bed \
../../data/mm10/mm9ToMm10.over.chain \
../../data/ENCODE/wgEncodeUwDnaseWbrainC57bl6MAdult8wksPkRep1.narrowPeak.mm10.bed \
../../data/ENCODE/wgEncodeUwDnaseWbrainC57bl6MAdult8wksPkRep1.narrowPeak.unMapped

####################
# lift Stroud et al cortex data
ls ../../data/Stroud_2017_Cell_data/8wk_WT_CTX_H3K*.bed |parallel --verbose 'liftOver {} ../../data/mm10/mm9ToMm10.over.chain {.}.mm10.bed {}.unMapped'
ls ../../data/Stroud_2017_Cell_data/8wk_WT_CTX_H3K*.mm10.bed |parallel --verbose 'grep -v 'chrM' {} >{.}noChrM.bed'
########################################################################################################################################################
# Convert bigwig to bedGraph
ls ../../data/ENCODE/wgEncodeLicrHistoneCortex*.bigWig | parallel --verbose  'bigWigToBedGraph {} {.}.mm9.bedGraph'

# Liftover from mm9 to mm10
ls ../../data/ENCODE/wgEncodeLicrHistoneCortex*.mm9.bedGraph | parallel --verbose  \
'liftOver {} ../../data/mm10/mm9ToMm10.over.chain {= s/.mm9.bedGraph// =}.mm10.bedGraph {.}.unMapped'


perl Download_ENCODE_data.pl \
../../data/ENCODE/Histone_brain_mm10_ENCODE_experiment_report.tsv  \
../../data/ENCODE/Histone_brain_mm10_ENCODE.txt Histone

perl Download_ENCODE_data.pl \
../../data/ENCODE/TFs_experiment_report_2024_4_10_17h_45m.tsv  \
../../data/ENCODE/TFs_brain_mm10_ENCODE.txt TFs
