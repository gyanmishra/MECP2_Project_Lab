#!/usr/bin/perl
use POSIX qw(strftime);

$SRAfile = $ARGV[0];
$GEO_accession = $ARGV[1];
$projectDir = "../data" . $GEO_accession;

# results directory
unless (-e "$projectDir/fastqc") {
    `mkdir -p $projectDir/fastqc`;
}
unless (-e "$projectDir/bam") {
    `mkdir -p $projectDir/bam`;
}
unless (-e "$projectDir/logs") {
    `mkdir -p $projectDir/logs`;
}
unless (-e "$projectDir/bigWig") {
    `mkdir -p $projectDir/bigWig`;
}
unless (-e "$projectDir/homer_tagDir") {
    `mkdir -p $projectDir/homer_tagDir`;
}

#print "$SRAfile\n" ;
#print "$GEO_accession\n" ; 
#print "$projectDir\n" ; 


$datestring = strftime "%F", localtime;
$datestring =~ s/-//g ; 

$Job_dir  = "../data/" . $GEO_accession . '/' . $GEO_accession .'_' .$datestring;

unless (-e $Job_dir) {
    `mkdir $Job_dir`;
}

@file  = `cat $SRAfile`;
foreach $line (@file) 
{
    if($line !~ 'Run'){
    chomp($line);
    @col = split(',',$line);
    chomp($col[0]);
    print "$col[0]\n";
    # provide a file having two column (1st column : Sample name and 2nd column : name of fastq file)
    $out_filename = $Job_dir . '/' . $col[0] . '.sh';

    print "$out_filename\n";
    open(FH, '>', $out_filename) or die $!;
    print FH "#!/bin/bash\n\n";
    print FH "#SBATCH --job-name MECP2_ChIPseq\n";
    print FH "#SBATCH -N 1\n";
    print FH "#SBATCH --partition=super\n";
    print FH "#SBATCH -t 0-12:0:0\n";
    print FH "#SBATCH -o $projectDir/logs/$col[0].out\n";
    print FH "#SBATCH -o $projectDir/logs/$col[0].err\n";
    print FH "#SBATCH --mail-type ALL\n";
    print FH "#SBATCH --mail-user GyanPrakash.Mishra@UTSouthwestern.edu\n\n\n";
    print FH "module load fastqc/0.11.8 parallel/20150122 perl/5.30.1 cutadapt/2.5 bowtie2/2.4.2 samtools/1.6 bedtools/2.29.2 homer/4.9 UCSC_userApps/v317 Trimmomatic/0.32 sra_toolkit/3.0.0 \n\n";
    #print FH "perl ../../../CUT_and_RUN/scripts/CNR_pipeline.pl ../../../CUT_and_RUN/data/CUTandRUN_seq_Data.txt\n\n";
    print FH "\n\n";
    print FH "# COMMAND GROUP 1\n\n";

    print FH "prefetch $col[0] -O $projectDir/\n";
    print FH "fastq-dump --split-3 $projectDir/$col[0]/$col[0].sra -O $projectDir/$col[0]/\n\n";
    # Merge samples fastq
    #print FH "cat $rawFile/$col[1]*_R1*fastq.gz >$projectDir/results/fastq/$col[0]_R1.fastq.gz\n";
    #print FH "cat $rawFile/$col[1]*_R2*fastq.gz >$projectDir/results/fastq/$col[0]_R2.fastq.gz\n";
    
    # Trim adapters
    #print FH "cutadapt -a  AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATC -A  AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATC -m 20 -o $projectDir/results/fastq/$col[4]_R1.cutadapt.fastq.gz -p $projectDir/results/fastq/$col[4]_R2.cutadapt.fastq.gz $rawFile/$col[5] $rawFile/$col[6]\n";
    #
    #print FH "java -jar \$Trimmomatic PE  $rawFile/$col[5] $rawFile/$col[6] $projectDir/results/fastq/$col[4]_R1.paired.fastq.gz $projectDir/results/fastq/$col[4]_R1.unpaired.fastq.gz $projectDir/results/fastq/$col[4]_R2.paired.fastq.gz $projectDir/results/fastq/$col[4]_R2.unpaired.fastq.gz ILLUMINACLIP:../../../tools/trimmomatic/adapters/TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 SLIDINGWINDOW:4:30 MINLEN:20\n\n";
    ## FASTQC 
    print FH "fastqc $projectDir/$col[0]/$col[0].fastq -o $projectDir/fastqc\n";
#
    ## Align to mm10 reference genome
    print FH "bowtie2 -p 4 --local --very-sensitive-local --no-unal --no-mixed --no-discordant --phred33 -I 10 -X 700 -x ../../../CUT_and_RUN/data/mm10/index/mm10 -U $projectDir/$col[0]/$col[0].fastq | samtools view -bS -  >$projectDir/bam/$col[0].bam\n";
    ##
    print FH "samtools sort -@ 8 -n -O BAM $projectDir/bam/$col[0].bam -o $projectDir/bam/$col[0].s.bam\n";
    print FH "samtools fixmate -m $projectDir/bam/$col[0].s.bam $projectDir/bam/$col[0].s.fm.bam\n";
    print FH "samtools sort -@ 8 $projectDir/bam/$col[0].s.fm.bam -o $projectDir/bam/$col[0].s.fm.sort.bam\n";
    ###
    ### Mark Duplicates 
    print FH "samtools markdup -r -s $projectDir/bam/$col[0].s.fm.sort.bam $projectDir/bam/$col[0].s.fm.sort.md.bam\n";
    ##
    ### Name sort bam file 
    print FH "samtools sort -@ 8 -n -O BAM $projectDir/bam/$col[0].s.fm.sort.md.bam -o $projectDir/bam/$col[0].final.bam\n";
    ##
    ### makeTagDirectory
    print FH "makeTagDirectory $projectDir/homer_tagDir/$col[0] $projectDir/bam/$col[0].final.bam\n";
    #
    ## makeUCSC bigWig
    print FH "makeUCSCfile $projectDir/homer_tagDir/$col[0] -bigWig ../../data/mm10/mm10.chrom.sizes -o $projectDir/bigWig/$col[0].bw\n\n";
    #print FH "makeUCSCfile $projectDir/results/$col[0].tagDir -o $projectDir/results/$col[0].bedGraph\n\n";

    print FH "#END OF SCRIPT\n";
    
    `sbatch $out_filename`;
    }
}

