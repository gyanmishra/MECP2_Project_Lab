#!/usr/bin/perl
use POSIX qw(strftime);

$projectDir = "../../../CUT_and_RUN";
$rawFile = "/archive/OBI/Neuroinformatics_Core/Stroud_lab/shared/Gyanstuff/";
$outputDir = "../../../CUT_and_RUN/results/";

# results directory
unless (-e "$projectDir/results") {
    `mkdir -p $projectDir/results`;
}
unless (-e "$outputDir/fastq") {
    `mkdir -p $outputDir/fastq`;
}
unless (-e "$outputDir/fastqc") {
    `mkdir -p $outputDir/fastqc`;
}
unless (-e "$outputDir/bam") {
    `mkdir -p $outputDir/bam`;
}
unless (-e "$outputDir/logs") {
    `mkdir -p $outputDir/logs`;
}
unless (-e "$outputDir/BigWig") {
    `mkdir -p $outputDir/BigWig`;
}
unless (-e "$outputDir/homer_tagDir") {
    `mkdir -p $outputDir/homer_tagDir`;
}


$datestring = strftime "%F", localtime;
$datestring =~ s/-//g ; 

$Job_dir  = '../../../CUT_and_RUN/CUT_and_RUN_' . $datestring;

unless (-e $Job_dir) {
    `mkdir $Job_dir`;
}


# provide a metadata file having two column 
# (1st column : Sample name and 2nd column : name of fastq file)
# or
# change index of 
# $col[3] = Sample Name 
# $col[4] = Fastq file 
# if sample name and fastq file name is in different column in sample metadata file

@file  = <> ;
my $header = shift (@file); 
foreach $line (@file) 
{
       if ($line !~ '#' || $line !~ 'Date'){
       
       chomp($line);
       @col = split('\t',$line);
       chomp($col[0]);      
        

       $col[5] =~ s/fastq.gz//g;
       $sampleName = $col[3];
       $sampleName =~ s/://g;

       $out_filename = $Job_dir . '/' . $sampleName . '.sh'; 

       print "$out_filename\n";
       open(FH, '>', $out_filename) or die $!;
       print FH "#!/bin/bash\n\n";
       print FH "#SBATCH --job-name CUT_and_RUN\n";
       print FH "#SBATCH -N 1\n";
       print FH "#SBATCH --partition=super\n";
       print FH "#SBATCH -t 0-12:0:0\n";
       print FH "#SBATCH -o $outputDir/logs/$sampleName.out\n";
       print FH "#SBATCH -o $outputDir/logs/$sampleName.err\n";
       print FH "#SBATCH --mail-type ALL\n";
       print FH "#SBATCH --mail-user GyanPrakash.Mishra@UTSouthwestern.edu\n\n\n";
       print FH "module load fastqc/0.11.8 parallel/20150122 perl/5.30.1 cutadapt/2.5 bowtie2/2.4.2 samtools/1.6 bedtools/2.29.2 homer/4.9 UCSC_userApps/v317 Trimmomatic/0.32\n\n";

       print FH "\n\n";
       print FH "# COMMAND GROUP 1\n\n";

       # Merge samples fastq
       print FH "cat $rawFile/$col[4]*_R1*fastq.gz >$projectDir/results/fastq/$sampleName.R1.fastq.gz\n\n";
       print FH "cat $rawFile/$col[4]*_R2*fastq.gz >$projectDir/results/fastq/$sampleName.R2.fastq.gz\n\n";     

       # Trim adapters
       print FH "cutadapt -a  AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATC \\\n-A  AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATC \\\n-m 20 \\\n-o $outputDir/results/fastq/$sampleName.R1.cutadapt.fastq.gz \\\n-p $outputDir/results/fastq/$sampleName.R2.cutadapt.fastq.gz \\\n$projectDir/results/fastq/$sampleName.R1.fastq.gz \\\n$projectDir/results/fastq/$sampleName.R2.fastq.gz\n\n";
       

       # FASTQC 
       print FH "fastqc $outputDir/fastq/$sampleName.R1.cutadapt.fastq.gz \\\n-o $outputDir/fastqc\n\n";
       print FH "fastqc $outputDir/fastq/$sampleName.R1.cutadapt.fastq.gz \\\n-o $outputDir/fastqc\n\n";

       # Align to mm10 reference genome
       print FH "bowtie2 -p 4 --local --very-sensitive-local --no-unal --no-mixed --no-discordant --phred33 -I 10 -X 700 \\\n-x ../../../CUT_and_RUN/data/mm10/index/mm10 \\\n-1 $outputDir/fastq/$sampleName.R1.cutadapt.fastq.gz \\\n-2 $outputDir/fastq/$sampleName.R2.cutadapt.fastq.gz | samtools view -bS -  \\\n>$outputDir/bam/$sampleName.bam\n\n";
       
       #
       print FH "samtools sort -@ 8 -n -O BAM \\\n$outputDir/bam/$sampleName.bam \\\n-o $outputDir/bam/$sampleName.s.bam\n\n";
       print FH "samtools fixmate -m \\\n$outputDir/bam/$sampleName.s.bam \\\n$outputDir/bam/$sampleName.s.fm.bam\n\n";
       print FH "samtools sort -@ 8 \\\n$outputDir/bam/$sampleName.s.fm.bam \\\n-o $outputDir/bam/$sampleName.s.fm.sort.bam\n\n";
       
       # Mark Duplicates 
       print FH "samtools markdup -r -s \\\n$outputDir/bam/$sampleName.s.fm.sort.bam \\\n$outputDir/bam/$sampleName.s.fm.sort.md.bam\n\n";
       
       # Name sort bam file 
       print FH "samtools sort -@ 8 -n -O BAM \\\n$outputDir/bam/$sampleName.s.fm.sort.md.bam \\\n-o $outputDir/bam/$sampleName.final.bam\n\n";
       
       # makeTagDirectory
       print FH "makeTagDirectory \\\n$outputDir/homer_tagDir/$sampleName \\\n$outputDir/bam/$sampleName.final.bam\n\n";
       
       # makeUCSC bedgraph
       print FH "makeUCSCfile \\\n$outputDir/homer_tagDir/$sampleName \\\n-o $outputDir/bigWig/$sampleName.bedGraph\n\n";

       # makeUCSC bigWig
       print FH "makeUCSCfile \\\n$outputDir/homer_tagDir/$sampleName \\\n-bigWig ../data/mm10/mm10.chrom.sizes \\\n-o $outputDir/bigWig/$sampleName.bw\n\n";
          
       # Unzip bedGraph 
       print FH "/usr/bin/gunzip $outputDir/bigWig/$sampleName.bedGraph.gz\n\n";
       
       # Convert bedGraph to bed 
       print FH "grep -v track \\\n$outputDir/bigWig/$sampleName.bedGraph \\\n>$outputDir/bigWig/$sampleName.bed\n";

       #`sbatch $out_filename`;
       }
    }




