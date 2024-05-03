#!/usr/bin/perl
use POSIX qw(strftime);

$projectDir = "../../../CUT_and_RUN";

# results directory
unless (-e "$projectDir") {
    `mkdir $projectDir`;
}

unless (-e "$projectDir/results") {
    `mkdir $projectDir/results`;
}


$datestring = strftime "%F", localtime;
$datestring =~ s/-//g ; 

$Job_dir  = '../../../CUT_and_RUN/' . $datestring . '_macs2';

unless (-e $Job_dir) {
    `mkdir $Job_dir`;
}

@file  = <> ; # Enter file with two column (1st column = target bam file; 2nd column = control bam file)

foreach $line (@file) 
{
    if ($line !~ '#'){
        
    chomp($line);
    @col = split('\t',$line);
    chomp($col[0]);
    #print "$col[0]\n";

    $filename = $col[0];

    
    $filename =~s/.final.bam//g;

    #print $col[1];
    $col[1] =~ s/.fastq.gz//g;
    #print "$col[1]\n";

    $out_filename = $Job_dir . '/' . $filename . '.sh';

    
    open(FH, '>', $out_filename) or die $!;
    print FH "#!/bin/bash\n\n";
    print FH "#SBATCH --job-name macs2\n";
    print FH "#SBATCH -N 1\n";
    print FH "#SBATCH --partition=super\n";
    print FH "#SBATCH -t 0-12:0:0\n";
    print FH "#SBATCH -o $Job_dir/$filename.out\n";
    print FH "#SBATCH -o $Job_dir/$filename.err\n";
    print FH "#SBATCH --mail-type ALL\n";
    print FH "#SBATCH --mail-user GyanPrakash.Mishra@UTSouthwestern.edu\n\n\n";
    print FH "module load parallel/20150122 perl/5.30.1 cutadapt/2.5 bowtie2/2.4.2 samtools/1.6 bedtools/2.29.2 homer/4.9 UCSC_userApps/v317 macs/2.1.0-20151222 \n\n";

        print FH "macs2 callpeak -t $projectDir/results/$col[0] -c $projectDir/results/$col[1] -f BAM -g mm --outdir $filename --call-summits -n $filename\n";

    `sbatch $out_filename`;
    }
}

