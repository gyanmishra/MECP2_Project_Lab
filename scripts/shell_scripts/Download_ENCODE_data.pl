#!/usr/bin/perl

use File::Basename;

#@file = `cat ../../data/ENCODE/Histone_brain_mm10_ENCODE_experiment_report.tsv`; # ENCODE metadata file
#@file2 = `cat ../../data/ENCODE/Histone_brain_mm10_ENCODE.txt`; # file with link for downloading ENCODE data
$file_1 = $ARGV[0];
$file_2 = $ARGV[1];
$flag = $ARGV[2]; 
$tissue = "";

if($flag == 'Histone'){
    $tissue = "forebrain"
}
if($flag == 'TFs'){
    $tissue = "brain"
}

@file = `cat $ARGV[0]`; # ENCODE metadata file
@file2 = `cat $ARGV[1]`; # file with link for downloading ENCODE data

# results directory
unless (-e "../../data/ENCODE") {
    `mkdir -p ../../data/ENCODE`;
}
foreach $line2 (@file2){
    chomp($line2);
    foreach $line (@file){
        chomp($line);
        @col = split("\t",$line);
        if($line2 =~ 'bed.gz'){
            #$link = $line2;
            $link = basename($line2);
            $link=~s/.bed.gz//g;
            #print "$link\n";
            if($col[9] =~ $tissue and $line =~ $link){
                #print "$col[1]\t$col[3]\t$col[4]\t$col[6]\t$col[9]\t$link\n";
                $fileName = $col[6] . '_' . $col[2] . '_' . $col[9] . '_' . $link . '.bed.gz';
                #print "$fileName\n";
                `wget $line2 -O ../../data/ENCODE/$fileName`;
                last;
            }
        }
        if($line2 =~ 'bigWig'){
            #$link = $line2;
            $link = basename($line2);
            $link=~s/bigWig//g;
            #print "$link\n";
            if($col[9] =~ 'brain' and $line =~ $link){
                #print "$col[1]\t$col[3]\t$col[4]\t$col[6]\t$col[9]\t$link\n";
                $fileName = $col[6] . '_' . $col[2] . '_' . $col[9] . '_' . $link . 'bigWig';
                #print "$fileName\n";
                `wget $line2 -O ../../data/ENCODE/$fileName`;
                last; 
            }
        }
    }
}
