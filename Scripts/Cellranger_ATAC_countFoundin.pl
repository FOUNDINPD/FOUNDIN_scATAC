#!/usr/bin/perl -w


##=============================================================================================
## USAGE (not using custom reference genome as the web summary files will show warnings)
## perl Cellranger_ATAC_countFoundin.pl --sampleFile Sample_test.txt --refDirName /media/root/dataE/FOUNDIN_vikas/Database/refdata-cellranger-atac-GRCh38-1.2.0 --resultDir /media/root/dataE/FOUNDIN_vikas/scATACseq/CountsATAC_cellRanger/ --localMem 50
## 
## 
##
## sampleFile should be the tab sep 3 column file. First column contains the sample ID (not the full fastq file name but the first part) example D17-8753_S1_L001_I1_001.fastq.gz should be only written D17-8753 (notice part before first underscore) and second column is the directory with fastq files for that sample and the third column is the ID for output.
## refDirName is the directory contaning genome reference build by cellranger mkref. Example command
##
## cellranger-atac mkref GRCh38_foundinGTF_atacseq --config GRCh38_Foundin_atac.config
##
##
## Author:   Vikas Bansal
## Date:     07.06.2020
##============================================================================================
#
# 
# 




use warnings;
use strict;

use Getopt::Long;


my $sample_file;
my $result_dir;
my $ref_dir;
my $local_mem;

GetOptions(
    "sampleFile=s"   => \$sample_file,
    "resultDir=s"    => \$result_dir,
    "refDirName=s"    => \$ref_dir,
    "localMem=s"    => \$local_mem,	
);

open FILE, "$sample_file";

while(<FILE>){
	if(!/^\#/){
	   chomp;
	   (my $sampleID, my $sample_folder, my $out_file)=split();
	    

		 
		  print "cd $result_dir && /media/root/dataE/FOUNDIN_vikas/Tools/cellranger-atac-1.2.0/cellranger-atac count --id=$out_file --reference=$ref_dir --fastqs=$sample_folder --sample=$sampleID --localmem=$local_mem > $out_file.err\n" ;

	} 
}

