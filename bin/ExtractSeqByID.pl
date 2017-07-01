#!/usr/bin/env perl

#################################################################################
##                                                                             ##
##                       Quasispecies Analysis Package                         ##
##                                                                             ##
#################################################################################
##                                                                             ##
##  A software suite designed for virus quasispecies analysis                  ##
##  See our website: <http://bioinfo.rjh.com.cn/labs/jhuang/tools/gap/>        ##
##                                                                             ##
##  Version 1.0                                                                ##
##                                                                             ##
##  Copyright (C) 2017 by Mingjie Wang, All rights reserved.                   ##
##  Contact:  huzai@sjtu.edu.cn                                                ##
##  Organization: Research Laboratory of Clinical Virology, Rui-jin Hospital,  ##
##  Shanghai Jiao Tong University, School of Medicine                          ##
##                                                                             ##
##  This file is a subprogram of GAP suite.                                    ##
##                                                                             ##
##  QAP is a free software; you can redistribute it and/or                     ##
##  modify it under the terms of the GNU General Public License                ##
##  as published by the Free Software Foundation; either version               ##
##  3 of the License, or (at your option) any later version.                   ##
##                                                                             ##
##  QAP is distributed in the hope that it will be useful,                     ##
##  but WITHOUT ANY WARRANTY; without even the implied warranty of             ##
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              ##
##  GNU General Public License for more details.                               ##
##                                                                             ##
##  You should have received a copy of the GNU General Public                  ##
##  License along with ViralFusionSeq; if not, see                             ##
##  <http://www.gnu.org/licenses/>.                                            ##
##                                                                             ##
#################################################################################

use diagnostics;
use strict;
use warnings;
use FindBin qw/$RealBin/;
use File::Spec;
use Getopt::Long;
use Pod::Usage;
use lib "$FindBin::Bin/../lib";
use Cwd qw/getcwd abs_path/;
use File::Basename;


####Use modules in this program####
use General;

####---------------------------####
####The program begins here
####---------------------------####

##Show welcome
print "You are now running subprogram: ";
printcol ("ExtractSeqByID","green");
print "\n";

##get workding directory
my $wk_dir = getcwd;
my $mainBin;
if ($RealBin =~ /(.*)\/bin/){
	$mainBin = $1;
}

####define command line arguments
my $help;
my $fq1;
my $fq2;
my $outputDir;
my $pairEnded = 0;
my $idFile;
my $format;

my $DateNow = `date +"%Y%m%d_%Hh%Mm%Ss"`;
chomp $DateNow;

GetOptions(
'1|fastq1|=s'      => \$fq1,
'2|fastq2|=s'      => \$fq2,
'o|outputDir=s'    => \$outputDir,
'h|help|'          => \$help,
'i|readID|=s'       => \$idFile,
'f|outFormat|=s'    => \$format
);



##check command line arguments
if (defined $help){
	pod2usage(-verbose=>2,-exitval=>1);
}


if (defined $outputDir){
	$outputDir = abs_path($outputDir) . "/";
	if (not -e $outputDir){
 		InfoWarn("The output directory $outputDir does NOT exist.",'yellow');
 		InfoWarn("Will mkdir $outputDir and use it as the output directory.",'yellow');
		#pod2usage(-verbose=>0,-exitval=>1);
		#exit;
		if (!-e $outputDir){
			my $cmd = "mkdir -p $outputDir";
			system($cmd);
		}else{
			InfoError("Mkdir Failed! Folder $outputDir already exists!","red");
			InfoError("Please specify another output directory using option -o/--outputDir");
			pod2usage(-verbose=>0,-exitval=>1);
			exit;
		}
	}
}else{
	$outputDir = File::Spec -> catfile($wk_dir,"qap_Results_for_ExtractSeqByID_$DateNow");
	InfoWarn("The output directory is not provided!",'yellow');
	InfoWarn("Will mkdir \"$outputDir\" and use it as the output directory.",'yellow');
	
	if (!-e "$outputDir"){
		my $cmd = "mkdir -p $outputDir";
		system($cmd);
	}else{
		InfoError("Mkdir Failed! $outputDir already exists!","red");
		InfoError("Please specify another output directory using option -o/--outputDir\n");
		pod2usage(-verbose=>0,-exitval=>1);
		exit;
	}

}

if (defined $format){
	if ($format eq 'fastq' || $format eq 'fasta'){
		#nothing
	}else{
		InfoError("The output file format should be one of \"fastq\" or \"fasta\".","red");
		exit;
	}
}else{
	$format = "fastq";
	Info("Output format is not provided. Using \"fastq\" as default.");
}

if (defined $idFile){
	if (not -e $idFile){
		InfoError("Input read ID file $idFile does NOT exist!",'red');
		pod2usage(-verbose=>0,-exitval=>1);
		exit;
	}
	$idFile = abs_path($idFile);
}else{
	InfoError("Input read ID file must be provided!",'red');
	pod2usage(-verbose=>1,-exitval=>1);
	exit;
}

if (defined $fq1){
	if (not -e $fq1){
		InfoError("Input fastq1 file $fq1 does NOT exist!",'red');
		pod2usage(-verbose=>0,-exitval=>1);
		exit;
	}
	$fq1 = abs_path($fq1);
}else{
	InfoError("Input fastq1 file must be provided!",'red');
	pod2usage(-verbose=>1,-exitval=>1);
	exit;
}

if (defined $fq2){
	if (not -e $fq2){
		InfoError("Input fastq2 file $fq2 does NOT exist!",'red');
		pod2usage(-verbose=>0,-exitval=>1);
		exit;
	}
	$fq2 = abs_path($fq2);
	$pairEnded = 1;
}

##start the analysis
my %id;
open T,"$idFile" or die "Can NOT open read ID file:$idFile:$!";
while(<T>){
	chomp;
	$id{$_} = 1;
}
close T;

##core program
&extractSeq($fq1,\%id,$outputDir,$format);
&extractSeq($fq2,\%id,$outputDir,$format) if $pairEnded;

##subprograms goes here
sub extractSeq {
	my $fq = shift;
	my $id = shift;
	my $outdir = shift;
	my $format = shift;
	
	my $fqName = basename($fq1);
	my %idWanted = %$id;
	
	Info("Extracting sequences from $fq.");
	
	##handle gzipped files
	if (isGzipped($fq)){
		my $tmpName = addTagForRawData($fqName, 'ungzipped', 1) . "_" . time();
		my $tmpout = File::Spec -> catfile($outdir,$tmpName);
		
		#uncompress
		Info("Uncompressing gzipped file $fqName");
		my $cmd = "gunzip -c $fq > $tmpout";
		system($cmd);
		#runcmd($cmd);
		$fq = $tmpout;
	}
	
	#output file name
	my $outName = addTagForRawData($fqName, 'extractedSeq', 1);
	$outName =~ s/\.gz$//; #remove .gz suffix is input fastq file is gzipped.
	
	my $outFile;
	if ($format eq 'fastq'){
		$outFile = File::Spec -> catfile($outdir, $outName);
	}elsif ($format eq 'fasta'){
		$outFile = File::Spec -> catfile($outdir, changeFastqSuffix2Fasta($outName));
	}else{
		InfoError("The output file format should be one of \"fastq\" or \"fasta\".","red");
		exit;
	}
	open RES,">$outFile" or die "Can NOT output extracted sequences data to file $outFile:$!\n";
	
	#get total line number
	my $wc = `wc -l $fq`;
	$wc =~ /^(\d+) /;
	my $totalLineNum = $1;
	my $totalSeqNum = int($totalLineNum / 4);
	my $hashNum = scalar(values(%idWanted));
	my $estimateTime = log($hashNum) * 10 * ($totalLineNum / 4000000);
	print(log($hashNum));
	if ($estimateTime > 5){
		#nothing
	}else{
		$estimateTime = int(rand() * 5) + 2;
	}
	my $estimateTimeExp = formatMinutes($estimateTime);
	
	#read in the input fastq file
	Info("Start to extract, this might take a long time. Please wait patiently. [ETA $estimateTimeExp]");
	if ($format eq 'fastq'){
		open T,$fq or die "Can NOT open the input fastq file:$!\n";
		my $i = 1;
		while(my $line1 = <T>){
			chomp $line1;
			chomp (my $line2 = <T>);
			chomp (my $line3 = <T>);
			chomp (my $line4 = <T>);
		
			$line1 =~ /(.*?)\s+/;
			$id = $1;
			
			if (exists $idWanted{$id}){
					print RES "$line1\n$line2\n$line3\n$line4\n";
			}
			InfoProcessBar($i, $totalSeqNum);
			$i++;
		}
		close T;
	}else{
		open T,$fq or die "Can NOT open the input fastq file:$!\n";
		my $i = 1;
		while(my $line1 = <T>){
			chomp $line1;
			chomp (my $line2 = <T>);
			chomp (my $line3 = <T>);
			chomp (my $line4 = <T>);
		
			$line1 =~ /(.*?)\s+/;
			$id = $1;
			
			if (exists $idWanted{$id}){
					$line1 =~ s/^\@//;
					print RES ">$line1\n$line2\n";
			}
			InfoProcessBar($i, $totalSeqNum);
			$i++;
		}
		close T;
	}
}

##run success
Info("Program completed!",'green');


####---------------------------####
####The program ends here
####---------------------------####


=pod 

=head1 NAME

qap -- Quasispecies analysis package

=head1 SYNOPSIS


       ______       ______       ______
      / ___  |     / ____ \     |  ___ \
     / /   | |    / |    | |    | |   \ \
    | |    | |    | |    |_|    | |    | |
     \ \___| |    \ |____\ \    | |___/ /
      \____  |     \_____/\_\   | |____/
           | |                  | |
           | |                  | |
           | |                  | |
           |_|                  |_|         v1.0
           



gap ExtractSeqByID [options]

Use --help to see more information.

gap is still in development. If you have encounted any problem in usage, please feel no hesitation to cotact us.

=head1 DESCRIPTION

This script implements a function for extracting sequences from fastq files with read IDs. The script has B<several> mandatory options that MUST appear last. 

=head1 OPTIONS

=over 5

=item --fastq1,-1 F<FILE> [Required]

Path to next generation sequencing raw data. REQUIRED for both single-end or paired-end reads. Both compressed files and uncompressed files are allowed. 

=item --fastq2,-2 F<FILE> [Optional]

Path to next generation sequencing raw data. REQUIRED for both paired-end reads. Both compressed files and uncompressed files are allowed. 

=item --idFile, -i F<FILE> [Required]

Path to the file containing the read IDs by using which the sequences would be extracted. Each ID should be listed per line.

=item --outFormat, -f F<STRING> [Optional]

The format of program output. The value should be one of 'fastq' or 'fasta'.

=item --outputDir,-o F<FILE> [Optional]

Path of the directory to storage result files. If NOT provided, the program will generate a folder automatically.

=item --help,-h

Display this detailed help information.

=back

=head1 EXAMPLE

=over 5

gap ExtractSeqByID -1 Data1_R1.fq.gz -2 Data1_R2.fq.gz -i readID.txt -f fastq -o ./seq

=back

=head1 AUTHOR

Mingjie Dr.Wang I<huzai@sjtu.edu.cn>

=head1 COPYRIGHT

Copyright (C) 2017, Mingjie Wang. All rights reserved.


