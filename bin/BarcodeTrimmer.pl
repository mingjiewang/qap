#!/usr/bin/env perl

#################################################################################
##                                                                             ##
##                       Quasispecies Analysis Package                         ##
##                                                                             ##
#################################################################################
##                                                                             ##
##  A software suite designed for virus quasispecies analysis                  ##
##  See our website: <http://bioinfo.rjh.com.cn/labs/jhuang/tools/qap/>        ##
##                                                                             ##
##  Version 1.0                                                                ##
##                                                                             ##
##  Copyright (C) 2017 by Mingjie Wang, All rights reserved.                   ##
##  Contact:  huzai@sjtu.edu.cn                                                ##
##  Organization: Research Laboratory of Clinical Virology, Rui-jin Hospital,  ##
##  Shanghai Jiao Tong University, School of Medicine                          ##
##                                                                             ##
##  This file is a subprogram of QAP suite.                                    ##
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
##  License along with QAP; if not, see                                        ##
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
printcol ("BarcodeTrimmer","green");
print "\n";

## check threads available or not
$| = 1;
sleep(1);
my $threads_usable = eval 'use threads; 1';
if ($threads_usable) {
	use threads;
	use threads::shared;
} else {
	Info("No threading is possible. Please install perl module: threads or recompile perl with option -Dusethreads","red");
}

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
my $matchPosition;
my $barcodeLen;
my $outputDir;
my $pairEnded = 0;
my $threads;
my $fasta;
my $format;

my $DateNow = `date +"%Y%m%d_%Hh%Mm%Ss"`;
chomp $DateNow;

GetOptions(
'1|fastq1|=s'        => \$fq1,
'2|fastq2|=s'        => \$fq2,
'f|fasta|=s'         => \$fasta,
'o|outputDir=s'      => \$outputDir,
'h|help|'            => \$help,
'l|barcodeLength|=s' => \$barcodeLen,
'p|matchPosition|=s' => \$matchPosition,
't|threads=s'        => \$threads
);


##check command line arguments
if (defined $help){
	pod2usage(-verbose=>2,-exitval=>1);
}

if(scalar(@ARGV) == 0){
	pod2usage(-verbose=>1,-exitval=>1);
}

if (defined $outputDir){
	$outputDir =~ s/\/$//;
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
	$outputDir = File::Spec -> catfile($wk_dir,"qap_Results_for_BarcodeTrimmer_$DateNow");
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

if (defined $threads){
	my $check_threads_positive = &CheckPositiveInt($threads);
	my $threads_max;
	if(CheckFile("/proc/cpuinfo")){
		$threads_max = `grep 'processor' /proc/cpuinfo | sort -u | wc -l`;
		chomp $threads_max;
		$threads_max =~ s/\s//g;
	}else{
		my $mac_threads = `sysctl hw.logicalcpu`;
		chomp $mac_threads;
		$mac_threads =~ s/.*\://;
		$mac_threads =~ s/\s//g;
		if($mac_threads >= 2){
			$threads_max = $mac_threads;
		}else{
			$threads_max = 2;
		}
	}

	if ($check_threads_positive && $threads <= $threads_max){
		#threads provided by user is ok, doing nothing
	}else{
		InfoError("Threads number wrong!",'red');
		InfoError("Please provide a threads number between 0 - $threads_max that this server could support.");

		pod2usage(-verbose=>2,-exitval=>1);
		exit;
	}
}else{
	$threads = 1;#if -t not provided, default is NOT use theads;
}

my @fq1files;
if (defined $fq1){
	$format = 'fastq';
	my @fq1 = split ",",$fq1;
	for my $tmp (@fq1){
		$tmp =~ s/\s//g;
		if (not -e $tmp){
			InfoError("Input fastq1 file $tmp does NOT exist!",'red');
			exit;
		}else{
			my $pathtmp = abs_path($tmp);
			push @fq1files,$pathtmp;
		}
	}
}

my @fq2files;
if (defined $fq2){
	$pairEnded = 1;
	my @fq2 = split ",",$fq2;
	for my $tmp (@fq2){
		$tmp =~ s/\s//g;
		if (not -e $tmp){
			InfoError("Input fastq2 file $tmp does NOT exist!",'red');
			exit;
		}else{
			my $pathtmp = abs_path($tmp);
			push @fq2files,$pathtmp;
		}
	}
}


my @fastaFiles;
if (defined $fasta){
	$format = 'fasta';
	my @fasta = split ",",$fasta;
	for my $tmp (@fasta){
		$tmp =~ s/\s//g;
		if (not -e $tmp){
			InfoError("Input fasta file $tmp does NOT exist!",'red');
			exit;
		}else{
			my $pathtmp = abs_path($tmp);
			push @fastaFiles,$pathtmp;
		}
	}
}


if (defined $fq1 xor defined $fasta){
	#nothing
}else{
	InfoError("Either --fastq1/-1 or --fastq should be provided!","red");
	pod2usage(-verbose=>1,-exitval=>1);
	exit;
}

if (defined $barcodeLen){
	if($barcodeLen >= 0){
		if($barcodeLen <= 4){
			InfoWarn("You barcode length is only ${barcodeLen}bps? You better double check it.","yellow");
		}else{
			#nothing
		}
	}else{
		InfoError("Allowed mismatches can NOT be negative.","red");
		exit;
	}
}else{
	InfoError("Barcode length --barcodeLength/-l MUST be provided!",'red');
	pod2usage(-verbose=>1,-exitval=>1);
	exit;
}

if (defined $matchPosition){
	if (uc($matchPosition) eq 'BOS'){
		#nothing
	}elsif(uc($matchPosition) eq 'EOS'){
		#nothing
	}else{
		InfoError("The position of barcodes MUST be specified either \"bos\" or \"eos\". ");
		exit;
	}
}else{
	InfoError("The position of barcodes MUST be specified. Please choose \"bos\" or \"eos\". ");
	exit;
}

##core program starts here
##get the fastx_trim program
my $DEBUG_MODE = 1;
my $fastx_trim = File::Spec -> catfile($RealBin, '3rdPartyTools', 'fastx_toolkit', 'fastx_trimmer');
if(CheckProgram($fastx_trim, __FILE__, __LINE__, $DEBUG_MODE)){
	#keep running
}else{
	InfoError("The program $fastx_trim does NOT exist. Exiting...");
	exit;
}


if ($format eq 'fastq'){
	if($threads == 1){
		for my $f (@fq1files){
			my $tmpName = basename($f);
			my $outName = addTagForRawData($tmpName, 'barcodeTrimmed', 1);
			my $outFile = File::Spec -> catfile($outputDir, $outName);
			
			&trimmer($f, $outFile, $matchPosition, $barcodeLen, $format);
		}
		if($pairEnded){
			for my $f (@fq2files){
				my $tmpName = basename($f);
				my $outName = addTagForRawData($tmpName, 'barcodeTrimmed', 1);
				my $outFile = File::Spec -> catfile($outputDir, $outName);
			
				&trimmer($f, $outFile, $matchPosition, $barcodeLen, $format);
			}
		}
	}else{
		my @outFiles1;
		my @matchPos;
		my @barcodeLen;
		my @format;
		for my $f (@fq1files){
			my $tmpName = basename($f);
			my $outName = addTagForRawData($tmpName, 'barcodeTrimmed', 1);
			my $outFile = File::Spec -> catfile($outputDir, $outName);
			push @outFiles1, $outFile;
			push @matchPos, $matchPosition;
			push @barcodeLen, $barcodeLen;
			push @format, $format;
		}
		runMultipleThreadsWith5Args(\&trimmer, \@fq1files, \@outFiles1, \@matchPos,\@barcodeLen, \@format, $threads);
		
		if($pairEnded){
			my @outFiles2;
			for my $f (@fq2files){
				my $tmpName = basename($f);
				my $outName = addTagForRawData($tmpName, 'barcodeTrimmed', 1);
				my $outFile = File::Spec -> catfile($outputDir, $outName);
				push @outFiles2, $outFile;
			}
			runMultipleThreadsWith5Args(\&trimmer, \@fq2files, \@outFiles2, \@matchPos,\@barcodeLen, \@format, $threads);
		}
	}
}else{
	#$format eq 'fasta'
	if ($threads == 1){
		for my $f (@fastaFiles){
			my $tmpName = basename($f);
			my $outName = addTagForFasta($tmpName, 'barcodeTrimmed', 1);
			my $outFile = File::Spec -> catfile($outputDir, $outName);
			
			&trimmer($f, $outFile, $matchPosition, $barcodeLen, $format);
		}
	}else{
		my @outFiles;
		my @matchPos;
		my @barcodeLen;
		my @format;
		for my $f (@fastaFiles){
			my $tmpName = basename($f);
			my $outName = addTagForRawData($tmpName, 'barcodeTrimmed', 1);
			my $outFile = File::Spec -> catfile($outputDir, $outName);
			push @outFiles, $outFile;
			push @matchPos, $matchPosition;
			push @barcodeLen, $barcodeLen;
			push @format, $format;
		}
		runMultipleThreadsWith5Args(\&trimmer, \@fastaFiles, \@outFiles, \@matchPos,\@barcodeLen, \@format, $threads);
		
	}
}


##sub-program goes here
sub trimmer {
	my $inputFile = shift;
	my $outputFile = shift;
	my $pos = shift;
	my $len = shift;
	my $format = shift;

	
	my $inputFileName = basename($inputFile);
	my $outdir = dirname($outputFile);

	##info
	Info("Start to trim $inputFileName.");
	if(isGzipped($inputFile)){
		my $tmpName;
		if($format eq 'fastq'){
			$tmpName = addTagForRawData($inputFileName, 'ungzipped', 1) . "_" . time();
		}else{
			#$format eq 'fasta
			$tmpName = addTagForFasta($inputFileName, 'ungzipped', 1) . "_" . time();
		}

		my $tmpout = File::Spec -> catfile($outdir,$tmpName);

		#uncompress
		Info("Uncompressing gzipped file $inputFileName");
		my $cmd = "gunzip -c $inputFile > $tmpout";
		system($cmd);
		#runcmd($cmd);
		$inputFile = $tmpout;

		#start to trim
		if(uc($pos) eq 'BOS'){
			my $startpos = $len + 1;
			$cmd = "$fastx_trim -f $startpos -i $inputFile -o $outputFile";
			runcmd($cmd);
		}else{
			$cmd = "$fastx_trim -t $len -i $inputFile -o $outputFile";
			runcmd($cmd);
		}

		#remove tmp file
		$cmd = "rm -rf $inputFile";
		system($cmd);
	}else{
		#start to trim
		if(uc($pos) eq 'BOS'){
			my $startpos = $len + 1;
			my $cmd = "$fastx_trim -f $startpos -i $inputFile -o $outputFile";
			runcmd($cmd);
		}else{
			my $cmd = "$fastx_trim -t $len -i $inputFile -o $outputFile";
			runcmd($cmd);
		}
	}
	Info("Trim completed for $inputFileName.");

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
           



qap BarcodeTrimmer [options]

Use --help to see more information.

qap is still in development. If you have encounted any problem in usage, please feel no hesitation to cotact us.

=head1 DESCRIPTION

This script implements a function to trim barcode from fastq/fasta sequencing data. The script has B<several> mandatory options that MUST appear last. 

=head1 OPTIONS

=over 5

=item --fastq1,-1 F<FILE> [Required]

Path to next generation sequencing raw data. REQUIRED for both single-end or paired-end reads. Both compressed files and uncompressed files are allowed. 

Several files are allowed and should be seperated by comma, e.g. -1 Data1_R1.fq.gz,Data2_R1.fq.gz,Data3_R1.fq.gz . 

=item --fastq2,-2 F<FILE> [Optional]

Path to next generation sequencing raw data. REQUIRED for both paired-end reads. Both compressed files and uncompressed files are allowed. 

Several files are allowed and should be seperated by comma, e.g. -2 Data1_R2.fq.gz,Data2_R2.fq.gz,Data3_R2.fq.gz. The number of fastq2 files should equal the number of fastq1 files.

=item --fasta,-f F<File> [Optional]

Path of fasta files to be trimmed. One of --fastq1 or --fasta MUST be provided. Both compressed files and uncompressed files are allowed. 

Several files are allowed and should be seperated by comma, e.g. -f Data1.fasta.gz,Data2.fasta.gz,Data3.fasta.gz. 

=item --barcodeLength,-b F<INTEGER> [Required]

The length of barcode to be trimmed. MUST be provided.

=item --matchPosition,-p F<STRING> [Required]

This argument indicates where are the barcodes located. Must be one of "BOS"(Beginning Of Sequence) or "EOS"(End Of Sequence)[Case insensitive]. The program will try to 

trim barcodes at the begnning or end of sequences.

=item --outputDir,-o F<FILE> [Optional]

Path of the directory to storage result files. If NOT provided, the program will generate a folder automatically.

=item --threads,-t F<INTEGER> [Optional]

Number of threads this program will use when computing. A positive integer should be provided. The default value is 1.

=item --help,-h

Display this detailed help information.

=back

=head1 EXAMPLE

=over 5

qap BarcodeTrimmer -1 Data1_R1.fq.gz,Data2_R1.fq.gz,Data3_R1.fq.gz -2 Data1_R2.fq.gz,Data2_R2.fq.gz,Data3_R2.fq.gz -l 8 -p bos -t 2 -o ./timmedSep

=back

=head1 AUTHOR

Mingjie Dr.Wang I<huzai@sjtu.edu.cn>

=head1 COPYRIGHT

Copyright (C) 2017, Mingjie Wang. All rights reserved.




