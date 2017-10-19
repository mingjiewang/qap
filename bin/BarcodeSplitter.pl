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
printcol ("BarcodeSplitter","green");
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
my $mismatch;
my $partial;
my $pairEnded = 0;
my $barcodeFile;
my $matchPosition;

my $DateNow = `date +"%Y%m%d_%Hh%Mm%Ss"`;
chomp $DateNow;

GetOptions(
'b|barcodeFile|=s' => \$barcodeFile,
'1|fastq1|=s'      => \$fq1,
'2|fastq2|=s'      => \$fq2,
'm|mismatch|=s'     => \$mismatch,
'r|partial|=s'      => \$partial,
'o|outputDir|=s'    => \$outputDir,
'p|matchPosition|=s'  => \$matchPosition,
'h|help|'          => \$help
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
	$outputDir = File::Spec -> catfile($wk_dir,"qap_Results_for_BarcodeSplitter_$DateNow");
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


if (defined $barcodeFile){
	if (not -e $barcodeFile){
		InfoError("Input barcode file $barcodeFile does NOT exist!",'red');
		pod2usage(-verbose=>0,-exitval=>1);
		exit;
	}
	$barcodeFile = abs_path($barcodeFile);
}else{
	InfoError("Input barcode file must be provided!",'red');
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

if (defined $mismatch){
	if($mismatch >= 0){
		if($mismatch >= 5){
			InfoWarn("You are allowing $mismatch mismatches when matching barcodes. It is too high, maybe a lower value would be better.","yellow");
		}else{
			#nothing
		}
	}else{
		InfoError("Allowed mismatches can NOT be negative.","red");
		exit;
	}
}else{
	$mismatch = 0;
}

if (defined $partial){
	if($partial >= 0){
		if($partial >= 5){
			InfoWarn("You are allowing $partial partial overlaps when matching barcodes. It is too high, maybe a lower value would be better.","yellow");
		}else{
			#nothing
		}
	}else{
		InfoError("Allowed overlaps can NOT be negative.","red");
		exit;
	}
}else{
	$partial = 0;
}

my $matchPositionExpression;
if (defined $matchPosition){
	if (uc($matchPosition) eq 'BOS'){
		$matchPositionExpression = "--bol";
	}elsif(uc($matchPosition) eq 'EOS'){
		$matchPositionExpression = "--eol";
	}else{
		InfoError("The position of barcodes MUST be specified either \"bos\" and \"eos\". ");
		exit;
	}
}else{
	InfoError("The position of barcodes MUST be specified. Please choose \"bos\" or \"eos\". ");
	exit;
}

##start to run fastx toolkit
my $fastxBarcodeSplitter = File::Spec -> catfile($RealBin,"3rdPartyTools","fastx_toolkit","fastx_barcode_splitter.pl");
if (isGzipped($fq1)){
	my $cmd = "zcat $fq1 | perl $fastxBarcodeSplitter --bcfile $barcodeFile $matchPositionExpression --mismatches $mismatch --partial $partial --prefix $outputDir --suffix \"_R1.fq\" ";
	runcmd($cmd);
}else{
	my $cmd = "cat $fq1 | perl $fastxBarcodeSplitter --bcfile $barcodeFile $matchPositionExpression --mismatches $mismatch --partial $partial --prefix $outputDir --suffix \"_R1.fq\" ";
	runcmd($cmd);
}

if ($pairEnded){
	if (isGzipped($fq2)){
		my $cmd = "zcat $fq2 | perl $fastxBarcodeSplitter --bcfile $barcodeFile $matchPositionExpression --mismatches $mismatch --partial $partial --prefix $outputDir --suffix \"_R2.fq\" ";
		runcmd($cmd);
	}else{
		my $cmd = "cat $fq2 | perl $fastxBarcodeSplitter --bcfile $barcodeFile $matchPositionExpression --mismatches $mismatch --partial $partial --prefix $outputDir --suffix \"_R2.fq\" ";
		runcmd($cmd);
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
           



qap BarcodeSplitter [options]

Use --help to see more information.

qap is still in development. If you have encounted any problem in usage, please feel no hesitation to cotact us.

=head1 DESCRIPTION

This script implements a function for splitting sequencing data into several files according to sample barcoding. The script has B<several> mandatory options that MUST appear last. 

=head1 OPTIONS

=over 5

=item --barcode,-b F<FILE> [Required]

Path to the file which contains all the barcodes of each sample. This file MUST be formatted accordingly: 1. NO headline; 2. Two columns all allowed. First column is the name of each sample

and the second column contains the barcode sequence for each sample. Columns are seprated by spaces.

=item --fastq1,-1 F<FILE> [Required]

Path to next generation sequencing raw data. REQUIRED for both single-end or paired-end reads. Both compressed files and uncompressed files are allowed.

=item --fastq2,-2 F<FILE> [Optional]

Path to next generation sequencing raw data. REQUIRED for both paired-end reads. Both compressed files and uncompressed files are allowed.

=item --mismatch,-m F<INTEGER> [Optional]

Max number of mismatches allowed. Default is 0.

=item --parital,-r F<INTEGER> [Optional]

Allow parital overlap of barcodes. Default is 0, which means no partial overlap allowed.

=item --matchPosition,-p F<STRING> [Required]

This argument indicates where are the barcodes located. Must be one of "BOS"(Beginning Of Sequence) or "EOS"(End Of Sequence)[Case insensitive]. The program will try to 

match barcodes at the begnning or end of sequences.

=item --outputDir,-o F<FILE> [Optional]

Path of the directory to storage the splitted reads data. If NOT provided, the program will generate a folder automatically.

=item --help,-h

Display this detailed help information.

=back

=head1 EXAMPLE

=over 5

qap BarcodeSplitter -1 ngs_R1.fq.gz -2 ngs_R2.fq.gz -m 0 -p 0 -b barcode.txt --matchPosition bos -o ./split

=back

=head1 AUTHOR

Mingjie Dr.Wang I<huzai@sjtu.edu.cn>

=head1 COPYRIGHT

Copyright (C) 2017, Mingjie Wang. All rights reserved.


