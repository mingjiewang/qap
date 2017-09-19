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
##  License along with QAP; if not, see                             ##
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
printcol ("ShannonEntropy","green");
print "\n";

##get workding directory
my $wk_dir = getcwd;
my $mainBin;
if ($RealBin =~ /(.*)\/bin/){
	$mainBin = $1;
}

####define command line arguments
my $help;
my $inputDir;
my $outputDir;
my $format;
my $suffix;
my $graphs;

my $DateNow = `date +"%Y%m%d_%Hh%Mm%Ss"`;
chomp $DateNow;

GetOptions(
'i|inputDir|=s'     => \$inputDir,
's|suffix|=s'       => \$suffix,
'f|format|=s'       => \$format,     
'o|outputDir|=s'    => \$outputDir,
'g|graphic|=s'      => \$graphs,
'h|help|'           => \$help
);


##check command line arguments
if (defined $help){
	pod2usage(-verbose=>2,-exitval=>1);
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
	$outputDir = File::Spec -> catfile($wk_dir,"qap_Results_for_ShannonEntropy_$DateNow");
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

if(defined $inputDir){
	$inputDir = abs_path($inputDir) . "/";
	if (not -e $inputDir){
		InfoError("Input directory $inputDir does NOT exist! Please check again.");
		exit;
	}
}else{
	InfoError("Input directory MUST be specified with -i/--inputDir\n");
	pod2usage(-verbose=>0,-exitval=>1);
	exit;
}

my $drawgraph = 0;
if(defined $graphs){
	if($graphs =~ /t/i){
		$drawgraph = 1;
	}elsif($graphs =~ /f/i){
		$drawgraph = 0
	}else{
		InfoError("Please specify whether draw graphs or not using -g/--graphs with \'T(true)\' or \'F(false)\'.");
	}
}else{
	$drawgraph = 0;
	InfoWarn("-g/--graphs not specified, using false as default.");
}


my $numberOfFiles = 0;
my @inputfiles;
if(defined $suffix){
	@inputfiles = glob ("${inputDir}/*.${suffix}");
	
	$numberOfFiles = scalar(@inputfiles);
	
	if ($numberOfFiles == 0){
		InfoError("There are NOT any files in $inputDir with suffix \'.${suffix}\'. Please check again.");
		exit;
	}
	
	Info("Find $numberOfFiles files.");
	my $i = 1;
	for my $f (@inputfiles){
		printf "[%02d] $f\n",$i;
		$i++;
	}
	
}else{
	InfoWarn("The suffix is not provided. The program will try to read in every file in $inputDir");
	
	@inputfiles = glob ("${inputDir}/*.*");
	
	$numberOfFiles = scalar(@inputfiles);
	
	if ($numberOfFiles == 0){
		InfoError("There are NOT any files in $inputDir with suffix \'.${suffix}\'. Please check again.");
		exit;
	}
	
	Info("Find $numberOfFiles files.");
	my $i = 1;
	for my $f (@inputfiles){
		printf "[%02d] $f\n",$i;
		$i++;
	}
	
}

if (defined $format){
	if(uc($format) eq 'FASTA' or uc($format) eq 'FASTQ'){
		#nothing
	}else{
		InfoError("The file format MUST be one of \'fasta\' or \'fastq\'.");
		pod2usage(-verbose=>0,-exitval=>1);
		exit;
	}
}else{
	InfoError("The file format MUST be specified using -f/--format.");
	pod2usage(-verbose=>0,-exitval=>1);
	exit;
}


##the core program starts here
Info("Start calculating...");
sleep(1);

#check r script 
my $rscript = File::Spec -> catfile($mainBin, 'bin', 'Rscripts', 'CalculateShannonEntropy.R');
if (-e $rscript){
	#nothing
}else{
	InfoError("Rscript $rscript is missing. Abortting...",'red');
	exit;
}

#output file
my $inputDirName = basename($inputDir);
my $resfileName = "ShannonEntropy" . "_" . $inputDirName . ".txt";
my $resfile = File::Spec -> catfile($outputDir, $resfileName);
open RES,">>$resfile" or die "Cannot output to $resfile:$!";
print RES "Sample\tValue\n";

#run sub program
my $i = 1;
for my $f (@inputfiles){
	&shannon($f, $format, $rscript, $resfile, $drawgraph);
	
	#process bar
	InfoProcessBar($i, $numberOfFiles);
	$i++;
}





sub shannon {
	my $file = shift;
	my $format = shift;
	my $rscript = shift;
	my $resfile = shift;
	my $drawgraph = shift;
	
	my $fileName = basename($file);
	
	my $rInputfile = $file . "." . time() . ".RInput";
	
	my $rRunLogFile = File::Spec -> catfile($outputDir, "${fileName}.RLog");
	
	if (uc($format) eq 'FASTA'){
		$fileName = removeFastaSuffix($fileName);
		
		#extract seq
		extractSeqFromFasta($file, $rInputfile);
		
		#calculate
		my $cmd = "Rscript $rscript -i $rInputfile -o $resfile -l $fileName -p $drawgraph > $rRunLogFile";
		system($cmd);
	}elsif(uc($format) eq 'FASTQ'){
		$fileName = removeFastqSuffix($fileName);
		
		#extract seq
		extractSeqFromFastq($file, $rInputfile);
		
		#calculate
		my $cmd = "Rscript $rscript -i $rInputfile -o $resfile -l $fileName -p $drawgraph > $rRunLogFile";
		system($cmd);
	}else{
		InfoError("The formar MUST be \'fasta\' or \'fastq\'.");
		exit;
	}
	
	#remove tmp data file
	my $inputFolder = dirname($file);
	my $cmd = "rm -rf $inputFolder/${fileName}*.RInput";
	system($cmd);
	
	$cmd = "rm -rf $rRunLogFile";
	system($cmd);
	
	#handle figure output
	if($wk_dir eq $outputDir){
		#nothing
	}else{
		#mv tif and pdf to outputDir
		my $png = File::Spec -> catfile($wk_dir, "${fileName}.png");
		$cmd = "mv $png $outputDir";
		system($cmd);
		
		my $pdf = File::Spec -> catfile($wk_dir, "${fileName}.pdf");
		$cmd = "mv $pdf $outputDir";
		system($cmd);
	}
	
}


##run success
print("\n");
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
           



qap ShannonEntropy [options]

Use --help to see more information.

qap is still in development. If you have encounted any problem in usage, please feel no hesitation to cotact us.

=head1 DESCRIPTION

This script implements a function for calculating Shannon entroy in batch. The script has B<several> mandatory options that MUST appear last. 

=head1 OPTIONS

=over 5

=item --inputDir,-i F<FILE> [Required]

Path to directory contaning all the files to be read in.

=item --suffix,-s F<STRING> [Optional]

Suffix of the files to be read in. If suffix is not provided, all the files in input directory will be read in.

=item --format,-f F<STRING> [Required]

The format of the files to be calculated. Should be one of 'fasta' or 'fastq'.

=item --graphs,-g F<BOOLEAN> [Optional]

Whether draw graphs or not for illustrating virus quansispecies population structure. T for true and F for false. Using F as default.

=item --outputDir,-o F<FILE> [Optional]

Path of the directory to storage result files. If NOT provided, the program will generate a folder automatically.

=item --help,-h

Display this detailed help information.

=back

=head1 EXAMPLE

=over 5

qap ShannonEntropy -i ./seq -f fasta -s fas -g T -o ./shannon

=back

=head1 AUTHOR

Mingjie Dr.Wang I<huzai@sjtu.edu.cn>

=head1 COPYRIGHT

Copyright (C) 2017, Mingjie Wang. All rights reserved.


