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
printcol ("ExtractSeqInR","green");
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
my $sampleLabel;
my $idFile;
my $format;
my $threads;

my $DateNow = `date +"%Y%m%d_%Hh%Mm%Ss"`;
chomp $DateNow;

GetOptions(
'1|fastq1|=s'       => \$fq1,
'2|fastq2|=s'       => \$fq2,
'o|outputDir|=s'    => \$outputDir,
'i|idFile|=s'       => \$idFile,
'f|outFormat|=s'    => \$format,
'l|sampleLabel|=s'  => \$sampleLabel,
'h|help|'           => \$help,
't|threads|=s'      => \$threads
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
	$outputDir = File::Spec -> catfile($wk_dir,"qap_Results_for_ExtractSeqInR_$DateNow");
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
#find the r script
my $Rscript = File::Spec -> catfile($RealBin, 'Rscripts', 'ExtractSeq.R');
if (-e $Rscript){
	#nothing
}else{
	InfoError("Rscript $Rscript is missing. Abortting...",'red');
	exit;
}

#handle gzipped files
if (isGzipped($fq1)){
	my $fqName = basename($fq1);
	my $tmpName = addTagForRawData($fqName, 'ungzipped', 1) . "_" . time();
	my $tmpout = File::Spec -> catfile($outputDir,$tmpName);
	
	#uncompress
	Info("Uncompressing gzipped file $fqName");
	my $cmd = "gunzip -c $fq1 > $tmpout";
	system($cmd);
	#runcmd($cmd);
	$fq1 = $tmpout;
	
	#handle the output file name and file path
	my $outName = addTagForRawData($fqName, 'extractedSeq', 1);
	$outName = $sampleLabel . "_" .$outName;
	$outName =~ s/\.gz$//; #remove .gz suffix is input fastq file is gzipped.
	
	#core program
	my $tmpRInputFile = File::Spec -> catfile($outputDir, $tmpName . ".RInput");
	if ($format eq 'fastq'){
		my $outFile = File::Spec -> catfile($outputDir, $outName);
		
		#run extraction
		formatFqIntoColumns($fq1, $tmpRInputFile, 'fastq');
		$cmd = "Rscript $Rscript -i $tmpRInputFile -d $idFile -o $outFile";
		runcmd($cmd);
		Info("Extraction completed for $fqName");
		
		#remove tmp files
		$cmd = "rm -rf $tmpRInputFile $fq1";
		system($cmd);
	}else{
		my $outFile = File::Spec -> catfile($outputDir, changeFastqSuffix2Fasta($outName));
		
		#extraction
		formatFqIntoColumns($fq1, $tmpRInputFile, 'fasta');
		my $cmd = "Rscript $Rscript -i $tmpRInputFile -d $idFile -o $outFile";
		runcmd($cmd);
		Info("Extraction completed for $fqName");
		
		#remove tmp files
		$cmd = "rm -rf $tmpRInputFile $fq1";
		system($cmd);
	}
	
}else{
	my $fqName = basename($fq1);
	my $outName = addTagForRawData($fqName, 'extractedSeq', 1);
	$outName = $sampleLabel . "_" .$outName;
	
	my $tmpRInputFile = File::Spec -> catfile($outputDir, $fqName . ".RInput");
	if ($format eq 'fastq'){
		my $outFile = File::Spec -> catfile($outputDir, $outName);
		
		#extraction
		formatFqIntoColumns($fq1, $tmpRInputFile, 'fastq');
		my $cmd = "Rscript $Rscript -i $tmpRInputFile -d $idFile -o $outFile";
		runcmd($cmd);
		Info("Extraction completed for $fqName");
		
		#remove tmp files
		$cmd = "rm -rf $tmpRInputFile";
		system($cmd);
	}else{
		my $outFile = File::Spec -> catfile($outputDir, changeFastqSuffix2Fasta($outName));
		
		#extraction
		formatFqIntoColumns($fq1, $tmpRInputFile, 'fasta');
		my $cmd = "Rscript $Rscript -i $tmpRInputFile -d $idFile -o $outFile";
		runcmd($cmd);
		Info("Extraction completed for $fqName");
		
		#remove tmp files
		$cmd = "rm -rf $tmpRInputFile";
		system($cmd);
	}
}


if($pairEnded){
	#handle gzipped files
	if (isGzipped($fq2)){
		my $fqName = basename($fq2);
		my $tmpName = addTagForRawData($fqName, 'ungzipped', 1) . "_" . time();
		my $tmpout = File::Spec -> catfile($outputDir,$tmpName);

		#uncompress
		Info("Uncompressing gzipped file $fqName");
		my $cmd = "gunzip -c $fq2 > $tmpout";
		system($cmd);
		#runcmd($cmd);
		$fq2 = $tmpout;

		#handle the output file name and file path
		my $outName = addTagForRawData($fqName, 'extractedSeq', 1);
		$outName =~ s/\.gz$//; #remove .gz suffix is input fastq file is gzipped.
		$outName = $sampleLabel . "_" .$outName;
		
		#core program
		my $tmpRInputFile = File::Spec -> catfile($outputDir, $tmpName . ".RInput");
		if ($format eq 'fastq'){
			my $outFile = File::Spec -> catfile($outputDir, $outName);

			#run extraction
			formatFqIntoColumns($fq2, $tmpRInputFile, 'fastq');
			$cmd = "Rscript $Rscript -i $tmpRInputFile -d $idFile -o $outFile";
			runcmd($cmd);
			Info("Extraction completed for $fqName");
			
			#remove tmp files
			$cmd = "rm -rf $tmpRInputFile $fq2";
			system($cmd);
		}else{
			my $outFile = File::Spec -> catfile($outputDir, changeFastqSuffix2Fasta($outName));
			
			#run extraction
			formatFqIntoColumns($fq2, $tmpRInputFile, 'fasta');
			$cmd = "Rscript $Rscript -i $tmpRInputFile -d $idFile -o $outFile";
			runcmd($cmd);
			Info("Extraction completed for $fqName");
			
			#remove tmp files
			$cmd = "rm -rf $tmpRInputFile $fq2";
			system($cmd);
		}

	}else{
		my $fqName = basename($fq2);
		my $outName = addTagForRawData($fqName, 'extractedSeq', 1);
		$outName = $sampleLabel . "_" .$outName;
		
		my $tmpRInputFile = File::Spec -> catfile($outputDir, $fqName . ".RInput");
		if ($format eq 'fastq'){
			my $outFile = File::Spec -> catfile($outputDir, $outName);
			
			#run extraction
			formatFqIntoColumns($fq2, $tmpRInputFile, 'fastq');
			my $cmd = "Rscript $Rscript -i $tmpRInputFile -d $idFile -o $outFile";
			runcmd($cmd);
			Info("Extraction completed for $fqName");
			
			#remove tmp files
			$cmd = "rm -rf $tmpRInputFile";
			system($cmd);
		}else{
			my $outFile = File::Spec -> catfile($outputDir, changeFastqSuffix2Fasta($outName));
			
			#run extraction
			formatFqIntoColumns($fq2, $tmpRInputFile, 'fasta');
			my $cmd = "Rscript $Rscript -i $tmpRInputFile -d $idFile -o $outFile";
			runcmd($cmd);
			Info("Extraction completed for $fqName");
			
			#remove tmp files
			$cmd = "rm -rf $tmpRInputFile";
			system($cmd);
		}
	}

}else{
	#nothing
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
           



qap ExtractSeqInR [options]

Use --help to see more information.

qap is still in development. If you have encounted any problem in usage, please feel no hesitation to cotact us.

=head1 DESCRIPTION

This script implements a function for extracting sequences from fastq files with read IDs in a fast way by using R. The script has B<several> mandatory options that MUST appear last. 

=head1 OPTIONS

=over 5

=item --fastq1,-1 F<FILE> [Required]

Path to next generation sequencing raw data. REQUIRED for both single-end or paired-end reads. Both compressed files and uncompressed files are allowed. 

=item --fastq2,-2 F<FILE> [Optional]

Path to next generation sequencing raw data. REQUIRED for both paired-end reads. Both compressed files and uncompressed files are allowed. 

=item --idFile,-i F<FILE> [Required]

Path to the file containing the read IDs by using which the sequences would be extracted. Each ID should be listed per line.

=item --sampleLabel,-s F<STRING> [Required]

The name of the sample. This argument MUST be provided. The output file will be named using sample label.

=item --outFormat, -f F<STRING> [Optional]

The format of program output. The value should be one of 'fastq' or 'fasta'.

=item --outputDir,-o F<FILE> [Optional]

Path of the directory to storage the result data. If NOT provided, the directory of raw data files would be used instead.

=item --threads,-t F<INTEGER> [Optional]

Threads used to run this program. A positive integer is required. Default value is 1.

=item --help,-h

Display this detailed help information.

=back

=head1 EXAMPLE

=over 5

qap ExtractSeqInR -1 Data1_R1.fq.gz -2 Data1_R2.fq.gz -i readID.txt -f fastq -o ./seq

=back

=head1 AUTHOR

Mingjie Dr.Wang I<huzai@sjtu.edu.cn>

=head1 COPYRIGHT

Copyright (C) 2017, Mingjie Wang. All rights reserved.


