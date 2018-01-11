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
printcol ("ExtractReadID","green");
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
my $outputDir;
my $pairEnded = 0;
my $sampleLabel;
my $threads;

my $DateNow = `date +"%Y%m%d_%Hh%Mm%Ss"`;
chomp $DateNow;

GetOptions(
'l|sampleLabel=s'  => \$sampleLabel,
'1|fastq1|=s'      => \$fq1,
'2|fastq2|=s'      => \$fq2,
'o|outputDir=s'    => \$outputDir,
'h|help|'          => \$help,
't|threads|=s'     => \$threads
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
	$outputDir = File::Spec -> catfile($wk_dir,"qap_Results_for_ExtractReadID_$DateNow");
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
	if(-e ("/proc/cpuinfo")){
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
}else{
	InfoError("Input fastq1 file must be provided!",'red');
	pod2usage(-verbose=>1,-exitval=>1);
	exit;
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

my @sampleLabel;
my $numberOfSamples;
if (defined $sampleLabel){
	@sampleLabel = split ",",$sampleLabel;
	my $num1 = scalar(@sampleLabel);
	my $num2 = scalar(@fq1files);
	my $num3 = scalar(@fq2files);
	my $equal = check3EqualNums($num1, $num2, $num3);
	
	if ($equal){
		#nothing
		$numberOfSamples = $num1;
	}else{
		InfoError("There are $num1 sample labels, $num2 fastq1 files and $num3 fastq2 files provided! Please check again.", "red");
		pod2usage(-verbose=>1,-exitval=>1);
		exit;
	}
}else{
	InfoError("Label of each sample must be provided!",'red');
	pod2usage(-verbose=>1,-exitval=>1);
	exit;
}


##start to extract
if ($threads == 1){
	for my $i (1..$numberOfSamples){
		my $fq1file = $fq1files[$i - 1];
		my $fq2file = $fq2files[$i - 1] if $pairEnded;
		my $sampleLabel = $sampleLabel[$i - 1];
		my $outputFileName = $sampleLabel . '.readID.txt';
		my $outputFile = File::Spec -> catfile($outputDir,$outputFileName);
		
		&getID($fq1file, $fq2file, $outputFile);
	}

}else{
	my @outputFiles;
	for my $i (1..$numberOfSamples){
		my $sampleLabel = $sampleLabel[$i - 1];
		my $outputFileName = $sampleLabel . '.readID.txt';
		my $outputFile = File::Spec -> catfile($outputDir,$outputFileName);
		push @outputFiles, $outputFile;
	}
	
	#run with multiple threads
	runMultipleThreadsWith3Args(\&getID, \@fq1files, \@fq2files, \@outputFiles, $threads);
	
}



##subprograms goes here
sub getID {
	my ($fq1,$fq2,$out) = @_;
	
	my $fq1Name = basename($fq1);
	my $fq2Name = basename($fq2);
	
	Info("Extracting IDs from $fq1Name and $fq2Name.");
	
	##handle gzipped files
	my $outdir = dirname($out);
	if (isGzipped($fq1Name)){
		my $tmpName = addTagForRawData($fq1Name, 'ungzipped', 1) . "_" . time();
		my $tmpout = File::Spec -> catfile($outdir,$tmpName);
		
		#uncompress
		Info("Uncompressing gzipped file $fq1Name");
		my $cmd = "gunzip -c $fq1 > $tmpout";
		system($cmd);
		#runcmd($cmd);
		$fq1 = $tmpout;
	}
	if (isGzipped($fq2Name)){
		my $tmpName = addTagForRawData($fq2Name, 'ungzipped', 1) . "_" . time();
		my $tmpout = File::Spec -> catfile($outdir,$tmpName);
		
		#uncompress
		Info("Uncompressing gzipped file $fq2Name");
		my $cmd = "gunzip -c $fq2 > $tmpout";
		system($cmd);
		#runcmd($cmd);
		$fq2 = $tmpout;
	}
	
	my @id;	
	open T1,$fq1 or die "Cannot open file $fq1:$!\n";
	while (my $line1 = <T1>){
		chomp $line1;
		chomp (my $line2 = <T1>);
		chomp (my $line3 = <T1>);
		chomp (my $line4 = <T1>);
		
		if($line1 =~ /(.*?)\s+/){
			push @id,$1;
		}elsif($line1 =~ /(.*?)[\/\\]1/){
			push @id,$1;
		}else{
			InfoError("Can NOT recognize read ID in $fq1Name.");
			exit;
		}
		
	}
	close T1;
	
	open T2,$fq2 or die "Cannot open file $fq2:$!\n";
	while (my $line1 = <T2>){
		chomp $line1;
		chomp (my $line2 = <T2>);
		chomp (my $line3 = <T2>);
		chomp (my $line4 = <T2>);
		
		if($line1 =~ /(.*?)\s+/){
			push @id,$1;
		}elsif($line1 =~ /(.*?)[\/\\]2/){
			push @id,$1;
		}else{
			InfoError("Can NOT recognize read ID in $fq2Name.");
			exit;
		}
	}
	close T2;
	
	my $idRmDup = DelDup(\@id);
	
	open OUT,">$out" or die "Cannot output to file $out:$!\n";
	map {print OUT "$_\n";} @$idRmDup;
	close OUT;
	
	my $idNum = scalar @$idRmDup;
	my $totalNum = scalar @id;
	Info("$idNum IDs from $totalNum reads are extracted and output to $out.");
	
	##rm tmp uncompressed files
	if (isGzipped($fq1Name)){
		my $cmd = "rm -rf $fq1";
		system($cmd);
	}
	if (isGzipped($fq2Name)){
		my $cmd = "rm -rf $fq2";
		system($cmd);
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
           



qap ExtractReadID [options]

Use --help to see more information.

qap is still in development. If you have encounted any problem in usage, please feel no hesitation to cotact us.

=head1 DESCRIPTION

This script implements a function for extracting read ID from fastq sequencing data. The script has B<several> mandatory options that MUST appear last. 

=head1 OPTIONS

=over 5

=item --fastq1,-1 F<FILE> [Required]

Path to next generation sequencing raw data. REQUIRED for both single-end or paired-end reads. Both compressed files and uncompressed files are allowed. 

Several files are allowed and should be seperated by comma, e.g. -1 Data1_R1.fq.gz,Data2_R1.fq.gz,Data3_R1.fq.gz . 

=item --fastq2,-2 F<FILE> [Optional]

Path to next generation sequencing raw data. REQUIRED for both paired-end reads. Both compressed files and uncompressed files are allowed. 

Several files are allowed and should be seperated by comma, e.g. -2 Data1_R2.fq.gz,Data2_R2.fq.gz,Data3_R2.fq.gz. The number of fastq2 files should equal the number of fastq1 files.

=item --sampleLabel, -l F<STRING> [Required]

The name of all samples. Multiple samples should be seperated by comma.

=item --outputDir,-o F<FILE> [Optional]

Path of the directory to storage read ID files. If NOT provided, the program will generate a folder automatically.

=item --threads, -t F<INTEGER> [Optional]

Number of threads this program will use when computing. A positive integer should be provided. The default value is 1.

=item --help,-h

Display this detailed help information.

=back

=head1 EXAMPLE

=over 5

qap ExtractReadID -1 Data1_R1.fq.gz,Data2_R1.fq.gz,Data3_R1.fq.gz -2 Data1_R2.fq.gz,Data2_R2.fq.gz,Data3_R2.fq.gz -l Data1,Data2,Data3  -o ./readID

=back

=head1 AUTHOR

Mingjie Dr.Wang I<huzai@sjtu.edu.cn>

=head1 COPYRIGHT

Copyright (C) 2017, Mingjie Wang. All rights reserved.
