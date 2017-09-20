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
use Mapper;
use QSR;

####---------------------------####
####The program begins here
####---------------------------####

##Show welcome
print "You are now running subprogram: ";
printcol ("ECnQSR","green");
print "\n";


## check threads available or not
$| = 1;
InfoPlain("Checking threading status");
sleep(1);
my $threads_usable = eval 'use threads; 1';
if ($threads_usable) {
	use threads;
	use threads::shared;
	InfoPlain("Perl threading enabled");
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
my $threads;
my $ref;
my $program;

my $DateNow = `date +"%Y%m%d_%Hh%Mm%Ss"`;
chomp $DateNow;

GetOptions(
'1|fastq1|=s'        => \$fq1,
'2|fastq2|=s'        => \$fq2,
'r|reference|=s'     => \$ref,
'o|outputDir|=s'     => \$outputDir,
'p|program|=s'       => \$program,
'h|help|'            => \$help,
't|threads|=s'       => \$threads
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
		my $cmd = "mkdir -p $outputDir";
		system($cmd);
	}else{
		InfoError("Folder $outputDir already exists. Please specify another output directory using option -o/--outputDir");
		exit;
	}
}else{
	$outputDir = File::Spec -> catfile($wk_dir,"qap_Results_for_ECnQSR_$DateNow");
	InfoWarn("The output directory is not provided!",'yellow');
	InfoWarn("Will mkdir \"$outputDir\" and use it as the output directory.",'yellow');
	
	if (!-e "$outputDir"){
		my $cmd = "mkdir -p $outputDir";
		system($cmd);
	}else{
		InfoError("Mkdir Failed! $outputDir already exists!","red");
		InfoError("Please specify another output directory using option -o/--outputDir\n");
		pod2usage(-verbose=>2,-exitval=>1);
		exit;
	}

}

if (defined $threads){
	my $check_threads_positive = &CheckPositiveInt($threads);
	my $threads_max;
	if(existFile("/proc/cpuinfo")){
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

if (defined $fq1){
	if(not existFile($fq1)){
		InfoError("Fastq1 file $fq1 does NOT exits.");
		exit(0);
	}
}else{
	InfoError("Fastq1 file MUST be provided using --fastq1/-1.");
	pod2usage(-verbose=>1,-exitval=>1);
	exit;
}

if (defined $fq2){
	if(not existFile($fq2)){
		InfoError("Fastq2 file $fq2 does NOT exits.");
		exit(0);
	}
}else{
	InfoError("Fastq2 file MUST be provided using --fastq2/-2.");
	pod2usage(-verbose=>1,-exitval=>1);
	exit;
}

if(defined $ref){
	if(existFile($ref)){
		#nothing
	}else{
		InfoError("$ref doesn't exist. Please check.");
		exit(0);
	}
}else{
	InfoError("Input reference file MUST be provided.");
	pod2usage(-verbose=>2,-exitval=>1);
	exit;
}

if (defined $threads){
	my $check_threads_positive = &CheckPositiveInt($threads);
	my $threads_max = `grep 'processor' /proc/cpuinfo | sort -u | wc -l`;
	chomp $threads_max;

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

my @program;
if (defined $program){
	my @tmp = split ',',$program;
	for my $t (@tmp){
		if(uc($t) eq 'SHORAH' || uc($t) eq 'QURE' || uc($t) eq 'PREDICTHAPLO' || uc($t) eq 'VIQUAS' ){
			push @program,$t;
		}else{
			InfoError("The mapping program MUST be among \"Shorah\", \"QuRe\", \"PredictHaplo\" or \"ViQuaS\" .","red");
			pod2usage(-verbose=>2,-exitval=>1);
			exit;
		}
	}
	
}else{
	InfoWarn("The mapping program --program/-p is not provided. Using \"ViQuaS\" as default.");
	push @program,'viquas';
}


##core program starts here
my $DEBUG_MODE = 1;
#check bwa
my $bwa_excu = File::Spec -> catfile($RealBin,'3rdPartyTools','bwa','bwa');
if(CheckProgram($bwa_excu, __FILE__, __LINE__, $DEBUG_MODE)){
	#keep running
}else{
	InfoError("The program $bwa_excu does NOT exist. Exiting...");
	exit;
}
#check fq2fa
my $fq2fa_excu = File::Spec -> catfile($RealBin,'3rdPartyTools','fastx_toolkit','fq2fa');
if(CheckProgram($fq2fa_excu, __FILE__, __LINE__, $DEBUG_MODE)){
	#keep running
}else{
	InfoError("The program $fq2fa_excu does NOT exist. Exiting...");
	exit;
}
#check samtools
my $samtools_excu = File::Spec -> catfile($RealBin,'3rdPartyTools','samtools','samtools');
if(CheckProgram($samtools_excu, __FILE__, __LINE__, $DEBUG_MODE)){
	#keep running
}else{
	InfoError("The program $samtools_excu does NOT exist. Exiting...");
	exit;
}
#check Shorah
my $shorah_excu = File::Spec -> catfile($RealBin,"3rdPartyTools","qsr","Shorah","shorah.py");
if(not existFile($shorah_excu)){
	InfoError("$shorah_excu does NOT exist. Please check.");
	exit(0);
}
#check PredictHaplo
my $predicthaplo_excu = File::Spec -> catfile($RealBin,"3rdPartyTools","qsr","PredictHaplo","PredictHaplo-Paired");
if(not existFile($predicthaplo_excu)){
	InfoError("$predicthaplo_excu does NOT exist. Please check.");
	exit(0);
}
#check QuRe
my $qure_excu_dir = File::Spec -> catfile($RealBin,"3rdPartyTools","qsr","QuRe");
if(not CheckFolder($qure_excu_dir)){
	InfoError("$qure_excu_dir does NOT exist. Please check.");
	exit(0);
}
#check ViQuaS
my $viquas_excu_dir = File::Spec -> catfile($RealBin,"3rdPartyTools","qsr","ViQuaS");
if(not CheckFolder($viquas_excu_dir)){
	InfoError("$viquas_excu_dir does NOT exist. Please check.");
	exit(0);
}

##start to run
for my $p (@program){
	qsr_pipeline($fq1,$fq2,$ref,$outputDir,$p,$threads,$fq2fa_excu,$samtools_excu,$bwa_excu,$shorah_excu,$predicthaplo_excu,$qure_excu_dir,$viquas_excu_dir);
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
           



qap ECnQSR [options]

Use --help to see more information.

qap is still in development. If you have encounted any problem in usage, please feel no hesitation to cotact us.

=head1 DESCRIPTION

This script implements a function for sequencing error correction (EC) and quasispecies reconstruction (QSR). The script has B<several> mandatory options that MUST appear last. 

=head1 OPTIONS

=over 5

=item --fastq1,-1 F<FILE> [Required]

Path to next generation sequencing raw data. REQUIRED for both paired-end reads and single-end reads. Both compressed files and uncompressed files are allowed. 

=item --fastq2,-2 F<FILE> [Required]

Path to next generation sequencing raw data. REQUIRED for paired-end reads. Both compressed files and uncompressed files are allowed. 

=item --refSeq,-r F<File> [Required]

Path to the reference fasta file. 

=item --outputDir,-o F<FILE> [Optional]

Path of the directory to storage result files. If NOT provided, the program will generate a folder automatically.

=item --program,-p F<STRING> [Optional]

The program used for error correction and quasispecies reconstruction. Choose among 'Shorah', 'QuRe', 'PredictHaplo' and 'ViQuaS' (case insensitive). The default value is 'ViQuaS' which is faster. 

You can choose multiple programs and separate your choice by comma, e.g. --program/-p Shorah,QuRe,PredictHaplo,ViQuaS. What you should know is that running EC and QSR using multiple programs will 

take quite long time.

=item --threads,-t F<INTEGER> [Optional]

Number of threads this program will use when computing. A positive integer should be provided. The default value is 1.

=item --help,-h

Display this detailed help information.

=back

=head1 EXAMPLE

=over 5

qap ECnQSR -1 test_R1.fastq -2 test_R2.fastq -r HBV.fasta -p Shorah,ViQuaS -t 5 -o ./qsr

=back

=head1 AUTHOR

Mingjie Dr.Wang I<huzai@sjtu.edu.cn>

=head1 COPYRIGHT

Copyright (C) 2017, Mingjie Wang. All rights reserved.
