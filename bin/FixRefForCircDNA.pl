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
use Mapper;
use Caller;

####---------------------------####
####The program begins here
####---------------------------####

##Show welcome
print "You are now running subprogram: ";
printcol ("FixRefForCircDNA","green");
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
my $inputfile;
my $outputDir;
my $threads;
my $ref;
my $spanLength;
my $mode;
my $graphic;
my $cutoff;

my $DateNow = `date +"%Y%m%d_%Hh%Mm%Ss"`;
chomp $DateNow;

GetOptions(
'1|fastq1|=s'        => \$fastq1,
'2|fastq2|=s'        => \$fastq2,
'r|reference|=s'     => \$ref,
'o|outputDir|=s'     => \$outputDir,
'h|help|'            => \$help,
't|threads|=s'       => \$threads
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
		my $cmd = "mkdir -p $outputDir";
		system($cmd);
	}else{
		InfoError("Folder $outputDir already exists. Please specify another output directory using option -o/--outputDir");
		exit;
	}
}else{
	$outputDir = File::Spec -> catfile($wk_dir,"qap_Results_for_FixRefForCircDNA_$DateNow");
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


##main program starts here



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
           



gap FixRefForCircDNA [options]

Use --help to see more information.

gap is still in development. If you have encounted any problem in usage, please feel no hesitation to cotact us.

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

=item --threads,-t F<INTEGER> [Optional]

Number of threads this program will use when computing. A positive integer should be provided. The default value is 1.

=item --help,-h

Display this detailed help information.

=back

=head1 EXAMPLE

=over 5

gap FixRefForCircDNA -1 test_R1.fq -2 test_R2.fq -r HBV.fasta -t 5 -o ./newref

=back

=head1 AUTHOR

Mingjie Dr.Wang I<huzai@sjtu.edu.cn>

=head1 COPYRIGHT

Copyright (C) 2017, Mingjie Wang. All rights reserved.




