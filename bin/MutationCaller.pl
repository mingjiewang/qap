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
use Mapper;
use Caller;

####---------------------------####
####The program begins here
####---------------------------####

##Show welcome
print "You are now running subprogram: ";
printcol ("MutationCaller","green");
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
my $inputfile;
my $outputDir;
my $threads;
my $ref;
my $program;

my $DateNow = `date +"%Y%m%d_%Hh%Mm%Ss"`;
chomp $DateNow;

GetOptions(
'i|inputFile|=s'     => \$inputfile,
'r|refSeq|=s'        => \$ref,
'o|outputDir=s'      => \$outputDir,
'h|help|'            => \$help,
't|threads|=s'       => \$threads,
'p|program|=s'       => \$program,
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
			pod2usage(-verbose=>2,-exitval=>1);
			exit;
		}
	}else{
		InfoError("Folder $outputDir already exists. Please specify another output directory using option -o/--outputDir");
		pod2usage(-verbose=>2,-exitval=>1);
		exit;
	}
}else{
	$outputDir = File::Spec -> catfile($wk_dir,"qap_Results_for_MutationCaller_$DateNow");
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


if(defined $inputfile){
	if(existFile($inputfile)){
		#nothing
	}else{
		InfoError("$inputfile doesn't exist. Please check.");
		exit(0);
	}
}else{
	InfoError("Input bam/sam file MUST be provided.");
	pod2usage(-verbose=>2,-exitval=>1);
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

my @program;
my @allowedProgram = qw/gatk lofreq varscan/;
if(defined $program){
	my @tmpprogram = split ',',$program;
	for my $p (@tmpprogram){
		if(isInArray($p,\@allowedProgram)){
			push @program,$p;
		}else{
			InfoError("$p is not allowed. Please choose a mutation calling program from \'gatk\' or \'lofreq\' or \'varscan\'.");
			pod2usage(-verbose=>2,-exitval=>1);
			exit;
		}
	}
}else{
	InfoWarn("Mutation calling program is not provided. Will use \'-p/--program gatk,lofreq,varscan\' as default.");
	@program = qw/gatk lofreq varscan/;
}


##the core program starts here
#check all the programs
my $DEBUG_MODE = 1;
#check samtools
my $samtools_excu = File::Spec -> catfile($mainBin,'bin','3rdPartyTools','samtools','samtools');
if(CheckProgram($samtools_excu, __FILE__, __LINE__, $DEBUG_MODE)){
	#keep running
}else{
	InfoError("The program $samtools_excu does NOT exist. Exiting...");
	exit;
}
#check picard
my $picard_excu = File::Spec -> catfile($mainBin,'bin','3rdPartyTools','caller','picard.jar');
if(CheckProgram($picard_excu, __FILE__, __LINE__, $DEBUG_MODE)){
	#keep running
}else{
	InfoError("The program $picard_excu does NOT exist. Exiting...");
	exit;
}
#check gatk
my $gatk_excu = File::Spec -> catfile($mainBin,'bin','3rdPartyTools','caller','GenomeAnalysisTK.jar');
if(CheckProgram($gatk_excu, __FILE__, __LINE__, $DEBUG_MODE)){
	#keep running
}else{
	InfoError("The program $gatk_excu does NOT exist. Exiting...");
	exit;
}
#check lofreq
my $lofreq_excu = File::Spec -> catfile($mainBin,'bin','3rdPartyTools','caller','lofreq');
if(CheckProgram($lofreq_excu, __FILE__, __LINE__, $DEBUG_MODE)){
	#keep running
}else{
	InfoError("The program $lofreq_excu does NOT exist. Exiting...");
	exit;
}
#check varscan
my $varscan_excu = File::Spec -> catfile($mainBin,'bin','3rdPartyTools','caller','VarScan.jar');
if(CheckProgram($varscan_excu, __FILE__, __LINE__, $DEBUG_MODE)){
	#keep running
}else{
	InfoError("The program $varscan_excu does NOT exist. Exiting...");
	exit;
}
#check fmtOutputPerlScript
my $fmtOutputPerlScript = File::Spec -> catfile($mainBin,'bin','PerlScripts','fmtVarscan2GATK.pl');
if(CheckFile($fmtOutputPerlScript, __FILE__, __LINE__, $DEBUG_MODE)){
	#keep running
}else{
	InfoError("The program $fmtOutputPerlScript does NOT exist. Exiting...");
	exit;
}
#check getMutationFrequency Perl Script
my $getMutFreqPerlScript = File::Spec -> catfile($mainBin,'bin','PerlScripts','getMutationFrequency.pl');
if(CheckFile($getMutFreqPerlScript, __FILE__, __LINE__, $DEBUG_MODE)){
	#keep running
}else{
	InfoError("The program $getMutFreqPerlScript does NOT exist. Exiting...");
	exit;
}
#check merge mutfreq R script
my $mergeMutFreqRScript = File::Spec -> catfile($mainBin, 'bin', 'Rscripts', 'MergeMutationFrequency.R');
if(CheckFile($mergeMutFreqRScript, __FILE__, __LINE__, $DEBUG_MODE)){
	#keep running
}else{
	InfoError("The program $mergeMutFreqRScript does NOT exist. Exiting...");
	exit;
}

#generate .dict and .fai file for $ref
prepareReference($ref,$picard_excu,$samtools_excu);

#start to run caller
mutationCallerPipeline($inputfile, $outputDir, $ref, \@program, $samtools_excu, $picard_excu, $gatk_excu, $lofreq_excu, $varscan_excu, $fmtOutputPerlScript, $getMutFreqPerlScript, $mergeMutFreqRScript, $threads);

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
           



qap MutationCaller [options]

Use --help to see more information.

qap is still in development. If you have encounted any problem in usage, please feel no hesitation to cotact us.

=head1 DESCRIPTION

This script implements a function to call variants from bam/sam files. The script has B<several> mandatory options that MUST appear last. 

=head1 OPTIONS

=over 5

=item --inputFile,-i F<FILE> [Required]

Path to input bam/sam file which is required.

=item --refSeq,-r F<File> [Required]

Path to the reference fasta file. 

=item --program,-p F<STRING> [Optional]

The program used for mutation calling. Choose between 'gatk', 'lofreq' and 'varscan'. Multiple callers are allowed, which should be seperated by comma, e.g. -p gatk,lofreq.

If multiple calling programs are selected, the program will generate a merged vcf file.

=item --outputDir,-o F<FILE> [Optional]

Path of the directory to storage result files. If NOT provided, the program will generate a folder automatically.

=item --threads,-t F<INTEGER> [Optional]

Number of threads this program will use when computing. A positive integer should be provided. The default value is 1.

=item --help,-h

Display this detailed help information.

=back

=head1 EXAMPLE

=over 5

qap MutationCaller -i test.bam -r HBV.fasta -p gatk,lofreq -t 5 -o ./mutation

=back

=head1 AUTHOR

Mingjie Dr.Wang I<huzai@sjtu.edu.cn>

=head1 COPYRIGHT

Copyright (C) 2017, Mingjie Wang. All rights reserved.

