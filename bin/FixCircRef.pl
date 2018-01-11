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
use File::Copy;
use List::Util qw/min max/;

####Use modules in this program####
use General;
use Mapper;
use REF;

####---------------------------####
####The program begins here
####---------------------------####

##Show welcome
print "You are now running subprogram: ";
printcol ("FixCircRef","green");
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
my $threads;
my $ref;
my $trim;

my $DateNow = `date +"%Y%m%d_%Hh%Mm%Ss"`;
chomp $DateNow;

GetOptions(
'1|fastq1|=s'        => \$fq1,
'2|fastq2|=s'        => \$fq2,
'r|refSeq|=s'        => \$ref,
'o|outputDir|=s'     => \$outputDir,
'h|help|'            => \$help,
't|threads|=s'       => \$threads,
'trim|=s'            => \$trim
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
		my $cmd = "mkdir -p $outputDir";
		system($cmd);
	}else{
		InfoError("Mkdir Failed! Folder $outputDir already exists!","red");
		InfoError("Please specify another output directory using option -o/--outputDir");
		pod2usage(-verbose=>2,-exitval=>1);
		exit;
	}
}else{
	$outputDir = File::Spec -> catfile($wk_dir,"qap_Results_for_FixCircRef_$DateNow");
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

if(defined $trim){
	if(uc($trim) eq 'Y' or uc($trim) eq 'N'){
		#nothing
	}else{
		InfoError("--trim $trim is not allowed. Please choose between \'Y\' or \'N\'.");
		pod2usage(-verbose=>1,-exitval=>1);
		exit;
	}
}else{
	$trim = 'Y';
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
	if(scalar(@fq1files) < 3){
		InfoError("At least 3 fastq1 files are required.");
		pod2usage(-verbose=>1,-exitval=>1);
		exit;
	}
}else{
	InfoError("Input fastq1 files MUST be provided!",'red');
	pod2usage(-verbose=>1,-exitval=>1);
	exit;
}

my @fq2files;
if (defined $fq2){
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
	if(scalar(@fq2files) < 3){
		InfoError("At least 3 fastq2 files are required.");
		pod2usage(-verbose=>1,-exitval=>1);
		exit;
	}
}else{
	InfoError("Input fastq2 files MUST be provided!",'red');
	pod2usage(-verbose=>1,-exitval=>1);
	exit;
}

##main program starts here
my $tmpdir = File::Spec -> catfile($outputDir,'tmp');
makedir($tmpdir);

##check ref file 
copy($ref,$tmpdir);
my $refname = basename($ref);
$ref = File::Spec -> catfile($tmpdir,$refname);
my $newref = removeFastaSuffix($ref) . '.2line' . ".fasta";
formatFastaToTwoLineMode($ref,$newref);

my $seqnum = CheckSeqNum($newref);
if($seqnum != 1){
	InfoError("There are $seqnum sequences in reference file $refname, which should be only one.");
	exit;
}

##check tools
my $DEBUG_MODE = 1;
#check bwa
my $bwaProgram = File::Spec -> catfile($RealBin,'3rdPartyTools','bwa','bwa');
if(CheckProgram($bwaProgram, __FILE__, __LINE__, $DEBUG_MODE)){
	#keep running
}else{
	InfoError("The program $bwaProgram does NOT exist. Exiting...");
	exit;
}
#check bowtie2
my $bowtie2Program = File::Spec -> catfile($RealBin,'3rdPartyTools','bowtie2','bowtie2');
if(CheckProgram($bowtie2Program, __FILE__, __LINE__, $DEBUG_MODE)){
	#keep running
}else{
	InfoError("The program $bowtie2Program does NOT exist. Exiting...");
	exit;
}
#check samtools
my $samtoolsProgram = File::Spec -> catfile($RealBin,'3rdPartyTools','samtools','samtools');
if(CheckProgram($samtoolsProgram, __FILE__, __LINE__, $DEBUG_MODE)){
	#keep running
}else{
	InfoError("The program $samtoolsProgram does NOT exist. Exiting...");
	exit;
}
#check rscript
my $rscript = File::Spec -> catfile($RealBin,'Rscripts','FindDepthInterval.R');
if(not existFile($rscript)){
	InfoError("The R script $rscript used for reference correction is missing.");
	exit;
}

##copy head to tail 500bp
my $fixref0 = $newref =~ s/.fasta$/.fix0.fasta/r;
copyRegion($newref,$fixref0,'head2tail',500);

##sam
my $samdir = File::Spec -> catfile($tmpdir,'bwa');
makedir($samdir);
my $depthdir = File::Spec -> catfile($tmpdir,'depth');
makedir($depthdir);

my @intervalfiles;
for my $i (0..scalar(@fq1files) - 1){
	my $sam = File::Spec -> catfile($samdir,"${i}.sam");
	#BWA_pipeline($bwaProgram,$fixref0,$fq1files[$i],$fq2files[$i],$sam,$threads);
	Bowtie2_pipeline($bowtie2Program,$fixref0,$fq1files[$i],$fq2files[$i],$sam,$threads);
	
	my $bam = File::Spec -> catfile($samdir,"${i}.PosSorted.bam");
	sam2SortedAndIndexedBam($samtoolsProgram,$sam,'POS',0,$threads);
	
	my $depth = File::Spec -> catfile($depthdir,"${i}.depth");
	my $cmd = "$samtoolsProgram depth -aa $bam > $depth";
	runcmd($cmd);
	
	my $interval = File::Spec -> catfile($depthdir,"${i}.interval");
	$cmd = "Rscript $rscript -i $depth -o $interval";
	runcmd($cmd);
	push @intervalfiles,$interval;
}

my %pos;
for my $interval (@intervalfiles){
	open I,$interval or die "Can NOT open interval file $interval:$!";
	while(<I>){
		chomp;
		$pos{$_}++;
	}
	close I;
}

my $position = File::Spec -> catfile($tmpdir,"FixPosition.txt");
open POS,">$position" or die "Can NOT output to $position:$!";
for my $k (sort {$pos{$a} <=> $pos{$b}} keys %pos){
	if($pos{$k} > 1){
		print POS "$k\n";
	}
}
close POS;

##manipulate ref file
open T,$position or die "Can NOT open $position:$!";
#chomp(my $finalStart = <T>); #the first line of $position should be the real start
my @pos;
while(<T>){
	chomp;
	push @pos,$_;
}
close T;

my $start = 1;
$start = min(@pos);
my $end = checkSeqLen($ref,$tmpdir);
if($trim =~ /y/i){
	$end = $start + checkSeqLen($ref,$tmpdir) - 1;
}else{
	$end = max(@pos);
}

my $finalRef = File::Spec -> catfile($outputDir, removeFastaSuffix(basename($ref)) . ".FixedRef.fasta");
#moveRegion($newref,$finalRef,'head2tail',$finalStart - 1);
subRegion($fixref0,$finalRef,$start,$end);
my $cmd = "sed -i \'s\/\\r\/\/g\' $finalRef";
system($cmd);

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
           



qap FixCircRef [options]

Use --help to see more information.

qap is still in development. If you have encounted any problem in usage, please feel no hesitation to cotact us.

=head1 DESCRIPTION

This script implements a function to fix reference files when studying virus with circular genome by using amplicon sequencing. The script has B<several> mandatory options that MUST appear last. 

=head1 OPTIONS

=over 5

=item --fastq1,-1 F<FILE> [Required]

Path to next generation sequencing raw data. REQUIRED for both paired-end reads and single-end reads. Both compressed files and uncompressed files are allowed. To get a more precise reference, multiple read files are required [at least 3 files]. Multiple files should be seperated by comma, e.g. test1_R1.fq,test2_R1.fq,test3_R1.fq 

=item --fastq2,-2 F<FILE> [Required]

Path to next generation sequencing raw data. REQUIRED for paired-end reads. Both compressed files and uncompressed files are allowed. Just like fastq1 files, multiple read files are required [at least 3 files]. When providing multiple read files, please make sure fastq2 files match the order of fastq1 files. 

=item --refSeq,-r F<File> [Required]

Path to the reference fasta file. 

=item --trim F<STRING> [Optional]

If the rebuilded reference sequence is longer than the real reference, choose how to trim the additional bases. 'Y' for trim and 'N' for keep. The default value is --trim Y.

=item --outputDir,-o F<FILE> [Optional]

Path of the directory to storage result files. If NOT provided, the program will generate a folder automatically.

=item --threads,-t F<INTEGER> [Optional]

Number of threads this program will use when computing. A positive integer should be provided. The default value is 1.

=item --help,-h

Display this detailed help information.

=back

=head1 EXAMPLE

=over 5

qap FixCircRef -1 test_R1.fq -2 test_R2.fq -r HBV.fasta -t 5 -o ./newref

=back

=head1 AUTHOR

Mingjie Dr.Wang I<huzai@sjtu.edu.cn>

=head1 COPYRIGHT

Copyright (C) 2017, Mingjie Wang. All rights reserved.




