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
use File::Copy;

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
my $bamfile;
my $samfile;

my $DateNow = `date +"%Y%m%d_%Hh%Mm%Ss"`;
chomp $DateNow;

GetOptions(
'1|fastq1|=s'        => \$fq1,
'2|fastq2|=s'        => \$fq2,
'b|bamFile|=s'       => \$bamfile,
's|samFile|=s'       => \$samfile,
'r|refSeq|=s'        => \$ref,
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

if (defined $fq1){
	if(not existFile($fq1)){
		InfoError("Fastq1 file $fq1 does NOT exits.");
		exit(0);
	}
}

if (defined $fq2){
	if(not existFile($fq2)){
		InfoError("Fastq2 file $fq2 does NOT exits.");
		exit(0);
	}
}else{
	if(defined $fq1){
		InfoError("Fastq2 file MUST be provided using --fastq2/-2.");
		pod2usage(-verbose=>1,-exitval=>1);
		exit;
	}
}

if(defined $bamfile){
	if(defined $fq1 or defined $samfile){
		InfoError("Only one of --fastq1/2 or --bamFile/--samFile should be defined.");
		pod2usage(-verbose=>1,-exitval=>1);
		exit;	
	}else{
		if(existFile($bamfile) ){
			if(isBamFile($bamfile)){
				#nothing
			}else{
				InfoError("Input BAM file $bamfile is incorrect formatted.");
				exit(0);
			}
		}else{
			InfoError("BAM file $bamfile does not exist or NOT accessible.");
			exit(0);
		}
	}
}

if(defined $samfile){
	if(defined $fq1 or defined $bamfile){
		InfoError("Only one of --fastq1/2, --bamFile, and --samFile should be defined.");
		pod2usage(-verbose=>1,-exitval=>1);
		exit;	
	}else{
		if(existFile($samfile) ){
			if(isSamFile($samfile)){
				#nothing
			}else{
				InfoError("Input SAM file $samfile is incorrect formatted.");
				exit(0);
			}
		}else{
			InfoError("SAM file $samfile does not exist or NOT accessible.");
			exit(0);
		}
	}
}

my $isFileInput = (defined $fq1 && defined $fq2) || defined $bamfile || defined $samfile;
if($isFileInput){
	#nothing
}else{
	InfoError("One of --fastq1/2, --bamFile, and --samFile MUST be defined.");
	pod2usage(-verbose=>1,-exitval=>1);
	exit;
}

if(!defined $fq1){
	$fq1 = "null";
}
if(!defined $fq2){
	$fq2 = "null";
}
if(!defined $bamfile){
	$bamfile = "null";
}
if(!defined $samfile){
	$samfile = "null";
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
#check bowtie2
my $bowtie2_excu = File::Spec -> catfile($RealBin,'3rdPartyTools','bowtie2','bowtie2');
if(CheckProgram($bowtie2_excu, __FILE__, __LINE__, $DEBUG_MODE)){
	#keep running
}else{
	InfoError("The program $bowtie2_excu does NOT exist. Exiting...");
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

#manipulate files
$ref = abs_path($ref);
$fq1 = abs_path($fq1) if $fq1 ne "null";
$fq2 = abs_path($fq2) if $fq2 ne "null";
$bamfile = abs_path($bamfile) if $bamfile ne "null";
$samfile = abs_path($samfile) if $samfile ne "null";

##manipulate data input
my $tmpdir = File::Spec -> catfile($outputDir,'tmp');
makedir($tmpdir);

if($fq1 ne "null"){
	my $fq1new = File::Spec -> catfile($tmpdir,removeFastqSuffix(basename($fq1)) . ".fastq");
	copy($fq1,$fq1new) if (not -e $fq1new);
	my $fq2new = File::Spec -> catfile($tmpdir,removeFastqSuffix(basename($fq2)) . ".fastq");
	copy($fq2,$fq2new) if (not -e $fq2new);
	
	my $sampleName = getCommonString(basename($fq1),basename($fq2));
	$sampleName =~ s/[_\.R]+$//i;
	
	#generate fasta file
	my $fasta1 = File::Spec -> catfile($tmpdir, removeFastqSuffix(basename($fq1)) . ".fasta");
	system("gunzip $fq1new") if (isGzipped($fq1new));
	my $cmd = "$fq2fa_excu -in $fq1 -out $fasta1";
	runcmd($cmd) if (not existFile($fasta1));
	
	my $fasta2 = File::Spec -> catfile($tmpdir, removeFastqSuffix(basename($fq2)) . ".fasta");
	system("gunzip $fq2new") if (isGzipped($fq2new));
	$cmd = "$fq2fa_excu -in $fq2 -out $fasta2";
	runcmd($cmd) if (not existFile($fasta2));
	
	my $fastafile = File::Spec -> catfile($tmpdir, $sampleName . ".merge.fasta");
	system("cat $fasta1 $fasta2 > $fastafile") if (not existFile($fastafile));
	
	#generate sam/bam file
	my $sam = File::Spec -> catfile($tmpdir,$sampleName . ".sam");
	Bowtie2_pipeline($bowtie2_excu,$ref,$fq1,$fq2,$sam,$threads) if (not existFile($sam));	
	
	my $bamfile = $sam =~ s/\.sam$/.PosSorted.bam/r;
	sam2SortedAndIndexedBam($samtools_excu,$sam,'Pos',0,$threads) if (not existFile($bamfile));
	
	my $bam_namesorted = $sam =~ s/\.sam$/.NameSorted.bam/r;
	sam2SortedAndIndexedBam($samtools_excu,$sam,'Name',0,$threads) if (not existFile($bam_namesorted));
	
	my $samfile = $bam_namesorted =~ s/\.bam$/.sam/r;
	bam2sam($samtools_excu,$bam_namesorted,$threads) if (not existFile($samfile));
	
	##format fasta file
	my $newref = File::Spec -> catfile($tmpdir, removeFastaSuffix(basename($ref)) . ".2line.fasta");
	formatFastaToTwoLineMode($ref,$newref);
	
	##start to run
	for my $p (@program){
		qsr_pipeline($fq1new,$fq2new,$fastafile,$bamfile, $samfile, $newref,$outputDir,$p,$threads,$shorah_excu,$predicthaplo_excu,$qure_excu_dir,$viquas_excu_dir);
	}
	
}elsif($bamfile ne "null"){
	my $bamfile_new = File::Spec -> catfile($tmpdir, removeBamSuffix(basename($bamfile)) . ".bam");
	copy($bamfile, $bamfile_new) if (not existFile($bamfile_new));
	$bamfile = $bamfile_new;
	
	#generate sam/bam file
	my $bamfile_PosSort = $bamfile =~ s/\.bam$/.PosSorted.bam/r;
	sortAndIndexBam($samtools_excu,$bamfile,'Pos',0,$threads) if (not existFile($bamfile_PosSort));
	
	my $bamfile_NameSorted = $bamfile =~ s/\.bam$/.NameSorted.bam/r;
	sortAndIndexBam($samtools_excu,$bamfile,'Name',0,$threads) if (not existFile($bamfile_NameSorted));
	
	my $samfile = $bamfile_NameSorted =~ s/\.bam$/.sam/r;
	bam2sam($samtools_excu,$bamfile_NameSorted,$threads) if (not existFile($samfile));
	
	#generate fastq files
	my $fq1new = File::Spec -> catfile($tmpdir,removeAllSuffix(basename($bamfile)) . "_R1.fq");
	my $fq2new = File::Spec -> catfile($tmpdir,removeAllSuffix(basename($bamfile)) . "_R2.fq");
	#check picard
	my $picard_excu = File::Spec -> catfile($RealBin,'3rdPartyTools','caller','picard.jar');
	if(existFile($picard_excu)){
		#keep running
	}else{
		InfoError("The program $picard_excu does NOT exist. Exiting...");
		exit;
	}
	my $cmd = "java -jar $picard_excu SamToFastq INPUT=$bamfile FASTQ=$fq1new SECOND_END_FASTQ=$fq2new";
	runcmd($cmd);
	
	my $sampleName = getCommonString(basename($fq1new),basename($fq2new));
	$sampleName =~ s/[_\.R]+$//i;
	
	#generate fasta file
	my $fasta1 = File::Spec -> catfile($tmpdir, removeFastqSuffix(basename($fq1new)) . ".fasta");
	system("gunzip $fq1new") if (isGzipped($fq1new));
	$cmd = "$fq2fa_excu -in $fq1new -out $fasta1";
	runcmd($cmd) if (not existFile($fasta1));
	
	my $fasta2 = File::Spec -> catfile($tmpdir, removeFastqSuffix(basename($fq2new)) . ".fasta");
	system("gunzip $fq2new") if (isGzipped($fq2new));
	$cmd = "$fq2fa_excu -in $fq2new -out $fasta2";
	runcmd($cmd) if (not existFile($fasta2));
	
	my $fastafile = File::Spec -> catfile($tmpdir, $sampleName . ".merge.fasta");
	system("cat $fasta1 $fasta2 > $fastafile") if (not existFile($fastafile));
	
	##format fasta file
	my $newref = File::Spec -> catfile($tmpdir, removeFastaSuffix(basename($ref)) . ".2line.fasta");
	formatFastaToTwoLineMode($ref,$newref);
	
	##start to run
	for my $p (@program){
		qsr_pipeline($fq1new,$fq2new,$fastafile,$bamfile_PosSort, $samfile, $newref,$outputDir,$p,$threads,$shorah_excu,$predicthaplo_excu,$qure_excu_dir,$viquas_excu_dir);
	}
}elsif($samfile ne "null"){
	my $samfile_new = File::Spec -> catfile($tmpdir, removeSamSuffix(basename($samfile)) . ".sam");
	copy($samfile, $samfile_new) if (not existFile($samfile_new));
	$samfile = $samfile_new;
	
	#generate sam/bam file
	my $bamfile_PosSort = $samfile =~ s/\.sam$/.PosSorted.bam/r;
	sam2SortedAndIndexedBam($samtools_excu,$samfile,'Pos',0,$threads) if (not existFile($bamfile_PosSort));
	
	my $bamfile_NameSorted = $bamfile_PosSort =~ s/\.bam$/.NameSorted.bam/r;
	sortAndIndexBam($samtools_excu,$bamfile_PosSort,'Name',0,$threads) if (not existFile($bamfile_NameSorted));
	
	my $samfile = $bamfile_NameSorted =~ s/\.bam$/.sam/r;
	bam2sam($samtools_excu,$bamfile_NameSorted,$threads) if (not existFile($samfile));
	
	#generate fastq files
	my $fq1new = File::Spec -> catfile($tmpdir,removeAllSuffix(basename($bamfile_PosSort)) . "_R1.fq");
	my $fq2new = File::Spec -> catfile($tmpdir,removeAllSuffix(basename($bamfile_PosSort)) . "_R2.fq");
	#check picard
	my $picard_excu = File::Spec -> catfile($RealBin,'3rdPartyTools','caller','picard.jar');
	if(existFile($picard_excu)){
		#keep running
	}else{
		InfoError("The program $picard_excu does NOT exist. Exiting...");
		exit;
	}
	my $cmd = "java -jar $picard_excu SamToFastq INPUT=$bamfile_PosSort FASTQ=$fq1new SECOND_END_FASTQ=$fq2new";
	runcmd($cmd);
	
	my $sampleName = getCommonString(basename($fq1new),basename($fq2new));
	$sampleName =~ s/[_\.R]+$//i;
	
	#generate fasta file
	my $fasta1 = File::Spec -> catfile($tmpdir, removeFastqSuffix(basename($fq1new)) . ".fasta");
	system("gunzip $fq1new") if (isGzipped($fq1new));
	$cmd = "$fq2fa_excu -in $fq1new -out $fasta1";
	runcmd($cmd) if (not existFile($fasta1));
	
	my $fasta2 = File::Spec -> catfile($tmpdir, removeFastqSuffix(basename($fq2new)) . ".fasta");
	system("gunzip $fq2new") if (isGzipped($fq2new));
	$cmd = "$fq2fa_excu -in $fq2new -out $fasta2";
	runcmd($cmd) if (not existFile($fasta2));
	
	my $fastafile = File::Spec -> catfile($tmpdir, $sampleName . ".merge.fasta");
	system("cat $fasta1 $fasta2 > $fastafile") if (not existFile($fastafile));
	
	##format fasta file
	my $newref = File::Spec -> catfile($tmpdir, removeFastaSuffix(basename($ref)) . ".2line.fasta");
	formatFastaToTwoLineMode($ref,$newref);
	
	##start to run
	for my $p (@program){
		qsr_pipeline($fq1new,$fq2new,$fastafile,$bamfile_PosSort, $samfile, $newref,$outputDir,$p,$threads,$shorah_excu,$predicthaplo_excu,$qure_excu_dir,$viquas_excu_dir);
	}
}else{
	InfoError("Something wrong with input arguments. Exiting...");
	exit(0);
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

=item --fastq1,-1 F<FILE> [Optional]

Path to next generation sequencing raw data. REQUIRED for both paired-end reads and single-end reads. Both compressed files and uncompressed files are allowed. 

=item --fastq2,-2 F<FILE> [Optional]

Path to next generation sequencing raw data. REQUIRED for paired-end reads. Both compressed files and uncompressed files are allowed. 

=item --bamFile,-b F<FILE> [Optional]

Path to the input BAM file. 

=item --samFile,-s F<FILE> [Optional]

Path to the input SAM file. One of input fastq or BAM or SAM file must be provided.

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
