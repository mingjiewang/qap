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

####---------------------------####
####The program begins here
####---------------------------####

##Show welcome
print "You are now running subprogram: ";
printcol ("MapReadsToRef","green");
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
my $pairEnded = 0;
my $threads;
my $ref;
my $outFormat;
my $sortMethod;
my $program;
my $sampleLabel;

my $DateNow = `date +"%Y%m%d_%Hh%Mm%Ss"`;
chomp $DateNow;

GetOptions(
'1|fastq1|=s'        => \$fq1,
'2|fastq2|=s'        => \$fq2,
'r|refSeq|=s'        => \$ref,
'o|outputDir=s'      => \$outputDir,
'h|help|'            => \$help,
't|threads|=s'       => \$threads,
'f|outFormat|=s'     => \$outFormat,
's|sort|=s'          => \$sortMethod,
'p|program|=s'       => \$program,
'l|sampleLabel|=s'   => \$sampleLabel
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
			pod2usage(-verbose=>2,-exitval=>1);
			exit;
		}
	}
}else{
	$outputDir = File::Spec -> catfile($wk_dir,"qap_Results_for_MapReadsToRef_$DateNow");
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

my @sampleLabel;
my @samFiles;
my $numberOfSamples;
if (defined $sampleLabel){
	@sampleLabel = split ",",$sampleLabel;
	for my $s (@sampleLabel){
		my $sam = File::Spec -> catfile($outputDir, $s . ".sam");
		push @samFiles,$sam;
	}
	$numberOfSamples = scalar(@sampleLabel);
}else{
	InfoError("Label of each sample must be provided!",'red');
	pod2usage(-verbose=>1,-exitval=>1);
	exit;
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
	InfoError("Input fastq1 files MUST be provided!",'red');
	pod2usage(-verbose=>2,-exitval=>1);
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
	
	#check whether equal
	my $num1 = scalar(@sampleLabel);
	my $num2 = scalar(@fq1files);
	my $num3 = scalar(@fq2files);
	
	if (check3EqualNums($num1, $num2, $num3)){
		#nothing
	}else{
		InfoError("There are $num1 sample labels, $num2 fastq1 files and $num3 fastq2 files provided! Please check again.", "red");
		pod2usage(-verbose=>1,-exitval=>1);
		exit;
	}
	
}else{
	for my $i (@fq1files){
		push @fq2files,'';
	}
}

my @ref;
my $isMultiRef = 0;
if (defined $ref){
	if ($ref !~ /\,/){
		if (not -e $ref){
			InfoError("Input reference file $ref does NOT exist!",'red');
			exit;
		}else{
			if(isGzipped($ref)){
				InfoError("The reference file you provided is gzipped. Please uncompress it using \"gunzip $ref\" and re-run the program. ");
				exit;
			}
			$ref = abs_path($ref);
		}
	}else{
		my @tmp = split ',',$ref;
		
		$isMultiRef = 1;
		
		#check the number of files provided are equal or not 
		my $num1 = scalar(@sampleLabel);
		my $num2 = scalar(@fq1files);
		my $num3 = scalar(@fq2files);
		my $num4 = scalar(@tmp);
		my @numToCheck = ($num1, $num2, $num3, $num4);
		
		if (checkMultipleEqualNums(\@numToCheck)){
			#nothing
		}else{
			InfoError("There are $num1 sample labels, $num2 fastq1 files, $num3 fastq2 files and $num4 reference files provided! Please check again.", "red");
			pod2usage(-verbose=>1,-exitval=>1);
			exit;
		}
		
		#check ref file status
		for my $r (@tmp){
			if (not -e $r){
				InfoError("Input reference file $r does NOT exist!",'red');
				exit;
			}else{
				if(isGzipped($r)){
					InfoError("The reference file you provided is gzipped. Please uncompress it using \"gunzip $r\" and re-run the program. ");
					exit;
				}
				my $rpath = abs_path($r);
				push @ref,$rpath;
			}
		}
	}
	
}else{
	InfoError("Input reference file MUST be provided!",'red');
	pod2usage(-verbose=>2,-exitval=>1);
	exit;
}

##get mapping program
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

if (defined $program){
	if(uc($program) eq 'BWA' || uc($program) eq 'BOWTIE2'){
		#nothing
	}else{
		InfoError("The mapping program MUST be one of \"bwa\" or \"bowtie2\".","red");
		pod2usage(-verbose=>2,-exitval=>1);
		exit;
	}
}else{
	InfoWarn("The mapping program --program/-p is not provided. Using \"bwa\" as default.");
	$program = 'bwa';
}

if (defined $sortMethod){
	if(uc($sortMethod) eq 'NAME' || uc($sortMethod) eq 'POS'){
		#nothing
	}else{
		InfoError("The output bam/sam MUST be sorted by \"name\" or \"pos\".","red");
		pod2usage(-verbose=>2,-exitval=>1);
		exit;
	}
}else{
	InfoWarn("Using \"--sort/-s pos\" as default. The output bam/sam file will be sorted by position.");
	$sortMethod = 'pos';
}

if (defined $outFormat){
	if(uc($outFormat) eq 'BAM' || uc($outFormat) eq 'SAM'){
		#nothing
	}else{
		InfoError("The output format MUST be one of \"bam\" or \"sam\".","red");
		pod2usage(-verbose=>2,-exitval=>1);
		exit;
	}
}else{
	InfoWarn("The output format --outFormat/-f is not provided. Using \"bam\" as default.");
	$outFormat = 'bam';
}

##core program starts here
if($threads > 1){
	my $threadsPerSample = int($threads / $numberOfSamples);
	if ($threadsPerSample > 1){
		#get input args ready
		my @threadsForEachSample;
		my @bowtie2Program;
		my @bwaProgram;
		my @samtoolsProgram;
		my @newref = ();
		if ($isMultiRef){
			@newref = @ref;
		}
		my @sortMethod;
		my @removeSam;
		my @bamfiles;
		for my $s (@samFiles){
			push @threadsForEachSample, $threadsPerSample;
			push @bowtie2Program, $bowtie2Program;
			push @bwaProgram, $bwaProgram;
			push @samtoolsProgram, $samtoolsProgram;
			push @newref, $ref if not $isMultiRef;
			push @sortMethod, $sortMethod;
			push @removeSam,0;
			if(uc($sortMethod) eq 'POS'){
				my $outbam = $s =~ s/\.sam$/.PosSorted.bam/r if isSamFile($s);
				push @bamfiles,$outbam;
			}else{
				my $outbam = $s =~ s/\.sam$/.NameSorted.bam/r if isSamFile($s);
				push @bamfiles,$outbam;
			}
		}
		

		# even though BWA_pipeline or Bowtie2_pipeline includes the index step, however index step MUST not go parallel, 
		# espically different samples are using the same index file. Because it will directly jump to the mapping step
		# before indexing. So index first.
		for my $f (@newref){
			if(uc($program) eq 'BWA'){
				BWA_index($bwaProgram, $f);
			}else{
				Bowtie2_index($bowtie2Program, $f);
			}
		}
		
		
		if (uc($program) eq 'BWA'){
			# my ($BWA_excu,$ref,$fq1,$fq2,$sam,$threads) = @_;
			#BWA_pipeline($bwaProgram, $ref, $fq1, $fq2, $sam, $threads);
			runMultipleThreadsWith6Args(\&BWA_pipeline, \@bwaProgram, \@newref, \@fq1files, \@fq2files, \@samFiles, \@threadsForEachSample, $numberOfSamples);
		}else{
			# my($Bowtie2_excu,$ref,$fq1,$fq2,$sam,$threads) = @_;
			#Bowtie2_pipeline($bowtie2Program, $ref, $fq1, $fq2, $sam, $threads);
			runMultipleThreadsWith6Args(\&Bowtie2_pipeline, \@bowtie2Program, \@newref, \@fq1files, \@fq2files, \@samFiles, \@threadsForEachSample, $numberOfSamples);
		}
		
		#runMultipleThreadsWith4Args(\&runMapping, \@fq1files, \@fq2files, \@samFiles, \@threadsForEachSample, $numberOfSamples);
		#sam2SortedAndIndexedBam($samtoolsProgram,$samFile,$sortMethod,0,$threads);
		#my ($samtoolsProgram,$samFile,$sortMethod,$removeSam,$threads) = @_;
		runMultipleThreadsWith5Args(\&handleSam, \@samtoolsProgram, \@samFiles, \@sortMethod, \@removeSam,\@threadsForEachSample, $numberOfSamples);
		
		if(uc($outFormat) eq 'SAM'){
			runMultipleThreadsWith3Args(\&bam2sam, \@samtoolsProgram, \@bamfiles, \@threadsForEachSample, $numberOfSamples);
		}
	}else{
		for my $i(0..scalar(@fq1files) - 1){
			my $newref;
			if($isMultiRef){
				$newref = $ref[$i];
			}else{
				$newref = $ref;
			}
			my $fq1 = $fq1files[$i];
			my $fq2 = $fq2files[$i];
			my $samfile = $samFiles[$i];
			my $bamfile;
			if(uc($sortMethod) eq 'POS'){
				$bamfile = $samfile =~ s/\.sam$/.PosSorted.bam/r if isSamFile($samfile);
			}else{
				$bamfile = $samfile =~ s/\.sam$/.NameSorted.bam/r if isSamFile($samfile);
			}
			
			
			&runMapping($newref, $fq1, $fq2, $samfile, $threads);
			
			&handleSam($samtoolsProgram,$samfile,$sortMethod,0,$threads);
			
			if(uc($outFormat) eq 'SAM'){
				bam2sam($samtoolsProgram, $bamfile, $threads);
			}
		}
	}
}else{
	for my $i(0..scalar(@fq1files) - 1){
		my $newref;
		if($isMultiRef){
			$newref = $ref[$i];
		}else{
			$newref = $ref;
		}
		my $fq1 = $fq1files[$i];
		my $fq2 = $fq2files[$i];
		my $samfile = $samFiles[$i];
		my $bamfile;
		if(uc($sortMethod) eq 'POS'){
			$bamfile = $samfile =~ s/\.sam$/.PosSorted.bam/r if isSamFile($samfile);
		}else{
			$bamfile = $samfile =~ s/\.sam$/.NameSorted.bam/r if isSamFile($samfile);
		}
		
		
		&runMapping($newref, $fq1, $fq2, $samfile, $threads);
		
		&handleSam($samtoolsProgram,$samfile,$sortMethod,0,$threads);
		
		if(uc($outFormat) eq 'SAM'){
			bam2sam($samtoolsProgram, $bamfile, $threads);
		}
	}
}


##sub program goes here
sub runMapping {
	my $ref = shift;
	my $fq1 = shift;
	my $fq2 = shift;
	my $sam = shift;
	my $threads = shift;
	
	if (uc($program) eq 'BWA'){
		# my ($BWA_excu,$ref,$fq1,$fq2,$sam,$threads) = @_;
		BWA_pipeline($bwaProgram, $ref, $fq1, $fq2, $sam, $threads);
	}else{
		# my($Bowtie2_excu,$ref,$fq1,$fq2,$sam,$threads) = @_;
		Bowtie2_pipeline($bowtie2Program, $ref, $fq1, $fq2, $sam, $threads);
	}
}

sub handleSam {
	my ($samtoolsProgram,$samFile,$sortMethod,$removeSam,$threads) = @_;
	
	sam2SortedAndIndexedBam($samtoolsProgram,$samFile,$sortMethod,$removeSam,$threads);
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
           



qap MapReadsToRef [options]

Use --help to see more information.

qap is still in development. If you have encounted any problem in usage, please feel no hesitation to cotact us.

=head1 DESCRIPTION

This script implements a function to map fastq/fasta sequencing data to the reference fasta file and generate a sam/bam file. The script has B<several> mandatory options that MUST appear last. 

=head1 OPTIONS

=over 5

=item --fastq1,-1 F<FILE> [Required]

Path to next generation sequencing raw data. REQUIRED for both paired-end reads and single-end reads. Both compressed files and uncompressed files are allowed. Several files are allowed and should be seperated by comma, e.g. -1 Data1_R1.fq.gz,Data2_R1.fq.gz,Data3_R1.fq.gz . 

=item --fastq2,-2 F<FILE> [Optional]

Path to next generation sequencing raw data. REQUIRED for paired-end reads. Both compressed files and uncompressed files are allowed. Several files are allowed and should be seperated by comma, e.g. -2 Data1_R2.fq.gz,Data2_R2.fq.gz,Data3_R2.fq.gz. The number of fastq2 files should equal the number of fastq1 files.

=item --sampleLabel, -l F<STRING> [Required]

The name of all samples. Multiple samples should be seperated by comma. The number of sample labels provided should equal to the number of fastq1/fastq2.

=item --refSeq,-r F<File> [Required]

Path to the reference fasta file. If all the sample use the same reference file, only one reference file should be provided. If the samples used different references, please provide all the reference files for each sequence data, and separated them by comma. e.g. -r hbv.fasta or -r hbv1.fasta,hbv2.fasta,hbv3.fasta. 

=item --program,-p F<STRING> [Optional]

The program used for mapping. Choose one between 'bwa' and 'bowtie2'. The default value is 'bwa'.

=item --outFormat,-f F<STRING> [Optional]

The format of the output result. Choose one between 'bam' or 'sam'. The default value is 'bam'.

=item --sort,-s F<STRING> [Optional]

The way to sort sam/bam output. Choose one between 'pos' and 'name' by with 'pos' means 'sorted by position' and 'name' means 'sorted by read name'. The default value is 'pos'. 

=item --outputDir,-o F<FILE> [Optional]

Path of the directory to storage result files. If NOT provided, the program will generate a folder automatically.

=item --threads,-t F<INTEGER> [Optional]

Number of threads this program will use when computing. A positive integer should be provided. The default value is 1.

=item --help,-h

Display this detailed help information.

=back

=head1 EXAMPLE

=over 5

qap MapReadsToRef -1 Data1_R1.fq.gz,Data2_R1.fq.gz,Data3_R1.fq.gz -2 Data1_R2.fq.gz,Data2_R2.fq.gz,Data3_R2.fq.gz -r HBV.fasta -l D1,D2,D3 -s pos -f bam -p bwa -t 5 -o ./timmedSep

=back

=head1 AUTHOR

Mingjie Dr.Wang I<huzai@sjtu.edu.cn>

=head1 COPYRIGHT

Copyright (C) 2017, Mingjie Wang. All rights reserved.




