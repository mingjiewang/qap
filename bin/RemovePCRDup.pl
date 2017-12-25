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

####Use modules in this program####
use General;
use Mapper;
use Caller;

####---------------------------####
####The program begins here
####---------------------------####

##Show welcome
print "You are now running subprogram: ";
printcol ("RemovePCRDup","green");
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
my $programFlag;
my $outformat;


my $DateNow = `date +"%Y%m%d_%Hh%Mm%Ss"`;
chomp $DateNow;

GetOptions(
'i|inputFile|=s'     => \$inputfile,
'o|outputDir|=s'     => \$outputDir,
'h|help|'            => \$help,
't|threads|=s'       => \$threads,
'p|program|=s'       => \$programFlag,
'f|outFormat|=s'     => \$outformat
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
	}
}else{
	$outputDir = File::Spec -> catfile($wk_dir,"qap_Results_for_RemovePCRDup_$DateNow");
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

my @inputfiles;
if(defined $inputfile){
	my @tmp = split ',',$inputfile;
	for my $tmp(@tmp){
		if(existFile($tmp)){
			my $pathtmp = abs_path($tmp);
			push @inputfiles,$pathtmp;
		}else{
			InfoError("Input file $tmp doesn't exist. Please check.");
			exit(0);
		}
	}
	
}else{
	InfoError("Input file MUST be provided.");
	pod2usage(-verbose=>2,-exitval=>1);
	exit;
}

my @program;
if(defined $programFlag){
	if ($programFlag =~ /^s/i){
		push @program,'samtools';
		
		if($programFlag =~ /sp/){
			push @program,'picard';
		}
	}
	if ($programFlag =~ /^p/i){
		push @program,'picard';
		
		if($programFlag =~ /ps/i){
			push @program,'samtools';
		}
	}
	if ($programFlag !~ /[ps]/i){
		InfoError("Only two values are allowed. \'s\' for \'samtools\', or \'p\' for \'picard\'.");
		pod2usage(-verbose=>2,-exitval=>1);
		exit;
	}
}else{
	InfoWarn{"The program used for PCR replicates removal is NOT defined. Will use -p \'s\' as default."};
	push @program,'samtools';
}

my $format;
if(defined $outformat){
	if (uc($outformat) eq 'BAM' or $outformat =~ /b/i){
		$format = 'bam';
	}elsif(uc($outformat) eq 'SAM' or $outformat =~ /s/i){
		$format = 'sam';
	}else{
		InfoError("The output format can ONLY be \'sam\' or \'bam\' format.");
		pod2usage(-verbose=>2,-exitval=>1);
		exit;
	}
}else{
	InfoWarn("The output format is NOT defined. Will use \'-f bam\' as default.");
	$format = 'bam';
}

##main program starts here
#check samtools
my $samtools_excu = File::Spec -> catfile($RealBin,'3rdPartyTools','samtools','samtools');
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

#preprocess input files
my $numberOfFiles = scalar(@inputfiles);
Info("Checking input $numberOfFiles files.");

my @newinputfiles;
for my $f (@inputfiles){
	copy($f,$outputDir);
	my $newf = File::Spec -> catfile($outputDir, basename($f)); 
	if (isSamFile($newf)){
		sam2SortedAndIndexedBam($samtools_excu, $newf, 'POS', 1, $threads);
		my $sortbam = $newf =~ s/\.sam$/.PosSorted.bam/r;
		push @newinputfiles,$sortbam;
	}elsif(isBamFile($newf)){
		sortAndIndexBam($samtools_excu, $newf, 'POS', 1, $threads);
		my $sortbam = $newf =~ s/\.bam$/.PosSorted.bam/r;
		push @newinputfiles,$sortbam;
	}else{
		InfoError("The input file $f is not a valid BAM/SAM file.");
		exit(0);
	}
}

#sort input bam files based on file size
my %tmp;
for my $f (@newinputfiles){
	my $filesize = getFileSize($f);
	$tmp{$f} = $filesize;
}
my @newinputfilessort;
for my $f (sort {$tmp{$a} <=> $tmp{$b}} keys %tmp){
	push @newinputfilessort,$f;
}


#core program
if ($threads > 1){
	my @programtouse;
	my @outdir;
	my @samtools_excu;
	my @picard_excu;
	my @format;
	for my $f (@newinputfilessort){
		my $programexp = join '-',@program;
		push @programtouse,$programexp;
		push @outdir,$outputDir;
		push @samtools_excu,$samtools_excu;
		push @picard_excu,$picard_excu;
		push @format,$format;
	}
	
	runMultipleThreadsWith6Args(\&rmdup, \@programtouse, \@newinputfilessort, \@outdir, \@samtools_excu, \@picard_excu, \@format, $threads);
}else{
	for my $f (@newinputfilessort){
		my $programexp = join '-',@program;
		
		&rmdup($programexp, $f, $outputDir, $samtools_excu, $picard_excu, $format);
	}
}

#rm copied raw data files
for my $f (@newinputfiles){
	my $findex = $f . ".bai";
	
	unlink $f if -e $f;
	unlink $findex if -e $findex;	
}


##run success
Info("Program completed!",'green');


##sub program starts here
sub rmdup {
	my $program = shift;
	my $inputfile = shift;
	my $outdir = shift;
	my $samtools_excu = shift;
	my $picard_excu = shift;
	my $format = shift;
	
	my $inputfilename = basename($inputfile);
	
	if($program eq 'samtools'){
		#run samtools
		my $outbam = File::Spec -> catfile($outdir, removeBamSuffix(basename($inputfile), 1) . ".samtoolsRmDup.bam");
		
		Info("Removing PCR duplicates for $inputfilename using samtools.");
		my $cmd = "$samtools_excu rmdup $inputfile $outbam";
		runcmd($cmd);

		if($format eq 'sam'){
			bam2sam($samtools_excu, $outbam, 1);
			unlink($outbam);
		}
		
	}elsif($program eq 'picard'){
		#run picard
		my $outbam = File::Spec -> catfile($outdir, removeBamSuffix(basename($inputfile), 1) . ".picardRmDup.bam");
		my $outbamindex = $outbam =~ s/\.bam$/.bai/r;
		my $metricfile = File::Spec -> catfile($outdir, removeBamSuffix(basename($inputfile), 1) . ".metrics");
		
		Info("Removing PCR duplicates for $inputfilename using picard.");
		my $cmd = "java -Xmx8g -jar $picard_excu MarkDuplicates I=$inputfile O=$outbam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true M=$metricfile";
		runcmd($cmd);
		
		if($format eq 'sam'){
			bam2sam($samtools_excu, $outbam, 1);
			unlink($outbam);
		}
		
		unlink $metricfile;
		unlink $outbamindex;
		
	}elsif($program eq 'samtools-picard'){
		#run samtools
		my $outbam = File::Spec -> catfile($outdir, removeBamSuffix(basename($inputfile), 1) . ".samtoolsRmDup.bam");
		
		Info("Removing PCR duplicates for $inputfilename using samtools.");
		my $cmd = "$samtools_excu rmdup $inputfile $outbam";
		runcmd($cmd);
		
		#run picard
		my $outbam2 = File::Spec -> catfile($outdir, removeBamSuffix(basename($inputfile), 1) . ".samtoolsPicardRmDup.bam");
		my $outbam2index = $outbam2 =~ s/\.bam$/.bai/r;
		my $metricfile = File::Spec -> catfile($outdir, removeBamSuffix(basename($inputfile), 1) . ".metrics");
		
		Info("Removing PCR duplicates for $inputfilename using picard.");
		$cmd = "java -Xmx8g -jar $picard_excu MarkDuplicates I=$outbam O=$outbam2 CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true M=$metricfile";
		runcmd($cmd);
		unlink $metricfile;
		unlink $outbam2index;
		
		if($format eq 'sam'){
			bam2sam($samtools_excu, $outbam, 1);
			bam2sam($samtools_excu, $outbam2, 1);
			unlink($outbam);
			unlink($outbam2);
		}
		
	}elsif($program eq 'picard-samtools'){
		#run picard
		my $outbam = File::Spec -> catfile($outdir, removeBamSuffix(basename($inputfile), 1) . ".picardRmDup.bam");
		my $outbamindex = $outbam =~ s/\.bam$/.bai/r;
		my $metricfile = File::Spec -> catfile($outdir, removeBamSuffix(basename($inputfile), 1) . ".metrics");
		
		Info("Removing PCR duplicates for $inputfilename using picard.");
		my $cmd = "java -Xmx8g -jar $picard_excu MarkDuplicates I=$inputfile O=$outbam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true M=$metricfile";
		runcmd($cmd);
		
		#run samtools
		my $outbam2 = File::Spec -> catfile($outdir, removeBamSuffix(basename($inputfile), 1) . ".picardSamtoolsRmDup.bam");
		
		Info("Removing PCR duplicates for $inputfilename using samtools.");
		$cmd = "$samtools_excu rmdup $outbam $outbam2";
		runcmd($cmd);
		
		unlink $metricfile;
		unlink $outbamindex;
		
		if($format eq 'sam'){
			bam2sam($samtools_excu, $outbam, 1);
			bam2sam($samtools_excu, $outbam2, 1);
			unlink($outbam);
			unlink($outbam2);
		}
		
	}else{
		InfoError("The program used to remove PCR duplicates should be one of \'Samtools\' or \'Picard\'.");
		exit(0);
	}
	
	
}



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
           



qap RemovePCRDup [options]

Use --help to see more information.

qap is still in development. If you have encounted any problem in usage, please feel no hesitation to cotact us.

=head1 DESCRIPTION

This script implements a function to remove PCR duplicates in BAM/SAM file. The script has B<several> mandatory options that MUST appear last. 

=head1 OPTIONS

=over 5

=item --inputFile,-i F<FILE> [Required]

Path to input BAM/SAM file which is required. You can also provide multiple files and seperate them by comma, e.g. -i test1.bam,test2.bam,test3.bam 

=item --outputDir,-o F<FILE> [Optional]

Path of the directory to storage result files. If NOT provided, the program will generate a folder automatically.

=item --program,-p F<STRING> [Optional]

Program used to remove PCR duplicates. Two values are available, 's' for 'samtools', and 'p' for 'picard'. You can also use 'sp' or 'ps' to use both softwares.

=item --outFormat,-f F<STRING> [Optional]

Format of output files. Two values are allowd: 'bam' or 'sam'. If --outFormat/-f is not specified, -f bam will be used as default.

=item --threads,-t F<INTEGER> [Optional]

Number of threads this program will use when computing. A positive integer should be provided. The default value is 1.

=item --help,-h

Display this detailed help information.

=back

=head1 EXAMPLE

=over 5

qap RemovePCRDup -i test1.bam,test2.bam -t 5 -p s -f sam -o ./rmdup

=back

=head1 AUTHOR

Mingjie Dr.Wang I<huzai@sjtu.edu.cn>

=head1 COPYRIGHT

Copyright (C) 2017, Mingjie Wang. All rights reserved.
