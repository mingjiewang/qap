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
printcol ("RawDataFiltraion","green");
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
my $input;
my $outputDir;
my $threads;
my $graphic;
my $DateNow = `date +"%Y%m%d_%Hh%Mm%Ss"`;
chomp $DateNow;

GetOptions(
'i|inputFile|=s'   => \$input,
'o|outputDir|=s'   => \$outputDir,
't|threads|=s'     => \$threads,
'h|help|'          => \$help,
'g|graphic|'       => \$graphic
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
	$outputDir = File::Spec -> catfile($wk_dir,"qap_Results_for_PickRobustOTU_$DateNow");
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


my @inputfiles;
if (defined $input){
	my @tmp = split ",",$input;
	for my $tmp (@tmp){
		$tmp =~ s/\s//g;
		if (not existFile($tmp)){
			InfoError("Input fastq file $tmp does NOT exist!",'red');
			exit;
		}else{
			my $pathtmp = abs_path($tmp);
			push @inputfiles,$pathtmp;
		}
	}
}else{
	InfoError("Input fastq files MUST be provided!",'red');
	pod2usage(-verbose=>1,-exitval=>1);
	exit;
}
my $numberOfFiles = scalar(@inputfiles);

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


##start to run fastqc
my $DEBUG_MODE = 1;
my $fastqc_excu = File::Spec -> catfile($RealBin,'3rdPartyTools','fastqc','fastqc');
if(CheckProgram($fastqc_excu, __FILE__, __LINE__, $DEBUG_MODE)){
	#keep running
}else{
	InfoError("The program $fastqc_excu does NOT exist. Exiting...");
	exit;
}

my $fastqc_res_dir = File::Spec -> catfile($outputDir, 'fastqcResults');
makedir($fastqc_res_dir);
	
if($threads > 1){
	my @outdir;
	my @fastqc_excu;
	my @threadsForEach;
	
	my $threadsForEachSample = int($threads / $numberOfFiles);
	
	for my $f (@inputfiles){
		push @outdir,$fastqc_res_dir;
		push @fastqc_excu,$fastqc_excu;
		
		if ($threadsForEachSample > 1){
			push @threadsForEach,$threadsForEachSample;
		}else{
			push @threadsForEach,1
		}
	}
	
	Info("Running Fastqc with multiple threads.");
	runMultipleThreadsWith4Args(\&qc, \@inputfiles, \@outdir, \@fastqc_excu, \@threadsForEach, $threads);
	
}else{
	Info("Running Fastqc with single thread.");
	for my $f (@inputfiles){
		&qc($f, $fastqc_res_dir, $fastqc_excu, 1);
	}
}


##merge fastqc results to a table for easy view
Info("Merging Fastqc results into single file.\n");
my $outfile = File::Spec -> catfile($outputDir, "QualityControlResult.tsv");
open OUT,">$outfile" or die "Can NOT output to $outfile:$!";
print OUT "File.Name\tBasic.Statistics\tPer.base.sequence.quality\tPer.tile.sequence.quality\tPer.sequence.quality.scores\tPer.base.sequence.content\tPer.sequence.GC.content\tPer.base.N.content\tSequence.Length.Distribution\tSequence.Duplication.Levels\tOverrepresented.sequences\tAdapter.Content\tKmer.Content\n";

for my $dir (glob "$fastqc_res_dir/*_fastqc"){
	if(CheckFolder($dir)){
		#nothing
	}else{
		InfoError("Something is wrong with Fastqc result folder $dir. Exiting...");
		exit;
	}
	
	my @result;
	my $checkfile;
	my $resfile = File::Spec -> catfile($dir, 'summary.txt');
	
	open R,"$resfile" or die "Can NOT open $resfile:$!\n";
	while(<R>){
		chomp;
		my @tmp = split "\t",$_;
		my $checkres = $tmp[0];
		push @result,$checkres;
		$checkfile = $tmp[2];
	}
	
	my $result = join "\t",@result;
	print OUT "$checkfile\t$result\n";
}

close OUT;

##draw plots
my $rscript = File::Spec -> catfile($RealBin,'Rscripts','QC2Heatmap.R');
if(not existFile($rscript)){
	InfoError("The R script $rscript is missing. Please check. Exiting...");
	exit;
}

if(defined $graphic){
	Info("Drawing graphs for QC results visualization.");
	
	my $fig = File::Spec -> catfile($outputDir,"QualityControlResult");
	my $cmd = "Rscript $rscript -i $outfile -o $fig";
	runcmd($cmd);
}


##run success
sleep(1);
Info("Program completed!",'green');


##sub program goes here
sub qc {
	my $file = shift;
	my $outdir = shift;
	my $fastqc_excu = shift;
	my $threads = shift;
	
	my $cmd = "$fastqc_excu -o $outdir --extract -f fastq -q -t $threads $file";
	runcmd($cmd);
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
           



qap RawDataQC [options]

Use --help to see more information.

qap is still in development. If you have encounted any problem in usage, please feel no hesitation to cotact us.

=head1 DESCRIPTION

This script implements a function for raw sequencin data quality control. The script has B<several> mandatory options that MUST appear last. 

=head1 OPTIONS

=over 5

=item --inputFile,-i F<FILE> [Required]

Path to next generation sequencing raw data. REQUIRED for both single-end or paired-end reads. Both compressed files and uncompressed files are allowed. Multiple files are allowed, which should be seperated by comma, e.g. --inputFile test1.fq.gz,test2.fq.gz.

=item --outputDir,-o F<FILE> [Optional]

Path of the directory to storage result files. If NOT provided, the program will generate a folder automatically.

=item --threads,-t F<INTEGER> [Optional]

Number of threads this program will use when computing. A positive integer should be provided. The default value is 1.

=item --graphic,-g [Optional]

Whether to show graphs about merged quality control results or not.

=item --help,-h

Display this detailed help information.

=back

=head1 EXAMPLE

=over 5

qap RawDataQC -i test1.fq,test2.fq -t 5 -o ./qc

=back

=head1 AUTHOR

Mingjie Dr.Wang I<huzai@sjtu.edu.cn>

=head1 COPYRIGHT

Copyright (C) 2017, Mingjie Wang. All rights reserved.


