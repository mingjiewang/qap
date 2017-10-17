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
##  This file is a subprogram of qap suite.                                    ##
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
use List::Util qw/max min/;

####Use modules in this program####
use General;
use NT;

####Flush cache
$| = 1;

####---------------------------####
####The program begins here
####---------------------------####

##Show welcome
print "You are now running subprogram: ";
printcol ("DoimantStrain","green");
print "\n";

## check threads available or not
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
my $inputDir;
my $outputDir;
my $suffix;
my $graphic;
my $threads;
my $cutoff;

my $DateNow = `date +"%Y%m%d_%Hh%Mm%Ss"`;
chomp $DateNow;

GetOptions(
'i|inputDir|=s'     => \$inputDir,
's|suffix|=s'       => \$suffix,
'o|outputDir|=s'    => \$outputDir,
't|threads|=s'      => \$threads,
'g|graphic|=s'      => \$graphic,
'h|help|'           => \$help,
'c|cutoff|=s'       => \$cutoff
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
		}
	}
}else{
	$outputDir = File::Spec -> catfile($wk_dir,"qap_Results_for_DoimantStrain_$DateNow");
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

if(defined $cutoff){
	if(CheckDouble($cutoff) and $cutoff >=0 and $cutoff<= 1){
		#Nothing
	}else{
		InfoError("Cutoff value should be a positive decimal number. Provide again with --cutoff/-c.");
		pod2usage(-verbose=>0,-exitval=>1);
		exit;
	}
}else{
	$cutoff = 0.001; # use 0.1% as default
}


if(defined $inputDir){
	$inputDir =~ s/\/$//;
	$inputDir = abs_path($inputDir) . "/";
	if (not -e $inputDir){
		InfoError("Input directory $inputDir does NOT exist! Please check again.");
		exit;
	}
}else{
	InfoError("Input directory MUST be specified with -i/--inputDir\n");
	pod2usage(-verbose=>0,-exitval=>1);
	exit;
}

my $tmpDir = File::Spec -> catfile($outputDir,'tmp');
makedir($tmpDir);

my $numberOfFiles = 0;
my @inputfiles;
if(defined $suffix){
	@inputfiles = glob ("${inputDir}/*.${suffix}");
	
	$numberOfFiles = scalar(@inputfiles);
	
	if ($numberOfFiles == 0){
		InfoError("There are NOT any files in $inputDir with suffix \'.${suffix}\'. Please check again.");
		exit;
	}
	
	Info("Find $numberOfFiles files.");
	my $i = 1;
	for my $f (@inputfiles){
		printf "[%02d] $f\n",$i;
		my $len = checkSeqLen($f, $tmpDir);
		$i++;
	}
	
}else{
	InfoWarn("The suffix is not provided. The program will try to read in every file in $inputDir");
	
	@inputfiles = glob ("${inputDir}/*.*");
	
	$numberOfFiles = scalar(@inputfiles);
	
	if ($numberOfFiles == 0){
		InfoError("There are NOT any files in $inputDir with suffix \'.${suffix}\'. Please check again.");
		exit;
	}
	
	Info("Find $numberOfFiles files.");
	my $i = 1;
	for my $f (@inputfiles){
		printf "[%02d] $f\n",$i;
		my $len = checkSeqLen($f, $tmpDir);
		$i++;
	}
	
}

sleep(1);

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

my $pieChart = 1;
if(defined $graphic){
	if($graphic =~ /^y/i){
		$pieChart = 1;
	}elsif($graphic =~ /^n/i){
		$pieChart = 0;		
	}else{
		InfoError("Argument error.Choose between --graphic Y or --graphic N.");
		pod2usage(-verbose=>1,-exitval=>1);
		exit;
	}
	
}else{
	InfoWarn("--graphic is NOT defined. Using --graphic Y as default.");
	$pieChart = 1;
}

##core program starts here
my @outfile;
my @tmpdir;
my @pieplot;
my @plotdir;
my @cutoff;
my $plotdir = File::Spec -> catfile($outputDir, 'plots');
makedir($plotdir);
my @rscript;
my $rscript = File::Spec -> catfile($mainBin,'bin','Rscripts','PieChart.R');
if($pieChart and !existFile($rscript)){
	InfoError("R script used to draw graphs is missing. Please check.");
	exit;
}

for my $f (@inputfiles){
	my $outfile = File::Spec -> catfile($outputDir, removeFastaSuffix(basename($f)) . "_dominantStrain.fasta");
	push @outfile,$outfile;
	push @tmpdir,$tmpDir;
	push @pieplot,$pieChart;
	push @plotdir,$plotdir;
	push @rscript,$rscript;
	push @cutoff,$cutoff;
	
	#getDominantStrain($f,$outfile,$tmpDir,$pieChart,$plotdir,$rscript,$cutoff);
}
runMultipleThreadsWith7Args(\&getDominantStrain, \@inputfiles, \@outfile, \@tmpdir, \@pieplot, \@plotdir, \@rscript, \@cutoff,$threads);

##run success
print("\n\n");
Info("Program completed!",'green');

##sub program starts here
sub getDominantStrain {
	my $file = shift;
	my $outfile = shift;
	my $tmpdir = shift;
	my $pieplot = shift;
	my $plotdir = shift;
	my $rscript = shift;
	my $cutoff = shift;
	
	Info("Start to idenfity dominant strain for $file.");
	
	#convert fasta file to two line format
	my $file2Name = removeFastaSuffix(basename($file)) . ".2line.fasta";
	my $file2 = File::Spec -> catfile($tmpdir,$file2Name);
	formatFastaToTwoLineMode($file,$file2);
	
	open T,$file2 or die "Can not open file $file2:$!\n";
	my %hash;
	my $totalSeq;
	while(my $line1 = <T>){
		chomp $line1;
		chomp (my $line2 = <T>);
		
		$hash{$line2}++;
		$totalSeq++;
	}
	
	open RES,">$outfile" or die "Can not output to $outfile:$!";
	my @sortkeys = sort {$hash{$b} <=> $hash{$a}} keys %hash;
	my $ratio = $hash{$sortkeys[0]} / $totalSeq;
	print RES ">Dominant Strain[Ratio:$ratio $hash{$sortkeys[0]}\/$totalSeq]\n$sortkeys[0]\n";
	close RES;
	
	if($pieplot){
		Info("Drawing pie chart for $file.");
		my $rinput = File::Spec -> catfile($tmpdir, basename($outfile) . ".rInput");
		open RINPUT,">$rinput" or die "Can not output to $rinput:$!";
		print RINPUT "Strain\tCount\n";
		for my $k(keys %hash){
			print RINPUT "$k\t$hash{$k}\n";
		}
		close RINPUT;
		
		## run rscript
		my $outplot = File::Spec -> catfile($plotdir,removeFastaSuffix(basename($outfile)));
		my $cmd = "Rscript $rscript --cutoff $cutoff --inputFile $rinput --outputFile $outplot";
		runcmd($cmd);
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
           



qap DoimantStrain [options]

Use --help to see more information.

qap is still in development. If you have encounted any problem in usage, please feel no hesitation to cotact us.

=head1 DESCRIPTION

This script implements a function for cutting fasta sequences with base intervals in batch. The script has B<several> mandatory options that MUST appear last. 

=head1 OPTIONS

=over 5

=item --inputDir,-i F<FILE> [Required]

Path to directory contaning all the files to be read in.

=item --suffix,-s F<STRING> [Optional]

Suffix of the files to be read in. If suffix is not provided, all the files in input directory will be read in.

=item --cutoff,-c F<FLOAT> [Optional]

Cutoff value of strain ratio. Strains with ratio larget than the cutoff vlaue will be shown in the pie chart. The default value is 0.001.

=item --outputDir,-o F<FILE> [Optional]

Path of the directory to storage result files. If NOT provided, the program will generate a folder automatically.

=item --graphic,-g [Optional] 

Whether draw graphs or not. Choose between 'Y'(Yes) and 'N'(No). -g Y is specified as default. 

=item --threads,-t F<INTEGER> [Optional]

Number of threads this program will use when computing. A positive integer should be provided. The default value is 1.

=item --help,-h

Display this detailed help information.

=back

=head1 EXAMPLE

=over 5

qap DoimantStrain -i ./seq -t 10 -g -o ./result

=back

=head1 AUTHOR

Mingjie Dr.Wang I<huzai@sjtu.edu.cn>

=head1 COPYRIGHT

Copyright (C) 2017, Mingjie Wang. All rights reserved.
