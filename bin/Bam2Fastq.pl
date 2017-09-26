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

####---------------------------####
####The program begins here
####---------------------------####

##Show welcome
print "You are now running subprogram: ";
printcol ("Bam2Fastq","green");
print "\n";

##get workding directory
my $wk_dir = getcwd;
my $mainBin;
if ($RealBin =~ /(.*)\/bin/){
	$mainBin = $1;
}

####define command line arguments
my $help;
my $inputFile;
my $outFq1;
my $outFq2;
my $threads;

my $DateNow = `date +"%Y%m%d_%Hh%Mm%Ss"`;
chomp $DateNow;

GetOptions(
'i|inputFile|=s'      => \$inputFile,
'1|outputFastq1|=s'   => \$outFq1,
'2|outputFastq2|=s'   => \$outFq2,
't|threads|=s'        => \$threads 
);


##check command line arguments
if (defined $help){
	pod2usage(-verbose=>2,-exitval=>1);
}

if(defined $inputFile){
	if(existFile($inputFile)){
		if(isBamFile($inputFile) or isSamFile($inputFile)){
			$inputFile = abs_path($inputFile);
		}else{
			InfoError("Input file $inputFile is incorrect formatted. Whether a BAM or SAM file is required.");
			pod2usage(-verbose=>1,-exitval=>1);
			exit(0);
		}
		
	}else{
		InfoError("Input file $inputFile does NOT exist.");
		exit;
	}
}else{
	InfoError("Input file $inputFile MUST be defined.");
	pod2usage(-verbose=>1,-exitval=>1);
	exit;
}

if(defined $outFq1){
	if($outFq1 !~ /^\//){
		$outFq1 = File::Spec -> catfile($wk_dir, $outFq1);
	}
	my $outdir = dirname($outFq1);
	makedir($outdir);
	
	if(defined $outFq2){
		if($outFq2 !~ /^\//){
			$outFq2 = File::Spec -> catfile($wk_dir, $outFq2);
		}
		my $outdir2 = dirname($outFq2);
		makedir($outdir2);
		if($outdir2 ne $outdir){
			InfoWarn("Output fastq1 file and fastq2 file are output to different directory.Please check.");
		}
	}else{
		InfoWarn("Only --outputFastq1 is defined. Is it a single-ended unpaired fastq file? Go on processing anyway...");
		$outFq2 = "null";
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


##start convert
my $DEBUG_MODE = 1;
#check picard
my $picard_excu = File::Spec -> catfile($RealBin,'3rdPartyTools','caller','picard.jar');
if(existFile($picard_excu)){
	#keep running
}else{
	InfoError("The program $picard_excu does NOT exist. Exiting...");
	exit;
}

&bam2fq($picard_excu,$inputFile,$outFq1,$outFq2,$threads);

##sub programs start here
sub bam2fq{
	my $picard_excu = shift;
	my $inputFile = shift;
	my $outFq1 = shift;
	my $outFq2 = shift;
	my $threads = shift;
	
	my $unpaired_fq = File::Spec -> catfile(dirname($outFq1), removeFastqSuffix(basename($outFq1)) . "_Unpaired.fastq");
	if($outFq2 ne "null"){
		my $cmd = "java -jar $picard_excu SamToFastq INPUT=$inputFile FASTQ=$outFq1 SECOND_END_FASTQ=$outFq2 UNPAIRED_FASTQ=$unpaired_fq";
		runcmd($cmd);
		
		if(not -s $unpaired_fq){
			unlink $unpaired_fq;
		}
	}else{
		my $cmd = "java -jar $picard_excu SamToFastq INPUT=$inputFile FASTQ=$outFq1";
		runcmd($cmd);
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
           



qap Bam2Fastq [options]

Use --help to see more information.

qap is still in development. If you have encounted any problem in usage, please feel no hesitation to cotact us.

=head1 DESCRIPTION

This script implements a function to convert bam/sam file to fastq files. The script has B<several> mandatory options that MUST appear last. 

=head1 OPTIONS

=over 5

=item --inputFile,-i F<FILE> [Required]

Path to input BAM/SAM file which is required. 

=item --outputFastq1,-1 F<FILE> [Required]

Path to output fastq1 file or fastq file for single-ended sequencing data.

=item --outputFastq2,-2 F<FILE> [Optional]

Path to output fastq2 file. Must be specified in pair-ended sequencing data.

=item --threads,-t F<INTEGER> [Optional]

Number of threads this program will use when computing. A positive integer should be provided. The default value is 1.

=item --help,-h

Display this detailed help information.

=back

=head1 EXAMPLE

=over 5

qap Bam2Fastq -i test.bam -1 test_R1.fq -2 test_R2.fq -t 10

=back

=head1 AUTHOR

Mingjie Dr.Wang I<huzai@sjtu.edu.cn>

=head1 COPYRIGHT

Copyright (C) 2017, Mingjie Wang. All rights reserved.


