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

####---------------------------####
####The program begins here
####---------------------------####

##Show welcome
print "You are now running subprogram: ";
printcol ("BarcodeSplitter","green");
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
my $outputFile;
my $threads;

my $DateNow = `date +"%Y%m%d_%Hh%Mm%Ss"`;
chomp $DateNow;

GetOptions(
'i|inputFile|=s'      => \$inputFile,
'o|outputFile|=s'     => \$outputFile,
't|threads|=s'        => \$threads 
);


##check command line arguments
if (defined $help){
	pod2usage(-verbose=>2,-exitval=>1);
}

my $tmpdir = File::Spec -> catfile("/home/webusers/tmp/","qap_fq2fa_" . time());
makedir($tmpdir);

if(defined $inputFile){
	if(existFile($inputFile)){
		my $newfile = File::Spec -> catfile($tmpdir,basename($inputFile). ".fastq");
		copy($inputFile,$newfile);
		$inputFile = $newfile;
	}else{
		InfoError("Input file $inputFile does NOT exist.");
		exit;
	}
}else{
	InfoError("Input file $inputFile MUST be defined.");
	pod2usage(-verbose=>1,-exitval=>1);
	exit;
}

if(not defined($outputFile)){
	InfoError("Output file MUST be defined.");
	pod2usage(-verbose=>1,-exitval=>1);
	exit;
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
##check fq2fa
my $fq2fa_excu = File::Spec -> catfile($RealBin,'3rdPartyTools','fastx_toolkit','fq2fa');
if(CheckProgram($fq2fa_excu, __FILE__, __LINE__, $DEBUG_MODE)){
	#keep running
}else{
	InfoError("The program $fq2fa_excu does NOT exist. Exiting...");
	exit;
}

my $outfa = $inputFile =~ s/\.fastq$/.fasta/r;
my $cmd = "$fq2fa_excu -in $inputFile -out $outfa";
runcmd($cmd);

copy($outfa,$outputFile);

##run success
Info("Program completed!",'green');

##remove output dir
system("rm -rf $tmpdir");

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
           



qap Fastq2Fasta [options]

Use --help to see more information.

qap is still in development. If you have encounted any problem in usage, please feel no hesitation to cotact us.

=head1 DESCRIPTION

This script implements a function to convert fastq file to fasta file. The script has B<several> mandatory options that MUST appear last. 

=head1 OPTIONS

=over 5

=item --inputFile,-i F<FILE> [Required]

Path to input fastq file which is required. 

=item --outputFile,-o F<FILE> [Required]

Path to output fasta file.

=item --threads,-t F<INTEGER> [Optional]

Number of threads this program will use when computing. A positive integer should be provided. The default value is 1.

=item --help,-h

Display this detailed help information.

=back

=head1 EXAMPLE

=over 5

qap Bam2Sam -i test.fq -o test.fa -t 10

=back

=head1 AUTHOR

Mingjie Dr.Wang I<huzai@sjtu.edu.cn>

=head1 COPYRIGHT

Copyright (C) 2017, Mingjie Wang. All rights reserved.

