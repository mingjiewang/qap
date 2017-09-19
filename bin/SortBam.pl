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
my $sortMethod;

my $DateNow = `date +"%Y%m%d_%Hh%Mm%Ss"`;
chomp $DateNow;

GetOptions(
'i|inputFile|=s'      => \$inputFile,
'o|outputFile|=s'     => \$outputFile,
't|threads|=s'        => \$threads,
's|sortMethod|=s'     => \$sortMethod
);


##check command line arguments
if (defined $help){
	pod2usage(-verbose=>2,-exitval=>1);
}

my $tmpdir = File::Spec -> catfile("/home/webusers/tmp/","qap_sortbam_" . time());
makedir($tmpdir);

if(defined $inputFile){
	if(existFile($inputFile)){
		my $newfile = File::Spec -> catfile($tmpdir,basename($inputFile). ".bam");
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
	my $threads_max = `grep 'processor' /proc/cpuinfo | sort -u | wc -l`;
	chomp $threads_max;

	if ($check_threads_positive && $threads <= $threads_max){
		#threads provided by user is ok, doing nothing
	}else{
		InfoError("Threads number wrong!",'red');
		InfoError("Please provide a threads number between 0 - $threads_max that this server could support.");
		
		pod2usage(-verbose=>1,-exitval=>1);
		exit;
	}
}else{
	$threads = 1;#if -t not provided, default is NOT use theads;
}
$threads = 1;

if(defined $sortMethod){
	if($sortMethod =~ /pos/i){
		$sortMethod = 'pos';
	}elsif($sortMethod =~ /name/i){
		$sortMethod = 'name';
	}else{
		InfoError("Sort method should be one of \'--sortMethod position\' or \'--sortMethod name\'.");
		pod2usage(-verbose=>1,-exitval=>1);
		exit;
	}
}else{
	InfoError("BAM sort method must be defined.");
	pod2usage(-verbose=>1,-exitval=>1);
	exit;
}

##start convert
my $DEBUG_MODE = 1;
#check samtools
my $samtools_excu = File::Spec -> catfile($RealBin,'3rdPartyTools','samtools','samtools');
if(CheckProgram($samtools_excu, __FILE__, __LINE__, $DEBUG_MODE)){
	#keep running
}else{
	InfoError("The program $samtools_excu does NOT exist. Exiting...");
	exit;
}
sortAndIndexBam($samtools_excu,$inputFile,$sortMethod,1,$threads);
if (uc($sortMethod) eq 'POS'){
	my $bamfile_sort = $inputFile =~ s/\.bam$/.PosSorted.bam/r;
	copy($bamfile_sort,$outputFile);
}else{
	my $bamfile_sort = $inputFile =~ s/\.bam$/.NameSorted.bam/r;
	copy($bamfile_sort,$outputFile);
}

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
           



qap SortBam [options]

Use --help to see more information.

qap is still in development. If you have encounted any problem in usage, please feel no hesitation to cotact us.

=head1 DESCRIPTION

This script implements a function to sort bam file. The script has B<several> mandatory options that MUST appear last. 

=head1 OPTIONS

=over 5

=item --inputFile,-i F<FILE> [Required]

Path to input BAM file which is required. 

=item --outputFile,-o F<FILE> [Required]

Path to output sorted BAM file.

=item --sortMethod,-s F<STRING> [Required]

Sort method should be selected from 'position' or 'name'.

=item --threads,-t F<INTEGER> [Optional]

Number of threads this program will use when computing. A positive integer should be provided. The default value is 1.

=item --help,-h

Display this detailed help information.

=back

=head1 EXAMPLE

=over 5

qap Bam2Sam -i test.bam -s position -o test.bam -t 10

=back

=head1 AUTHOR

Mingjie Dr.Wang I<huzai@sjtu.edu.cn>

=head1 COPYRIGHT

Copyright (C) 2017, Mingjie Wang. All rights reserved.


