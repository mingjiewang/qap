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
printcol ("IGV","green");
print "\n";

printcol("########################CAUTION########################",'white');
printcol("##                                                   ##",'white');
printcol("##  IGV must be run under a GUI allowed environment  ##",'white');
printcol("##                                                   ##",'white');
printcol("#######################################################",'white');
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
my $ref;
my $gui;
my $threads;

my $DateNow = `date +"%Y%m%d_%Hh%Mm%Ss"`;
chomp $DateNow;

GetOptions(
'i|inputFile|=s'     => \$inputfile,
'r|refSeq|=s'        => \$ref,
'o|outputDir|=s'     => \$outputDir,
't|threads|=s'       => \$threads,
'h|help|'            => \$help,
'g|gui|'             => \$gui
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
	$outputDir = File::Spec -> catfile($wk_dir,"qap_Results_for_IGV_$DateNow");
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
	if(-e "/proc/cpuinfo"){
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
	InfoError("Input MSA fasta file MUST be provided.");
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

##core program starts here
#first check jgv program
my $igvProgram = File::Spec -> catfile($mainBin,'bin','3rdPartyTools','igv','igv.jar');
if (not existFile($igvProgram)){
	InfoError("The main program of IGV is missing. Please check.");
	exit(0);
}
#check samtools
my $DEBUG_MODE = 1;
my $samtoolsProgram = File::Spec -> catfile($RealBin,'3rdPartyTools','samtools','samtools');
if(CheckProgram($samtoolsProgram, __FILE__, __LINE__, $DEBUG_MODE)){
	#keep running
}else{
	InfoError("The program $samtoolsProgram does NOT exist. Exiting...");
	exit;
}

#first sort bam/sam file
copy($inputfile, $outputDir);
$inputfile = File::Spec -> catfile($outputDir,basename($inputfile));
 
my $inputfilename;
if(isSamFile($inputfile)){
	$inputfilename = removeSamSuffix(basename($inputfile),1);
	
	#my ($samtools_excu, $samfile, $sortMethod, $clear_sam, $threads) = @_;
	sam2SortedAndIndexedBam($samtoolsProgram,$inputfile,'Pos',1,$threads);
	
	#update $inputfile
	$inputfile =~ s/\.sam$/.PosSorted.bam/;
	
	
}elsif(isBamFile($inputfile)){
	$inputfilename = removeBamSuffix(basename($inputfile),1);
	
	#my ($samtools_excu, $bamfile, $sortMethod, $clear, $threads) = @_;
	sortAndIndexBam($samtoolsProgram,$inputfile,'Pos',1,$threads);
	
	#update $inputfile 
	$inputfile =~ s/\.bam$/.PosSorted.bam/;
	
}

#index genome file
my $indexfile = $ref . ".fai";
if (existFile($indexfile)){
	#nothing
}else{
	my $cmd = "$samtoolsProgram faidx $ref";
	runcmd($cmd);
}

##run igv
#write batch file
my $batchfile = File::Spec -> catfile($outputDir,$inputfilename . ".batch");
my $outpngname = $inputfilename . ".IGV.png";
my $outsvgname = $inputfilename . ".IGV.svg";

open B,">$batchfile" or die "Can NOT output to $batchfile:$!";
print B "new\n";
print B "genome $ref\n";
print B "load $inputfile\n";
print B "snapshotDirectory $outputDir\n";
print B "goto all\n";
print B "sort position\n";
print B "collapse\n";
print B "snapshot $outpngname\n";
print B "snapshot $outsvgname\n";
if (defined $gui){
	#nothing
}else{
	print B "exit\n";
}
close B;

#run batch file
my $cmd = "java -Xmx6g -jar $igvProgram -g $ref -b $batchfile";
runcmd($cmd);

#move to tmp dir
my $tmpdir = File::Spec -> catfile($outputDir,'tmp');
makedir($tmpdir);

move($inputfile,$tmpdir) if existFile($inputfile);
move($inputfile . ".bai",$tmpdir) if existFile($inputfile . ".bai");
move($batchfile,$tmpdir) if existFile($batchfile);

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
           



qap IGV [options]

Use --help to see more information.

qap is still in development. If you have encounted any problem in usage, please feel no hesitation to cotact us.

=head1 DESCRIPTION

This script implements a function to visualize the mapping details of BAM/SAM file. The script has B<several> mandatory options that MUST appear last. 

=head1 OPTIONS

=over 5

=item --inputFile,-i F<FILE> [Required]

Path to input bam or sam file which is required.

=item --refSeq,-r F<File> [Required]

Path to the reference fasta file. 

=item --outputDir,-o F<FILE> [Optional]

Path of the directory to storage result files. If NOT provided, the program will generate a folder automatically.

=item --gui,-g [Optional]

Whether to keep the GUI of IGV or not. If -g is specified, the program will leave IGV running after generating figures.

=item --threads,-t F<INTEGER> [Optional]

Number of threads this program will use when computing. A positive integer should be provided. The default value is 1.

=item --help,-h

Display this detailed help information.

=back

=head1 EXAMPLE

=over 5

qap IGV -i test.bam -r HBV.fasta -o ./igv

=back

=head1 AUTHOR

Mingjie Dr.Wang I<huzai@sjtu.edu.cn>

=head1 COPYRIGHT

Copyright (C) 2017, Mingjie Wang. All rights reserved.


