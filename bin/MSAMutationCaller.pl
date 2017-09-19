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
use Caller;

####---------------------------####
####The program begins here
####---------------------------####

##Show welcome
print "You are now running subprogram: ";
printcol ("MSAMutationCaller","green");
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
my $inputfile;
my $outputDir;
my $ref;
my $cutoff;

my $DateNow = `date +"%Y%m%d_%Hh%Mm%Ss"`;
chomp $DateNow;

GetOptions(
'i|inputFile|=s'     => \$inputfile,
'r|refSeq|=s'     => \$ref,
'o|outputDir=s'      => \$outputDir,
'c|cutoff|=s'        => \$cutoff,
'h|help|'            => \$help
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
	}else{
		InfoError("Folder $outputDir already exists. Please specify another output directory using option -o/--outputDir");
		pod2usage(-verbose=>2,-exitval=>1);
		exit;
	}
}else{
	$outputDir = File::Spec -> catfile($wk_dir,"qap_Results_for_MSAMutationCaller_$DateNow");
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

if(defined $cutoff){
	if ($cutoff > 0 and $cutoff < 1){
		#nothing
	}else{
		InfoError("The cutoff value should be >= 0 and <= 1.");
		exit(0);
	}
}else{
	InfoWarn("The cufoff value of variant frequency is not provided. Using 0.01 as default.");
	$cutoff = 0.01;
}


##the core program starts here
#first check the seq len of MSA file
my $seqlen = checkSeqLen($inputfile,$outputDir);
#conver to r input
my $rinputfile = convert2Rinput($inputfile,$outputDir);

#locate the amplicon in the reference genome
my $DEBUG_MODE = 1;

#check blat
my $blatProgram = File::Spec -> catfile($RealBin,'3rdPartyTools','blat','blat');
if(CheckProgram($blatProgram, __FILE__, __LINE__, $DEBUG_MODE)){
	#keep running
}else{
	InfoError("The program $blatProgram does NOT exist. Exiting...");
	exit;
}

#check r script
my $rscript = File::Spec -> catfile($RealBin,'Rscripts','CallVariantsFromMSA.R');
if(existFile($rscript)){
	#keep running
}else{
	InfoError("The script $rscript does NOT exist. Exiting...");
	exit;
}

#get location from blat output
my ($psl,$start,$end) = blatPipeline($blatProgram,$inputfile,$ref,$outputDir);

#get total ref seq
my $ref2 = File::Spec -> catfile($outputDir,removeFastaSuffix(basename($ref)) . ".2line.fasta");
formatFastaToTwoLineMode($ref,$ref2);
system("sed -i \'s/\r//g\' $ref2");
open REF,"$ref2" or die "Can NOT open $ref2:$!";
chomp (my $refname = <REF>);
$refname =~ s/^>//;
chomp (my $totalref = <REF>);
close REF;

#get ref
my $refseq = substr $totalref, $start - 1, $end - $start + 1;
my $refseqfile = File::Spec -> catfile($outputDir, ${refname} . "_" .time() . ".txt");
open  REFSEQ,">$refseqfile" or die "Can NOT output to $refseqfile:$!";
print REFSEQ "$refseq\n";
close REFSEQ;

#check ref and seq length
if (length($refseq) == $seqlen){
	#nothing
}else{
	InfoError("The length of MSA sequences and the reference fasta file do NOT equal! Please check.");
	exit;
}

#get mutation sites
Info("Calling variants from $inputfile");
getMutationFromMSA($rscript,$refseqfile,$rinputfile,$outputDir,$start,$cutoff);


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
           



qap MSAMutationCaller [options]

Use --help to see more information.

qap is still in development. If you have encounted any problem in usage, please feel no hesitation to cotact us.

=head1 DESCRIPTION

This script implements a function to call variants from multiple sequence alignment (MSA). The script has B<several> mandatory options that MUST appear last. 

=head1 OPTIONS

=over 5

=item --inputFile,-i F<FILE> [Required]

Path to input MSA file (fasta format) which is required.

=item --refSeq,-r F<File> [Required]

Path to the reference fasta file. 

=item --outputDir,-o F<FILE> [Optional]

Path of the directory to storage result files. If NOT provided, the program will generate a folder automatically.

=item --cutoff,-c F<DOUBLE> [Optional]

The cutoff value of variant frequency. The default value is 0.01.

=item --help,-h

Display this detailed help information.

=back

=head1 EXAMPLE

=over 5

qap MSAMutationCaller -i test.fasta -r HBV.fasta -o ./mutation

=back

=head1 AUTHOR

Mingjie Dr.Wang I<huzai@sjtu.edu.cn>

=head1 COPYRIGHT

Copyright (C) 2017, Mingjie Wang. All rights reserved.

