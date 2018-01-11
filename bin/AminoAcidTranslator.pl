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
printcol ("AminoAcidTranslator","green");
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
my $threads;

my $DateNow = `date +"%Y%m%d_%Hh%Mm%Ss"`;
chomp $DateNow;

GetOptions(
'i|inputDir|=s'     => \$inputDir,
's|suffix|=s'       => \$suffix,
'o|outputDir|=s'    => \$outputDir,
't|threads|=s'      => \$threads,
'h|help|'           => \$help
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
		}else{
			InfoError("Mkdir Failed! Folder $outputDir already exists!","red");
			InfoError("Please specify another output directory using option -o/--outputDir");
			pod2usage(-verbose=>0,-exitval=>1);
			exit;
		}
	}
}else{
	$outputDir = File::Spec -> catfile($wk_dir,"qap_Results_for_AminoAcidTranslator_$DateNow");
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
		$i++;
	}

}

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


##core program starts here
if($threads > 1){
	Info("Start translating using multipe threads. Please wait...");

	my @outfiles;
	for my $f (@inputfiles){
		my $outfile = File::Spec -> catfile($outputDir, removeFastaSuffix(basename($f)) . ".aa.fasta");
		push @outfiles,$outfile;
	}

	&runMultipleThreadsWith2Args(\&translateFileWithInfo, \@inputfiles, \@outfiles, $threads);
}else{
	Info("Start translating using single thread. Please wait...");

	# we have to sort the inputfiles by basename length, and translate them one by one.
	# Otherwise, the progress bar will hit a bug and leave a mark of the previous inputfile name.

	# sorting
	my %fastaNameLen;
	for my $f (@inputfiles){
		my $filename = basename($f);
		$fastaNameLen{$f} = length($filename);
	}

	my @inputfilesSortByNameLength = sort {$fastaNameLen{$a} <=> $fastaNameLen{$b}} keys %fastaNameLen;

	# translation
	my $i = 1;
	for my $f (@inputfilesSortByNameLength){
		my $outfile = File::Spec -> catfile($outputDir, removeFastaSuffix(basename($f)) . ".aa.fasta");

		my $seqID = basename($f);
		&ProcessBarForAATranslator($i, $numberOfFiles, $seqID);

		&translateFile($f, $outfile);

		$i++;
	}
}


##sub-program starts here
sub ProcessBarForAATranslator {
	local $| = 1;
	my $i = $_[0] || return 0;
	my $n = $_[1] || return 0;
	my $seqID = $_[2] || return 0;

	chomp (my $date = `date`);

    #output
	print "\rINFO    : Progress [ ".("=" x int(($i/$n)*50)).(" " x (50 - int(($i/$n)*50)))." ] ";
	printf("%4.2f%%",$i/$n*100);
	print " ($seqID)";

	local $| = 0;
}

sub translateFile {
	my $fastafile = shift;
	my $outfile = shift;

	# reformat input fasta file into 2 line format
	my $fastafile2 = File::Spec -> catfile($outputDir, removeFastaSuffix(basename($fastafile)) . ".refmt.fasta" );

	formatFastaToTwoLineMode($fastafile, $fastafile2);

	# output file handle
	open RES,">$outfile" or die "Can NOT output to file $outfile:$!\n";

	# read in file
	open T,$fastafile2 or die "Can NOT open fasta file $fastafile2:$!\n";
	while(my $line1 = <T>){
		chomp $line1;
		chomp (my $line2 = <T>);

		if ($line1 =~ /^>/){
			print RES "$line1\n";
		}else{
			InfoError("The input fasta file is NOT right formatted.");
			#exit;
			return 0; # Stop translating for current file and keep running the program
		}

		my $seqID = $line1 =~ s/^>//r;

		my $aa = translation($line2, $seqID);

		print RES "$aa\n";
	}
	close T;
	close RES;

	# remove refomat fasta file
	my $cmd = "rm -rf $fastafile2";
	system($cmd);
}

sub translateFileWithInfo {
	my $fastafile = shift;
	my $outfile = shift;

	# update status
	Info("Start to translate $fastafile.");

	# reformat input fasta file into 2 line format
	my $fastafile2 = File::Spec -> catfile($outputDir, removeFastaSuffix(basename($fastafile)) . ".refmt.fasta" );

	formatFastaToTwoLineMode($fastafile, $fastafile2);

	# output file handle
	open RES,">$outfile" or die "Can NOT output to file $outfile:$!\n";

	# read in file
	open T,$fastafile2 or die "Can NOT open fasta file $fastafile2:$!\n";
	while(my $line1 = <T>){
		chomp $line1;
		chomp (my $line2 = <T>);

		if ($line1 =~ /^>/){
			print RES "$line1\n";
		}else{
			InfoError("The input fasta file is NOT right formatted.");
			#exit;
			return 0; # Stop translating for current file and keep running the program
		}

		my $seqID = $line1 =~ s/^>//r;

		my $aa = translation($line2, $seqID);

		print RES "$aa\n";
	}
	close T;
	close RES;

	# remove refomat fasta file
	my $cmd = "rm -rf $fastafile2";
	system($cmd);
}



##run success
print("\n\n");
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




qap AminoAcidTranslator [options]

Use --help to see more information.

qap is still in development. If you have encounted any problem in usage, please feel no hesitation to cotact us.

=head1 DESCRIPTION

This script implements a function for translating nucleotide sequences into amino acid (protein) sequences in batch. The script has B<several> mandatory options that MUST appear last.

=head1 OPTIONS

=over 5

=item --inputDir,-i F<FILE> [Required]

Path to directory contaning all the files to be read in.

=item --suffix,-s F<STRING> [Optional]

Suffix of the files to be read in. If suffix is not provided, all the files in input directory will be read in.

=item --outputDir,-o F<FILE> [Optional]

Path of the directory to storage result files. If NOT provided, the program will generate a folder automatically.

=item --threads,-t F<INTEGER> [Optional]

Number of threads this program will use when computing. A positive integer should be provided. The default value is 1.

=item --help,-h

Display this detailed help information.

=back

=head1 EXAMPLE

=over 5

qap AminoAcidTranslator -i ./seq -t 10 -s fasta -o ./aa

=back

=head1 AUTHOR

Mingjie Dr.Wang I<huzai@sjtu.edu.cn>

=head1 COPYRIGHT

Copyright (C) 2017, Mingjie Wang. All rights reserved.
