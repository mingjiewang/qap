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
printcol ("Diversity","green");
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
my $inputDir;
my $outputDir;
my $suffix;
my $seqType;
my $threads;

my $DateNow = `date +"%Y%m%d_%Hh%Mm%Ss"`;
chomp $DateNow;

GetOptions(
'i|inputDir|=s'     => \$inputDir,
's|suffix|=s'       => \$suffix,
'o|outputDir|=s'    => \$outputDir,
'p|seqType|=s'      => \$seqType,
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
	$outputDir = File::Spec -> catfile($wk_dir,"qap_Results_for_Diversity_$DateNow");
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

if(defined $seqType){
	if(uc($seqType) eq 'NNT'){
		Info("The files you provided are regarded as non-coding nucleotide sequences.");
	}elsif(uc($seqType) eq 'CNT'){
		Info("The files you provided are regarded as protein coding nucleotide sequences and use the standard condon table.");
	}elsif(uc($seqType) eq 'AA'){
		Info("The files you provided are regarded as amino acid sequences.");
	}else{
		InfoError("The -p/-seqType MUST be one of \'nnt\' or \'cnt\' or \'aa\'.");
		
		pod2usage(-verbose=>2,-exitval=>1);
		exit;
	}
}else{
	InfoError("-p/--seqType MUST be defined.");
	pod2usage(-verbose=>2,-exitval=>1);
	exit;
}


##the core program starts here
Info("Start calculating...");

#check mao
my $distance_nt_mao = File::Spec -> catfile($mainBin, 'lib', 'mao', 'distance_nt.mao');
if(uc($seqType) eq 'NNT'){
	&checkMao($distance_nt_mao);
}

my $ds_mao = File::Spec -> catfile($mainBin, 'lib', 'mao', 'ds.mao');
my $dn_mao = File::Spec -> catfile($mainBin, 'lib', 'mao', 'dn.mao');
if(uc($seqType) eq 'CNT'){
	&checkMao($ds_mao);
	&checkMao($dn_mao);
}

my $distance_aa_mao = File::Spec -> catfile($mainBin, 'lib', 'mao', 'distance_aa.mao');
if(uc($seqType) eq 'AA'){
	&checkMao($distance_aa_mao);
}

#check megacc
my $megacc_excu = File::Spec -> catfile($mainBin, 'bin', '3rdPartyTools', 'megacc', 'megacc');
my $DEBUG_MODE = 1;
if(CheckProgram($megacc_excu, __FILE__, __LINE__, $DEBUG_MODE)){
	#keep running
}else{
	InfoError("The program $megacc_excu does NOT exist. Exiting...");
	exit;
}


#start to run
#mkdir to store result for each file
$outputDir .= "/DetailedResults";
if(-e $outputDir){
	#nothing
}else{
	mkdir $outputDir or die "Can not mkdir $outputDir:$!";
}

if ($threads > 1){
	#go parrallel
	
	my @megacc_excu;
	for my $f(@inputfiles){
		push @megacc_excu,$megacc_excu;
	}
	
	if(uc($seqType) eq 'NNT'){
		#mao
		my @mao;
		map {push @mao,$distance_nt_mao;} @inputfiles;
		
		#output
		my @outputfile;
		for my $f(@inputfiles){
			my $outname = removeFastaSuffix(basename($f));
			my $outfile = File::Spec -> catfile($outputDir, "${outname}_distanceNT.txt.txt");
			push @outputfile,$outfile;
		}
		
		
		#run
		runMultipleThreadsWith4Args(\&diversity_nt, \@megacc_excu, \@mao, \@inputfiles, \@outputfile, $threads);	
		
		#extract result
		my $resfolder = $outputDir =~ s/DetailedResults$//r;
		my $resfile = File::Spec -> catfile($resfolder, "distanceNT.txt");
		
		open RES,">$resfile" or die "Can not output to file $resfile:$!";
		
		print RES "Sample\tValue\n";
		
		for my $f(@inputfiles){
			my $outname = removeFastaSuffix(basename($f));
			my $outfile = File::Spec -> catfile($outputDir, "${outname}_distanceNT.txt");
			
			my $res = &extractRes($outfile);
			print RES "$outname\t$res\n";
		}
		
		close RES;
			
	}elsif(uc($seqType) eq 'CNT'){
		#mao
		my @maods;
		my @maodn;
		map {push @maods,$ds_mao; push @maodn,$dn_mao;} @inputfiles;
		
		#output
		my @dsoutputfile;
		my @dnoutputfile;
		for my $f(@inputfiles){
			my $outname = removeFastaSuffix(basename($f));
			my $dsoutfile = File::Spec -> catfile($outputDir, "${outname}_dS.txt.txt");
			push @dsoutputfile,$dsoutfile;
			my $dnoutfile = File::Spec -> catfile($outputDir, "${outname}_dN.txt.txt");
			push @dnoutputfile,$dnoutfile;
		}
		
		
		#run
		runMultipleThreadsWith6Args(\&dsdn, \@megacc_excu, \@maods, \@maodn, \@inputfiles, \@dsoutputfile, \@dnoutputfile, $threads);
		
		#extract result
		my $resfolder = $outputDir =~ s/DetailedResults$//r;
		my $resfileds = File::Spec -> catfile($resfolder, "dS.txt");
		my $resfiledn = File::Spec -> catfile($resfolder, "dN.txt");
		
		open RESDS,">$resfileds" or die "Can not output to file $resfileds:$!";
		open RESDN,">$resfiledn" or die "Can not output to file $resfiledn:$!";
		
		print RESDS "Sample\tValue\n";
		print RESDN "Sample\tValue\n";
		
		for my $f(@inputfiles){
			my $outname = removeFastaSuffix(basename($f));
			my $outfileds = File::Spec -> catfile($outputDir, "${outname}_dS.txt");
			my $outfiledn = File::Spec -> catfile($outputDir, "${outname}_dN.txt");
			
			my $resds = &extractRes($outfileds);
			my $resdn = &extractRes($outfiledn);
			print RESDS "$outname\t$resds\n";
			print RESDN "$outname\t$resdn\n";
		}
		
		close RESDS;
		close RESDN;

	}else{
		#$seqTpye eq 'aa'
		
		#mao
		my @mao;
		map {push @mao,$distance_aa_mao} @inputfiles;
		
		#output
		my @outputfile;
		for my $f(@inputfiles){
			my $outname = removeFastaSuffix(basename($f));
			my $outfile = File::Spec -> catfile($outputDir, "${outname}_distanceAA.txt.txt");
			push @outputfile,$outfile;
		}
		
		
		#run
		runMultipleThreadsWith4Args(\&diversity_aa, \@megacc_excu, \@mao, \@inputfiles, \@outputfile, $threads);
		
		#extract result
		my $resfolder = $outputDir =~ s/DetailedResults$//r;
		my $resfile = File::Spec -> catfile($resfolder, "distanceAA.txt");
		
		open RES,">$resfile" or die "Can not output to file $resfile:$!";
		
		print RES "Sample\tValue\n";
		
		for my $f(@inputfiles){
			my $outname = removeFastaSuffix(basename($f));
			my $outfile = File::Spec -> catfile($outputDir, "${outname}_distanceAA.txt");
			
			my $res = &extractRes($outfile);
			print RES "$outname\t$res\n";
		}
		
		close RES;
		
	}
}else{
	if(uc($seqType) eq 'NNT'){
		#run
		for my $f (@inputfiles){
			my $outname = removeFastaSuffix(basename($f));
			my $outfile = File::Spec -> catfile($outputDir, "${outname}_distanceNT.txt.txt");
		
			&diversity_nt($megacc_excu, $distance_nt_mao, $f, $outfile);
		}
		
		#extract result
		my $resfolder = $outputDir =~ s/DetailedResults$//r;
		my $resfile = File::Spec -> catfile($resfolder, "distanceNT.txt");
		
		open RES,">$resfile" or die "Can not output to file $resfile:$!";
		
		print RES "Sample\tValue\n";
		
		for my $f(@inputfiles){
			my $outname = removeFastaSuffix(basename($f));
			my $outfile = File::Spec -> catfile($outputDir, "${outname}_distanceNT.txt");
			
			my $res = &extractRes($outfile);
			print RES "$outname\t$res\n";
		}
		
		close RES;
	}elsif(uc($seqType) eq 'CNT'){
		#run
		for my $f (@inputfiles){
			my $outname = removeFastaSuffix(basename($f));
			my $dsoutfile = File::Spec -> catfile($outputDir, "${outname}_dS.txt.txt");
			my $dnoutfile = File::Spec -> catfile($outputDir, "${outname}_dN.txt.txt");
		
			&dsdn($megacc_excu, $ds_mao, $dn_mao, $f, $dsoutfile, $dnoutfile);
		}
		
		#extract result
		my $resfolder = $outputDir =~ s/DetailedResults$//r;
		my $resfileds = File::Spec -> catfile($resfolder, "dS.txt");
		my $resfiledn = File::Spec -> catfile($resfolder, "dN.txt");
		
		open RESDS,">$resfileds" or die "Can not output to file $resfileds:$!";
		open RESDN,">$resfiledn" or die "Can not output to file $resfiledn:$!";
		
		print RESDS "Sample\tValue\n";
		print RESDN "Sample\tValue\n";
		
		for my $f(@inputfiles){
			my $outname = removeFastaSuffix(basename($f));
			my $outfileds = File::Spec -> catfile($outputDir, "${outname}_dS.txt");
			my $outfiledn = File::Spec -> catfile($outputDir, "${outname}_dN.txt");
			
			my $resds = &extractRes($outfileds);
			my $resdn = &extractRes($outfiledn);
			print RESDS "$outname\t$resds\n";
			print RESDN "$outname\t$resdn\n";
		}
		
		close RESDS;
		close RESDN;
	}else{
		#$seqTpye eq 'aa'
		
		#run
		for my $f (@inputfiles){
			my $outname = removeFastaSuffix(basename($f));
			my $outfile = File::Spec -> catfile($outputDir, "${outname}_distanceAA.txt.txt");
		
			&diversity_aa($megacc_excu, $distance_aa_mao, $f, $outfile);
		}
		
		#extract result
		my $resfolder = $outputDir =~ s/DetailedResults$//r;
		my $resfile = File::Spec -> catfile($resfolder, "distanceAA.txt");
		
		open RES,">$resfile" or die "Can not output to file $resfile:$!";
		
		print RES "Sample\tValue\n";
		
		for my $f(@inputfiles){
			my $outname = removeFastaSuffix(basename($f));
			my $outfile = File::Spec -> catfile($outputDir, "${outname}_distanceAA.txt");
			
			my $res = &extractRes($outfile);
			print RES "$outname\t$res\n";
		}
		
		close RES;
	}
}


##the sub program starts here
sub diversity_nt {
	my $megacc_excu = shift;
	my $mao = shift;
	my $data = shift;
	my $outfile = shift;
	
	my $cmd = "$megacc_excu -f fasta -a $mao -d $data --noSumamry -o $outfile";
	runcmd($cmd);
}

sub dsdn {
	my $megacc_excu = shift;
	my $dsmao = shift;
	my $dnmao = shift;
	my $data = shift;
	my $dsoutfile = shift;
	my $dnoutfile = shift;
	
	#dS
	my $cmd = "$megacc_excu -f fasta -a $dsmao -d $data --noSumamry -o $dsoutfile";
	runcmd($cmd);
	
	#dN
	$cmd = "$megacc_excu -f fasta -a $dnmao -d $data --noSumamry -o $dnoutfile";
	runcmd($cmd);
}

sub diversity_aa {
	my $megacc_excu = shift;
	my $mao = shift;
	my $data = shift;
	my $outfile = shift;
	
	my $cmd = "$megacc_excu -f fasta -a $mao -d $data --noSumamry -o $outfile";
	runcmd($cmd);
}
 
sub checkMao {
	my $mao = shift;
	
	if(existFile($mao)){
		#nothing
	}else{
		InfoError("The configuration file $mao is missing. Please fix it.");
		InfoError("Exiting...");
		exit;
	}
	
}

sub extractRes {
	my $file = shift;
	
	my $num = 'NA';
	
	open T,$file or die "Cannot open $file:$!";
	while(<T>){
		chomp;
		if(/\s+([0-9\.]+)\s+[0-9\.]+/){
			$num = $1;
			last;
		}
	}
	close T;
	
	if ($num eq 'NA'){
		InfoError("The sequence diversity calculated for $file is ill-formated or unfinished. Please check again.");
	}
	return $num;
}


##run success
print("\n");
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
           



qap Diversity [options]

Use --help to see more information.

qap is still in development. If you have encounted any problem in usage, please feel no hesitation to cotact us.

=head1 DESCRIPTION

This script implements a function for calculating diversity in batch. The script has B<several> mandatory options that MUST appear last. 

=head1 OPTIONS

=over 5

=item --inputDir,-i F<FILE> [Required]

Path to directory contaning all the files to be read in.

=item --suffix,-s F<STRING> [Optional]

Suffix of the files to be read in. If suffix is not provided, all the files in input directory will be read in.

=item --setType,-p F<STRING> [Required]

The type of your sequence files. Should be one of 'nnt' for non-coding nucleotide sequences, 'cnt' for coding nucleotide sequences and 'aa' for amino acid sequences. 

If -p nnt is specified, the program will calculate nucleotide overall mean distance for input files. If -p cnt is provided, the program will calculate dS and dN. If -p aa is provided 

the program will calculate amino acide mean distance for input files.

=item --outputDir,-o F<FILE> [Optional]

Path of the directory to storage result files. If NOT provided, the program will generate a folder automatically.

=item --threads,-t F<INTEGER> [Optional]

Number of threads this program will use when computing. A positive integer should be provided. The default value is 1.

=item --help,-h

Display this detailed help information.

=back

=head1 EXAMPLE

=over 5

qap Diversity -i ./seq -s fas -p nnt -t 10 -o ./shannon

=back

=head1 AUTHOR

Mingjie Dr.Wang I<huzai@sjtu.edu.cn>

=head1 COPYRIGHT

Copyright (C) 2017, Mingjie Wang. All rights reserved.


