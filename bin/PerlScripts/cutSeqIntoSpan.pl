#!/usr/bin/env perl

use diagnostics;
use strict;
use warnings;
use FindBin qw/$RealBin/;
use File::Spec;
use Getopt::Long;
use Pod::Usage;
use lib "$FindBin::Bin/../../lib";
use Cwd qw/getcwd abs_path/;
use File::Basename;
use List::Util qw/max min/;

####Use modules in this program####
use General;

####Flush cache
$| = 1;

####---------------------------####
####The program begins here
####---------------------------####


##get workding directory
my $wk_dir = getcwd;
my $mainBin;
if ($RealBin =~ /(.*)\/bin/){
	$mainBin = $1;
}

####define command line arguments
my $help;
my $inputFile;
my $outputDir;
my $suffix;
my $tileLength;
my $stepSize;

my $DateNow = `date +"%Y%m%d_%Hh%Mm%Ss"`;
chomp $DateNow;

GetOptions(
'i|inputFile|=s'    => \$inputFile,
's|suffix|=s'       => \$suffix,
'o|outputDir|=s'    => \$outputDir,
'l|tileLength=s'    => \$tileLength,
'e|stepSize=s'      => \$stepSize,
'h|help|'           => \$help
);


##check command line arguments
if (defined $help){
	pod2usage(-verbose=>2,-exitval=>1);
}

if (defined $outputDir){
	$outputDir = abs_path($outputDir) . "/";
	if (not -e $outputDir){
		#pod2usage(-verbose=>0,-exitval=>1);
		#exit;
		my $cmd = "mkdir -p $outputDir";
		system($cmd);
	}
}else{
	$outputDir = File::Spec -> catfile($wk_dir,"qap_Results_for_cutSeqIntoSpan_$DateNow");
	
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

if(defined $inputFile){
	if (not -e $inputFile){
		InfoError("Input file $inputFile does NOT exist! Please check again.");
		exit;
	}
}else{
	InfoError("Input file MUST be specified with -i/--inputFile\n");
	pod2usage(-verbose=>0,-exitval=>1);
	exit;
}

my $numberOfFiles = 1;

sleep(1);


if (defined $tileLength){
	if (CheckPositiveInt($tileLength)){
		# nothing
	}else{
		InfoError("The tile length MUST be a positive integer. Please specify the right length with --tileLength/-l.");
		
		pod2usage(-verbose=>2,-exitval=>1);
		exit;
	}
}else{
	InfoError("The tile length MUST be provided. Please specify the tile length with --tileLength/-l.");
		
	pod2usage(-verbose=>2,-exitval=>1);
	exit;
}

if (defined $stepSize){
	if (CheckPositiveInt($stepSize)){
		if($stepSize < $tileLength){
			#nothing
		}else{
			InfoError("The step size MUST be smaller than tile length. Please specify the right step size with --stepSize/-e.");
		
			pod2usage(-verbose=>2,-exitval=>1);
			exit;
		}
	}else{
		InfoError("The step size MUST be a positive integer. Please specify the right step size with --stepSize/-e.");
		
		pod2usage(-verbose=>2,-exitval=>1);
		exit;
	}
}else{
	$stepSize = $tileLength; #If stepsize is not provided, then it equals tile length. tiles have no overlaps.
}

##core program starts here
#check the length of all sequences 
Info("Start checking sequences length.");
my $len = checkSeqLen($inputFile, $outputDir); #if len does not equal, the program will exit. if $len is output, then it passed the len check
#check tile length 
if ($tileLength < $len){
	#nothing
}else{
	InfoError("The tile length is larger than sequence length. Please provide another tile length using --tileLength/-l.");
	
	pod2usage(-verbose=>2,-exitval=>1);
	exit;
}

#start to cut sequence
Info("Start cutting sequences.");
&cutTileWithInfo($inputFile,$outputDir,$tileLength,$stepSize);



##sub program starts here
sub cutTileWithInfo {
	my $inputfile = shift;
	my $outputdir = shift;
	my $tileLength = shift;
	my $stepSize = shift;
	
	Info("Start cutting $inputfile.");
	
	#convert fasta file to two line format
	my $inputfile2Name = removeFastaSuffix(basename($inputfile)) . ".2line.fasta";
	my $inputfile2 = File::Spec -> catfile($outputDir,$inputfile2Name);
	formatFastaToTwoLineMode($inputfile,$inputfile2);
	
	open T,$inputfile2 or die "Can not open fasta file $inputfile2:$!";
	while(my $line1 = <T>){
		chomp $line1;
		chomp (my $line2 = <T>);
		
		$line2 = uc($line2);
		
		my $len = length($line2);
		
		my $start = 1;
		for (my $start = 1; $start - 1 + $tileLength <= $len - 1 + $stepSize; $start += $stepSize){
			my $cut = substr $line2,$start - 1, $tileLength;
			
			my $end = $start - 1 + $tileLength;
			if($end > $len){
				$end = $len;
			}
			
			#output
			my $outputfileName = $start . "_" . $end . ".fasta";
			my $outputsubdir = $outputdir;
			if (-e $outputsubdir){
				#nothing
			}else{
				mkdir $outputsubdir or die "Can not mkdir $outputsubdir:$!";
			}
			my $outputfile = File::Spec -> catfile($outputsubdir,$outputfileName);
			
			open RES,">>$outputfile" or die "Can NOT output to file $outputfile:$!";
			
			print RES "$line1\n$cut\n";
			
			close RES;
		}
		
	}
	close T;
	
	system("rm -rf $inputfile2");
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
           



perl cutSeqIntoSpan.pl [options]

Use --help to see more information.

gap is still in development. If you have encounted any problem in usage, please feel no hesitation to cotact us.

=head1 DESCRIPTION

This script implements a function for cutting fasta sequences into tiles. The script has B<several> mandatory options that MUST appear last. 

=head1 OPTIONS

=over 5

=item --inputFile,-i F<FILE> [Required]

Path to input file contaning MSA results.

=item --outputDir,-o F<FILE> [Optional]

Path of the directory to storage result files. If NOT provided, the program will generate a folder automatically.

=item --tileLength,-l F<INTEGER> [Required]

Length of the tile (short sequences) to cut. 

=item --stepSize,-e F<INTEGER> [Optional]

Size to move the step when cutting tiles. For example, if --tileLength 100 -stepSize 20, the sequence will be cutted from 1-100,20-120,40-140,60-160 etc. If --stepSize is not defined, the program will not cut the sequences into tiles with a stepwise manner.

=item --help,-h

Display this detailed help information.

=back

=head1 EXAMPLE

=over 5

perl cutSeqIntoSpan -i ./seq.fasta -l 100 -e 20 -o ./cutTile

=back

=head1 AUTHOR

Mingjie Dr.Wang I<huzai@sjtu.edu.cn>

=head1 COPYRIGHT

Copyright (C) 2017, Mingjie Wang. All rights reserved.

