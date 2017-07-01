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


####Use modules in this program####
use General;
use Mapper;
use Caller;

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
my $inputfile;
my $outputDir;
my $ref;
my $cutoff;

my $DateNow = `date +"%Y%m%d_%Hh%Mm%Ss"`;
chomp $DateNow;

GetOptions(
'i|inputFile|=s'     => \$inputfile,
'r|reference|=s'     => \$ref,
'o|outputDir=s'      => \$outputDir,
'c|cutoff|=s'        => \$cutoff,
'h|help|'            => \$help
);

if (defined $outputDir){
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
	$outputDir = File::Spec -> catfile($wk_dir,"qap_Results_for_getMutationCounts_$DateNow");
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
#conver to r input
my $rinputfile = convert2Rinput($inputfile,$outputDir);

#locate the amplicon in the reference genome
my $DEBUG_MODE = 1;


#check r script
my $rscript = File::Spec -> catfile($RealBin,'..','Rscripts','CallVariantsFromMSA.R');
if(existFile($rscript)){
	#keep running
}else{
	InfoError("The script $rscript does NOT exist. Exiting...");
	exit;
}


#get mutation sites
Info("Calling variants from $inputfile");
getMutationFromMSA($rscript,$ref,$rinputfile,$outputDir,1,$cutoff);



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
           



perl getMutationCounts.pl [options]

Use --help to see more information.

gap is still in development. If you have encounted any problem in usage, please feel no hesitation to cotact us.

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

perl getMutationCounts.pl -i test.bam -r HBV.fasta -o ./mutation

=back

=head1 AUTHOR

Mingjie Dr.Wang I<huzai@sjtu.edu.cn>

=head1 COPYRIGHT

Copyright (C) 2017, Mingjie Wang. All rights reserved.

