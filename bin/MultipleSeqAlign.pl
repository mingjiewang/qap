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

####Flush cache
$| = 1;

####---------------------------####
####The program begins here
####---------------------------####

##Show welcome
print "You are now running subprogram: ";
printcol ("MultipleSeqAlign","green");
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
my $program;
my $quiet;
my $verbose;
my $threads;

my $DateNow = `date +"%Y%m%d_%Hh%Mm%Ss"`;
chomp $DateNow;

GetOptions(
'i|inputDir|=s'     => \$inputDir,
's|suffix|=s'       => \$suffix,
'o|outputDir|=s'    => \$outputDir,
't|threads|=s'      => \$threads,
'p|program|=s'      => \$program,
'q|quiet|'          => \$quiet,
'v|verbose|'        => \$verbose,
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
	$outputDir = File::Spec -> catfile($wk_dir,"qap_Results_for_MultipleSeqAlign_$DateNow");
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

my $msaprogram;
if (defined $program){
	if($program == 1){
		$msaprogram = 'co'; # co for clustalo
	}elsif($program == 2){
		$msaprogram = 'cw'; # cw for clustalw
	}elsif($program == 3){
		$msaprogram = 'mu'; # mu for muscle
	}else{
		InfoError("Please provide a program for MSA among \'Clustal Omega\'(1) or \'ClustalW2\'(2) or \'MUSCLE\'(3).");
		
		pod2usage(-verbose=>2,-exitval=>1);
		exit;
	}
	
}else{
	$msaprogram = 'co'; #using clustalo as default
}

my $showinfo = 1;
if(defined $quiet and not defined $verbose){
	$showinfo = 0;
}elsif(not defined $quiet and defined $verbose){
	$showinfo = 1;
}elsif(not defined $quiet and not defined $verbose){
	$showinfo = 1;
}else{
	InfoError("Both -quiet/-q and --verbose/-v are defined. Only one of them chould be used.");
	
	pod2usage(-verbose=>2,-exitval=>1);
	exit;
}


##core program starts here

#first thing to do is to check the msa program
#co for clustalo
my $coExcu = File::Spec -> catfile($RealBin,'3rdPartyTools','msa','clustalo','clustalo');
#cw for clustalw
my $cwExcu = File::Spec -> catfile($RealBin,'3rdPartyTools','msa','clustalw','clustalw');
#mu for muscle
my $muExcu = File::Spec -> catfile($RealBin,'3rdPartyTools','msa','muscle','muscle');

if($msaprogram eq 'co'){
	if(CheckProgram($coExcu, __FILE__, __LINE__, $DEBUG_MODE)){
		#keep running
	}else{
		InfoError("The program $coExcu does NOT exist. Exiting...");
		exit;
	}
}elsif($msaprogram eq 'cw'){
	if(CheckProgram($cwExcu, __FILE__, __LINE__, $DEBUG_MODE)){
		#keep running
	}else{
		InfoError("The program $cwExcu does NOT exist. Exiting...");
		exit;
	}
}else{
	#$msaprogram eq 'mu'
	if(CheckProgram($muExcu, __FILE__, __LINE__, $DEBUG_MODE)){
		#keep running
	}else{
		InfoError("The program $muExcu does NOT exist. Exiting...");
		exit;
	}
}

## core program, start the alignment
if($threads > 1){
	Info("Start MSA using multiple threads. Please wait patiently.");
	my $threadsForEachSample = int($threads / $numberOfFiles);
	if($threadsForEachSample > 1){
		my @program;
		my @coExcu;
		my @cwExcu;
		my @muExcu;
		my @outputfiles;
		my @threadsForEachSample;
		my @showinfo;
		for my $f (@inputfiles){
			push @program,$msaprogram;
			push @coExcu,$coExcu;
			push @cwExcu,$cwExcu;
			push @muExcu,$muExcu;
			my $outfilename = removeFastaSuffix(basename($f)) . "_ALIGN.fasta";
			my $outfile = File::Spec -> catfile($outputDir, $outfilename);
			push @outputfiles,$outfile;
			push @threadsForEachSample,$threadsForEachSample;
			push @showinfo,$showinfo;
		}
		
		runMultipleThreadsWith8Args(\&MSAWithInfo,\@program,\@coExcu,\@cwExcu,\@muExcu,\@inputfiles,\@outputfiles,\@threadsForEachSample,\@showinfo,$numberOfFiles);
		
	}else{
		my @program;
		my @coExcu;
		my @cwExcu;
		my @muExcu;
		my @outputfiles;
		my @threadsForEachSample;
		my @showinfo;
		for my $f (@inputfiles){
			push @program,$msaprogram;
			push @coExcu,$coExcu;
			push @cwExcu,$cwExcu;
			push @muExcu,$muExcu;
			my $outfilename = removeFastaSuffix(basename($f)) . "_ALIGN.fasta";
			my $outfile = File::Spec -> catfile($outputDir, $outfilename);
			push @outputfiles,$outfile;
			push @threadsForEachSample,1;
			push @showinfo,$showinfo;
		}
		
		runMultipleThreadsWith8Args(\&MSAWithInfo,\@program,\@coExcu,\@cwExcu,\@muExcu,\@inputfiles,\@outputfiles,\@threadsForEachSample,\@showinfo,$threads);
		
	}
}else{
	Info("Start MSA using single thread. Please wait patiently if your input data is large.");
	if ($showinfo){
		for my $f(@inputfiles){
			my $outfilename = removeFastaSuffix(basename($f)) . "_ALIGN.fasta";
			my $outfile = File::Spec -> catfile($outputDir, $outfilename);
			
			&MSAWithInfo($msaprogram,$coExcu,$cwExcu,$muExcu,$f,$outfile,1,$showinfo);
		}
	}else{
		my $i = 1;
		for my $f(@inputfiles){
			my $outfilename = removeFastaSuffix(basename($f)) . "_ALIGN.fasta";
			my $outfile = File::Spec -> catfile($outputDir, $outfilename);
			
			InfoSTDERRProcessBar($i,$numberOfFiles);
			
			&MSAQuiet($msaprogram,$coExcu,$cwExcu,$muExcu,$f,$outfile,1);
			
			$i++;
		}
	}
	
}


##sub program starts here
#MSA is only used when -q is specified and single threads used. So only progress bar is shown. No other output to interupt the bar.
sub MSAQuiet {
	my $program = shift;
	my $coExcu = shift;
	my $cwExcu = shift;
	my $muExcu = shift;
	my $inputfile = shift;
	my $outputfile = shift;
	my $threads = shift;
	
	
	if ($program eq 'co'){
		#./clustalo -i ./test.fasta --infmt=fa -o ./align.fasta --outfmt=fa --threads=10 --wrap=1000000 --force --verbose
		my $cmd = "$coExcu -i $inputfile -o $outputfile --outfmt=fa --threads=$threads --wrap=1000000 --force";
		system($cmd);
	}elsif($program eq 'cw'){
		#even though -QUIET is specified, ClustalW2 still output running details. so change STDERR.
		my $stdoutlog = File::Spec -> catfile(dirname($outputfile),"runningDetail.log");
		open STDOUT,">>$stdoutlog" or die "Can NOT output log to $stdoutlog:$!";
	
		#./clustalw -INFILE=./test.fasta -OUTFILE=./align.fasta -OUTPUT=fasta -OUTORDER=INPUT -QUIET
		my $cmd = "$cwExcu -INFILE=$inputfile -OUTFILE=$outputfile -OUTPUT=fasta -OUTORDER=INPUT -QUIET";
		system($cmd);
		
		system("rm -rf $stdoutlog");
		
		#rm dnd tree
		my $dndtree = File::Spec -> catfile(dirname($inputfile),removeFastaSuffix(basename($inputfile)) . ".dnd");
		$cmd = "rm -rf $dndtree";
		system($cmd);
	}elsif($program eq 'mu'){
		#./muscle -in ./test.fasta -out ./align.fasta -maxhours 1 -quiet
		my $cmd = "$muExcu -in $inputfile -out $outputfile -maxhours 1 -quiet";
		system($cmd);
	}else{
		InfoError("Please provide a program for MSA among \'Clustal Omega\'(1) or \'ClustalW2\'(2) or \'MUSCLE\'(3).");
		
		pod2usage(-verbose=>2,-exitval=>1);
		exit;
	}
	
	#format the output fasta file to 2 line mode
	my $finalout = $outputfile =~ s/_ALIGN.fasta$/_MSA.fasta/r;
	formatFastaToTwoLineMode($outputfile,$finalout);
	system("rm -rf $outputfile");
}

sub MSAWithInfo {
	my $program = shift;
	my $coExcu = shift;
	my $cwExcu = shift;
	my $muExcu = shift;
	my $inputfile = shift;
	my $outputfile = shift;
	my $threads = shift;
	my $showinfo = shift;
	
	Info("Start MSA for $inputfile.");
	if ($program eq 'co'){
		if($showinfo){
			#./clustalo -i ./test.fasta --infmt=fa -o ./align.fasta --outfmt=fa --threads=10 --wrap=1000000 --force --verbose
			my $cmd = "$coExcu -i $inputfile -o $outputfile --outfmt=fa --threads=$threads --wrap=1000000 --force --verbose";
			runcmd($cmd);
		}else{
			#./clustalo -i ./test.fasta --infmt=fa -o ./align.fasta --outfmt=fa --threads=10 --wrap=1000000 --force
			my $cmd = "$coExcu -i $inputfile -o $outputfile --outfmt=fa --threads=$threads --wrap=1000000 --force";
			runcmd($cmd);
		}
	}elsif($program eq 'cw'){
		if($showinfo){
			#./clustalw -INFILE=./test.fasta -OUTFILE=./align.fasta -OUTPUT=fasta -OUTORDER=INPUT 
			my $cmd = "$cwExcu -INFILE=$inputfile -OUTFILE=$outputfile -OUTPUT=fasta -OUTORDER=INPUT";
			runcmd($cmd);
		}else{
			#./clustalw -INFILE=./test.fasta -OUTFILE=./align.fasta -OUTPUT=fasta -OUTORDER=INPUT 
			my $cmd = "$cwExcu -INFILE=$inputfile -OUTFILE=$outputfile -OUTPUT=fasta -OUTORDER=INPUT -QUIET";
			runcmd($cmd);
		}
		
		#rm dnd tree
		my $dndtree = File::Spec -> catfile(dirname($inputfile),removeFastaSuffix(basename($inputfile)) . ".dnd");
		my $cmd = "rm -rf $dndtree";
		system($cmd);
	}elsif($program eq 'mu'){
		if($showinfo){
			#./muscle -in ./test.fasta -out ./align.fasta -maxhours 1 
			my $cmd = "$muExcu -in $inputfile -out $outputfile -maxhours 1 ";
			runcmd($cmd);
		}else{
			#./muscle -in ./test.fasta -out ./align.fasta -maxhours 1 -quiet
			my $cmd = "$muExcu -in $inputfile -out $outputfile -maxhours 1 -quiet";
			runcmd($cmd);
		}
	}else{
		InfoError("Please provide a program for MSA among \'Clustal Omega\'(1) or \'ClustalW2\'(2) or \'MUSCLE\'(3).");
		
		pod2usage(-verbose=>2,-exitval=>1);
		exit;
	}
	
	#format the output fasta file to 2 line mode
	my $finalout = $outputfile =~ s/_ALIGN.fasta$/_MSA.fasta/r;
	formatFastaToTwoLineMode($outputfile,$finalout);
	system("rm -rf $outputfile");
}

sub InfoSTDERRProcessBar{
	local $| = 1;
	my $i = $_[0] || return 0;
	my $n = $_[1] || return 0;
	
	chomp (my $date = `date`);
    
    #output               
	print STDERR "\rINFO    : Progress [ ".("=" x int(($i/$n)*50)).(" " x (50 - int(($i/$n)*50)))." ] ";
	my $tmp = sprintf("%4.2f %%",$i/$n*100);
	print STDERR $tmp;
	local $| = 0;
}

sub InfoSTDERR {
	#change output from stdout to stderr, because stdout is changed to log file to omit clustalw output
    my ($info,$color) = @_;
    if (!defined $color){
        $color = "white";
    }

    chomp (my $date = `date`);
    
    #color ref
    my %ColorShell = ("black"     => 30,
                      "red"       => 31,
                      "green"     => 32,
                      "yellow"    => 33,
                      "blue"      => 34,
                      "purple"    => 35,
                      "turquoise" => 36,
                      "white"     => 37,);
    my $outcolor = 37; #If inputcolor is not in the ColorShell, use white as default
    if (defined $ColorShell{$color}){
        $outcolor = $ColorShell{$color};
    }
    #output               
    print STDERR "INFO    \@ \[${date}\]: ";
    system "echo -e \"\e[0;${outcolor};1m${info}\e[m\" >&2"; # echo the information to STDERR in shell by using '&2'
}


##run success
#change output from stdout to stderr, because stdout is changed to log file to omit clustalw output
print STDERR ("\n\n");
&InfoSTDERR("Program completed!",'green');


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
           



qap MultipleSeqAlign [options]

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

=item --outputDir ,-o F<FILE> [Optional]

Path of the directory to storage result files. If NOT provided, the program will generate a folder automatically.

=item --threads,-t F<INTEGER> [Optional]

Number of threads this program will use when computing. A positive integer should be provided. The default value is 1.

=item --program,-p F<STRING> [Optional]

Program to use for sequence alignment. '1' for 'Clustal Omega'; '2' for 'ClustalW2'; '3' for 'MUSCLE'. The program uses 'Clustal Omega' as default.

=item --quiet,-q [Optional]

Keep the running detail to minimum if -q/--quiet is specified.

=item --verbose,-v [Optional]

Output the running details to screen. Only one of -q and -v can be specified. -v is on as default.

=item --help,-h

Display this detailed help information.

=back

=head1 EXAMPLE

=over 5

qap MultipleSeqAlign -i ./seq -t 10 -s fasta -p 1 -q -o ./align

=back

=head1 AUTHOR

Mingjie Dr.Wang I<huzai@sjtu.edu.cn>

=head1 COPYRIGHT

Copyright (C) 2017, Mingjie Wang. All rights reserved.

