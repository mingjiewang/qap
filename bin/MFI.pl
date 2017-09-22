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
printcol ("MFI","green");
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
my $threads;
my $ref;
my $spanLength;
my $mode;
my $graphic;
my $cutoff;

my $DateNow = `date +"%Y%m%d_%Hh%Mm%Ss"`;
chomp $DateNow;

GetOptions(
'i|inputFile|=s'     => \$inputfile,
'r|reference|=s'     => \$ref,
'o|outputDir|=s'     => \$outputDir,
'm|mode|=s'          => \$mode,
'l|spanLength|=s'    => \$spanLength,
'h|help|'            => \$help,
't|threads|=s'       => \$threads,
'g|graphic|'         => \$graphic,
'c|cutoff|=s'        => \$cutoff
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
		my $cmd = "mkdir -p $outputDir";
		system($cmd);
	}else{
		InfoError("Folder $outputDir already exists. Please specify another output directory using option -o/--outputDir");
		exit;
	}
}else{
	$outputDir = File::Spec -> catfile($wk_dir,"qap_Results_for_MFI_$DateNow");
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

if(defined $mode){
	if(uc($mode) eq 'SEQ' or uc($mode) eq 'SPAN'){
		#nothing
	}else{
		InfoError("The mode should be selected between \'seq\' or \'span\'.");
		pod2usage(-verbose=>2,-exitval=>1);
		exit(0);
	}
}

if(defined $spanLength){
	if (uc($mode) eq 'SEQ'){
		InfoError("Running mode \'--mode/-m\' seq and \'--spanLength' can NOT be specified at the same time.");
		pod2usage(-verbose=>2,-exitval=>1);
		exit;
	}
	
	if (CheckPositiveInt($spanLength)){
		$mode = 'span';
	}else{
		InfoError("The length of span should be a positive number");
		exit(0);
	}
}else{
	if (uc($mode) eq 'SEQ'){
		#nothing
	}elsif(uc($mode) eq 'SPAN'){
		InfoError("--spanLength/-l should be provided when running mode \'--mode/-m span\' is selected.");
		pod2usage(-verbose=>2,-exitval=>1);
		exit;
	}else{
		Info("Both running mode and span length are NOT provided. Using \'--mode/-m seq\' as default.");
		$mode = 'seq';
	}
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

if(defined $graphic){
	if(uc($mode) eq 'SEQ'){
		InfoWarn("Ignoring \'-g\'. Graphs will NOT be drawn under 'seq' mode.");
	}
}

##core program starts here
#check rscript
my $rscript = File::Spec -> catfile($mainBin,'bin','Rscripts','CalculateMFI.R');
if (! existFile($rscript)){
	InfoError("The R script $rscript does NOT exist, please check.");
	exit(0);
}

#check blat
my $blatProgram = File::Spec -> catfile($RealBin,'3rdPartyTools','blat','blat');
if(CheckProgram($blatProgram, __FILE__, __LINE__, $DEBUG_MODE)){
	#keep running
}else{
	InfoError("The program $blatProgram does NOT exist. Exiting...");
	exit;
}

#first check the seq len of MSA file
my $seqlen = checkSeqLen($inputfile,$outputDir);

#tmp data file
my $tmpdir = File::Spec -> catfile($outputDir,"tmp");
makedir($tmpdir);

#start to run
if(uc($mode) eq 'SEQ'){
	#run MutationCallerFromMSA to quantify mutation counts
	my $mc = File::Spec -> catfile($mainBin,"bin","PerlScripts","getMutationCounts.pl");
	my $cmd = "perl $mc -i $inputfile -r $ref -o $tmpdir -c $cutoff";
	runcmd($cmd);
	
	#calculate MFI
	my $mutationFile = File::Spec -> catfile($tmpdir,'result',removeFastaSuffix(basename($inputfile)) . ".txt");
	my $outputfile = File::Spec -> catfile($outputDir,removeFastaSuffix(basename($inputfile)) . ".MFI.txt");
	$cmd = "Rscript $rscript -i $mutationFile -o $outputfile -l $seqlen";
	runcmd($cmd);
}elsif(uc($mode) eq 'SPAN'){
	#first cut seq into tiles
	my $cutspanscript = File::Spec -> catfile($mainBin,'bin','PerlScripts','cutSeqIntoSpan.pl');
	my $inputdir = dirname($inputfile);
	my $spandir = File::Spec -> catfile($tmpdir,'span');
	makedir($spandir);
	
	my $cmd = "perl $cutspanscript -i $inputfile -o $spandir -l $spanLength ";
	runcmd($cmd);
	
	#get and sort cutted span files
	my %tmpname;
	for my $f (glob "$spandir/*.fasta"){
		my $filename = basename($f);
		$filename =~ /(\d+)_\d+.*/;
		$tmpname{$f} = $1;
	}
	my @spanfiles;
	for my $n (sort {$tmpname{$a} <=> $tmpname{$b}} keys %tmpname){
		push @spanfiles,$n;
	}
	#map {print "$_\n"} @spanfiles; 
	
	## Get ref file for each span seq file. 
	## Because span seqs might be cutted very short. 
	## So ref seq for span seqs should better be directly cutted from ref file instead of mapping to the ref genome,
	## which might cause non unique mapping positions
	
	#get location from blat output
	my $refdir = File::Spec -> catfile($tmpdir,'ref');
	makedir($refdir);
	my ($psl,$start,$end) = blatPipeline($blatProgram,$inputfile,$ref,$refdir);
	
	#get total ref seq
	my $ref2 = File::Spec -> catfile($refdir,removeFastaSuffix(basename($ref)) . ".2line.fasta");
	formatFastaToTwoLineMode($ref,$ref2);
	system("sed -i \'s/\r//g\' $ref2");
	open REF,"$ref2" or die "Can NOT open $ref2:$!";
	chomp (my $refname = <REF>);
	$refname =~ s/^>//;
	chomp (my $totalref = <REF>);
	close REF;
	
	#get ref
	my $refseq = substr $totalref, $start - 1, $end - $start + 1;
	
	#check ref and seq length
	if (length($refseq) == $seqlen){
		#nothing
	}else{
		InfoError("The length of MSA sequences and the reference fasta file do NOT equal! Please check.");
		exit;
	}
	
	#cut ref into pieces
	my @subref;
	my @start;
	my @end;
	for (my $i = 0; $i <= $seqlen; $i += $spanLength){
		my $subseq = substr $refseq,$i,$spanLength;
		my $start = $i + 1;
		push @start,$start;
		my $end = $i + $spanLength;
		if ($end > $seqlen){
			$end = $seqlen;
		}
		push @end,$end;
		my $subRefSeqFile = File::Spec -> catfile($refdir, $start . "_" . $end . "_" . time() . ".fasta" );
		open T,">$subRefSeqFile" or die "Can NOT output to file $subRefSeqFile:$!";
		print T "$subseq\n";
		close T;
		
		push @subref,$subRefSeqFile;
	}
	
	#calculate mutation counts for each file
	my $mutationdir = File::Spec -> catfile($tmpdir,'mutation');
	makedir($mutationdir);
	
	#scripts used
	my $mc = File::Spec -> catfile($mainBin,"bin","PerlScripts","getMutationCountsV2.pl");
	
	#output mfi dir
	my $mfidir = File::Spec -> catfile($tmpdir,'MFI');
	makedir($mfidir);
	
	#final output file
	my $resfile = File::Spec -> catfile($outputDir,removeFastaSuffix(basename($inputfile)) . ".MFI.txt");
	if (-e $resfile){
		system("rm -rf $resfile");
	}
	open RES,">>$resfile" or die "Can NOT output to file $resfile:$!\n";
	print RES "Start\tEnd\tMFI\n";
	
	if($threads > 1){
		my @mc;
		my @outdir;
		my @cutoff;
		my @subseqlen;
		my @rscript;
		my @mutationfiles;
		my @mfifiles;
		for my $f (@spanfiles){
			push @mc,$mc;
			push @outdir,$mutationdir;
			push @cutoff,$cutoff;
			my $subseqlen = checkSeqLen($f,$spandir);
			push @subseqlen,$subseqlen;
			push @rscript,$rscript;
			my $mutationfile = File::Spec -> catfile($mutationdir, 'result',removeSuffix(basename($f)) . ".txt");
			push @mutationfiles,$mutationfile;
			my $mfifile = File::Spec -> catfile($mfidir,removeSuffix(basename($f)) . ".MFI.txt");
			push @mfifiles,$mfifile;
		}
		
		#get mut count
		runMultipleThreadsWith5Args(\&getMutationCount,\@mc,\@spanfiles,\@subref,\@outdir,\@cutoff,$threads);
		
		#get mfi
		runMultipleThreadsWith4Args(\&calculateMFI,\@rscript,\@mutationfiles,\@mfifiles,\@subseqlen,$threads);
		
		#merge mfi results
		for my $i (0..scalar(@spanfiles) - 1){
			my $mfifile = File::Spec -> catfile($mfidir,removeSuffix(basename($spanfiles[$i])) . ".MFI.txt");
			open T,$mfifile or die "Can not open $mfifile:$!\n";
			chomp (my $mfi = <T>);
			close T;
			if($mfi eq 'NA'){
				$mfi = 0;
			}
			
			my $start = $start[$i];
			my $end = $end[$i];
			
			print RES "$start\t$end\t$mfi\n";
		}
	}else{
		for my $i (0..scalar(@spanfiles) - 1){
			#get mut count
			&getMutationCount($mc,$spanfiles[$i],$subref[$i],$mutationdir,$cutoff);
			
			#get mfi
			my $mutationfile = File::Spec -> catfile($mutationdir, 'result',removeSuffix(basename($spanfiles[$i])) . ".txt");
			my $subseqlen = checkSeqLen($spanfiles[$i],$spandir);
			my $mfifile = File::Spec -> catfile($mfidir,removeSuffix(basename($spanfiles[$i])) . ".MFI.txt");
			&calculateMFI($rscript,$mutationfile,$mfifile,$subseqlen);
			
		}
		
		#merge mfi results
		for my $i (0..scalar(@spanfiles) - 1){
			my $mfifile = File::Spec -> catfile($mfidir,removeSuffix(basename($spanfiles[$i])) . ".MFI.txt");
			open T,$mfifile or die "Can not open $mfifile:$!\n";
			chomp (my $mfi = <T>);
			close T;
			if($mfi eq 'NA'){
				$mfi = 0;
			}
			
			my $start = $start[$i];
			my $end = $end[$i];
			
			print RES "$start\t$end\t$mfi\n";
		}
	}
	
	close RES;
	
	##draw graphs
	if(defined $graphic){
		Info("Drawing plots for MFI results");
		
		my $plotRscript = File::Spec -> catfile($mainBin,'bin','Rscripts','MFIPlots.R');
		if(not existFile($plotRscript)){
			InfoError("The R script for MFI visualization is missing.Please check.");
			exit;
		}
		
		my $outputplotdir = File::Spec -> catfile($outputDir,'plots');
		makedir($outputplotdir);
		
		my $outputplot = File::Spec -> catfile($outputplotdir, removeFastaSuffix(basename($inputfile)));
		my $cmd = "Rscript $plotRscript -i $resfile -o $outputplot";
		runcmd($cmd);
	}
}


##run success
Info("Program completed!",'green');

##sub-program starts here
sub getMutationCount {
	my $perlscript = shift;
	my $inputfile = shift;
	my $ref = shift;
	my $outputdir = shift;
	my $cutoff = shift;
	
	my $cmd = "perl $perlscript -i $inputfile -r $ref -o $outputdir -c $cutoff";
	runcmd($cmd);
}

sub calculateMFI {
	my $rscript = shift;
	my $mutationfile = shift;
	my $outputfile = shift;
	my $seqlen = shift;
	
	my $cmd = "Rscript $rscript -i $mutationfile -o $outputfile -l $seqlen";
	runcmd($cmd);
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
           



qap MFI [options]

Use --help to see more information.

qap is still in development. If you have encounted any problem in usage, please feel no hesitation to cotact us.

=head1 DESCRIPTION

This script implements a function to calculate mutation frequency index (MFI). The script has B<several> mandatory options that MUST appear last. 

=head1 OPTIONS

=over 5

=item --inputFile,-i F<FILE> [Required]

Path to input MSA file (fasta format) which is required.

=item --refSeq,-r F<File> [Required]

Path to the reference fasta file. 

=item --outputDir,-o F<FILE> [Optional]

Path of the directory to storage result files. If NOT provided, the program will generate a folder automatically.

=item --mode,-m F<STRING> [Optional]

While mode to use for MFI calculation. Please select between 'seq' and 'span' by which 'seq' means calculate MFI for whole length sequence and 'span' means cutting the sequences into spans and then calculate MFI. If both --mode and --spanLength are not specified, the program will use '--mode/-m seq' as default. 

=item --spanLength,-l F<INTEGER> [Optional]

If --mode/-m span is specified, the span length MUST be provided which should be a positive integer. If --spanLenth/-l is specified, the running

mode will be set to '--mode/-m span'.

=item --graphic,-g [Optional]

Whether to draw graphs for MFI visualization if '--mode/-m span' is selected. If '--mode/-m seq' is specified, --graphic/-g will not be used.

=item --threads,-t F<INTEGER> [Optional]

Number of threads this program will use when computing. A positive integer should be provided. The default value is 1.

=item --cutoff,-c F<DOUBLE> [Optional]

The cutoff value of variant frequency. The default value is 0.01.

=item --help,-h

Display this detailed help information.

=back

=head1 EXAMPLE

=over 5

qap MFI -i test.fasta -r HBV.fasta -m span -l 50 -t 5 -o ./mfi

=back

=head1 AUTHOR

Mingjie Dr.Wang I<huzai@sjtu.edu.cn>

=head1 COPYRIGHT

Copyright (C) 2017, Mingjie Wang. All rights reserved.

