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

####---------------------------####
####The program begins here
####---------------------------####

##Show welcome
print "You are now running subprogram: ";
printcol ("PhyloTree","green");
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
my $inputDir;
my $outputDir;
my $suffix;
my $seqType;
my $threads;
my $treeMethod;
my $treeModel;
my $bootstrap;
my $bootstrapNumber;
my $treeRates;
my $treeGaps;
my $treeStyle;
my $nodeLabel;
my $edgeLabel;
my $colorScale;
my $fontScale;


my $DateNow = `date +"%Y%m%d_%Hh%Mm%Ss"`;
chomp $DateNow;

GetOptions(
'i|inputDir|=s'        => \$inputDir,
's|suffix|=s'          => \$suffix,
'o|outputDir|=s'       => \$outputDir,
'p|seqType|=s'         => \$seqType,
't|threads|=s'         => \$threads,
'e|method|=s'          => \$treeMethod,
'm|model|=s'           => \$treeModel,
'b|bootstrap|=s'       => \$bootstrap,
'n|bootstrapNumber|=s' => \$bootstrapNumber,
'r|variantRates|=s'    => \$treeRates,
'g|gapTreatment|=s'    => \$treeGaps,
'l|treeStyle|=s'       => \$treeStyle,
'h|help|'              => \$help,
'nodeLabel|=s'         => \$nodeLabel,
'edgeLabel|=s'         => \$edgeLabel,
'colorScale|=s'        => \$colorScale,
'fontScale|=s'         => \$fontScale
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
		}
	}else{
		InfoWarn("The output directory $outputDir already exist.",'yellow');
	}
}else{
	$outputDir = File::Spec -> catfile($wk_dir,"qap_Results_for_PhyloTree_$DateNow");
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

#mkdir to store result for each file
my $tmpDir = File::Spec -> catfile($outputDir, "tmp");
makedir($tmpDir);

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

if(defined $seqType){
	if(uc($seqType) eq 'NNT'){
		Info("The files you provided are regarded as non-coding nucleotide sequences.");
		$seqType = 'nnt';
	}elsif(uc($seqType) eq 'CNT'){
		Info("The files you provided are regarded as protein coding nucleotide sequences and use the standard condon table.");
		$seqType = 'cnt';
	}elsif(uc($seqType) eq 'AA'){
		Info("The files you provided are regarded as amino acid sequences.");
		$seqType = 'aa';
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

if(defined $treeMethod){
	if($treeMethod =~ /me/i){
		$treeMethod = "ME";
	}elsif($treeMethod =~ /ml/i){
		$treeMethod = 'ML';
	}elsif($treeMethod =~ /mp/i){
		$treeMethod = 'MP';
	}elsif($treeMethod =~ /nj/i){
		$treeMethod = 'NJ';
	}elsif($treeMethod =~ /upgma/i){
		$treeMethod = 'UPGMA';
	}else{
		InfoError("The statistical method used for phylogenetic analysis should be one of \'ME\' or \'ML\' or \'MP\' or \'NJ\' or \'UPGMA\'.");
		pod2usage(-verbose=>2,-exitval=>1);
		exit;
	}
}else{
	$treeMethod = "NJ";
}

my %hash_method = ( "ML" => "Maximum Likelihood",
					"NJ" => "Neighbor-joining",
					"ME" => "Minimum Evolution method",
					"UPGMA" => "UPGMA",
					"MP" => "Maximum Parsimony");
Info("[ARGUEMENT] Using $hash_method{$treeMethod} to conduct the phylogenetic tree.");
sleep 1;

if(defined $treeModel){
	my @treeModel;
	if($seqType eq 'nnt' or $seqType eq 'cnt'){
		if($treeMethod eq 'ML'){
			@treeModel = qw/Jukes-Cantor_model Kimura_2-parameter_model Tamura_3-parameter_model Hasegawa-Kishino-Yano_model Tamura-Nei_model General_Time_Reversible_model/;
		}elsif($treeMethod eq "ME" or $treeMethod eq "NJ" or $treeMethod eq "UPGMA"){
			@treeModel = qw/No._of_differences p-distance Jukes-Cantor_model Kimura_2-parameter_model Tajima-Nei_model Tamura_3-parameter_model  Tamura-Nei_model Maximum_Composite_Likelihood LogDet/;
		}elsif($treeMethod eq 'MP'){
			InfoError("--model is not availabe if \'--method MP\' is specified.");
			pod2usage(-verbose=>2,-exitval=>1);
			exit;
		}else{
			InfoError("The statistical method used for phylogenetic analysis should be one of \'ME\' or \'ML\' or \'MP\' or \'NJ\' or \'UPGMA\'.");
			pod2usage(-verbose=>2,-exitval=>1);
			exit;
		}
	}elsif($seqType eq 'aa'){
		if($treeMethod eq 'ML'){
			@treeModel = qw/Poisson_model Equal_input_model Dayhoff_model Dayhoff_model_with_Freqs._(+F) Jones-Taylor-Thornton_(JTT)_model JTT_with_Freqs._(+F)_model WAG_model WAG_with_Freqs._(+F)_model LG_model LG_with_Freqs._(+F)_model General_Reversible_Mitochondrial_(mtREV) mtREV_with_Freqs._(+F)_model General_Reversible_Chloroplast_(cpREV) cpREV_with_Freqs._(+F)_model General_Reverse_Transcriptase_model_(rtREV) rtREV_with_Freqs._(+F)_model/; 
		}elsif($treeMethod eq "ME" or $treeMethod eq "NJ" or $treeMethod eq "UPGMA"){
			@treeModel = qw/No._of_differences p-distance Poisson_model Equal_input_model Dayhoff_model Jones-Taylor-Thornton_(JTT)_model/;
		}elsif($treeMethod eq 'MP'){
			InfoError("--model is not availabe if \'--method MP\' is specified.");
			pod2usage(-verbose=>2,-exitval=>1);
			exit;
		}else{
			InfoError("The statistical method used for phylogenetic analysis should be one of \'ME\' or \'ML\' or \'MP\' or \'NJ\' or \'UPGMA\'.");
			pod2usage(-verbose=>2,-exitval=>1);
			exit;
		}
	}else{
		InfoError("The -p/-seqType MUST be one of \'nnt\' or \'cnt\' or \'aa\'.");
		
		pod2usage(-verbose=>2,-exitval=>1);
		exit;
	}
	if(CheckPositiveInt($treeModel) and $treeModel <= scalar(@treeModel)){
		$treeModel = $treeModel[$treeModel - 1];
		Info("[ARGUEMENT] Using model $treeModel to conduct phylogenetic tree.");
	}else{
		InfoError("Some thing is wrong with \"--model $treeModel --method $treeMethod\". Please select correct statistical model for tree modelling.")
	}
}else{
	if($seqType eq 'nnt' or $seqType eq 'cnt'){
		if($treeMethod ne "MP"){
			$treeModel = 'Kimura_2-parameter_model';
			Info("[ARGUEMENT] Using model $treeModel to conduct phylogenetic tree.");
		}
	}elsif($seqType eq 'aa'){
		if($treeMethod ne "MP"){
			$treeModel = 'Jones-Taylor-Thornton_(JTT)_model';
			Info("[ARGUEMENT] Using model $treeModel to conduct phylogenetic tree.");
		}
	}
}
sleep 1;

if(defined ($bootstrap)){
	if($bootstrap =~ /^y/i){
		$bootstrap = "Bootstrap method";
	}elsif($bootstrap =~ /^n/i){
		$bootstrap = "None";
	}else{
		InfoError("--bootstrap should be one of \'Y\' or \'N\'.");
		pod2usage(-verbose=>2,-exitval=>1);
		exit;
	}
}else{
	$bootstrap = "None";
}

if(defined $bootstrapNumber){
	if($bootstrap eq 'Bootstrap method'){
		if(CheckPositiveInt($bootstrapNumber)){
			if($bootstrapNumber < 50){
				InfoWarn("The minimum number of bootstrap number is 50. Set --bootstrapNumber to 50.");
				$bootstrapNumber = 50;
			}
		}else{
			InfoError("--bootstrapNumber should be a positive integer.");
			pod2usage(-verbose=>2,-exitval=>1);
			exit;
		}
	}else{
		InfoError("--bootstrapNumber is available ONLY if \'--bootstrap Y\' is defined.");
		pod2usage(-verbose=>2,-exitval=>1);
		exit;
	}
}else{
	if($bootstrap eq "Bootstrap method"){
		$bootstrapNumber = 50;
		Info("[ARGUEMENT] Using --bootstrapNumber 50 for bootstrap test as default.");
	}else{
		$bootstrapNumber = "Not Applicable";
	}
}
sleep 1;

my $treeGapsCutoff = "Not Applicable";
if(defined $treeGaps){
	if($treeGaps =~ /^c/i){
		$treeGaps = 'Complete Deletion';
	}elsif($treeGaps =~ /^u/i){
		$treeGaps = 'Use all sites';
	}elsif($treeGaps =~ /^p/i){
		$treeGaps = 'Partial deletion';
		$treeGapsCutoff = 95;
	}else{
		InfoError("Gaps/missing sites treatment should be one of \'--gapTreatment C/U/P\'.");
		pod2usage(-verbose=>2,-exitval=>1);
		exit;
	}
}else{
	$treeGaps = 'Complete Deletion';
}
Info("[ARGUEMENT] Using $treeGaps for gaps/missing sites treatment.");
sleep 1;

my $gammaParam = "Not Applicable";
my $gammaCategory = "Not Applicable";
if(defined $treeRates){
	if($treeMethod eq 'ML'){
		if($treeRates =~ /^u/i){
			$treeRates = 'Uniform Rates';
		}elsif(uc($treeRates) eq 'G'){
			$treeRates = 'Gamma Distributed (G)';
			$gammaCategory = '5';
		}elsif($treeRates =~ /^i/i){
			$treeRates = 'Has Invariant Sites (I)';
		}elsif(uc($treeRates) eq 'GI'){
			$treeRates = 'Gamma Distributed With Invariant Sites (G+I)';
			$gammaCategory = '5';
		}else{
			InfoError("Please select among --variantRates G/I/GI.");
			pod2usage(-verbose=>2,-exitval=>1);
			exit;
		}
		Info("[ARGUEMENT] Using $treeRates to calculate rates among variant sites and number of Discrete Gamma Categories is set to 5.");
	}elsif($treeMethod ne 'MP' and $treeMethod ne 'ML'){
		if($treeRates =~ /^u/i){
			$treeRates = 'Uniform Rates';
		}elsif(uc($treeRates) eq 'G'){
			$treeRates = 'Gamma Distributed (G)';
			$gammaParam = '1.00';
		}elsif($treeRates =~ /^i/i){
			$treeRates = 'Has Invariant Sites (I)';
		}elsif(uc($treeRates) eq 'GI'){
			$treeRates = 'Gamma Distributed With Invariant Sites (G+I)';
			$gammaParam = '1.00';
		}else{
			InfoError("Please select among --variantRates G/I/GI.");
			pod2usage(-verbose=>2,-exitval=>1);
			exit;
		}
		Info("[ARGUEMENT] Using $treeRates to calculate rates among variant sites and Gamma parameter is set to 5.");
	}else{
		InfoError("--variantRates is not available when \'--method MP\' is specified.");
		pod2usage(-verbose=>2,-exitval=>1);
		exit;
	}
}else{
	if($treeMethod ne 'MP'){
		$treeRates = 'Uniform Rates';
		Info("[ARGUEMENT] Using $treeRates to calculate rates among variant sites.");
	}
}
sleep 1;

if(defined $treeStyle){
	if($treeStyle =~ /^p/i){
		$treeStyle = 'phylogram';
	}elsif($treeStyle =~ /^c/i){
		$treeStyle = 'cladogram';
	}elsif($treeStyle =~ /^f/i){
		$treeStyle = 'fan';
	}elsif($treeStyle =~ /^u/i){
		$treeStyle = 'unrooted';
	}elsif($treeStyle =~ /^r/i){
		$treeStyle = 'radial';
	}else{
		InfoError("Choose tree style among \'--treeStyle p/c/f/u/r\'.");
		pod2usage(-verbose=>2,-exitval=>1);
		exit;
	}
}else{
	$treeStyle = 'phylogram';
}
Info("[ARGUEMENT] Using $treeStyle to plot the phylogenetic tree.");
sleep 1;

if(defined $colorScale){
	if($colorScale =~ /^y/i){
		$colorScale = 'Y';
	}elsif($colorScale =~ /^n/i){
		$colorScale = 'N';
	}else{
		InfoError("--colorScale should be \'Y\' or \'N\'.");
		pod2usage(-verbose=>2,-exitval=>1);
		exit;
	}
}else{
	$colorScale = 'N';
}

if(defined $fontScale){
	if($fontScale =~ /^y/i){
		$fontScale = 'Y';
	}elsif($fontScale =~ /^n/i){
		$fontScale = 'N';
	}else{
		InfoError("--colorScale should be \'Y\' or \'N\'.");
		pod2usage(-verbose=>2,-exitval=>1);
		exit;
	}
}else{
	$fontScale = 'N';
}

if(defined $edgeLabel){
	if($edgeLabel =~ /^y/i){
		$edgeLabel = 'Y';
	}elsif($edgeLabel =~ /^n/i){
		$edgeLabel = 'N';
	}else{
		InfoError("--colorScale should be \'Y\' or \'N\'.");
		pod2usage(-verbose=>2,-exitval=>1);
		exit;
	}
}else{
	$edgeLabel = 'N';
}

if(defined $nodeLabel){
	if($nodeLabel =~ /^y/i){
		$nodeLabel = 'Y';
	}elsif($nodeLabel =~ /^n/i){
		$nodeLabel = 'N';
	}else{
		InfoError("--colorScale should be \'Y\' or \'N\'.");
		pod2usage(-verbose=>2,-exitval=>1);
		exit;
	}
}else{
	$nodeLabel = 'N';
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

#check the length of all sequences 
Info("Start checking sequences length.");
my $i = 1;
for my $f (@inputfiles){
	my $len = checkSeqLen($f, $outputDir); #if len does not equal, the program will exit. if $len is output, then it passed the len check
	
	my $name = basename($f);
	printf "[%02d] $name:$len bp\n",$i;
		
	sleep(1);
	$i++;
}

#convert to 2 line and check seq abundance
Info("Checking sequence abundances and converting file format.");
my @inputfiles2lineMod;
my @abundFile;
$i = 1;
for my $f (@inputfiles){
	my $inputfile2line = File::Spec -> catfile($tmpDir, removeFastaSuffix(basename($f)) . ".2line.fasta");
	formatFastaToTwoLineMode($f,$inputfile2line);
	
	my $inputfile2lineMod = File::Spec -> catfile($tmpDir, removeFastaSuffix(basename($f)) . ".2line.mod.fasta");
	push @inputfiles2lineMod, $inputfile2lineMod;
	my ($res,$abundfile) = checkAbundance($inputfile2line, $inputfile2lineMod);
	push @abundFile,$abundfile;
	
	my $name = basename($f);
	printf "[%02d] $name:$res\n",$i;
		
	sleep(1);
	$i++;
}


##the core program starts here
Info("Start calculating...");

#check mao
my @maoFiles;
for my $i (qw/ME ML MP NJ UPGMA/){
	for my $k (qw/nnt cnt aa/){
		my $mao = File::Spec -> catfile($mainBin, 'lib', 'mao', "${i}_${k}.mao");
		&checkMao($mao);
		push @maoFiles,$mao;
	}
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
#run megacc
##get mao files ready
Info("Preparing configuration files.");
sleep(1);
my $maoFileBak = File::Spec -> catfile($mainBin, "lib", 'mao', $treeMethod . "_" . $seqType . ".mao");
my $maoFileToUse = File::Spec -> catfile($tmpDir, $treeMethod . "_" . $seqType . ".mao");
copy($maoFileBak, $maoFileToUse);

if($treeMethod ne 'MP'){
	$treeModel =~ s/_/ /g;
	&changeMao($maoFileToUse, "Model\/Method", $treeModel);
	&changeMao($maoFileToUse, "Test of Phylogeny", $bootstrap);
	&changeMao($maoFileToUse, "No. of Bootstrap Replications", $bootstrapNumber);
	&changeMao($maoFileToUse, "Gaps\/Missing Data Treatment", $treeGaps);
	&changeMao($maoFileToUse, "Site Coverage Cutoff", $treeGapsCutoff);
	&changeMao($maoFileToUse, "Rates among Sites", $treeRates);
	&changeMao($maoFileToUse, "Gamma Parameter", $gammaParam);
	&changeMao($maoFileToUse, "No of Discrete Gamma Categories", '5');
}else{
	&changeMao($maoFileToUse, "Test of Phylogeny", $bootstrap);
	&changeMao($maoFileToUse, "No. of Bootstrap Replications", $bootstrapNumber);
	&changeMao($maoFileToUse, "Gaps\/Missing Data Treatment", $treeGaps);
	&changeMao($maoFileToUse, "Site Coverage Cutoff", $treeGapsCutoff);
}



#draw trees with nwk file
Info("Formatting fasta files.");
my @megacc_excu;
my @maoFileToUse;
my @outfile;
my @nwkfiles;
my @nwkfilesConsensus;
for my $f (@inputfiles){
	push @megacc_excu,$megacc_excu;
	push @maoFileToUse, $maoFileToUse;
	
	my $outfile = $tmpDir . "/" . removeFastaSuffix(basename($f));
	push @outfile,$outfile;
	
	my $nwkfile = removeAllSuffix($outfile) . ".nwk";
	my $nwkfileConsensus = removeAllSuffix($outfile) . "_consensus.nwk";
	push @nwkfiles,$nwkfile;
	push @nwkfilesConsensus, $nwkfileConsensus;
}

runMultipleThreadsWith4Args(\&runMega, \@megacc_excu, \@inputfiles2lineMod, \@maoFileToUse, \@outfile, $threads);

## draw trees
Info("Plotting phylogenetic trees");
my $rscript = File::Spec -> catfile($mainBin, 'bin', 'Rscripts','PhylogeneticTree.R');
if(!existFile($rscript)){
	InfoError("R script $rscript is missing. Please check. Exiting...");
	exit(0);
}

my $plotDir = File::Spec -> catfile($outputDir, 'treePlots');
makedir($plotDir);

my @rscript;
my @fontScale;
my @colorScale;
my @nodeLabel;
my @edgeLabel;
my @treeStyle;
my @outputImage;
for my $f (@inputfiles2lineMod){
	push @rscript,$rscript;
	push @fontScale,$fontScale;
	push @colorScale,$colorScale;
	push @nodeLabel,$nodeLabel;
	push @edgeLabel,$edgeLabel;
	push @treeStyle,$treeStyle;
	my $outputImage = File::Spec -> catfile($plotDir, removeAllSuffix(basename($f)) . "_" . $treeMethod . "_tree");
	push @outputImage,$outputImage;
}

runMultipleThreadsWith10Args(\&drawTree,\@rscript,\@nwkfiles,\@nwkfilesConsensus,\@abundFile,\@colorScale,\@fontScale,\@edgeLabel,\@nodeLabel,\@treeStyle,\@outputImage,$threads);

#move tree files
my $resDir = File::Spec -> catfile($outputDir, "treeFiles");
makedir($resDir);

for my $f (@nwkfiles){
	my $newf = File::Spec -> catfile($resDir, basename($f));
	copy($f, $newf);
}

for my $f (@nwkfilesConsensus){
	my $newf = File::Spec -> catfile($resDir, basename($f));
	copy($f, $newf);
}


##run success
Info("Program completed!",'green');

##the sub program starts here
sub drawTree{
	my $rscript = shift;
	my $nwkfile = shift;
	my $nwkfile_cs = shift;
	my $abundfile = shift;
	my $colorScale = shift;
	my $fontScale = shift;
	my $edgeLabel = shift;
	my $nodeLabel = shift;
	my $treeStyle = shift;
	my $outputFile = shift;
	
	my $cmd = "Rscript $rscript --inputFile1 $nwkfile --inputFile2 $nwkfile_cs --abundFile $abundfile --colorScale $colorScale --fontScale $fontScale --edgeLabel $edgeLabel --nodeLabel $nodeLabel --treeStyle $treeStyle --outputFile $outputFile";
	runcmd($cmd);
}

sub runMega {
	my $megacc_excu = shift;
	my $file = shift;
	my $mao = shift;
	my $outfile = shift;
	
	my $cmd = "$megacc_excu -a $mao -d $file -f fasta -o $outfile";
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

sub changeMao(){
	my $maofile = shift;
	my $arg = shift;
	my $value = shift;
	
	my $maofile2 = $maofile . ".tmp";
	move($maofile,$maofile2);
	
	open T,$maofile2 or die "Can not open mao file $maofile2:$!";
	open OUT,">$maofile" or die "Can not output to $maofile:$!";
	
	while(<T>){
		chomp;
		if(/^$arg/){
			$_ =~ s/\= .*/= $value/;
			print OUT "$_\n";
		}else{
			print OUT "$_\n";
		}
	}
	close T;
	close OUT;
	
	unlink($maofile2) or die "Can not delete file $maofile2:$!";
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
           



qap PhyloTree [options]

Use --help to see more information.

qap is still in development. If you have encounted any problem in usage, please feel no hesitation to cotact us.

=head1 DESCRIPTION

This script implements a function for calculating PhyloTree in batch. The script has B<several> mandatory options that MUST appear last. 

=head1 OPTIONS

=over 5

=item --inputDir,-i F<FILE> [Required]

Path to directory contaning all the files to be read in.

=item --suffix,-s F<STRING> [Optional]

Suffix of the files to be read in. If suffix is not provided, all the files in input directory will be read in.

=item --seqType,-p F<STRING> [Required]

The type of your sequence files. Should be one of 'nnt' for non-coding nucleotide sequences, 'cnt' for coding nucleotide sequences and 'aa' for amino acid sequences. 

=item --method,-e F<STRING> [Optional]

Method used to construct phylogenetic tree. There are 5 methods available. 'ML' for 'Maximum Likelihood', 'NJ' for 'Neighbor-joining', 'ME' for 'Minimum Evolution method', 'UPGMA' for 'Unweighted Pair Group Method with Arithmetic Mean', and 'MP' for 'Maximum Parsimony'. Default value is '--method NJ'.

=item --model,-m F<STRING> [Optional]

Statistical model used for calculation. This parameter is a little complicated. 

--treeModel is only available when --method MP is not defined. Choose model with index. 

For --seqtype nnt/cnt and --method ML, following models are available: 

[1]Jukes-Cantor model; [2]Kimura 2-parameter model; [3]Tamura 3-parameter model; [4]Hasegawa-Kishino-Yano model; [5]Tamura-Nei model; and [6]General Time Reversible model. 

For --seqType nnt/cnt and --method ME/NJ/UPGMA, the following methods are available: 

[1]No. of differences; [2]p-distance; [3]Jukes-Cantor model; [4]Kimura 2-parameter model; [5]Tajima-Nei model; [6]Tamura 3-parameter model; [7]Tamura-Nei model; [8]Maximum Composite Likelihood; [9]LogDet. 

For --seqType aa and --treeMethod ML, the following methods are available: 

[1]Poisson model; [2]Equal input model; [3]Dayhoff model; [4]Dayhoff model with Freqs. (+F); [5]Jones-Taylor-Thornton (JTT) model; [6]JTT with Freqs. (+F) model; [7]WAG model; [8]WAG with Freqs. (+F) model; [9]LG model; [10]LG with Freqs. (+F) model; [11]General Reversible Mitochondrial (mtREV); [12]mtREV with Freqs. (+F) model; [13]General Reversible Chloroplast (cpREV); [14]cpREV with Freqs. (+F) model; [15]General Reverse Transcriptase model (rtREV); [16]rtREV with Freqs. (+F) model. 

For --seqType aa and --treeMethod ME/NJ/UPGMA, the following methods are available:

[1]No. of differences; [2]p-distance; [3]Poisson model; [4]Equal input model; [5]Dayhoff model; [6]Jones-Taylor-Thornton (JTT) model.

If --model is not specified, a default model of 'Kimura 2-parameter model' is used for --seqType nnt/cnt and 'Jones-Taylor-Thornton (JTT) model' is used for --seqType aa. 

=item --bootstrap,-b F<STRING> [Optional]

Whether use bootstrap for phylogeny test. Choose between 'Y' and 'N'. To shorten running time, --bootstrap N is used as default.

=item --bootstrapNumber,-n F<INTEGER> [Optional]

When bootstrap is specified, use --bootstrapNumber to define the number of bootstraps used. Default value is 500.

=item --variantRates,-s F<STRING> [Optional]

Method to treat variants rates. Choose among 'U', 'G' ,'I' and 'GI'. 'U' for 'Uniform Rates', 'G' for 'Gamma Distributed', 'I' for 'Has Invariant Sites' and 'GI' for Gamma Distributed With Invariant Sites'. The default value is 'U'.

=item --gapTreatment,-g F<STRING> [Optional]

Method to treat gaps or missing sites. Choose among [C]omplete Deletion, [U]se all sites, and [P]artial deletion. 'C' for 'Complete Deletion' is used as default.

=item --treeStyle,-l F<STRING> [Optional]

The style of output tree. Choose among [P]hylogram, [C]ladogram, [F]an, [U]nrooted, and [R]adial. Default value is 'P'.

=item --nodeLabel F<STRING> [Optional]

Whether add node labels on phylogenetic tree. Choose between 'Y' or 'N'. Default value is 'N'(No).

=item --edgeLabel F<STRING> [Optional]

Whether add edge labels on phylogenetic tree. Choose between 'Y' or 'N'. Default value is 'N'(No).

=item --colorScale F<STRING> [Optional]

Whether use color scale to show sequence names according to sequence abundance. Choose between 'Y' or 'N'. Default value is 'N'(No).

=item --fontScale F<STRING> [Optional]

Whether use font size to show sequence names according to sequence abundance. Choose between 'Y' and 'N'. Default value is 'N'(No).

=item --outputDir,-o F<FILE> [Optional]

Path of the directory to storage result files. If NOT provided, the program will generate a folder automatically.

=item --threads,-t F<INTEGER> [Optional]

Number of threads this program will use when computing. A positive integer should be provided. The default value is 1.

=item --help,-h

Display this detailed help information.

=back

=head1 EXAMPLE

=over 5

qap PhyloTree -i ./seq -s fas -p nnt -t 10 -o ./shannon

=back

=head1 AUTHOR

Mingjie Dr.Wang I<huzai@sjtu.edu.cn>

=head1 COPYRIGHT

Copyright (C) 2017, Mingjie Wang. All rights reserved.


