#!/usr/bin/env perl

package QSR;

use diagnostics;
use strict;
use warnings;
use FindBin qw/$RealBin/;
use File::Spec;
use Getopt::Long;
use Pod::Usage;
use lib "$FindBin::Bin";
use Cwd qw/getcwd abs_path/;
use File::Basename;
use File::Copy;

####Use modules in this program####
use General;
use Mapper;

use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

our $VERSION     = 1.00;
our @ISA         = qw(Exporter);
our @EXPORT      = qw( shorah_pipeline predicthaplo_pipeline qure_pipeline viquas_pipeline qsr_pipeline );
our @EXPORT_OK   = qw( shorah_pipeline predicthaplo_pipeline qure_pipeline viquas_pipeline qsr_pipeline );
our %EXPORT_TAGS = ( DEFAULT => [qw( &shorah_pipeline &predicthaplo_pipeline &qure_pipeline &viquas_pipeline &qsr_pipeline )]);


##global variables
##get workding directory
my $wk_dir = getcwd;
my $mainBin;
if ($RealBin =~ /(.*)\/bin/){
	$mainBin = $1;
}

##subprogram goes here

sub shorah_pipeline {
	my $amplicon = shift;
	my $shorah_excu = shift;
	my $bamfile = shift;
	my $ref = shift;
	my $outdir = shift;
	
	chdir $outdir or die "Can NOT chdir to $outdir:$!";
	my $cmd;
	if($amplicon eq 'Y'){
		$cmd = "python $shorah_excu amplicon -b $bamfile -f $ref -d"
	}else{
		$cmd = "python $shorah_excu -b $bamfile -f $ref";
	}
	$cmd = "bash -c 'source activate python3 && $cmd && source deactivate'";

	runcmd($cmd,1);
	
	chdir $RealBin or die "Can NOT chdir to $RealBin:$!";
}

sub predicthaplo_pipeline {
	my $predicthaplo_excu = shift;
	my $ref = shift;
	my $samfile = shift;
	my $outdir = shift;
	my $samplename = shift;
	
	chdir $outdir or die "Can NOT chdir to $outdir:$!";
	
	#get seq len
	my $seqlen = checkSeqLen($ref,$outdir);
	
	#prepare config file
	my $configfile = File::Spec -> catfile($outdir,"config". "_$samplename");
	open T,">$configfile" or die "Can NOT output to configure file $configfile:$!";
	
	print T "% configuration file for the Virushaplotyper\n";
	print T "% prefix\n";
	print T "${samplename}_\n";
	print T "% filename of reference sequence (FASTA)\n";
	print T "$ref\n";
	print T "% do_visualize (1 = true, 0 = false)\n";
	print T "1\n";
	print T "% filname of the aligned reads (sam format)\n";
	print T "$samfile\n";
	print T "% have_true_haplotypes  (1 = true, 0 = false)\n";
	print T "0\n";
	print T "% filname of the true haplotypes (MSA in FASTA format) (fill in any dummy filename if there is no \"true\" haplotypes)\n";
	print T "somedummy.fasta\n";
	print T "% do_local_analysis  (1 = true, 0 = false) (must be 1 in the first run)\n";
	print T "1\n";
	print T "% max_reads_in_window;\n";
	print T "10000\n";
	print T "% entropy_threshold\n";
	print T "4e-2\n";
	print T "%reconstruction_start\n";
	print T "1\n";
	print T "%reconstruction_stop\n";
	print T "$seqlen\n";
	print T "%min_mapping_qual\n";
	print T "30\n";
	print T "%min_readlength\n";
	print T "100\n";
	print T "%max_gap_fraction (relative to alignment length)\n";
	print T "0.05\n";
	print T "%min_align_score_fraction (relative to read length)\n";
	print T "0.35\n";
	print T "%alpha_MN_local (prior parameter for multinomial tables over the nucleotides)\n";
	print T "25\n";
	print T "%min_overlap_factor (reads must have an overlap with the local reconstruction window of at least this factor times the window size)\n";
	print T "0.85\n";
	print T "%local_window_size_factor (size of  local reconstruction window relative to the median of the read lengths)\n";
	print T "0.7\n";
	print T "% max number of clusters (in the truncated Dirichlet process)\n";
	print T "25\n";
	print T "% MCMC iterations\n";
	print T "501\n";
	print T "% include deletions (0 = no, 1 = yes)\n";
	print T "1\n";
	
	my $cmd = "$predicthaplo_excu $configfile";
	runcmd($cmd,1);

	chdir $RealBin or die "Can NOT chdir to $RealBin:$!";
}

sub qure_pipeline {
	my $qure_excu = shift;
	my $fastafile = shift;
	my $reffile = shift;
	my $outdir = shift;
	
	chdir $outdir or die "Can NOT chdir to $outdir:$!";
	
	my $cmd = "java -classpath $qure_excu -Xmx10g QuRe $fastafile $reffile";
	runcmd($cmd,1);
	
	chdir $RealBin or die "Can NOT chdir to $RealBin:$!";
}

sub viquas_pipeline {
	my $viquas_dir = shift;
	my $ref = shift;
	my $bamfile = shift;
	
	chdir $viquas_dir or die "Can NOT chdir to $viquas_dir:$!";
	
	my $cmd = "Rscript ViQuaS.R $ref $bamfile";
	runcmd($cmd,1);
	
	chdir $RealBin or die "Can NOT chdir to $RealBin:$!";
}

sub qsr_pipeline {
	my $amplicon = shift;
	my $fq1 = shift;
	my $fq2 = shift;
	my $fasta = shift;
	my $bamfile = shift;
	my $samfile = shift;
	my $ref = shift;
	my $outdir = shift;
	my $program = shift;
	my $threads = shift;
	my $shorah_excu = shift;
	my $predicthaplo_excu = shift;
	my $qure_excu_dir = shift;
	my $viquas_excu_dir = shift;
	
	$ref = abs_path($ref);
	$fq1 = abs_path($fq1) if $fq1 ne "null";
	$fq2 = abs_path($fq2) if $fq2 ne "null";
	$fasta = abs_path($fasta) if $fasta ne "null";
	$bamfile = abs_path($bamfile) if $bamfile ne "null";
	$samfile = abs_path($samfile) if $samfile ne "null";
	
	my $sampleName = getCommonString(basename($fq1),basename($fq2));
	$sampleName =~ s/[_\.R]+$//i;
	
	#start to run qsr
	if (uc($program) eq 'SHORAH'){
		my $shorah_dir = File::Spec -> catfile($outdir,'Shorah');
		makedir($shorah_dir);
		
		Info("Running Shorah for ECnQSR");
		shorah_pipeline($amplicon, $shorah_excu,$bamfile,$ref,$shorah_dir);
		
	}elsif(uc($program) eq 'PREDICTHAPLO'){
		my $predicthaplo_dir = File::Spec -> catfile($outdir,'PredictHaplo');
		makedir($predicthaplo_dir);
		
		Info("Running PredictHaplo for ECnQSR");
		predicthaplo_pipeline($predicthaplo_excu,$ref,$samfile,$predicthaplo_dir,$sampleName);
		
	}elsif(uc($program) eq 'QURE'){
		my $qure_dir = File::Spec -> catfile($outdir,'QuRe');
		makedir($qure_dir);
		
		##copy input file to out dir
		my $qurefastafile = File::Spec -> catfile($qure_dir,basename($fasta));
		copy($fasta,$qure_dir);
		my $qurereffile = File::Spec -> catfile($qure_dir,basename($ref));
		copy($ref,$qure_dir);
		
		Info("Running QuRe for ECnQSR");
		qure_pipeline($qure_excu_dir, $qurefastafile, $qurereffile, $qure_dir);
		
	}elsif(uc($program) eq 'VIQUAS'){
		my $viquas_dir = File::Spec -> catfile($outdir,'ViQuaS');
		makedir($viquas_dir);
		
		#copy input file to out dir 
		my $viquas_bam = File::Spec -> catfile($viquas_dir,basename($bamfile));
		copy($bamfile,$viquas_dir);
		my $viquas_ref = File::Spec -> catfile($viquas_dir,basename($ref));
		copy($ref,$viquas_dir);
		
		#copy source code to out dir
		system("cp -r $viquas_excu_dir/* $viquas_dir");
		Info("Running ViQuaS for ECnQSR");
		viquas_pipeline($viquas_dir,basename($viquas_ref),basename($viquas_bam));
		
	}else{
		InfoError("The mapping program MUST be among \"Shorah\", \"QuRe\", \"PredictHaplo\" or \"ViQuaS\" .","red");
		exit;
	}
}


1;
