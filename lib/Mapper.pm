#!/usr/bin/env perl

package Mapper;

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
use List::Util qw/min max/;

####Use modules in this program####
use General;

use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

our $VERSION     = 1.00;
our @ISA         = qw(Exporter);
our @EXPORT      = qw(Bowtie2_index Bowtie2_pipeline Bowtie2_pipeline_fasta BWA_index BWA_pipeline BWA_pipeline_fasta sam2SortedAndIndexedBam sortAndIndexBam 
isSamFile isBamFile bam2sam bam2newsam existBowtieIndex existBowtie2Index existBWAIndex blatPipeline blastPipeline);
our @EXPORT_OK   = qw(Bowtie2_index Bowtie2_pipeline Bowtie2_pipeline_fasta BWA_index BWA_pipeline BWA_pipeline_fasta sam2SortedAndIndexedBam sortAndIndexBam 
isSamFile isBamFile bam2sam bam2newsam existBowtieIndex existBowtie2Index existBWAIndex blatPipeline blastPipeline);
our %EXPORT_TAGS = ( DEFAULT => [qw(&Bowtie2_index &Bowtie2_pipeline &Bowtie2_pipeline_fasta &BWA_index &BWA_pipeline_fasta &BWA_pipeline &sam2SortedAndIndexedBam &sortAndIndexBam 
&isSamFile &isBamFile &bam2sam &bam2newsam &existBowtieIndex &existBowtie2Index &existBWAIndex &blatPipeline &blastPipeline)]);


##global variables
##get workding directory
my $wk_dir = getcwd;
my $mainBin;
if ($RealBin =~ /(.*)\/bin/){
	$mainBin = $1;
}

##subprogram goes here
sub Bowtie2_index {
	my($Bowtie2_excu,$ref) = @_;
	
	my $Bowtie2_dir = dirname($Bowtie2_excu);
	my $Bowtie2_build_excu = File::Spec -> catfile($Bowtie2_dir,'bowtie2-build');
	my $DEBUG_MODE = 1;
	#check bowtie2-build
	if(CheckProgram($Bowtie2_build_excu, __FILE__, __LINE__, $DEBUG_MODE)){
		#keep running
	}else{
		InfoError("The program $Bowtie2_build_excu does NOT exist. Exiting...");
		exit;
	}
	
	##fix \r in windows os
	my $sedcmd = "sed -i \'s\/\\r\/\/g\' $ref";
	system($sedcmd);
	
	#build index
	my $Bowtie2_index = removeFastaSuffix($ref);
	my $cmd = "$Bowtie2_build_excu $ref $Bowtie2_index";
	if(existBowtie2Index($Bowtie2_index)){
		#nothing
	}else{
		InfoPlain("Building Bowtie2 index for $ref.");
		runcmd($cmd);
	}
	
}

sub Bowtie2_pipeline{
	my($Bowtie2_excu,$ref,$fq1,$fq2,$sam,$threads) = @_;
	
	my $Bowtie2_dir = dirname($Bowtie2_excu);
	my $Bowtie2_build_excu = File::Spec -> catfile($Bowtie2_dir,'bowtie2-build');
	my $DEBUG_MODE = 1;
	#check bowtie2-build
	if(CheckProgram($Bowtie2_build_excu, __FILE__, __LINE__, $DEBUG_MODE)){
		#keep running
	}else{
		InfoError("The program $Bowtie2_build_excu does NOT exist. Exiting...");
		exit;
	}
	
	##fix \r in windows os
	my $sedcmd = "sed -i \'s\/\\r\/\/g\' $ref";
	system($sedcmd);
	
	#build index
	my $Bowtie2_index = removeFastaSuffix($ref);
	my $cmd = "$Bowtie2_build_excu $ref $Bowtie2_index";
	if(existBowtie2Index($Bowtie2_index)){
		#nothing
	}else{
		InfoPlain("Building Bowtie2 index for $ref.");
		runcmd($cmd);
	}
	
	#run bowtie2 mapping
	InfoPlain("Start to map reads to $ref using Bowtie2.");
	if (existFile($sam)){
		InfoWarn("$sam already exist,skip Bowtie2 mapping.");
	}else{
		if ($fq2){
			InfoPlain("Running Bowtie2 for paired-end read mapping");
			my $cmd = "$Bowtie2_excu --very-sensitive-local --threads $threads -x $Bowtie2_index -1 $fq1 -2 $fq2 -S $sam";
			runcmd($cmd);
		}else{
			InfoPlain("Running Bowtie2 for single-end read mapping");
			my $cmd = "$Bowtie2_excu --very-sensitive-local --threads $threads -x $Bowtie2_index -U $fq1 -S $sam";
			runcmd($cmd);
		}
	}
}

sub Bowtie2_pipeline_fasta{
	my($Bowtie2_excu,$ref,$fasta,$sam,$threads) = @_;
	
	my $Bowtie2_dir = dirname($Bowtie2_excu);
	my $Bowtie2_build_excu = File::Spec -> catfile($Bowtie2_dir,'bowtie2-build');
	my $DEBUG_MODE = 1;
	#check bowtie2-build
	if(CheckProgram($Bowtie2_build_excu, __FILE__, __LINE__, $DEBUG_MODE)){
		#keep running
	}else{
		InfoError("The program $Bowtie2_build_excu does NOT exist. Exiting...");
		exit;
	}
	
	##fix \r in windows os
	my $sedcmd = "sed -i \'s\/\\r\/\/g\' $ref";
	system($sedcmd);
	
	#build index
	my $Bowtie2_index = removeFastaSuffix($ref);
	my $cmd = "$Bowtie2_build_excu $ref $Bowtie2_index";
	if(existBowtie2Index($Bowtie2_index)){
		#nothing
	}else{
		InfoPlain("Building Bowtie2 index for $ref.");
		runcmd($cmd);
	}
	
	#run bowtie2 mapping
	InfoPlain("Start to map reads to $ref using Bowtie2.");
	if (existFile($sam)){
		InfoWarn("$sam already exist,skip Bowtie2 mapping.");
	}else{
		InfoPlain("Running Bowtie2 for mapping");
		my $cmd = "$Bowtie2_excu --very-sensitive-local -x $Bowtie2_index -f $fasta -S $sam";
		runcmd($cmd);
	}
}

sub BWA_index {
	my ($BWA_excu,$ref) = @_;
	
	##fix \r in windows os
	my $sedcmd = "sed -i \'s\/\\r\/\/g\' $ref";
	system($sedcmd);
	
	##build index file
	my $BWA_index = removeFastaSuffix($ref);
	my $refSize = getFileSize($ref);
	if ($refSize < 10000000){
		my $cmd = "$BWA_excu index -a is -p $BWA_index $ref";
		if (existBWAIndex($BWA_index)){
			#nothing
		}else{
			InfoPlain("Building BWA index for $ref.");
			runcmd($cmd);
		}
	}else{
		my $cmd = "$BWA_excu index -a bwtsw -p $BWA_index $ref";
		if (existBWAIndex($BWA_index)){
			#nothing
		}else{
			runcmd($cmd);
		}
	}
}

sub BWA_pipeline {
	my ($BWA_excu,$ref,$fq1,$fq2,$sam,$threads) = @_;
	
	##fix \r in windows os
	my $sedcmd = "sed -i \'s\/\\r\/\/g\' $ref";
	system($sedcmd);
	
	##build index file
	my $BWA_index = removeFastaSuffix($ref);
	my $refSize = getFileSize($ref);
	if ($refSize < 10000000){
		my $cmd = "$BWA_excu index -a is -p $BWA_index $ref";
		if (existBWAIndex($BWA_index)){
			#nothing
		}else{
			InfoPlain("Building BWA index for $ref.");
			runcmd($cmd);
		}
	}else{
		my $cmd = "$BWA_excu index -a bwtsw -p $BWA_index $ref";
		if (existBWAIndex($BWA_index)){
			#nothing
		}else{
			runcmd($cmd);
		}
	}
	
	##start to map reads to ref 
	InfoPlain("Start to map reads to $ref using BWA");
	if($fq2){
		InfoPlain("Running BWA for paired-end read mapping.");
		##mem
		if (existFile($sam)){
			InfoPlain("$sam already exist,skip BWA mem mapping.");
		}else{
			my $cmd = "$BWA_excu mem -t $threads $BWA_index $fq1 $fq2 > $sam";
			runcmd($cmd);
		}
	}else{
		InfoPlain("Running BWA for single-end read mapping");
		##mem
		if (existFile($sam)){
			InfoPlain("$sam already exist,skip BWA mem mapping.");
		}else{
			my $cmd = "$BWA_excu mem -t $threads $BWA_index $fq1 > $sam";
			runcmd($cmd);
		}
	}
}

sub BWA_pipeline_fasta {
	my ($BWA_excu,$ref,$fasta,$sam,$threads) = @_;
	
	##fix \r in windows os
	my $sedcmd = "sed -i \'s\/\\r\/\/g\' $ref";
	system($sedcmd);
	
	##build index file
	my $BWA_index = removeFastaSuffix($ref);
	my $refSize = getFileSize($ref);
	if ($refSize < 10000000){
		my $cmd = "$BWA_excu index -a is -p $BWA_index $ref";
		if (existBWAIndex($BWA_index)){
			#nothing
		}else{
			InfoPlain("Building BWA index for $ref.");
			runcmd($cmd);
		}
	}else{
		my $cmd = "$BWA_excu index -a bwtsw -p $BWA_index $ref";
		if (existBWAIndex($BWA_index)){
			#nothing
		}else{
			runcmd($cmd);
		}
	}
	
	##start to map reads to ref 
	InfoPlain("Start to map reads to $ref using BWA");
	##mem
	if (existFile($sam)){
		InfoPlain("$sam already exist,skip BWA mem mapping.");
	}else{
		my $cmd = "$BWA_excu mem -t $threads $BWA_index $fasta > $sam";
		runcmd($cmd);
	}
}


sub sam2SortedAndIndexedBam {
	my ($samtools_excu, $samfile, $sortMethod, $clear_sam, $threads) = @_;
	
	my $DEBUG_MODE = 1;
	#check samtools
	if(CheckProgram($samtools_excu, __FILE__, __LINE__, $DEBUG_MODE)){
		#keep running
	}else{
		InfoError("The program $samtools_excu does NOT exist. Exiting...");
		exit;
	}
	
	my $samfileName = basename($samfile);
	
	InfoPlain("Converting $samfile to bam file");
	my $bamfile_unsort = $samfile =~ s/\.sam$/.bam/r if isSamFile($samfile);
	if (existFile($bamfile_unsort)){
		InfoPlain("$bamfile_unsort already exist, skipping convert sam2bam");
	}else{
		my $cmd1 = "$samtools_excu view -h -@ $threads -Sb $samfile > $bamfile_unsort";
		runcmd($cmd1);
	}
	
	my $bamfile;
	if (uc($sortMethod) eq 'POS'){
		InfoPlain("Sorting bam file $bamfile_unsort with coordinates.");
		#samtools sort -@ 30 -O bam -T bwa -o bwa.sort.bam bwa.bam
		my $DateNow = time();
		
		my $prefix = $samfileName =~ s/\.sam$/.PosSorted.${DateNow}/r;
		$bamfile = $samfile =~ s/\.sam$/.PosSorted.bam/r;
		if (existFile($bamfile)){
			InfoPlain("$bamfile already exists, skipping samtools sort.");
		}else{
			my $cmd2 = "$samtools_excu sort -@ $threads -O bam -T $prefix -o $bamfile $bamfile_unsort";
			runcmd($cmd2);
		}
	}else{
		InfoPlain("Sorting bam file $bamfile_unsort with reads name.");
		#samtools sort -@ 30 -O bam -T bwa -o bwa.sort.bam bwa.bam
		my $DateNow = time();
		
		my $prefix = $samfileName =~ s/\.sam$/.NameSorted.${DateNow}/r;
		$bamfile = $samfile =~ s/\.sam$/.NameSorted.bam/r;
		if (existFile($bamfile)){
			InfoPlain("$bamfile already exists, skipping samtools sort.");
		}else{
			my $cmd2 = "$samtools_excu sort -n -@ $threads -O bam -T $prefix -o $bamfile $bamfile_unsort";
			runcmd($cmd2);
		}
	}
	
	InfoPlain("Indexing sorted bam file.");
	if (existFile("${bamfile}.bai")){
		InfoPlain("The index of $bamfile already exists, skipping \"samtools index\".");
	}else{
		if (uc($sortMethod) eq 'POS'){
			my $cmd = "$samtools_excu index $bamfile";
			runcmd($cmd);
		}else{
			Info("$bamfile is sorted by name, so skipping \"samtools index\". ");
		}
	}
	
	my $cmd_clear = "rm -rf $samfile $bamfile_unsort";
	if ($clear_sam){
		system($cmd_clear);
	}else{
		#nothing;
	}
}


sub sortAndIndexBam {
	my ($samtools_excu, $bamfile, $sortMethod, $clear, $threads) = @_;
	
	my $DEBUG_MODE = 1;
	#check samtools
	if(CheckProgram($samtools_excu, __FILE__, __LINE__, $DEBUG_MODE)){
		#keep running
	}else{
		InfoError("The program $samtools_excu does NOT exist. Exiting...");
		exit;
	}
	
	my $bamfile_sort;
	my $bamfileName = basename($bamfile);
	
	if (uc($sortMethod) eq 'POS'){
		$bamfile_sort = $bamfile =~ s/\.bam$/.PosSorted.bam/r;
		if(existFile($bamfile_sort)){
			InfoPlain("Sorted bam file $bamfile_sort already exist, skipping \"samtools sort\".");
		}else{
			InfoPlain("Sorting bam file $bamfile with coordinates.");
			#samtools sort -@ 30 -O bam -T bwa -o bwa.sort.bam bwa.bam
			my $DateNow = time();
			
			my $prefix = $bamfileName =~ s/\.bam$/.PosSorted.${DateNow}/r;
			my $cmd = "$samtools_excu sort -@ $threads -O bam -T $prefix -o $bamfile_sort $bamfile";
			runcmd($cmd);
		}
		
	}else{
		$bamfile_sort = $bamfile =~ s/\.bam$/.NameSorted.bam/r;
		if(existFile($bamfile_sort)){
			InfoPlain("Sorted bam file $bamfile_sort already exist, skipping \"samtools sort\".");
		}else{
			InfoPlain("Sorting bam file $bamfile with read names.");
			#samtools sort -n -@ 30 -O bam -T bwa -o bwa.sort.bam bwa.bam
			my $DateNow = time();
			
			my $prefix = $bamfileName =~ s/\.bam$/.NameSorted.${DateNow}/r;
			my $cmd = "$samtools_excu sort -n -@ $threads -O bam -T $prefix -o $bamfile_sort $bamfile";
			runcmd($cmd);
		}
	}
	

	InfoPlain("Indexing sorted bam file.");
	if (existFile("${bamfile_sort}.bai")){
		InfoPlain("The index of $bamfile_sort already exists, skipping \"samtools index\".");
	}else{
		if (uc($sortMethod) eq 'POS'){
			my $cmd = "$samtools_excu index $bamfile_sort";
			runcmd($cmd);
		}else{
			Info("$bamfile_sort is sorted by name, so skipping \"samtools index\". ");
		}
	}
	
	my $cmd_clear = "rm -rf $bamfile";
	if ($clear){
		system($cmd_clear);
	}else{
		#nothing;
	}
}

sub isSamFile {
	my $file = shift;
	
	my $filename = basename($file);
	#if ($filename =~ /\.sam$/){
	#	return 1;
	#}else{
	#	InfoWarn("The suffix of sam file $file should be \".sam\". ");
	#	return 0;
	#}
	
	my $fileInfo = `file $file`;
	chomp $fileInfo;
	
	if ($fileInfo =~ /ASCII text/i and $filename =~ /\.sam$/){
		return 1;
	}
	
	return 0;
	
}

sub isBamFile {
	my $file = shift;
	
	my $filename = basename($file);
	#if ($filename =~ /\.bam$/){
	#	return 1;
	#}else{
	#	InfoWarn("The suffix of bam file $file should be \".bam\". ");
	#	return 0;
	#}
	
	my $fileInfo = `file $file`;
	chomp $fileInfo;
	
	if ($fileInfo =~ /compressed data/i and $filename =~ /\.bam$/){
		return 1;
	}
	
	return 0;
	
	
}

sub bam2sam {
	my $samtools_excu = shift;
	my $bamfile = shift;
	my $threads = shift;
	
	my $DEBUG_MODE = 1;
	#check samtools
	if(CheckProgram($samtools_excu, __FILE__, __LINE__, $DEBUG_MODE)){
		#keep running
	}else{
		InfoError("The program $samtools_excu does NOT exist. Exiting...");
		exit;
	}
	
	my $samfile = $bamfile =~ s/\.bam$/.sam/r if isBamFile($bamfile);
	
	Info("Converting bam file $bamfile to $samfile");
	my $cmd = "$samtools_excu view -h -@ $threads $bamfile > $samfile";
	runcmd($cmd);
}

sub bam2newsam {
	my $samtools_excu = shift;
	my $bamfile = shift;
	my $samfile = shift;
	my $threads = shift;
	
	my $DEBUG_MODE = 1;
	#check samtools
	if(CheckProgram($samtools_excu, __FILE__, __LINE__, $DEBUG_MODE)){
		#keep running
	}else{
		InfoError("The program $samtools_excu does NOT exist. Exiting...");
		exit;
	}
	
	Info("Converting bam file $bamfile to $samfile");
	my $cmd = "$samtools_excu view -h -@ $threads $bamfile > $samfile";
	runcmd($cmd);
}

sub existBowtieIndex {
	my $index = shift;

	my $exist = 1;
	
	my $file1 = $index . ".1.ebwt";
	my $file2 = $index . ".2.ebwt";
	my $file3 = $index . ".3.ebwt";
	my $file4 = $index . ".4.ebwt";
	my $file5 = $index . ".rev.1.ebwt";
	my $file6 = $index . ".rev.2.ebwt";
	
	my @allfiles = ($file1, $file2, $file3, $file4, $file5, $file6);
	
	for my $f (@allfiles){
		if (!existFile($f)){
			$exist = 0;
			last;
		}
	}	

	return $exist;
}

sub existBowtie2Index {
	my $index = shift;

	my $exist = 1;
	
	my $file1 = $index . ".1.bt2";
	my $file2 = $index . ".2.bt2";
	my $file3 = $index . ".3.bt2";
	my $file4 = $index . ".4.bt2";
	my $file5 = $index . ".rev.1.bt2";
	my $file6 = $index . ".rev.2.bt2";
	
	my @allfiles = ($file1, $file2, $file3, $file4, $file5, $file6);
	
	for my $f (@allfiles){
		if (!existFile($f)){
			$exist = 0;
			last;
		}
	}	

	return $exist;
}

sub existBWAIndex {
	my $index = shift;
	
	my $exist = 1;
	
	my $file1 = $index . ".amb";
	my $file2 = $index . ".ann";
	my $file3 = $index . ".bwt";
	my $file4 = $index . ".pac";
	my $file5 = $index . ".sa";
	
	my @allfiles = ($file1, $file2, $file3, $file4, $file5);
	
	for my $f (@allfiles){
		if (! existFile($f)){
			$exist = 0;
			last;
		}
	}
	
	return $exist;
}

sub blatPipeline {
	my $blatProgram = shift;
	my $input = shift;
	my $ref = shift;
	my $outputdir = shift;
	
	#run blat
	Info("Locating the position of sequences using BLAT");
	my $psl = File::Spec -> catfile($outputdir,removeFastaSuffix(basename($input)) . ".psl");
	my $cmd = "$blatProgram $ref $input $psl";
	system($cmd);
	
	#parse blat output
	open T,"$psl" or die "Can NOT open blat output $psl:$!";
	my $dump = <T>;
	$dump = <T>;
	$dump = <T>;
	$dump = <T>;
	$dump = <T>;
	
	my @start;
	my @end;
	while(<T>){
		chomp;
		my @tmp = split "\t",$_;
		my $start = $tmp[15] + 1;
		my $query_size = $tmp[10];
		#my $end = $tmp[16];
		my $end = $start + $query_size - 1;
		
		push @start,$start;
		push @end,$end;
	}
	
	my $outstart = getEleWithMaxFreq(\@start);
	my $outend = getEleWithMaxFreq(\@end);
	
	return($psl,$outstart,$outend);
}

sub blastPipeline {
	my $makeblastdb_excu = shift;
	my $blastn_excu = shift;
	my $ref = shift;
	my $query = shift;
	my $output = shift;
	my $outfmt6 = shift;
	my $threads = shift;
	
	Info("Running Blast for $query");
	#makedb
	my $cmd = "$makeblastdb_excu -in $ref -dbtype nucl -parse_seqids";
	runcmd($cmd);
	
	#blastn
	if($outfmt6){
		$cmd = "$blastn_excu -db $ref -query $query -outfmt 6 -out $output -num_threads $threads";
	}else{
		$cmd = "$blastn_excu -db $ref -query $query -out $output -num_threads $threads";
	}
	runcmd($cmd);
	
	return 1;
}


1;

