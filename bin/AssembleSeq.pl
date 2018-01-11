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
use List::Util qw/max min/;
use File::Copy;

####Use modules in this program####
use General;
use Mapper;

####Flush cache
$| = 1;

####Set random seed
srand( time() ^ ( $$ + ( $$<<15 ) ) );

####---------------------------####
####The program begins here
####---------------------------####

##Show welcome
print "You are now running subprogram: ";
printcol ("AssembleSeq","green");
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
my $inputFile;
my $outputDir;
my $threads;
my $bedFile;
my $ampliconNumber;
my $insertion;

my $DateNow = `date +"%Y%m%d_%Hh%Mm%Ss"`;
chomp $DateNow;

GetOptions(
'i|inputFile|=s'         => \$inputFile,
'o|outputDir|=s'         => \$outputDir,
'b|bedFile|=s'           => \$bedFile,
't|threads|=s'           => \$threads,
'n|ampliconNumber|=s'    => \$ampliconNumber,
's|withIns|=s'           => \$insertion,
'h|help|'                => \$help
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
	}else{
		InfoError("The output directory $outputDir already exist.");
		InfoError("Please specify another output directory using option -o/--outputDir");
		pod2usage(-verbose=>0,-exitval=>1);
		exit;
	}
}else{
	$outputDir = File::Spec -> catfile($wk_dir,"qap_Results_for_AssembleSeq_$DateNow");
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

if(defined $inputFile){
	if (not -e $inputFile){
		InfoError("Input file $inputFile does NOT exist! Please check again.");
		exit;
	}
}else{
	InfoError("Input SAM/BAM file MUST be specified with -i/--inputDir\n");
	pod2usage(-verbose=>0,-exitval=>1);
	exit;
}

if(defined $bedFile){
	if (not -e $bedFile){
		InfoError("Input bed file $bedFile does NOT exist! Please check again.");
		exit;
	}
}else{
	if(not defined $ampliconNumber){
		InfoError("Either --bedFile or --ampliconNumber should be specified.");
		pod2usage(-verbose=>0,-exitval=>1);
		exit;
	}
}

if(defined $ampliconNumber){
	if(CheckPositiveInt($ampliconNumber)){
		#nothing
	}else{
		InfoError("The number of amplicons MUST be positive integer.");
		pod2usage(-verbose=>0,-exitval=>1);
		exit;
	}
}else{
	if(not defined $bedFile){
		InfoError("Either --bedFile or --ampliconNumber should be specified.");
		pod2usage(-verbose=>0,-exitval=>1);
		exit;
	}
}

if(defined $insertion){
	if ($insertion =~ /^1/){
		$insertion = 1;
	}elsif($insertion =~ /^2/){
		$insertion = 2;
	}else{
		InfoError("--withIns/-s MUST be specified with \'1\' or \'2\'.");
		pod2usage(-verbose=>0,-exitval=>1);
		exit;
	}
}else{
	$insertion = 2;
	InfoWarn("--withIns/-s is not provided. Choose \'-s 2\' as default.")
}

##core program starts here


#check samtools
my $samtoolsProgram = File::Spec -> catfile($RealBin,'3rdPartyTools','samtools','samtools');
if(CheckProgram($samtoolsProgram, __FILE__, __LINE__, $DEBUG_MODE)){
	#keep running
}else{
	InfoError("The program $samtoolsProgram does NOT exist. Exiting...");
	exit;
}

#check file format and handle files
my $samfile_sort;
system("cp -r $inputFile $outputDir");
$inputFile = File::Spec -> catfile($outputDir, basename($inputFile));
if(isBamFile($inputFile)){
	sortAndIndexBam($samtoolsProgram, $inputFile, 'NAME', 0, $threads);
	my $bamfile_sort = $inputFile =~ s/\.bam$/.NameSorted.bam/r;
	
	bam2sam($samtoolsProgram, $bamfile_sort, $threads);
	$samfile_sort = $bamfile_sort =~ s/\.bam$/.sam/r;

}elsif(isSamFile($inputFile)){
	sam2SortedAndIndexedBam($samtoolsProgram, $inputFile, 'NAME', 0, $threads);
	my $bamfile_sort = $inputFile =~ s/\.sam$/.NameSorted.bam/r;
	
	bam2sam($samtoolsProgram, $bamfile_sort, $threads);
	$samfile_sort = $bamfile_sort =~ s/\.bam$/.sam/r;
	
}else{
	InfoError("$inputFile MUST be sam or bam file, no matter sorted or not.");
	exit;
}


#parse sorted sam file and get reference length
open T,$samfile_sort or die "Can not read inputfile $samfile_sort:$!";
my $genomeLen;
my @genomeLen;
my $numberOfReferences = 0;
while(<T>){
	chomp;
	# @SQ SN:refC_AB014381    LN:3215
	if(/\@SQ\s+SN.*?LN\:(\d+)/){
		my $tmp = $1;
		push @genomeLen,$tmp;
		$numberOfReferences++;
	}
}
close T;

if(scalar(@genomeLen) > 1){
	#this is a virus with several genome segments
}elsif(scalar(@genomeLen) == 1){
	#this virus only have one segment or this is a virus with circular DNA/RNA
	$genomeLen = $genomeLen[0];
}else{
	#maybe there is no headline in this sam file
	InfoWarn("The head information of SAM file $samfile_sort is missing.");
}


#parse bed file
my @chr;
my @startPos;
my @endPos;
my @ampliconLens;
my @ampliconNames;
if(defined $bedFile){
	open BED,$bedFile or die "Can NOT open file $bedFile:$!";
	my $i = 1;
	while (<BED>){
		my ($chr,$start,$end,$name) = split(/\s+/,$_); #split both \t and space
		push @chr,$chr;
		push @startPos,$start;
		push @endPos,$end;
		push @ampliconNames,$name;
		if($start > $end){
			if(defined $genomeLen){
				my $len = ($genomeLen - $start + 1) + ($end - 1 + 1);
				push @ampliconLens,$len; 
			}else{
				push @ampliconLens,'NA';
			}
		}elsif($start < $end){
			my $len = $end - $start + 1;
			push @ampliconLens,$len;
		}else{
			InfoWarn("BED file LINE $i: Start position equal end position. Please Check.");
			push @ampliconLens,1;
		}
		
		$i++;
	}
	close BED;
	
	if(defined $ampliconNumber){
		if(scalar(@ampliconNames) == $ampliconNumber){
			#nothing
		}else{
			InfoError("The number of amplicons you provided in BED file and CMD argument are unequal. Please check.");
		}
	}else{
		$ampliconNumber = scalar(@ampliconNames);
	}
}


#start to handle the sam file
Info("Start to assemble sequences from $inputFile. Please wait...");
#cut head line
my $samfile_sort_nohead = $samfile_sort =~ s/\.sam$/.noHead.sam/r;
open NH,">$samfile_sort_nohead" or die "Can not output to file $samfile_sort_nohead:$!";

open T,$samfile_sort or die "Can not open file $samfile_sort:$!";
while(<T>){
	chomp;
	if (/^\@/){
		next;
	}else{
		print NH "$_\n";
	}
} 
close T;
close NH;

#start to assemble
open T,$samfile_sort_nohead or die "Can not read file $samfile_sort_nohead:$!";
my $i = 1;
my $k = 2;
my $numberOfWrongID = 0;
my $numberOfAssembledRP = 0;
my $numberOfTotalRP = 0;
my $numberOfPassFiltrationRP = 0;
while (my $line1 = <T>){
	chomp $line1;
	my $FH_pos1 = tell(T);
	chomp (my $line2 = <T>);
	my $FH_pos2 = tell(T);
	
	$numberOfTotalRP++;
	
	my @arr1 = split "\t",$line1;
	my @arr2 = split "\t",$line2;
	my ($id1,$ref1,$pos1,$matExp1,$interval1,$seq1,$qual1) = @arr1[0,2,3,5,8,9,10];
	my ($id2,$ref2,$pos2,$matExp2,$interval2,$seq2,$qual2) = @arr2[0,2,3,5,8,9,10];
	
	#check whether input is right
	if($numberOfWrongID >= $numberOfTotalRP * 0.1){
		InfoError("Something MUST be wrong with the SAM file $samfile_sort_nohead. Maybe the file is not sorted by read name. Please check.");
		exit(0);
	}
	
	#check if $id1 eq $id2
	if($id1 eq $id2){
		# add the LINE count, nothing else to do
		$i += 2;
		$k += 2;
	}else{
		InfoWarn("The sequence ID in LINE $i and LINE $k are different. Please Check.");
		$numberOfWrongID++;
		
		# set find handle read position to previous line
		seek(T,$FH_pos1,0);
		
		# add the LINE count
		$i++;
		$k++;
		
		next;
	}
	
	#filter
	#'N' means skipped region from the reference. As viruses do not have introns, so seq with 'N' is excluded.
	if($matExp1 =~ /N/i or $matExp2 =~ /N/i){
		next;
	}
	#'H' means hard clipped reads which should be excluded
	if($matExp1 =~ /H/i or $matExp2 =~ /H/i){
		next;
	}
	#interval TLEN filtration
	if ($interval1 != $interval2 * -1){
		next;
	}
	if ($interval1 == 0 or $interval2 == 0){
		next;
	}
	if(defined $bedFile){
		next if abs($interval1) > max(@ampliconLens); #no need to check $interval2, because $interval1 eq $interval2
		next if abs($interval1) < &getMin(\@ampliconLens); #there may be 'NA' in @ampliconLens, so get the real min number with &getmin
	}
	##map position filtration
	#if(abs($pos1 - $pos2) < 100){ #the minimum length of illumina reads is 75, so set cutoff to 100
	#	next;
	#}
	#ref filteration
	if($ref1 ne $ref2){
		next;
	}
	
	#if --withIns is false, discard sequences with insertions.
	if ($insertion == 1){
		next if $matExp1 =~ /I/ or $matExp2 =~ /I/;
	}
	
	$numberOfPassFiltrationRP++;
	
	#cut seq by using CIGAR
	my ($seq1Cut,$qual1Cut,$seq1Count) = processSeqUsingCIGAR($matExp1,$seq1,$qual1);
	my ($seq2Cut,$qual2Cut,$seq2Count) = processSeqUsingCIGAR($matExp2,$seq2,$qual2);
	
	#discard sequences without overlap
	next if abs($interval1) > length($seq1Cut) + length($seq2Cut);

	
	#set output file ready
	my $outfilename = removeSamBamSuffix(basename($inputFile));
	my $outfile;
	if($numberOfReferences == 1){
		$outfile = File::Spec -> catfile($outputDir, ${outfilename} . "_" . min($pos1,$pos2) .".fasta");
	}elsif($numberOfReferences > 1){
		$outfile = File::Spec -> catfile($outputDir, ${outfilename} . "_" . $ref1 . "_" . min($pos1,$pos2) .".fasta");
	}else{
		InfoWarn("The head information of SAM file $inputFile is missing.");
	}
	
	
	open RES,">>$outfile" or die "Can NOT output to file $outfile:$!";
		
	#start to assemle
	if($interval1 > 0 and $interval2 < 0){
		my $seqCountToCut = 0;
		if($matExp1 =~ /(\d+)S$/){
			my $tmp = $1;
			$seqCountToCut += $tmp;
		}
		if($matExp2 =~ /^(\d+)S/){
			my $tmp = $1;
			$seqCountToCut += $tmp;
		}
		my $overlapLength = $seq1Count + $seq2Count - $seqCountToCut - abs($interval1);
		my $contig1 = substr $seq1Cut, 0, length($seq1Cut) - $overlapLength; #unoverlapped region of seq1
		if($overlapLength > length($seq2Cut)){
			next;
		}
		my $contig2 = substr $seq2Cut, $overlapLength, length($seq2Cut) - $overlapLength; #unonverlapped region of seq2
		
		my $seq1Overlap = substr $seq1Cut, length($seq1Cut) - $overlapLength, $overlapLength;
		my $seq2Overlap = substr $seq2Cut, 0, $overlapLength;
		my $qual1Overlap = substr $qual1Cut, length($qual1Cut) - $overlapLength, $overlapLength;
		my $qual2Overlap = substr $qual2Cut, 0, $overlapLength;
		
		#get overlap region 
		my $overlapSeq = &assembleSeqBaseOnQual($id1, $seq1Overlap, $seq2Overlap, $qual1Overlap, $qual2Overlap);
		next if (not defined $overlapSeq);
		
		my $assembledSeq = $contig1 . $overlapSeq . $contig2;
		#print "$seq1Overlap\n$seq2Overlap\n\n";
		print RES ">$id1\n$assembledSeq\n";
		
		$numberOfAssembledRP++;
	}elsif($interval1 < 0 and $interval2 > 0){
		my $seqCountToCut = 0;
		if($matExp2 =~ /(\d+)S$/){
			my $tmp = $1;
			$seqCountToCut += $tmp;
		}
		if($matExp1 =~ /^(\d+)S/){
			my $tmp = $1;
			$seqCountToCut += $tmp;
		}
		my $overlapLength = $seq2Count + $seq1Count - $seqCountToCut - abs($interval1);
		my $contig2 = substr $seq2Cut, 0, length($seq2Cut) - $overlapLength; #unoverlapped region of seq1
		if($overlapLength > length($seq1Cut)){
			next;
		}
		my $contig1 = substr $seq1Cut, $overlapLength, length($seq1Cut) - $overlapLength; #unonverlapped region of seq2
		
		my $seq2Overlap = substr $seq2Cut, length($seq2Cut) - $overlapLength, $overlapLength;
		my $seq1Overlap = substr $seq1Cut, 0, $overlapLength;
		my $qual2Overlap = substr $qual2Cut, length($qual2Cut) - $overlapLength, $overlapLength;
		my $qual1Overlap = substr $qual1Cut, 0, $overlapLength;
		
		#get overlap region 
		my $overlapSeq = &assembleSeqBaseOnQual($id1, $seq1Overlap, $seq2Overlap, $qual1Overlap, $qual2Overlap);
		next if (not defined $overlapSeq);
		
		my $assembledSeq = $contig2 . $overlapSeq . $contig1;
		#print "$seq2Overlap\n$seq1Overlap\n\n";
		print RES ">$id1\n$assembledSeq\n";
		
		$numberOfAssembledRP++;
	}else{
		#InfoWarn("The absolute value of TLEN in SAM file (column 9) at line $i and line $k should be equal.");
		#print "$id1\t $interval1\t $id2\t $interval2\n";
	}
	
	
	#handle file handle
	close RES;
}
close T;


##output the final statistic information
Info("Total $numberOfTotalRP read paires were handled of which $numberOfPassFiltrationRP / $numberOfTotalRP passed filtration, and $numberOfAssembledRP / $numberOfTotalRP read paires were assembled.");

##select the right mapped fasta file
my %fileSize;
for my $f (glob "$outputDir/*.fasta"){
	my @args = stat($f);
	my $size = $args[7];
	$fileSize{$f} = $size;
}

my @filesSortedBySize = sort {$fileSize{$b} <=> $fileSize{$a}} keys %fileSize; #sort the fasta files by file size from large to small
my @wantedFiles = @filesSortedBySize[0..$ampliconNumber - 1];


#mkdir tmp data file
my $tmpDir = File::Spec -> catfile($outputDir, 'tmpData');
mkdir $tmpDir or die "Can NOT mkdir $tmpDir:$!";

#move all data to tmp folder
for my $f (glob "$outputDir/*.*"){
	#system("mv $f $tmpDir");
	move($f,$tmpDir);
}

#mkdir tmp2 data dir
my $tmpDir2 = File::Spec -> catfile($tmpDir, 'tmp');
makedir($tmpDir2);

#get wanted files to tmp dir
for my $f (@wantedFiles){
	my $tmpFilePath = File::Spec -> catfile(dirname($f), 'tmpData', basename($f));
	my $tmpFilePath2 = File::Spec -> catfile(dirname($f), 'tmpData', 'tmp', basename($f));
	#system("cp $tmpFilePath $outputDir");
	copy($tmpFilePath, $tmpFilePath2);
	
	my $len = getMostSeqLen($tmpFilePath2, $tmpDir2);
	my $tmpFilePath2WithLen = removeFastaSuffix($tmpFilePath2) . ".lenFilter.fasta";
	extractSeqWithLen($tmpFilePath2, $len, $tmpFilePath2WithLen);
	my $csFile = removeFastaSuffix($tmpFilePath2WithLen) . ".cs.fasta";
	my $cs = &getConsensusSeq($tmpFilePath2WithLen, $csFile);
	
	&formatToEqualLen($tmpFilePath2, $f, $len, $cs);
}


##sub program starts here
sub getConsensusSeq {
	my $inputfile = shift;
	my $outputfile = shift;
	
	my $rinputName = removeFastaSuffix(basename($inputfile)) . ".RInput";
	my $rinput = File::Spec -> catfile(dirname($outputfile), $rinputName);
	open R,">$rinput" or die "Can NOT output to $rinput:$!\n";
	
	open T,"$inputfile" or die "Can NOT open $inputfile:$!\n";
	while(my $line1 = <T>){
		chomp $line1;
		chomp (my $line2 = <T>);
		
		print R "$line2\n";
	}
	close T;
	close R;
	
	#get cs
	#check rscript
	my $rscript = File::Spec -> catfile($mainBin,'bin','Rscripts','CalculateCS.R');
	if (existFile($rscript)){
		#nothing
	}else{
		InfoError("R script $rscript is missing. Please check.");
		InfoError("Aborting...");
		exit(0);
	}
	my $fileLabel = removeFastaSuffix(basename($inputfile));
	my $cmd = "Rscript $rscript -i $rinput -o $outputfile -l $fileLabel";
	system($cmd);
	
	open C,"$outputfile" or die "Can not open $outputfile:$!";
	my $dump = <C>;
	chomp (my $cs = <C>);
	close C;
	
	return $cs;
}

sub formatToEqualLen{
	my $file = shift;
	my $outfile = shift;
	my $len = shift;
	my $cs = shift;
	
	open T,$file or die "Can not open file $file:$!";
	open OUT,">$outfile" or die "Can not output to file $outfile:$!";
	
	while(my $line1 = <T>){
		chomp $line1;
		chomp (my $line2 = <T>);
		
		if (length($line2) == $len){
			print OUT "$line1\n$line2\n";
		}elsif(length($line2) > $len){
			my $outseq = substr $line2,0,$len;
			print OUT "$line1\n$outseq\n";
		}elsif(length($line2) < $len){
			my $outseq = $line2 . substr($cs,length($line2),$len - length($line2));
			print OUT "$line1\n$outseq\n";
		}else{
			InfoError("Something wrong with sequence [$line1] in $file");
		}
		
	}
	close T;
	close OUT;
	
}

sub getMin {
	my $arr = shift;
	
	my @arr = @$arr;
	my @newArrWithoutChar;
	for my $i (@arr){
		if ($i =~ /[a-z]/i){
			next;
		}else{
			push @newArrWithoutChar,$i;
		}
	}
	
	my $min = min(@newArrWithoutChar);
	
	return $min;
}

sub processSeqUsingCIGAR {
	my $matExp = shift;
	my $seq = shift;
	my $qual = shift;
	
	my $numberOfUnRecognizeCIGAR = 0;
	my $seqOut;
	my $qualOut;
	my $countOut;
	while ($matExp =~ /^(\d+)([MSID])/){
		my $num = $1;
		my $mat = $2;
		
		if($mat eq 'M'){
			my $cutSeq = substr $seq,0,$num;
			$seqOut .= $cutSeq;
			
			my $cutQual = substr $qual,0,$num;
			$qualOut .= $cutQual;
			
			substr $seq,0,$num, '';
			substr $qual,0,$num, '';
			
			$countOut += $num;
		}elsif($mat eq 'S'){
			substr $seq,0,$num,'';
			substr $qual,0,$num,'';	
			
			$countOut += $num;
		}elsif($mat eq 'I'){
			my $cutSeq = substr $seq,0,$num;
			#$seqOut .= $cutSeq; #do NOT concentrate insertion sites to the output sequence
			
			my $cutQual = substr $qual,0,$num;
			#$qualOut .= $cutQual; #do NOT concentrate insertion sites to the output quality 
			
			substr $seq,0,$num, '';
			substr $qual,0,$num, '';
			
		}elsif($mat eq 'D'){
			$seqOut .= '-' x $num;
			
			$qualOut .= '+' x $num; #phred quality 10
			
			$countOut += $num;
		}else{
			$numberOfUnRecognizeCIGAR++;
		}
		
		$matExp =~ s/^${num}${mat}//;
		
		if ($numberOfUnRecognizeCIGAR > 30){
			InfoWarn("There may be something wrong with CIGAR expression of the SAM file. Please check.");
		}
	}
	
	#print "$seqOut\n$qualOut\n";
	return($seqOut,$qualOut,$countOut);
}

sub assembleSeqBaseOnQual {
	my $id = shift;
	my $seq1 = shift;
	my $seq2 = shift;
	my $qual1 = shift;
	my $qual2 = shift;
	
	#check length
	my @len = (length($seq1), length($seq2), length($qual1), length($qual2));
	if (checkMultipleEqualNums(\@len)){
		#nothing
	}else{
		#InfoError("The length of overlap region of assembled read pair \[${id}\] are different.");
		my $len1 = length($seq1);
		my $len2 = length($seq2);
		my $len3 = length($qual1);
		my $len4 = length($qual2);
		#print "$id\t$len1\t$len2\t$len3\t$len4\n";
		return;
	}
	
	#cut into array
	my @seq1 = split '',$seq1;
	my @seq2 = split '',$seq2;
	my @qual1 = split '',$qual1;
	my @qual2 = split '',$qual2;
	
	my $outSeq;
	for my $i (0..scalar(@seq1) - 1){
		my $seq1Char = $seq1[$i];
		my $seq2Char = $seq2[$i];
		
		my $qual1Char = $qual1[$i];
		my $qual2Char = $qual2[$i];
		
		#convert ASCII quality character to integer
		my $qual1Num = ord($qual1Char);
		my $qual2Num = ord($qual2Char);
		
		if ($qual1Num > $qual2Num){
			$outSeq .= $seq1Char;
		}elsif($qual1Num < $qual2Num){
			$outSeq .= $seq2Char;
		}else{
			#quality number equals
			if($seq1Char eq $seq2Char){
				$outSeq .= $seq1Char;
			}else{
				if(rand() > 0.5){
					$outSeq .= $seq1Char;
				}else{
					$outSeq .= $seq2Char;
				}
			}
		}
	}
	
	return($outSeq);
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
           



qap AssembleSeq [options]

Use --help to see more information.

qap is still in development. If you have encounted any problem in usage, please feel no hesitation to cotact us.

=head1 DESCRIPTION

This script implements a function to assemble and extract amplicon sequences from sam/bam files. The script has B<several> mandatory options that MUST appear last. 

=head1 OPTIONS

=over 5

=item --inputFile,-i F<FILE> [Required]

Path to the sam/bam file to be read in.

=item --outputDir,-o F<FILE> [Optional]

Path of the directory to storage the assembled sequences. If NOT provided, the program will generate a folder automatically.

=item --ampliconNumber,-n F<INTEGER> [Optional]

The number of amplicons in the sam/bam file. Either one of --ampliconNumber or -- bedFile MUST be provided.

=item --bedFile,-b F<FILE> [Optional]

The intervals of amplicons. Providing a bed file of amplified intervals is highly recommended. 4 columns are required,[1]chr; [2]start position; [3]end position; [4] interval name.

More information about bed formats, please refer to [http://www.ensembl.org/info/website/upload/bed.html]. Either one of --bedFile or --ampliconNumber MUST be provided.

=item --withIns,-s F<INTEGER> [Optional]

Whether the output fasta file should include sequences with insertions. Choose between '1'(discard sequences with insertions) or '2'(keep sequences with insertions but cut the insertion sites). Default value if '2'.

=item --threads,-t F<INTEGER> [Optional]

Number of threads this program will use when computing. A positive integer should be provided. The default value is 1.

=item --help,-h

Display this detailed help information.

=back

=head1 EXAMPLE

=over 5

qap AssembleSeq -i test.sam -n 10 -t 10 -o ./seq

=back

=head1 AUTHOR

Mingjie Dr.Wang I<huzai@sjtu.edu.cn>

=head1 COPYRIGHT

Copyright (C) 2017, Mingjie Wang. All rights reserved.



