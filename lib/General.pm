#!perl

package General;

use diagnostics;
use strict;
use warnings;
use Cwd qw/getcwd abs_path/;
use File::Basename;
use List::Util qw/min max/;
use FindBin qw/$RealBin/;
use File::Spec;
use POSIX qw/strftime/;

$| = 1; 

use vars qw ($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);
our $VERSION         = 1.00;
our @ISA             = qw (Exporter);
our @EXPORT          = qw ($DEBUG_MODE DelDup  GetSeqFromFasta ParseSamFlag RevCom ProcessBar InfoProcessBar windows2linux
CheckProgram CheckFile CheckFolder CheckPositiveInt CheckDouble Info InfoError InfoWarn InfoDie InfoPlain 
printcol runcmd SplitBigFile CutArray runMultipleThreads runMultipleThreadsWith10Args runMultipleThreadsWith9Args runMultipleThreadsWith8Args runMultipleThreadsWith7Args 
runMultipleThreadsWith6Args runMultipleThreadsWith5Args runMultipleThreadsWith4Args runMultipleThreadsWith3Args 
runMultipleThreadsWith2Args Fastq2Fasta CutLine mix existFile addTagForRawData addTagForFasta isGzipped 
check3EqualNums changeFastqSuffix2Fasta formatMinutes formatFqIntoColumns removeFastaSuffix removeFastqSuffix 
removeSamBamSuffix getFileSize checkMultipleEqualNums getSystemInfo extractSeqFromFasta extractSeqFromFastq 
formatFastaToTwoLineMode checkSeqLen formatFastaForRInput makedir isInArray removeBamSuffix removeSamSuffix 
removeSamBamSuffix2 convert2Rinput getCommonString CheckSeqNum getEleWithMaxFreq removeSuffix checkAbundance removeAllSuffix) ;
our @EXPORT_OK       = qw ($DEBUG_MODE DelDup  GetSeqFromFasta ParseSamFlag RevCom ProcessBar InfoProcessBar windows2linux
CheckProgram CheckFile CheckFolder CheckPositiveInt CheckDouble Info InfoError InfoWarn InfoDie InfoPlain 
printcol runcmd SplitBigFile CutArray runMultipleThreads runMultipleThreadsWith10Args runMultipleThreadsWith9Args runMultipleThreadsWith8Args runMultipleThreadsWith7Args 
runMultipleThreadsWith6Args runMultipleThreadsWith5Args runMultipleThreadsWith4Args runMultipleThreadsWith3Args 
runMultipleThreadsWith2Args Fastq2Fasta CutLine mix existFile addTagForRawData addTagForFasta isGzipped 
check3EqualNums changeFastqSuffix2Fasta formatMinutes formatFqIntoColumns removeFastaSuffix removeFastqSuffix 
removeSamBamSuffix getFileSize checkMultipleEqualNums getSystemInfo extractSeqFromFasta extractSeqFromFastq 
formatFastaToTwoLineMode checkSeqLen formatFastaForRInput makedir isInArray removeBamSuffix removeSamSuffix 
removeSamBamSuffix2 convert2Rinput getCommonString CheckSeqNum getEleWithMaxFreq removeSuffix checkAbundance removeAllSuffix);
our %EXPORT_TAGS     = ( DEFAULT => [qw (&DelDup  &GetSeqFromFasta &ParseSamFlag &RevCom &ProcessBar &InfoProcessBar &windows2linux
&CheckProgram &CheckFile &CheckFolder &CheckPositiveInt &CheckDouble &Info &InfoError &InfoWarn &InfoDie &InfoPlain 
&printcol &runcmd &SplitBigFile &CutArray &runMultipleThreads &runMultipleThreadsWith10Args &runMultipleThreadsWith9Args &runMultipleThreadsWith8Args &runMultipleThreadsWith7Args 
&runMultipleThreadsWith6Args &runMultipleThreadsWith5Args &runMultipleThreadsWith4Args &runMultipleThreadsWith3Args 
&runMultipleThreadsWith2Args &Fastq2Fasta &CutLine &mix &existFile &addTagForRawData &addTagForFasta &isGzipped 
&check3EqualNums &changeFastqSuffix2Fasta &formatMinutes &formatFqIntoColumns &removeFastaSuffix &removeFastqSuffix 
&removeSamBamSuffix &getFileSize &checkMultipleEqualNums &getSystemInfo &extractSeqFromFasta &extractSeqFromFastq 
&formatFastaToTwoLineMode &checkSeqLen &formatFastaForRInput &makedir &isInArray &removeBamSuffix &removeSamSuffix 
&removeSamBamSuffix2 &convert2Rinput &getCommonString &CheckSeqNum &getEleWithMaxFreq &removeSuffix &checkAbundance &removeAllSuffix) ] );

##Debug mode prints more information while running 
our $DEBUG_MODE = 1;

##subprograms begin here##
sub DelDup {
	##delete the duplicates. Input a ref of an array, output the ref of dup-deleted array.
	my $arr = shift;
	my @arr = @$arr;
	my %hash;
	map {$hash{$_}++;} @arr;
	my $res = [keys %hash];
	$res;
}

sub GetSeqFromFasta {
	my $start = shift;
	my $len = shift;
	#my $fa = shift;
	my $filename = shift;

	open (T,"$filename") || die "Can\'t open $filename:$!\n";
	my $str_len;
	my $str_old;
	my $str_want;
	my $str_len_want;
	my $i;
	while (<T>) {
		chomp;
		next if /^>/;

		if ($str_len > $start and $i < 50){
			$str_want .= $str_old;
			$str_len_want = $str_len - length $_ if $i == 0;
			$i++;
		}

		#print "Now:$_\nLst:$str_old\n";
		$str_len += length $_;
		$str_old = $_;

	}

	close (T);
	#print $str_len_want;
	substr $str_want,$start - $str_len_want - 1,$len;

}


sub ParseSamFlag {
	my $input = shift;

	my %hash = (1	=> 'read paired',
	2 => 'read mapped in proper pair',
	4 => 'read unmapped',
	8 => 'mate unmapped',
	16 => 'read reverse strand',
	32 => 'mate reverse strand',
	64 => 'first in pair',
	128 => 'second in pair',
	256 => 'not primary alignment',
	512 => 'read fails platform/vendor quality checks',
	1024 => 'read in PCR or optical duplicate',
	2048 => 'supplementary alignment',
	);

	my @flag = reverse sort {$a <=> $b} keys %hash;
	my $value;

	my $check_later = $input;
	my @want;
	for my $flag (@flag){

		if ($input >= $flag){
			$input -= $flag;
			push @want,$flag;
		}
	}
	my $sum;
	map {$sum += $_;} @want;
	die "wrong flag num!\n" if $sum != $check_later;
	my @res;
	map {push @res,$hash{$_};} sort {$a <=> $b} @want;
	my $res = join ";",@res;
	$res;
}



sub RevCom {
	my $str = shift;
#	$str = uc $str;
	my $rev = reverse $str;
	$rev =~ tr/NATGCagct/NTACGtcga/;
	$rev;
}

sub ProcessBar {
	local $| = 1;
	my $i = $_[0] || return 0;
	my $n = $_[1] || return 0;
	
    #output               
	print   "\r [ ".("=" x int(($i/$n)*50)).(" " x (50 - int(($i/$n)*50)))." ] ";
	printf("%2.1f %%",$i/$n*100);
	local $| = 0;
}

sub InfoProcessBar{
	local $| = 1;
	my $i = $_[0] || return 0;
	my $n = $_[1] || return 0;
	
	chomp (my $date = `date`);
    
    #output               
	print   "\rINFO    : Progress [ ".("=" x int(($i/$n)*50)).(" " x (50 - int(($i/$n)*50)))." ] ";
	printf("%4.2f %%",$i/$n*100);
	local $| = 0;
}

sub CheckProgram {
	my ($program, $srcname, $subline, $debug) = @_;
	if (!-e $program) {
		 if (1) {
			 printf STDERR "DEBUG:  Called from \'%s\' at line \'%u\'...\n", $srcname, $subline;
		}
		printf STDERR "ERROR:  The path \'%s\' does not exist!\n", $program;
		return 0;
	}
	if (!-x $program) {
		if ($debug) {
			printf STDERR "DEBUG:  Called from \'%s\' at line \'%u\'...\n", $srcname, $subline;
		}
		printf STDERR "ERROR:  The path \'%s\' is not an executable program!\n", $program;
		return 0;
	}
	return 1;
}

sub CheckFile {
	my ($filename, $srcname, $subline, $debug) = @_;
	if (!-e $filename) {
		if (1) {
			printf STDERR "DEBUG:  Called from \'%s\' at line \'%u\'...\n", $srcname, $subline;
		}
		printf STDERR "ERROR:  The path \'%s\' does not exist!\n", $filename;
		return 0;
	}
	if (!-f $filename) {
		if ($debug) {
			printf STDERR "DEBUG:  Called from \'%s\' at line \'%u\'...\n", $srcname, $subline;
		}
		printf STDERR "ERROR:  The file \'%s\' does not exist!\n", $filename;
		return 0;
	}
	return 1;
}

sub CheckFolder {
	my ($filename, $srcname, $subline, $debug) = @_;
	if (!-e $filename) {
		if (1) {
			printf STDERR "DEBUG:  Called from \'%s\' at line \'%u\'...\n", $srcname, $subline;
		}
		printf STDERR "ERROR:  The path \'%s\' does not exist!\n", $filename;
		return 0;
	}
	if (!-d $filename) {
		if ($debug) {
			printf STDERR "DEBUG:  Called from \'%s\' at line \'%u\'...\n", $srcname, $subline;
		}
		printf STDERR "ERROR:  The file \'%s\' does not exist!\n", $filename;
		return 0;
	}
	return 1;
}

sub CheckPositiveInt {
	my $value = @_;
	if ($value <= 0){
		return 0;
	}
	if ($value =~ /\./){
		return 0;
	}
	
	return 1;
}

sub CheckDouble {
	my $value = @_;
	if ($value =~ /[^\.0-9]/){
		return 0;
	}
	
	return 1;
}

sub Info {
    my ($info,$color) = @_;
    if (!defined $color){
        $color = "turquoise";
    }

    #chomp (my $date = `date`);
	my $date = strftime("%Y-%m-%d %H:%M:%S %Z", localtime);
    
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
    print "INFO    \@ \[${date}\]: ";
    system "echo \"\e[0;${outcolor};1m${info}\e[m\" ";
    #print "\\u001B[0;${outcolor}m${info}";
}

sub InfoError {
    my ($info,$color) = @_;
    if (!defined $color){
        $color = "red";
    }

    #chomp (my $date = `date`);
	my $date = strftime("%Y-%m-%d %H:%M:%S %Z", localtime);
    
    #color ref
    my %ColorShell = ("black"     => 30,
                      "red"       => 31,
                      "green"     => 32,
                      "yellow"    => 33,
                      "blue"      => 34,
                      "purple"    => 35,
                      "turquoise" => 36,
                      "white"     => 37,);
    my $outcolor = 31; #If inputcolor is not in the ColorShell, use red as default
    if (defined $ColorShell{$color}){
        $outcolor = $ColorShell{$color};
    }
    #output               
    print "ERROR   \@ \[${date}\]: ";
    system "echo \"\e[0;${outcolor};1m${info}\e[m\" ";
}

sub InfoWarn {
    my ($info,$color) = @_;
    if (!defined $color){
        $color = "yellow";
    }

    #chomp (my $date = `date`);
	my $date = strftime("%Y-%m-%d %H:%M:%S %Z", localtime);
    
    #color ref
    my %ColorShell = ("black"     => 30,
                      "red"       => 31,
                      "green"     => 32,
                      "yellow"    => 33,
                      "blue"      => 34,
                      "purple"    => 35,
                      "turquoise" => 36,
                      "white"     => 37,);
    my $outcolor = 33; #If inputcolor is not in the ColorShell, use yellow as default
    if (defined $ColorShell{$color}){
        $outcolor = $ColorShell{$color};
    }
    #output               
    print "WARNING \@ \[${date}\]: ";
    system "echo \"\e[0;${outcolor};1m${info}\e[m\" ";
}

sub InfoDie {
    my ($info,$color) = @_;
    if (!defined $color){
        $color = "white";
    }

    #chomp (my $date = `date`);
	my $date = strftime("%Y-%m-%d %H:%M:%S %Z", localtime);
    
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
    print "Error   \@ \[${date}\]: ";
    system "echo \"\e[0;${outcolor};1m${info}\e[m\" ";
    die "Aborting Misson$!\n";
}

sub InfoPlain {
	my $info = shift;
	#chomp (my $date = `date`);
	my $date = strftime("%Y-%m-%d %H:%M:%S %Z", localtime);
	
    print "INFO    \@ \[${date}\]: $info\n";
}

sub printcol {
    my ($info,$color) = @_;
    if (!defined $color){
        $color = "white";
    }

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
    print " ";    
    system "echo \"\e[0;${outcolor};1m${info}\e[m\" ";           
}

sub runcmd {
	my ($cmd,$flag) = @_;
	if (!defined $flag){
		$flag = 1;
	}

	if ($flag){
		if (ref $cmd){
			for my $tmpcmd (@$cmd){
				&Info($tmpcmd,"white");
				system $tmpcmd;
			}
		}else{
			&Info($cmd,"white");
			system ($cmd);
		}
	}else{
		if (ref $cmd){
			for my $tmpcmd (@$cmd){
				&Info($tmpcmd,"white");
			}
		}else{
			&Info($cmd,"white");
		}
	}
}

sub SplitBigFile {
	my $FileName = shift;
	my $LineNum = shift;
	
	$FileName = basename($FileName);
	system "split -l $LineNum --suffix-length=3 --additional-suffix=.${FileName}.tmp --numeric-suffixes=1 $FileName";
	
}

sub CutArray {
	#cut array into small arrays with specified length
	my $InputArr = shift;
	my $StepSize = shift; #$StepSize means the length of spliced small array

	my $arr_len = scalar @$InputArr;
	my $n = 0;
	my $arr_cut; ##ref of arrary
	
	for (my $i = 0;$i < $arr_len + $StepSize;$i += $StepSize){
		my $k = $i + $StepSize - 1;
		if ($k < $arr_len){
			$arr_cut->[$n] = [@$InputArr[$i..$k]];
		}else{
			$arr_cut->[$n] = [@$InputArr[$i..$arr_len - 1]];
			last;
		}
		$n++;
	}
	
	#use YAML;
	#print Dump($arr_cut);
	
	return $arr_cut;
}

sub runMultipleThreads {
	my ($program,$bigArray,$threadsNum) = @_; #$bigArray is the array of input args
	my $SplitArray = &CutArray($bigArray,$threadsNum);
	for my $SmallArray (@$SplitArray){
		my @threads = ();
		for my $args (@$SmallArray){
			my $thr = threads -> create($program,$args);
			push @threads,$thr;
		}

		for my $tmp (@threads){
			$tmp -> join();
		}
	}
}

sub runMultipleThreadsWith10Args{
	my ($program,$bigArray1,$bigArray2,$bigArray3,$bigArray4,$bigArray5,$bigArray6,$bigArray7,$bigArray8,$bigArray9,$bigArray10,$threadsNum) = @_; #$bigArray is the array of input args
	
	my $index = [0 .. scalar(@$bigArray1) - 1];
	my $SplitIndex = &CutArray($index, $threadsNum);
	
	for my $SmallArray (@$SplitIndex){
		my @threads = ();
		for my $j (@$SmallArray){
			my $arg1 = $bigArray1 -> [$j];
			my $arg2 = $bigArray2 -> [$j];
			my $arg3 = $bigArray3 -> [$j];
			my $arg4 = $bigArray4 -> [$j];
			my $arg5 = $bigArray5 -> [$j];
			my $arg6 = $bigArray6 -> [$j];
			my $arg7 = $bigArray7 -> [$j];
			my $arg8 = $bigArray8 -> [$j];
			my $arg9 = $bigArray9-> [$j];
			my $arg10 = $bigArray10 -> [$j];
			
			my $thr = threads -> create($program,$arg1,$arg2,$arg3,$arg4,$arg5,$arg6,$arg7,$arg8,$arg9,$arg10);
			push @threads,$thr;
		}

		for my $tmp (@threads){
			$tmp -> join();
		}
	}
}

sub runMultipleThreadsWith9Args{
	my ($program,$bigArray1,$bigArray2,$bigArray3,$bigArray4,$bigArray5,$bigArray6,$bigArray7,$bigArray8,$bigArray9,$threadsNum) = @_; #$bigArray is the array of input args
	
	my $index = [0 .. scalar(@$bigArray1) - 1];
	my $SplitIndex = &CutArray($index, $threadsNum);
	
	for my $SmallArray (@$SplitIndex){
		my @threads = ();
		for my $j (@$SmallArray){
			my $arg1 = $bigArray1 -> [$j];
			my $arg2 = $bigArray2 -> [$j];
			my $arg3 = $bigArray3 -> [$j];
			my $arg4 = $bigArray4 -> [$j];
			my $arg5 = $bigArray5 -> [$j];
			my $arg6 = $bigArray6 -> [$j];
			my $arg7 = $bigArray7 -> [$j];
			my $arg8 = $bigArray8 -> [$j];
			my $arg9 = $bigArray9 -> [$j];
			
			my $thr = threads -> create($program,$arg1,$arg2,$arg3,$arg4,$arg5,$arg6,$arg7,$arg8,$arg9);
			push @threads,$thr;
		}

		for my $tmp (@threads){
			$tmp -> join();
		}
	}
}

sub runMultipleThreadsWith8Args {
	my ($program,$bigArray1,$bigArray2,$bigArray3,$bigArray4,$bigArray5,$bigArray6,$bigArray7,$bigArray8,$threadsNum) = @_; #$bigArray is the array of input args
	
	my $index = [0 .. scalar(@$bigArray1) - 1];
	my $SplitIndex = &CutArray($index, $threadsNum);
	
	for my $SmallArray (@$SplitIndex){
		my @threads = ();
		for my $j (@$SmallArray){
			my $arg1 = $bigArray1 -> [$j];
			my $arg2 = $bigArray2 -> [$j];
			my $arg3 = $bigArray3 -> [$j];
			my $arg4 = $bigArray4 -> [$j];
			my $arg5 = $bigArray5 -> [$j];
			my $arg6 = $bigArray6 -> [$j];
			my $arg7 = $bigArray7 -> [$j];
			my $arg8 = $bigArray8 -> [$j];
			
			my $thr = threads -> create($program,$arg1,$arg2,$arg3,$arg4,$arg5,$arg6,$arg7,$arg8);
			push @threads,$thr;
		}

		for my $tmp (@threads){
			$tmp -> join();
		}
	}
}

sub runMultipleThreadsWith7Args {
	my ($program,$bigArray1,$bigArray2,$bigArray3,$bigArray4,$bigArray5,$bigArray6,$bigArray7,$threadsNum) = @_; #$bigArray is the array of input args
	
	my $index = [0 .. scalar(@$bigArray1) - 1];
	my $SplitIndex = &CutArray($index, $threadsNum);
	
	for my $SmallArray (@$SplitIndex){
		my @threads = ();
		for my $j (@$SmallArray){
			my $arg1 = $bigArray1 -> [$j];
			my $arg2 = $bigArray2 -> [$j];
			my $arg3 = $bigArray3 -> [$j];
			my $arg4 = $bigArray4 -> [$j];
			my $arg5 = $bigArray5 -> [$j];
			my $arg6 = $bigArray6 -> [$j];
			my $arg7 = $bigArray7 -> [$j];
			
			my $thr = threads -> create($program,$arg1,$arg2,$arg3,$arg4,$arg5,$arg6,$arg7);
			push @threads,$thr;
		}

		for my $tmp (@threads){
			$tmp -> join();
		}
	}
}

sub runMultipleThreadsWith6Args {
	my ($program,$bigArray1,$bigArray2,$bigArray3,$bigArray4,$bigArray5,$bigArray6,$threadsNum) = @_; #$bigArray is the array of input args
	
	my $index = [0 .. scalar(@$bigArray1) - 1];
	my $SplitIndex = &CutArray($index, $threadsNum);
	
	for my $SmallArray (@$SplitIndex){
		my @threads = ();
		for my $j (@$SmallArray){
			my $arg1 = $bigArray1 -> [$j];
			my $arg2 = $bigArray2 -> [$j];
			my $arg3 = $bigArray3 -> [$j];
			my $arg4 = $bigArray4 -> [$j];
			my $arg5 = $bigArray5 -> [$j];
			my $arg6 = $bigArray6 -> [$j];
			
			my $thr = threads -> create($program,$arg1,$arg2,$arg3,$arg4,$arg5,$arg6);
			push @threads,$thr;
		}

		for my $tmp (@threads){
			$tmp -> join();
		}
	}
}


sub runMultipleThreadsWith5Args {
	my ($program,$bigArray1,$bigArray2,$bigArray3,$bigArray4,$bigArray5,$threadsNum) = @_; #$bigArray is the array of input args
	
	my $index = [0 .. scalar(@$bigArray1) - 1];
	my $SplitIndex = &CutArray($index, $threadsNum);
	
	for my $SmallArray (@$SplitIndex){
		my @threads = ();
		for my $j (@$SmallArray){
			my $arg1 = $bigArray1 -> [$j];
			my $arg2 = $bigArray2 -> [$j];
			my $arg3 = $bigArray3 -> [$j];
			my $arg4 = $bigArray4 -> [$j];
			my $arg5 = $bigArray5 -> [$j];
			
			my $thr = threads -> create($program,$arg1,$arg2,$arg3,$arg4,$arg5);
			push @threads,$thr;
		}

		for my $tmp (@threads){
			$tmp -> join();
		}
	}
}

sub runMultipleThreadsWith4Args {
	my ($program,$bigArray1,$bigArray2,$bigArray3,$bigArray4,$threadsNum) = @_; #$bigArray is the array of input args
	
	my $index = [0 .. scalar(@$bigArray1) - 1];
	my $SplitIndex = &CutArray($index, $threadsNum);
	
	for my $SmallArray (@$SplitIndex){
		my @threads = ();
		for my $j (@$SmallArray){
			my $arg1 = $bigArray1 -> [$j];
			my $arg2 = $bigArray2 -> [$j];
			my $arg3 = $bigArray3 -> [$j];
			my $arg4 = $bigArray4 -> [$j];
			
			my $thr = threads -> create($program,$arg1,$arg2,$arg3,$arg4);
			push @threads,$thr;
		}

		for my $tmp (@threads){
			$tmp -> join();
		}
	}
}

sub runMultipleThreadsWith3Args {
	my ($program,$bigArray1,$bigArray2,$bigArray3,$threadsNum) = @_; #$bigArray is the array of input args
	
	my $index = [0 .. scalar(@$bigArray1) - 1];
	my $SplitIndex = &CutArray($index, $threadsNum);
	
	for my $SmallArray (@$SplitIndex){
		my @threads = ();
		for my $j (@$SmallArray){
			my $arg1 = $bigArray1 -> [$j];
			my $arg2 = $bigArray2 -> [$j];
			my $arg3 = $bigArray3 -> [$j];
			
			my $thr = threads -> create($program,$arg1,$arg2,$arg3);
			push @threads,$thr;
		}

		for my $tmp (@threads){
			$tmp -> join();
		}
	}
}

sub runMultipleThreadsWith2Args {
	my ($program,$bigArray1,$bigArray2,$threadsNum) = @_; #$bigArray is the array of input args
	
	my $index = [0 .. scalar(@$bigArray1) - 1];
	my $SplitIndex = &CutArray($index, $threadsNum);
	
	for my $SmallArray (@$SplitIndex){
		my @threads = ();
		for my $j (@$SmallArray){
			my $arg1 = $bigArray1 -> [$j];
			my $arg2 = $bigArray2 -> [$j];
			
			my $thr = threads -> create($program,$arg1,$arg2);
			push @threads,$thr;
		}

		for my $tmp (@threads){
			$tmp -> join();
		}
	}
}

sub Fastq2Fasta {
	my ($fq,$fa) = @_;
	my $cmd = "awk \'BEGIN\{P=1\}\{if(P==1\|\|P==2)\{gsub(\/^\[@\]/,\">\");print\}; if(P==4)P=0; P++\}\' $fq \> $fa";
	system($cmd);
}

sub CutLine {
	my ($line,$n) = @_;

	my $len = length $line;
	my @pos_rand;
	my $min_seq_range = $len / ($n * 3);

	for my $i (1..$n-1){
		my $tmp = int(rand($len));
		if ((grep {abs($tmp - $_) < $min_seq_range} @pos_rand) || $tmp == 0 || $tmp == $len){
			redo;
		}
		push @pos_rand,$tmp;
	}
	my @pos = (0,@pos_rand,$len);
	@pos = sort {$a <=> $b} @pos;
	#print "@pos\n";

	my @cutlines;
	for(my $k = 0; $k < scalar(@pos) - 1; $k++){
		my $start = $pos[$k];
		my $end = $pos[$k + 1];
		my $len = $end - $start;
		my $lineCut = substr $line,$start,$len;
		push @cutlines,$lineCut;
	}

	return (\@cutlines,\@pos);
}

sub mix {
	my($smallArr_ref,$bigArr_ref) = @_;
	
	my @smallArr = @$smallArr_ref;
	my @bigArr = @$bigArr_ref;
	my @mix;
	my $mixLen = scalar(@smallArr) + scalar(@bigArr);
	for my $k (0..$mixLen - 1){
		if($k % 2 == 0){
			my $index = $k / 2;
			$mix[$k] = $bigArr[$index];
		}else{
			my $index = ($k - 1) / 2;
			$mix[$k] = $smallArr[$index];
		}
	}
	return \@mix;
}

sub existFile {
	my $file = shift;
	
	my $exist_flag = 0;
	if (-e $file and -s $file){
		$exist_flag = 1;
	}
	
	return $exist_flag;
}

sub gtf2genelist {
	my ($gtf,$genelist) = @_;
	
	open OUT,">$genelist" or InfoDie("Can NOT output to file $genelist:$!",'red');
	print OUT "id\tname\n";
	open GTF,"$gtf" or InfoDie("Can NOT open file $gtf:$!",'red');
	while (<GTF>){
		chomp;
		#gene_id "ENSG00000232023"; gene_name "AC009410.1";
		if($_ =~ /gene_id \"(ENSG\d+)\"; gene_name \"(.*?)\";/){
			my $gene_id = $1;
			my $gene_name = $2;
			print OUT "$gene_id\t$gene_name\n";
		}
	}
}

sub addTagForRawData {
	my $rawDataName = shift;
	my $tag = shift;
	my $removegz = shift;
	
	my $newName = '';
	
	if($rawDataName =~ /\.fastq\.gz$/){
		$newName = $rawDataName =~ s/\.fastq\.gz$/.${tag}.fastq.gz/r;
	}elsif($rawDataName =~ /\.fq\.gz$/){
		$newName = $rawDataName =~ s/\.fq\.gz$/.${tag}.fq.gz/r;
	}elsif($rawDataName =~ /\.fastq$/){
		$newName = $rawDataName =~ s/\.fastq$/.${tag}.fastq/r;
	}elsif($rawDataName =~ /\.fq$/){
		$newName = $rawDataName =~ s/\.fq$/.${tag}.fq/r;
	}else{
		InfoError("The name of $rawDataName is NOT compliant with the naming conventions, which should be ended with \".fastq/.fq/.fastq.gz/.fa.gz\".");
		exit;
	}
	
	if ($removegz and $newName =~ /\.gz$/){
		$newName =~ s/\.gz$//;
	}
	return $newName;
}

sub addTagForFasta {
	my $oldName = shift;
	my $tag = shift;
	my $rmgz = shift;
	
	my $newName;
	if ($oldName =~ /\.fasta$/){
		$newName = $oldName =~ s/\.fasta$/.${tag}.fasta/r;
	}elsif($oldName =~ /\.fasta\.gz$/){
		$newName = $oldName =~ s/\.fasta\.gz$/.${tag}.fasta.gz/r;
	}elsif($oldName =~ /\.fa$/){
		$newName = $oldName =~ s/\.fa$/.${tag}.fa/r;
	}elsif($oldName =~ /\.fa\.gz/){
		$newName = $oldName =~ s/\.fa\.gz$/.${tag}.fa.gz/r;
	}elsif($oldName =~ /\.fas$/){
		$newName = $oldName =~ s/\.fas$/.${tag}.fas/r;
	}elsif($oldName =~ /\.fas\.gz$/){
		$newName = $oldName =~ s/\.fas\.gz$/.${tag}.fas.gz/r;
	}else{
		InfoError("The suffix of $oldName should be one of \".fasta/.fasta.gz/.fas/.fas.gz/.fa/.fa.gz\".");
		exit;
	}
	
	if ($rmgz){
		$newName =~ s/\.gz//;
	}
	
	return $newName;
}

sub isGzipped {
	my $file = shift;
	
	my $isGzipped = 0;
	if ($file =~ /\.gz$/){
		$isGzipped = 1;
	}
	
	return $isGzipped;
}

sub check3EqualNums {
	my $n1 = shift;
	my $n2 = shift;
	my $n3 = shift;
	
	if ($n1 == $n2 && $n1 == $n3 && $n2 == $n3){
		return 1;
	}else{
		return 0;
	}

}

sub changeFastqSuffix2Fasta{
	my $name = shift;
	
	if ($name =~ /\.fastq$/){
		$name =~ s/\.fastq$/.fasta/;
	}elsif ($name =~ /\.fq$/){
		$name =~ s/\.fq$/.fasta/;
	}else{
		InfoError("The suffix of $name should be .fastq or .fq");
	}
	
	return $name;
}

sub formatMinutes {
	my $n = shift;
	
	my $h = int($n / 60);
	my $m = $n % 60;
	
	my $out;
	if ($h > 0){
		$out = $h . " Hour(s) " . $m . " Minute(s)";
	}else{
		$out = $m . " Minute(s)";
	}
	return $out;
}

sub formatFqIntoColumns {
	my $file = shift;
	my $out = shift;
	my $format = shift;
	
	my $fileName = basename($file);
	Info("Formatting fastq file $fileName for R input.");
	
	open OUT,">$out" or die "Can NOT output to file $out:$!\n";
	if ($format eq 'fastq'){
		print OUT "ID\tSeq\tAnno\tQual\n";
	}else{
		print OUT "ID\tSeq\n";
	}
	
	if ($format eq 'fastq'){
		open T,$file or die "Can NOT open file $file:$!\n";
		while (my $line1 = <T>){
			chomp $line1;
			chomp (my $line2 = <T>);
			chomp (my $line3 = <T>);
			chomp (my $line4 = <T>);
		
			print OUT "$line1\t$line2\t$line3\t$line4\n";
		}
		close T;
	}else{
		open T,$file or die "Can NOT open file $file:$!\n";
		while (my $line1 = <T>){
			chomp $line1;
			chomp (my $line2 = <T>);
			chomp (my $line3 = <T>);
			chomp (my $line4 = <T>);
		
			print OUT "$line1\t$line2\n";
		}
		close T;
	}
}

sub removeFastaSuffix {
	my $oldName = shift;
	
	my $newName;
	if ($oldName =~ /\.fasta$/){
		$newName = $oldName =~ s/\.fasta$//r;
	}elsif($oldName =~ /\.fasta\.gz$/){
		$newName = $oldName =~ s/\.fasta\.gz$//;
	}elsif($oldName =~ /\.fa$/){
		$newName = $oldName =~ s/\.fa$//r;
	}elsif($oldName =~ /\.fa\.gz/){
		$newName = $oldName =~ s/\.fa\.gz$//r;
	}elsif($oldName =~ /\.fas$/){
		$newName = $oldName =~ s/\.fas$//r;
	}elsif($oldName =~ /\.fas\.gz$/){
		$newName = $oldName =~ s/\.fas\.gz$//r;
	}elsif($oldName =~ /\.fsa$/){
		$newName = $oldName =~ s/\.fsa$//r;
	}elsif($oldName =~ /\.fsa\.gz$/){
		$newName = $oldName =~ s/\.fsa\.gz$//r;
	}else{
		InfoError("The suffix of $oldName should be one of \".fasta/.fasta.gz/.fas/.fas.gz/.fa/.fa.gz/.fsa/.fsa.gz\".");
		#return $oldName;
		exit;
	}
	
	return $newName;
}

sub removeFastqSuffix {
	my $oldName = shift;
	
	my $newName;
	if ($oldName =~ /\.fastq$/){
		$newName = $oldName =~ s/\.fastq$//r;
	}elsif($oldName =~ /\.fastq\.gz$/){
		$newName = $oldName =~ s/\.fastq\.gz$//;
	}elsif($oldName =~ /\.fq$/){
		$newName = $oldName =~ s/\.fq$//r;
	}elsif($oldName =~ /\.fq\.gz/){
		$newName = $oldName =~ s/\.fq\.gz$//r;
	}else{
		InfoError("The suffix of $oldName should be one of \".fastq/.fastq.gz/.fq/.fq.gz/\".");
		exit;
	}
	
	return $newName;
}

sub removeSamBamSuffix {
	my $oldName = shift;
	
	my $newName;
	if ($oldName =~ /\.bam$/){
		$newName = $oldName =~ s/\.bam$//r;
	}elsif($oldName =~ /\.sam$/){
		$newName = $oldName =~ s/\.sam$//;
	}elsif($oldName =~ /\.NameSorted\.bam$/){
		$newName = $oldName =~ s/\.NameSorted\.bam$//r;
	}elsif($oldName =~ /\.NameSorted\.sam$/){
		$newName = $oldName =~ s/\.NameSorted\.sam$//r;
	}elsif($oldName =~ /\.PosSorted\.bam$/){
		$newName = $oldName =~ s/\.PosSorted\.bam$//r;
	}elsif($oldName =~ /\.PosSorted\.sam$/){
		$newName = $oldName =~ s/\.PosSorted\.sam$//r;
	}else{
		InfoError("The suffix of $oldName should be one of \'.bam\' or \'.sam\'.");
		exit;
	}
	
	return $newName;
}

sub getFileSize {
	my $file = shift;
	
	my $size;
	
	if(-e $file){
		my @args = stat ($file);
		$size = $args[7];
	}else{
		InfoError("$file does NOT exist! Please check again.");
		exit;
	}
	
	return $size;
}

sub checkMultipleEqualNums {
	my $ref = shift;
	
	my @num = @$ref;
	
	if(min(@num) == max(@num)){
		return 1;
	}else{
		return 0;
	}
}

sub getSystemInfo {
	my $system = 'linux';
	my $cmd = `cat /etc/*-release`;
	chomp $cmd;
	
	if ($cmd =~ /centos/i){
		$system = 'centos';
	}elsif($cmd =~ /ubuntu/i){
		$system = 'ubuntu';
	}elsif($cmd =~ /red hat/i or $cmd =~ /redhat/i){
		$system = 'redhat';
	}elsif($cmd =~ /fedora/i){
		$system = 'fedora';
		if ($cmd =~ /release (\d+)/i){
			my $release = $1;
			$system .= $release;
		}
	}elsif($cmd =~ /openSUSE/i){
		$system = 'opensuse';
	}elsif($cmd =~ /debian/i){
		$system = 'debian';
	}else{
		$system = 'other linux';
	}
	
	return $system;
}

sub extractSeqFromFasta {
	my $file = shift;
	my $outfile = shift;
	
	my $outfile_tmp = $outfile . '.' . time() . '.tmp';
	open TMP,">$outfile_tmp" or die "Cannot output to file $outfile_tmp:$!";
	
	open T,$file or die "Cannot open file $file:$!";
	while (<T>){
		chomp;
		if (/^>/){
			print TMP "\n$_\n";
		}else{
			print TMP "$_";
		}
	}
	close T;
	close TMP;
	
	open RES,">$outfile" or die "Cannot output to file $outfile:$!";
	
	open R,$outfile_tmp or die "Cannot open file $outfile_tmp:$!";
	while(<R>){
		chomp;
		if($_ eq ''){
			next;
		}elsif(/^>/){
			next;
		}else{
			print RES "$_\n";
		}
	}
	close R;
	close RES;
	
	my $cmd = "rm -rf $outfile_tmp";
	system($cmd);
}

sub extractSeqFromFastq {
	my $file = shift;
	my $outfile = shift;
	
	open RES,">$outfile" or die "Can not output to file $outfile:$!";
	
	open T,$file or die "Can not open file $file:$!";
	while(my $line1 = <T>){
		chomp $line1;
		chomp (my $line2 = <T>);
		chomp (my $line3 = <T>);
		chomp (my $line4 = <T>);
		
		print RES "$line2\n";
	}
	close T;
	close RES;
	
}

sub formatFastaToTwoLineMode {
	my $file = shift;
	my $outfile = shift;
	
	my $outfile_tmp = $outfile . '.' . time() . '.tmp';
	open TMP,">$outfile_tmp" or die "Cannot output to file $outfile_tmp:$!";
	
	open T,$file or die "Cannot open file $file:$!";
	while (<T>){
		chomp;
		if (/^>/){
			print TMP "\n$_\n";
		}else{
			$_ = uc($_);
			$_ =~ s/[^A-Z\-\*]//g;
			print TMP "$_";
		}
	}
	close T;
	close TMP;
	
	open RES,">$outfile" or die "Cannot output to file $outfile:$!";
	
	open R,$outfile_tmp or die "Cannot open file $outfile_tmp:$!";
	while(<R>){
		chomp;
		if($_ eq ''){
			next;
		}elsif(/^>/){
			print RES "$_\n";
		}else{
			print RES "$_\n";
		}
	}
	#print RES "\n";
	close R;
	close RES;
	
	my $cmd = "rm -rf $outfile_tmp";
	system($cmd);
}

sub checkSeqLen {
	my $file = shift;
	my $outputDir = shift;
	
	#convert fasta file to two line format
	my $file2Name = removeFastaSuffix(basename($file)) . ".2line.fasta";
	my $file2 = File::Spec -> catfile($outputDir,$file2Name);
	formatFastaToTwoLineMode($file,$file2);
	
	open T,$file2 or die "Can not open fasta file $file2:$!";
	my @len;
	while(my $line1 = <T>){
		chomp $line1;
		chomp (my $line2 = <T>);
		
		next if($line1 eq '' or $line2 eq '');
		
		my $len = length $line2;
		push @len,$len;
	}
	close T;
	
	system("rm -rf $file2");	
	
	if(min(@len) == max(@len)){
		return $len[0];
	}else{
		#if sequences length does not equal, no need to return, directly exit
		InfoError("Sequences in $file does not have the same length. Please check the file.");
		InfoError("Aborting...");
		exit(0);
	}
	
}

sub formatFastaForRInput {
	#remove the sequence ID line, and output the sequence line for R input
	my $inputfile = shift;
	my $outputfile = shift;
	
	open RES,">$outputfile" or die "Can NOT output to file $outputfile:$!";
	
	open T,$inputfile or die "Can NOT open file $inputfile:$!";
	while(my $line1 = <T>){
		chomp $line1;
		chomp (my $line2 = <T>);
		
		print RES "$line2\n";
		
	}
	close T;
	close RES;
}

sub makedir {
	my $dir = shift;
	
	if (-e $dir){
		return 1;
	}else{
		mkdir $dir or die "Can NOT mkdir $dir:$!";
	}
	
}

sub isInArray {
	my $id = shift;
	my $arr = shift;
	my @arr = @$arr;
	
	for my $f (@arr){
		if ($id eq $f){
			return 1;
		}
	}
	
	return 0;
}

sub removeBamSuffix {
	my $name = shift;
	my $removeSortFlag = shift;
	
	$name =~ s/\.bam$//i;
	
	if ($removeSortFlag){
		$name =~ s/\.PosSorted//i;
		$name =~ s/\.NameSorted//i;
		$name =~ s/\.rmdup//i;
		$name =~ s/\.withRG//i;
	}
	
	return $name;
}

sub removeSamSuffix {
	my $name = shift;
	my $removeSortFlag = shift;
	
	$name =~ s/\.sam$//i;
	
	if ($removeSortFlag){
		$name =~ s/\.PosSorted//i;
		$name =~ s/\.NameSorted//i;
		$name =~ s/\.rmdup//i;
		$name =~ s/\.withRG//i;
	}
	
	return $name;
}

sub removeSamBamSuffix2 {
	my $name = shift;
	my $removeFlag = shift;
	
	if ($name =~ /\.bam$/){
		removeBamSuffix($name,1);
	}elsif($name =~ /\.sam$/){
		removeSamSuffix($name,1);
	}else{
		InfoError("The suffix of $name should be .bam or .sam");
		exit;
	}
}

sub convert2Rinput {
	my $inputfile = shift;
	my $outputDir = shift;
	
	#format fasta to 2 line mode
	my $inputfile2 = File::Spec -> catfile($outputDir,removeFastaSuffix(basename($inputfile)) . ".2line.fasta");
	formatFastaToTwoLineMode($inputfile,$inputfile2);
	
	my $rinput = File::Spec -> catfile($outputDir,removeFastaSuffix(basename($inputfile)) . ".RInput");
	open R,">$rinput" or die "Can NOT output to file $rinput:$!";
	
	open T,$inputfile2 or die "Can NOT open fasta file $inputfile2:$!";
	while(my $line1 = <T>){
		chomp $line1;
		chomp (my $line2 = <T>);
		
		print R "$line2\n";
	}
	close T;
	close R;
	
	return $rinput;

}

sub getCommonString {
	my $str1 = shift;
	my $str2 = shift;
	
	my @str1 = split '',$str1;
	my @str2 = split '',$str2;
	
	if (not scalar @str1 == scalar @str2){
		InfoWarn("\'$str1\' and \'$str2\' does NOT have the same length.");
	}
	
	my @out;
	my $out;
	for my $i (0..scalar(@str1) - 1){
		my $char1 = $str1[$i];
		my $char2 = $str2[$i];
		
		if ($char1 eq $char2){
			push @out,$char1;
		}else{
			$out = join '',@out;
			return $out;
		}
	}
	
	$out = join '',@out;
	return $out;
}

sub CheckSeqNum {
	my $file = shift;
	
	my $seqnum = 0;
	open T,$file or die "Can NOT open $file:$!";
	while(<T>){
		if (/^>/){
			$seqnum++;
		}
	}
	close T;
		
	return $seqnum;
}

sub windows2linux {
	my $file = shift;
	
	if(not existFile($file)){
		InfoError("$file does NOT exist.");
		return 0;
	}
	
	my $cmd = "sed -i \'s\/\\r\/\/g\' $file ";
	runcmd($cmd);
	
	return 1;
}

sub getEleWithMaxFreq {
	my $arr = shift;
	
	my @arr = @$arr;
	my %hash;
	for my $e (@arr){
		$hash{$e}++;
	}
	
	my @sort;
	for my $k (sort {$hash{$b} <=> $hash{$a}} keys %hash){
		push @sort,$k;
	}
	
	return $sort[0];
	
}

sub removeSuffix {
	my $name = shift;
	
	$name =~ s/(.*)\..*?//;
	
	return $1;
}

sub checkAbundance{
	my $file = shift;
	my $outfile = shift;
	
	my $outdir = dirname($outfile);
	
	my $outMsg = "Sequence abundance not detected";
	my $abundFlag = 0;
	
	my $outAbundanceFile = File::Spec -> catfile($outdir, removeFastaSuffix(basename($file)) . ".abund");
	open RES,">$outAbundanceFile" or die "Can not output to file $outAbundanceFile:$!";
	print RES "seqID\tcount\n";
	
	open OUTSEQ,">$outfile" or die "Can not output to file $outfile:$!";
	
	open T,$file or die "Cannot open file $file:$!";
	my $seqCount;
	while(my $line1 = <T>){
		chomp $line1;
		chomp (my $line2 = <T>);
		
		if($line1 =~ /^>/){
			if($line1 =~ /^>(.*?);size=(\d+)/i){
				print RES "$1\t$2\n";
				$abundFlag++;
			}else{
				my $id = $line1 =~ s/^>//r; 
				print RES "$id\t1\n";
			}
			$line1 =~ s/;size=.*//;
			print OUTSEQ "$line1\n$line2\n";
		}else{
			InfoError("Fasta file $file is wrong formatted. Please check.");
			exit(0);
		}
		
		$seqCount++;
	}
	
	if($abundFlag > 0 and $abundFlag < $seqCount){
		$outMsg = "Sequence abundance partly detected and saved";
	}elsif($abundFlag == $seqCount){
		$outMsg = "Sequence abundance fully detected and saved";
	}else{
		$outMsg = "Sequence abundance not detected and created";
	}
	
	return ($outMsg,$outAbundanceFile);
}

sub removeAllSuffix {
    my $name = shift;

    my $newname = basename($name);
    while($newname =~ /\..*$/){
        $newname =~ s/\..*$//;
    }   
    
    my $outname = File::Spec -> catfile(dirname($name), $newname);
    
    return $outname;
}


1;


