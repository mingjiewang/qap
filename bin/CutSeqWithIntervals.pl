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
printcol ("CutSeqWithIntervals","green");
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
my $intervalFile;
my $threads;

my $DateNow = `date +"%Y%m%d_%Hh%Mm%Ss"`;
chomp $DateNow;

GetOptions(
'i|inputDir|=s'     => \$inputDir,
's|suffix|=s'       => \$suffix,
'o|outputDir|=s'    => \$outputDir,
't|threads|=s'      => \$threads,
'v|interval|=s'     => \$intervalFile,
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
	$outputDir = File::Spec -> catfile($wk_dir,"qap_Results_for_CutSeqWithIntervals_$DateNow");
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

sleep(1);

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

my @start;
my @end;
my @name;
if(defined $intervalFile){
	if(existFile($intervalFile)){
		open T,$intervalFile or die "Can NOT open interval file $intervalFile:$!\n";
		
		my $i = 1;
		while(<T>){
			chomp;
			
			if ($_ eq ''){
				next;
			}
			
			# check interval file format
			if(/[\d+\.\.\d+,?]+\s*\w*/){
				# pass
				my @arr = split(/\s+/,$_);
				if(scalar(@arr) == 1){
					my $pos = $arr[0];
					if($pos =~ /,/){
						my @startEndPair = split ',',$pos;
						
						my @allstart;
						my @allend;
						for my $se (@startEndPair){
							my @tmp = split(/\.\./,$se);
							
							my $start = $tmp[0];
							push @allstart,$start;
							
							my $end = $tmp[1];
							push @allend,$end; 
						}
						
						if(scalar(@allend) == scalar(@allstart)){
							Info("Interval file LINE $i: Joined subintervals detected.");
						}else{
							InfoError("The format of interval file \(LINE\:$i\) MUST be in the right format. Please check.");
				
							pod2usage(-verbose=>2,-exitval=>1);
							exit;
						}
						
						my $start = join ',',@allstart;
						my $end = join ',',@allend;
						
						push @start,$start;
						push @end,$end;
						
					}else{
						my @tmp = split(/\.\./,$pos);
						my $start = $tmp[0];
						my $end = $tmp[1];
						
						push @start,$start;
						push @end,$end;
					}

					push @name,'';
				}elsif(scalar(@arr) == 2){
					my $pos = $arr[0];
					if($pos =~ /,/){
						my @startEndPair = split ',',$pos;
						
						my @allstart;
						my @allend;
						for my $se (@startEndPair){
							my @tmp = split(/\.\./,$se);
							
							my $start = $tmp[0];
							push @allstart,$start;
							
							my $end = $tmp[1];
							push @allend,$end; 
						}
						
						if(scalar(@allend) == scalar(@allstart)){
							Info("Interval file LINE $i: Joined subintervals detected.");
						}else{
							InfoError("The format of interval file \(LINE\:$i\) MUST be in the right format. Please check.");
				
							pod2usage(-verbose=>2,-exitval=>1);
							exit;
						}
						
						my $start = join ',',@allstart;
						my $end = join ',',@allend;
						
						push @start,$start;
						push @end,$end;
						
					}else{
						my @tmp = split(/\.\./,$pos);
						my $start = $tmp[0];
						my $end = $tmp[1];
						
						push @start,$start;
						push @end,$end;
					}

					push @name,$arr[1];
				}else{
					InfoError("The format of interval file \(LINE\:$i\) MUST be in the right format. Please check.");
				
					pod2usage(-verbose=>2,-exitval=>1);
					exit;
				}
			}else{
				InfoError("The format of interval file \(LINE\:$i\) MUST be in the right format. Please check.");
				
				pod2usage(-verbose=>2,-exitval=>1);
				exit;
			}
			
			$i++;
		}
		close T;
		
	}elsif(not existFile($intervalFile)){
		InfoError("Input interval file $intervalFile does NOT exist! Please check again.");
		exit;
	}
}else{
	InfoError("Interval file MUST be provided by using --interva/-v.");
	pod2usage(-verbose=>0,-exitval=>1);
	exit;
}


##core program starts here
##get output folder ready
for my $i(0..scalar(@start) - 1){
	my $foldername = "From" . $start[$i] . "To" . $end[$i];
	my $folder = File::Spec -> catfile($outputDir, $foldername);
	
	$folder =~ s/,/./g;
	
	mkdir $folder or die "Can NOT mkdir $folder:$!\n";
}


##start the cutting program
if($threads > 1){
	Info("Start cutting sequences using multipe threads. Please wait...");
	
	my @bigStart;
	my @bigEnd;
	my @bigName;
	for my $f (@inputfiles){
		push @bigStart,\@start;
		push @bigEnd,\@end;
		push @bigName,\@name;
	}
	
	&runMultipleThreadsWith4Args(\&cutSeqWithMultipleIntervalsWithInfo, \@inputfiles, \@bigStart, \@bigEnd, \@bigName, $threads);
	
}else{
	Info("Start cutting sequences using single thread. Please wait...");
	
	my $i = 1;
	for my $f (@inputfiles){
		
		&cutSeqWithMultipleIntervals($f, \@start, \@end, \@name);
		
		InfoProcessBar($i, $numberOfFiles);
		$i++;
	}
}


##sub program starts here
sub cutSeq{
	my $inputfile = shift;
	my $outfile = shift;
	my $start = shift;
	my $end = shift;
	
	# reformat input fasta file into 2 line format
	my $inputfile2 = File::Spec -> catfile($outputDir, removeFastaSuffix(basename($inputfile)) . ".refmt.fasta" );
	
	formatFastaToTwoLineMode($inputfile, $inputfile2);
	
	# output file handle
	open RES,">$outfile" or die "Can NOT output to file $outfile:$!\n";
	
	# read in file
	open T,$inputfile2 or die "Can NOT open fasta file $inputfile2:$!\n";
	while(my $line1 = <T>){
		chomp $line1;
		chomp (my $line2 = <T>);
		
		if ($line1 =~ /^>/){
			print RES "$line1\n";
		}else{
			InfoError("The input fasta file is NOT right formatted.");
			#exit;
			return; # Stop translating for current file and keep running the program
		}
		
		my $cutseq;
		if ($start =~ /,/){
			my @start = split ',',$start;
			my @end = split ',',$end;
			
			for my $k (0..scalar(@start) - 1){
				my $startPos = $start[$k];
				my $endPos = $end[$k];
				
				my $subcutseq = substr $line2, $startPos - 1, $endPos - $startPos + 1;
				
				$cutseq .= $subcutseq;
			}
		}else{
			$cutseq = substr $line2, $start - 1, $end - $start + 1;
		}
		
		print RES "$cutseq\n";
	}
	close T;
	close RES;
	
	# remove refomat fasta file
	my $cmd = "rm -rf $inputfile2";
	system($cmd);
	
	
}

sub cutSeqWithMultipleIntervals {
	my $inputfile = shift;
	my $start = shift;
	my $end = shift;
	my $name = shift;
	
	my @start = @$start;
	my @end = @$end;
	my @name = @$name;
	# check pos number
	my @numToCheck = (scalar(@start), scalar(@end), scalar(@name));
	if (checkMultipleEqualNums(\@numToCheck)){
		# nothing
	}else{
		InfoError("The number of start positions, end positions and interval names are not equal. Please check.");
		exit;
	}
	
	for my $i (0..scalar(@start) - 1){
		my $startPos = $start[$i];
		my $endPos = $end[$i];
		my $intervalName = $name[$i];
		
		# output subfolder
		my $outfoldername = "From" . $start[$i] . "To" . $end[$i];
		$outfoldername =~ s/\,/./g;
		my $outfolder = File::Spec -> catfile($outputDir, $outfoldername);
		
		my $outname;
		if ($intervalName ne ''){
			if ($startPos =~ /,/){
				my $tmpStart = $startPos =~ s/,/./gr;
				my $tmpEnd = $endPos =~ s/,/./gr;
				$outname = removeFastaSuffix(basename($inputfile)) . "_From" . "$tmpStart" . "To" . "$tmpEnd" . "_" . "$intervalName" . ".fasta";
			}else{
				$outname = removeFastaSuffix(basename($inputfile)) . "_From" . "$startPos" . "To" . "$endPos" . "_" . "$intervalName" . ".fasta";
			}
		}else{
			if ($startPos =~ /,/){
				my $tmpStart = $startPos =~ s/,/./gr;
				my $tmpEnd = $endPos =~ s/,/./gr;
				$outname = removeFastaSuffix(basename($inputfile)) . "_From" . "$tmpStart" . "To" . "$tmpEnd" . "_" . "$intervalName" . ".fasta";
			}else{
				$outname = removeFastaSuffix(basename($inputfile)) . "_From" . "$startPos" . "To" . "$endPos" . ".fasta";
			}
			
		}
		
		my $outfile = File::Spec -> catfile($outfolder, $outname);
		
		&cutSeq($inputfile, $outfile, $startPos, $endPos);

	}
}


sub cutSeqWithMultipleIntervalsWithInfo {
	my $inputfile = shift;
	my $start = shift;
	my $end = shift;
	my $name = shift;
	
	Info("Start to cut $inputfile.");
	
	my @start = @$start;
	my @end = @$end;
	my @name = @$name;
	# check pos number
	my @numToCheck = (scalar(@start), scalar(@end), scalar(@name));
	if (checkMultipleEqualNums(\@numToCheck)){
		# nothing
	}else{
		InfoError("The number of start positions, end positions and interval names are not equal. Please check.");
		exit;
	}
	
	for my $i (0..scalar(@start) - 1){
		my $startPos = $start[$i];
		my $endPos = $end[$i];
		my $intervalName = $name[$i];
		
		# output subfolder
		my $outfoldername = "From" . $start[$i] . "To" . $end[$i];
		$outfoldername =~ s/\,/./g;
		my $outfolder = File::Spec -> catfile($outputDir, $outfoldername);
		
		my $outname;
		if ($intervalName ne ''){
			if ($startPos =~ /,/){
				my $tmpStart = $startPos =~ s/,/./gr;
				my $tmpEnd = $endPos =~ s/,/./gr;
				$outname = removeFastaSuffix(basename($inputfile)) . "_From" . "$tmpStart" . "To" . "$tmpEnd" . "_" . "$intervalName" . ".fasta";
			}else{
				$outname = removeFastaSuffix(basename($inputfile)) . "_From" . "$startPos" . "To" . "$endPos" . "_" . "$intervalName" . ".fasta";
			}
		}else{
			if ($startPos =~ /,/){
				my $tmpStart = $startPos =~ s/,/./gr;
				my $tmpEnd = $endPos =~ s/,/./gr;
				$outname = removeFastaSuffix(basename($inputfile)) . "_From" . "$tmpStart" . "To" . "$tmpEnd" . "_" . "$intervalName" . ".fasta";
			}else{
				$outname = removeFastaSuffix(basename($inputfile)) . "_From" . "$startPos" . "To" . "$endPos" . ".fasta";
			}
			
		}
		
		my $outfile = File::Spec -> catfile($outfolder, $outname);
		
		&cutSeq($inputfile, $outfile, $startPos, $endPos);

	}
}


##run success
print("\n\n");
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
           



qap CutSeqWithIntervals [options]

Use --help to see more information.

qap is still in development. If you have encounted any problem in usage, please feel no hesitation to cotact us.

=head1 DESCRIPTION

This script implements a function for cutting fasta sequences with base intervals in batch. The script has B<several> mandatory options that MUST appear last. 

=head1 OPTIONS

=over 5

=item --inputDir,-i F<FILE> [Required]

Path to directory contaning all the files to be read in.

=item --interval,-v F<FILE> [Required]

Path to the interval file which contains two columns: 1. Start Position..End Position; 2. Interval Name. The first column are mandatory and the second one is optional. Columns should be seperated by tab or space and interval names should not contain any special characters. 

If there are several sub-intervals within the same interval, the first column should be provided with format: Start Position1..End Postion1,Start Postion2..End Position2. e.g. 154..369,590..834 OrfName  

=item --suffix,-s F<STRING> [Optional]

Suffix of the files to be read in. If suffix is not provided, all the files in input directory will be read in.

=item --outputDir,-o F<FILE> [Optional]

Path of the directory to storage result files. If NOT provided, the program will generate a folder automatically.

=item --threads,-t F<INTEGER> [Optional]

Number of threads this program will use when computing. A positive integer should be provided. The default value is 1.

=item --help,-h

Display this detailed help information.

=back

=head1 EXAMPLE

=over 5

qap CutSeqWithIntervals -i ./seq -t 10 -s fasta -v ./intervals.txt -o ./cut

=back

=head1 AUTHOR

Mingjie Dr.Wang I<huzai@sjtu.edu.cn>

=head1 COPYRIGHT

Copyright (C) 2017, Mingjie Wang. All rights reserved.

