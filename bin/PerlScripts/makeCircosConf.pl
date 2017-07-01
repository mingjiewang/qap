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

####Flash caches####
$| = 1;

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
my $outputfile;

my $DateNow = `date +"%Y%m%d_%Hh%Mm%Ss"`;
chomp $DateNow;

GetOptions(
'o|outputFile=s'     => \$outputfile,
'h|help|'            => \$help
);


##check command line argments
if (defined $outputfile){
	if(existFile($outputfile)){
		InfoError("The output configuration file $outputfile already exists. Please specify another output file using -o/--outputFile.");
		pod2usage(-verbose=>1,-exitval=>1);
		exit;
	}
}else{
	InfoError("Output configure file is not defined. Please specify the output file using option -o/--outputFile\n");
	pod2usage(-verbose=>1,-exitval=>1);
	exit;
}

##output file
open CONF,">$outputfile" or die "Can not output to file $outputfile:$!";
print CONF "#########################################\n";
print CONF "####                                 ####\n";
print CONF "##    CONFIGURATION FILE FOR CIRCOS    ##\n";
print CONF "####                                 ####\n";
print CONF "#########################################\n\n";


##start to ask questions
Info("Before running circos, there are several questions to answer.\n");
sleep(1);

my $qNum = 1;

my $bedfile = File::Spec -> catfile(dirname($outputfile), basename($outputfile) . "." . time() . ".bed");
open BED,">$bedfile" or die "Can NOT output to bed file:$!";

my $circular = 0;
my @chrname;
my @chrlength;

while($circular != 1 or $circular != 2){
	print "${qNum}. Is the virus circular[1] or linear[2]?\n";
	print "Your answer: ";
	chomp (my $circular = <STDIN>);
	
	if ($circular !~ /[12]/){
		print("ERROR! Please select between \'1\' or \'2\'.\n");
		redo;
	}
	
	if($circular == 1){
		$qNum++;
		print "\n${qNum}. Name of your virus?\n";
		print "Your answer: ";
		chomp (my $name = <STDIN>);
		push @chrname,$name;
		
		my $length = 0;
		while($length <= 1 || $length =~ /\./){
			$qNum++;
			print "\n${qNum}. Length of virus genome?\n";
			print "Your answer: ";
			chomp ($length = <STDIN>);
			
			if($length >= 1 and $length !~ /\./){
				push @chrlength,$length;
				last;
			}
			
			print "ERROR! Genome length should be a positive integer greater than 1\n";
		}
		
		print BED "$name\t1\t$length\n";
		
		last;
	}elsif($circular == 2){
		my $numberofchr = 0;
		while($numberofchr < 1 || $numberofchr =~ /\./){
			$qNum++;
			print "\n${qNum}. How many fragments(chromosomes) does the virus have?\n";
			print "Your answer: ";
			
			chomp($numberofchr = <STDIN>);
			
			last if $numberofchr >= 1 and $numberofchr !~ /\./;
			
			print "ERROR! Number of fragments should be a positive integer.\n";
		}
		
		$qNum++;
		print "\n${qNum}. Information about the fragment(s)/chromosome(s) of virus genome.\n";
		
		for my $i (1..$numberofchr){
			my $ith;
			if ($i == 1){
				$ith = '1st';
			}elsif($i == 2){
				$ith = '2nd';
			}elsif($i == 3){
				$ith = '3rd';
			}else{
				$ith = $i . "th";
			}
			
			print "--${qNum}.$i Name and length of the $ith fragment/chromosome? (name,length)\n";
			print "--Your answer: ";
			chomp (my $answer = <STDIN>);
			
			if($answer !~ /,/){
				print("ERROR! Wrong input. Input should be like \'name,length\'.\n\n");
				redo;
			}
			
			my @tmp = split ",",$answer;
			my $name = $tmp[0];
			if(not $name){
				print("ERROR! Invalid name provided.\n\n");
				redo;
			}
			
			my $length = $tmp[1];
			if(not $length){
				print("ERROR! Invalid length provided.\n\n");
				redo;
			}
			
			if ($length !~ /[0-9]/){
				print("ERROR! Number of fragments should be a positive integer.\n\n");
				redo;
			}
			
			if ($length <= 1 or $length =~ /[\.A-Za-z]/){
				print("ERROR! Number of fragments should be a positive integer greater than 1.\n\n");
				redo;
			}
			
			push @chrname,$name;
			push @chrlength,$length;
			
			print BED "$name\t1\t$length\n";
		}
		
		print CONF "\nbedfile     = $bedfile\n";
		
		last;
	}else{
		print("ERROR! Please select between \'1\' or \'2\'.\n");
	}
	
	
}
close BED;

#number of tracks
$qNum++;
my $numberoftrack = 0;
while ($numberoftrack <= 0 || $numberoftrack > 6){
	print "\n${qNum}. How many data tracks do you want to show?\n";
	print "Your answer: ";
	chomp ($numberoftrack = <STDIN>);
	
	if(not CheckPositiveInt($numberoftrack)){
		print "ERROR! Number of tracks should be a positive integer.\n";
		redo;
	}
	
	if($numberoftrack > 6){
		print("ERROR! Tracks should not be more than 6.\n");
		redo;
	}
}


#track data
$qNum++;
print "\n${qNum}. Path to the data track files.\n";

my $pathToData;
my @trackData;
my $dataType;
my @dataType;
my $numberOfLinkData = 0;
TRACKDATA: for my $i (1..$numberoftrack){
	my $ith;
	if ($i == 1){
		$ith = '1st';
	}elsif($i == 2){
		$ith = '2nd';
	}elsif($i == 3){
		$ith = '3rd';
	}else{
		$ith = $i . "th";
	}
	
	print "--${qNum}.$i Path to the $ith data file?\n";
	print "--Your answer: ";
	chomp (my $answer = <STDIN>);
	
	#check track data exist
	if (not existFile($answer)){
		print("ERROR! The data file you provided does NOT exist.\n\n");
		redo;
	}
	
	#check track data format
	my $feature = 0;
	FEATURE: while($feature != 1 or $feature != 2){
		print "----${qNum}.$i.1 Does this data show features of links[1] or regions[2]?\n";
		print "----Your answer: ";
		chomp ($feature = <STDIN>);
		
		if(defined $feature){
			if ($feature =~ /[^12]/){
				print("ERROR! Please select between \'1\' or \'2\'.\n");
				redo FEATURE;
			}
		}else{
			print("ERROR! Please select between \'1\' or \'2\'.\n");
			redo FEATURE;
		}
		
		if($feature == 1){
			if($numberOfLinkData > 1){
				print("ERROR! You have already provided a data file for links. Only one linkage data is allowed.\n\n");
				redo TRACKDATA;
			}
			
			if($i < $numberoftrack){
				print("ERROR! Please make sure the data file for links goes last.\n");
				redo TRACKDATA;
			}
			
			print "[+]Checking link data format...";
	
			sleep(1);
			my $res = &checkLinkData($answer, \@chrname, \@chrlength);
			if ($res eq 'ok'){
				print "PASS\n\n";
				$dataType .= 'link,';
				push @dataType,'link';
				$numberOfLinkData++;
				
				sleep 1;
				last FEATURE;
			}else{
				print "FAILED\n";
				sleep 1;
				print "ERROR DETAIL: $res\n\n";
				sleep 1;
				redo TRACKDATA;
			}
		}elsif($feature == 2){
			print "[+]Checking region data format...";
			
			sleep(1);
			my $res = &checkRegionData($answer, \@chrname, \@chrlength);
			if ($res eq 'ok'){
				print "PASS\n\n";
				$dataType .= 'region,';
				push @dataType,'region';
				
				sleep 1;
				last FEATURE;
			}else{
				print "FAILED\n";
				sleep 1;
				print "ERROR DETAIL: $res\n\n";
				sleep 1;
				redo TRACKDATA;
			}
		}else{
			print("ERROR! Please select between \'1\' or \'2\'.\n");
			redo FEATURE;
		}
	}
	
	
	
	$answer = abs_path($answer);
	push @trackData,$answer;
	$answer .= ",";
	$pathToData .= $answer; 
}

$pathToData =~ s/,$//;
$dataType =~ s/,$//;

print CONF "trackData   = $pathToData\n";
print CONF "dataType    = $dataType\n";

#how to show track data
$qNum++;
print "${qNum}. Plot type to show track data.\n";

my @plotType;
for my $i (1..$numberoftrack){
	my $ith;
	if ($i == 1){
		$ith = '1st';
	}elsif($i == 2){
		$ith = '2nd';
	}elsif($i == 3){
		$ith = '3rd';
	}else{
		$ith = $i . "th";
	}
	
	my $filename = basename($trackData[$i - 1]);
	
	if($dataType[$i - 1] eq 'link'){
		print "----${qNum}.$i Data \<$filename\> in the $ith track is link data and will be shown as links or ribbons.\n";
		push @plotType,'link';
		sleep 1;
		next;
	}
	
	my $typechoice = 0;
	while($typechoice != 1 or $typechoice != 2 or $typechoice != 3){
		print "----${qNum}.$i To show \<$filename\> in the $ith track, choose a plot type: [1]Heatmap; [2]Histogram; [3]Line plot.\n";
		print "----Your answer: ";
		chomp ($typechoice = <STDIN>);
		
		if($typechoice =~ /[^123]/){
			print ("Error! Please choose plots among \'1\' or \'2\' or \'3\'.\n\n");
			redo;
		}
		
		my $plot;
		if($typechoice == 1){
			$plot = 'heatmap';
		}elsif($typechoice == 2){
			$plot = 'histogram';
		}elsif($typechoice == 3){
			$plot = 'line'
		}else{
			print ("Error! Please choose plots among \'1\' or \'2\' or \'3\'.\n\n");
			redo;
		}
		
		push @plotType,$plot;
		last;
	}

}

my $plotType = join ",",@plotType;
print CONF "plotType    = $plotType\n";
close CONF;

print "\n";
Info("Completed! Configurations written into $outputfile.");
print "\n";

##sub programs starts here
sub checkLinkData {
	my $file = shift;
	my $name = shift;
	my $len = shift;
	
	my @name = @$name;
	my @len = @$len;
	if (scalar(@name) != scalar(@len)){
		return 'Number of chromosome(s)/fragment(s) names and lengths does NOT equal!';
	}
	
	my %len;
	for my $i (0..scalar(@len) - 1){
		$len{$name[$i]} = $len[$i];
	}
	
	open T,"$file" or die "Can NOT open data track file $file:$!";
	my $line;
	while(<T>){
		chomp;
		$line++;
		
		next if $_ =~ /^$/;
		
		my @tmp = split /\s+/,$_;
		my $numberofele = scalar(@tmp);
		
		if ($numberofele == 6){
			#nothing
		}else{
			return "$numberofele data columns found at LINE $line, which should be 6.";
		}
		
		my $chr1 = $tmp[0];
		my $start1 = $tmp[1];
		my $end1 = $tmp[2];
		my $chr2 = $tmp[3];
		my $start2 = $tmp[4];
		my $end2 = $tmp[5];
		
		if (isInArray($chr1, $name) == 0 or isInArray($chr2, $name) == 0){
			return "Name of fragment/chromosome is not valid at LINE $line.";
		}
		
		if(not(CheckPositiveInt($start1) &&  CheckPositiveInt($end1) && CheckPositiveInt($start2) && CheckPositiveInt($end2))){
			return "Start/End position at LINE $line contains non-positive integer.";
		}
		
		if($start1 > $len{$chr1} or $start2 > $len{$chr2} or $end1 > $len{$chr1} or $end2 > $len{$chr2}){
			return "Start/End position at LINE $line is greater than fragment/chromosome length.";
		}
	}
	close T;
	
	return 'ok';
}

sub checkRegionData {
	my $file = shift;
	my $name = shift;
	my $len = shift;
	
	my @name = @$name;
	my @len = @$len;
	if (scalar(@name) != scalar(@len)){
		return 'Number of chromosome(s)/fragment(s) names and lengths does NOT equal!';
	}
	
	my %len;
	for my $i (0..scalar(@len) - 1){
		$len{$name[$i]} = $len[$i];
	}
	
	open T,"$file" or die "Can NOT open data track file $file:$!";
	my $line;
	while(<T>){
		chomp;
		$line++;
		
		my @tmp = split /\s+/,$_;
		my $numberofele = scalar(@tmp);
		
		if ($numberofele == 4){
			#nothing
		}else{
			return "$numberofele data columns found at LINE $line, which should be 4.";
		}
		
		my $chr1 = $tmp[0];
		my $start1 = $tmp[1];
		my $end1 = $tmp[2];
		my $value = $tmp[3];
		
		if (isInArray($chr1, $name) == 0){
			return "Name of fragment/chromosome is not valid at LINE $line.";
		}
		
		if(CheckPositiveInt($start1) == 0 or CheckPositiveInt($end1) == 0){
			return "Start/End position at LINE $line is not positive integer.";
		}
		
		if($start1 > $len{$chr1} or $end1 > $len{$chr1}){
			return "Start/End position at LINE $line is greater than fragment/chromosome length.";
		}
	}
	close T;
	
	return 'ok';
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
           



perl makeCircosConf.pl [options]

Use --help to see more information.

gap is still in development. If you have encounted any problem in usage, please feel no hesitation to cotact us.

=head1 DESCRIPTION

This script implements a function to make configuration files for subprogram Circos. The script has B<several> mandatory options that MUST appear last. 

=head1 OPTIONS

=over 5

Path to the reference fasta file. 

=item --outputFile,-o F<FILE> [Optional]

Path of the output configuration file. If NOT provided, the program will generate a file automatically.

=item --help,-h

Display this detailed help information.

=back

=head1 EXAMPLE

=over 5

perl makeCircosConf.pl -o circos.conf

=back

=head1 AUTHOR

Mingjie Dr.Wang I<huzai@sjtu.edu.cn>

=head1 COPYRIGHT

Copyright (C) 2017, Mingjie Wang. All rights reserved.

