#!/usr/bin/env perl

##Use third party modules
use diagnostics;
use strict;
use warnings;
use Getopt::Long;

## check threads available or not
$| = 1;

my $inputfile;
my $outputfile;

GetOptions(
'i|inputFile|=s'  => \$inputfile,
'o|outputFile|=s' => \$outputfile,
);

##core program
open OUT,">$outputfile" or die "Can NOT output to file $outputfile:$!";

open VCF,$inputfile or die "Can NOT open vcf file $inputfile:$!";
while(<VCF>){
	chomp;
	if(/^\#/){
		next;
	}
	my @arr = split "\t",$_;
	
	my ($chr,$pos,$refbase,$mutbase,$info,$genotype,$genotypeInfo) = @arr[0,1,3,4,7,8,9];
	
	my $freq = "NA";
	my $depth = "NA";
	if($info =~ /AF=([0-9\.]+).*?DP=(\d+)/){
		#lofreq and GATK output both have AF field, however varscan does not. 
		#GATK format
		$freq = $1;
		$depth = $2;
	}elsif($info =~ /DP=(\d+).*?AF=([0-9\.]+)/){
		#Lofreq format
		$depth = $1;
		$freq = $2;
	}else{
		#varscan have the FREQ written in the genotype field.
		#VarScan format
		if($genotype =~ /DP.*?FREQ/){
			my @tmp = split ':',$genotypeInfo;
			$depth = $tmp[3];
			$freq = $tmp[6];
			
			$freq =~ s/\%//;
			$freq /= 100;
		}
	}
	
	print OUT "$chr\t$pos\t$refbase\t$mutbase\t$freq\t$depth\n";
	
}
close VCF;
close OUT;


