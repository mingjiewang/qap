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

####---------------------------####
####The program begins here
####---------------------------####

##Show welcome
print "You are now running subprogram: ";
printcol ("OTUBarplot","green");
print "\n";

## check threads available or not
$| = 1;
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
my $legend;
my $barWidth;
my $relative;
my $sort;

my $DateNow = `date +"%Y%m%d_%Hh%Mm%Ss"`;
chomp $DateNow;

GetOptions(
'i|inputFile|=s'     => \$inputFile,
'o|outputDir|=s'     => \$outputDir,
'h|help|'            => \$help,
'w|barWidth|=s'      => \$barWidth,
'l|figLegend|=s'     => \$legend,
'r|relative|=s'      => \$relative,
's|sort|=s'          => \$sort
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
		my $cmd = "mkdir -p $outputDir";
		system($cmd);
	}
}else{
	$outputDir = File::Spec -> catfile($wk_dir,"qap_Results_for_OTUBarplot_$DateNow");
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

if(defined $inputFile){
	if(not existFile($inputFile)){
		InfoError("Input file $inputFile does not exist.");
		exit;
	}
	$inputFile = abs_path($inputFile);
}else{
	InfoError("Input OTU table MUST be provided.");
	pod2usage(-verbose=>0,-exitval=>1);
	exit;
}

if(defined $barWidth){
	if(CheckDouble($barWidth) and $barWidth >= 0 and $barWidth <= 1){
		#nothing
	}else{
		InfoError("--barWidth/-w should be a decimal between 0 and 1.");
		pod2usage(-verbose=>0,-exitval=>1);
		exit;
	}
}else{
	$barWidth = 1;
}

if(defined $legend){
	if($legend =~ /^y/i){
		$legend = "Y";
	}elsif($legend =~ /^n/i){
		$legend = "N";
	}else{
		InfoError("--figLegend/-l should be \'Y\' or \'N\'.");
		pod2usage(-verbose=>0,-exitval=>1);
		exit;
	}
}else{
	$legend = 'N';
}

if(defined $relative){
	if($relative =~ /^y/i){
		$relative = "Y";
	}elsif($relative =~ /^n/i){
		$relative = "N";
	}else{
		InfoError("--relative should be \'Y\' or \'N\'.");
		pod2usage(-verbose=>0,-exitval=>1);
		exit;
	}
}else{
	$relative = 'N';
}

if(defined $sort){
	if($sort =~ /^y/i){
		$sort = "Y";
	}elsif($sort =~ /^n/i){
		$sort = "N";
	}else{
		InfoError("--sort should be \'Y\' or \'N\'.");
		pod2usage(-verbose=>0,-exitval=>1);
		exit;
	}
}else{
	$sort = 'N';
}

sleep(1);

##start to calculate
#rscript
my $rscript = File::Spec -> catfile($RealBin,"Rscripts","OTUBarplot.R");
if(not existFile($rscript)){
	InfoError("R script $rscript is missing. Please check. Exiting...");
	exit(0);
}

#run r script
Info("Calculating relative abundance and drawing OTU barplot.");
my $outfig = File::Spec -> catfile($outputDir,"OTUBarplot");
my $cmd = "Rscript $rscript --inputFile $inputFile --legend $legend --barWidth $barWidth --outputFile $outfig --relative $relative --sort $sort";
sleep(1);
runcmd($cmd);
sleep(2);

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
           



qap OTUBarplot [options]

Use --help to see more information.

qap is still in development. If you have encounted any problem in usage, please feel no hesitation to cotact us.

=head1 DESCRIPTION

This script implements a function to draw bar plot based on QS abundance. The script has B<several> mandatory options that MUST appear last. 

=head1 OPTIONS

=over 5

=item --inputFile,-i F<FILE> [Required]

Path to OTU table to be read in. A normalized OTU count table is preferred.

=item --barWidth,-w F<DOUBLE> [Optional]

A decimal between 0 and 1 to specify the width of bars in the barplot, which 0 refers to a narrow bar and 1 refers to a wide bar. Default value is 1.

=item --figLegend,-l F<BOOLEAN> [Optional]

Whether to draw figure legends or not. Choose between 'Y'(Yes) and 'N'(No). Default value is 'N'.

=item --relative,-r F<BOOLEAN> [Optional]

Whether draw bars with relative height (Total height normalized to 1.0), or using abunance value as true height. Choose between 'Y'(Yes) and 'N'(No). Default value is 'N'.

=item --sort,-s F<BOOLEAN> [Optional]

Whether sort bars with heights. Choose between 'Y'(Yes) and 'N'(No). If set to 'Y', bars will be sorted from height to low, and OTUs within bars will also be sorted by abudances.

=item --outputDir,-o F<FILE> [Optional]

Path of the directory to storage result files. If NOT provided, the program will generate a folder automatically.

=item --help,-h

Display this detailed help information.

=back

=head1 EXAMPLE

=over 5

qap OTUBarplot -i ./OTUTable.txt -w 1 -l T -o ./result

=back

=head1 AUTHOR

Mingjie Dr.Wang I<huzai@sjtu.edu.cn>

=head1 COPYRIGHT

Copyright (C) 2017, Mingjie Wang. All rights reserved.

