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
printcol ("SampleCorrelation","green");
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
my $heatColor;
my $withLabel;
my $corTable;
my $triangle;
my $sort;
my $dendrogram;
my $line;

my $DateNow = `date +"%Y%m%d_%Hh%Mm%Ss"`;
chomp $DateNow;

GetOptions(
'i|inputFile|=s'     => \$inputFile,
'o|outputDir|=s'     => \$outputDir,
'h|help|'            => \$help,
'c|color|=s'         => \$heatColor,
'l|withLabel|=s'     => \$withLabel,
'b|corTable|=s'      => \$corTable,
'a|triangle|=s'      => \$triangle,
's|sort|=s'          => \$sort,
'd|dendrogram|=s'    => \$dendrogram,
'e|line|=s'          => \$line
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
	$outputDir = File::Spec -> catfile($wk_dir,"qap_Results_for_SampleCorrelation_$DateNow");
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

if(defined $heatColor){
	if($heatColor =~ /bluered/i){
		$heatColor = "blue";
	}elsif($heatColor =~ /greenred/i){
		$heatColor = "green";
	}elsif($heatColor =~ /pinkred/i){
		$heatColor = "pink";
	}else{
		InfoError("Input heatmap color -c/--heatColor should be one of \'BlueRed\' or \'GreenRed\' or \'PinkRed\'.");
		pod2usage(-verbose=>0,-exitval=>1);
		exit;
	}
}else{
	$heatColor = "blue";
}

if(defined $withLabel){
	if($withLabel =~ /^y/i){
		$withLabel = "Y";
	}elsif($withLabel =~ /^n/i){
		$withLabel = "N";
	}else{
		InfoError("--withLabel/-l should be \'Y\' or \'N\'.");
		pod2usage(-verbose=>0,-exitval=>1);
		exit;
	}
}else{
	$withLabel = "Y";
}

if(defined $corTable){
	if($corTable =~ /^y/i){
		$corTable = "Y";
	}elsif($corTable =~ /^n/i){
		$corTable = "N";
	}else{
		InfoError("--$corTable/-b should be \'Y\' or \'N\'.");
		pod2usage(-verbose=>0,-exitval=>1);
		exit;
	}
}else{
	$corTable = "Y";
}

if(defined $triangle){
	if($triangle =~ /^y/i){
		$triangle = "Y";
	}elsif($triangle =~ /^n/i){
		$triangle = "N";
	}else{
		InfoError("--triangle/-a should be \'Y\' or \'N\'.");
		pod2usage(-verbose=>0,-exitval=>1);
		exit;
	}
}else{
	$triangle = "Y";
}

if(defined $sort){
	if($sort =~ /^y/i){
		$sort = "Y";
	}elsif($sort =~ /^n/i){
		$sort = "N";
	}else{
		InfoError("--sort/-s should be \'Y\' or \'N\'.");
		pod2usage(-verbose=>0,-exitval=>1);
		exit;
	}
}else{
	$sort = "Y";
}

if(defined $dendrogram){
	if($dendrogram =~ /^y/i){
		$dendrogram = "Y";
	}elsif($dendrogram =~ /^n/i){
		$dendrogram = "N";
	}else{
		InfoError("--dendrogram/-d should be \'Y\' or \'N\'.");
		pod2usage(-verbose=>0,-exitval=>1);
		exit;
	}
}else{
	$dendrogram = "Y";
}

if(defined $line){
	if($line =~ /^y/i){
		$line = "Y";
	}elsif($line =~ /^n/i){
		$line = "N";
	}else{
		InfoError("--line/-e should be \'Y\' or \'N\'.");
		pod2usage(-verbose=>0,-exitval=>1);
		exit;
	}
}else{
	$line = "Y";
}

sleep(1);

##start to calculate
#rscript
my $rscript = File::Spec -> catfile($RealBin,"Rscripts","SampleCorrelation.R");
if(not existFile($rscript)){
	InfoError("R script $rscript is missing. Please check. Exiting...");
	exit(0);
}

#run r script
Info("Calculating sample correlations and drawing pairwise heatmap.");
my $outfig = File::Spec -> catfile($outputDir,"SampleCorrelation");
my $cmd = "Rscript $rscript --inputFile $inputFile --heatColor $heatColor --corTable $corTable --withLabel $withLabel --outputFile $outfig --triangle $triangle --sort $sort --line $line --dendrogram $dendrogram ";
runcmd($cmd);

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
           



qap SampleCorrelation [options]

Use --help to see more information.

qap is still in development. If you have encounted any problem in usage, please feel no hesitation to cotact us.

=head1 DESCRIPTION

This script implements a function to calculate correlations among sampels and draw heatmap. The script has B<several> mandatory options that MUST appear last. 

=head1 OPTIONS

=over 5

=item --inputFile,-i F<FILE> [Required]

Path to OTU table to be read in. A normalized OTU count table is preferred.

=item --color,-c F<STRING> [Optional]

Color used to draw heatmap. Choose between 'BlueRed', 'GreenRed', or 'PinkRed'. Default value is 'BlueRed'.

=item --withLabel,-l F<BOOLEAN> [Optional]

Whether add correlation labels to the figure or not. Choose bewtwen 'Y' for 'Yes' and 'N' for 'No'. Default value is 'Y'.

=item --corTable,-b F<BOOLEAN> [Optional]

Whether output correlation table or not. Choose between 'Y' for 'Yes' and 'N' for 'No'. Default value is 'Y'.

=item --triangle,-a F<BOOLEAN> [Optional]

Whether draw a triangle or a rectangle of sample correlation matrix. Choose between 'Y' for 'Yes' and 'N' for 'No'. Default value is 'Y'.

=item --sort,-s F<BOOLEAN> [Optional]

Whether reorder the samples by distance or not. Choose between 'Y' for 'Yes' and 'N' for 'No'. Default value is 'Y'.

=item --dendrogram,-d F<BOOLEAN> [Optional]

Whether draw sample dendrogram or not. Choose between 'Y' for 'Yes' and 'N' for 'No'. Default value is 'Y'.

=item --line,-e F<BOOLEAN> [Optional]

Whether draw lines to on colored correlation matrix or not. Choose between 'Y' for 'Yes' and 'N' for 'No'. Default value is 'N'.

=item --outputDir,-o F<FILE> [Optional]

Path of the directory to storage result files. If NOT provided, the program will generate a folder automatically.

=item --help,-h

Display this detailed help information.

=back

=head1 EXAMPLE

=over 5

qap SampleCorrelation -i ./OTUTable.txt -g G1,G1,G2,G2,G2,G3,G3 -o ./result

=back

=head1 AUTHOR

Mingjie Dr.Wang I<huzai@sjtu.edu.cn>

=head1 COPYRIGHT

Copyright (C) 2017, Mingjie Wang. All rights reserved.

