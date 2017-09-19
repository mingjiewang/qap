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
printcol ("SamplePCA","green");
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
my $inputFile;
my $outputDir;
my $group;
my $dimension;
my $plotStyle;

my $DateNow = `date +"%Y%m%d_%Hh%Mm%Ss"`;
chomp $DateNow;

GetOptions(
'i|inputFile|=s'     => \$inputFile,
'o|outputDir|=s'     => \$outputDir,
'h|help|'            => \$help,
'g|group|=s'         => \$group,
'd|dimension|=s'     => \$dimension,
's|plotStyle|=s'     => \$plotStyle
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
		my $cmd = "mkdir -p $outputDir";
		system($cmd);
	}
}else{
	$outputDir = File::Spec -> catfile($wk_dir,"qap_Results_for_SamplePCA_$DateNow");
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

if(defined $plotStyle){
	if($dimension == 3){
		InfoError("--plotStyle/-s is ONLY available when --dimension/-d 2");
		pod2usage(-verbose=>0,-exitval=>1);
		exit;
	}
	
	if($plotStyle == 1 or $plotStyle == 2 or $plotStyle == 3){
		#nothing
	}else{
		InfoError("Specify plot style with --plotStyle/-s between \'1\' or \'2\' or \'3\'.");
		pod2usage(-verbose=>0,-exitval=>1);
		exit;
	}
}else{
	$plotStyle = 1;
}

if(defined $dimension){
	if($dimension == 2 or $dimension == 3){
		#nothing
	}else{
		InfoError("Specify plot dimesions with --dimension/-d between \'2\' or \'3\'.");
		pod2usage(-verbose=>0,-exitval=>1);
		exit;
	}
}else{
	$dimension = 3;
}



if(not defined $group){
	InfoError("Sample group information MUST be provided.");
	pod2usage(-verbose=>0,-exitval=>1);
	exit;
}

sleep(1);

##start to calculate
#rscript
my $rscript = File::Spec -> catfile($RealBin,"Rscripts","SamplePCA.R");
if(not existFile($rscript)){
	InfoError("R script $rscript is missing. Please check. Exiting...");
	exit(0);
}

#run r script
Info("Running sample principle component analysis and drawing PCA scatter plot.");
my $outfig = File::Spec -> catfile($outputDir,"SamplePCA");
my $cmd = "Rscript $rscript --inputFile $inputFile --group $group --figStyle $plotStyle --dimension $dimension --outputFile $outfig";
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
           



qap SamplePCA [options]

Use --help to see more information.

qap is still in development. If you have encounted any problem in usage, please feel no hesitation to cotact us.

=head1 DESCRIPTION

This script implements a function to perform PCA on multiple samples and draw PCA scatter plot. The script has B<several> mandatory options that MUST appear last. 

=head1 OPTIONS

=over 5

=item --inputFile,-i F<FILE> [Required]

Path to OTU table to be read in. A normalized OTU count table is preferred.

=item --group,-g F<STRING> [Required]

Classification information for each sample. The order of given groups should be same with samples in OTU table. Please seperate each sample by comma, e.g. -g case,case,case,control,control,control.

=item --dimension,-d F<INTEGER> [Optional]

Whether draw 3-dimension PCA plot or 2-dimension PCA biplot. Choose between '2' or '3'. Default value is '3'.

=item --plotStyle,-s F<INTEGER> [Optional]

If --dimension/-d 3 is specified, use this parameter to specify the style of output 2-dimension PCA biplot. Choose among '1', or '2' or '3'. '1' means no sample label and no group circle; '2' means no sample label but with group circle; '3' means with sample label but no group circle. These plot style might be quite different, choose a suitable one. Default value is 1.

=item --outputDir,-o F<FILE> [Optional]

Path of the directory to storage result files. If NOT provided, the program will generate a folder automatically.

=item --help,-h

Display this detailed help information.

=back

=head1 EXAMPLE

=over 5

qap SamplePCA -i ./OTUTable.txt -l T -s 2 -d 2 -g G1,G1,G2,G2,G2,G3,G3 -o ./result

=back

=head1 AUTHOR

Mingjie Dr.Wang I<huzai@sjtu.edu.cn>

=head1 COPYRIGHT

Copyright (C) 2017, Mingjie Wang. All rights reserved.

