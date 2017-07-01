#!/usr/bin/env perl

#################################################################################
##                                                                             ##
##                       Quasispecies Analysis Package                         ##
##                                                                             ##
#################################################################################
##                                                                             ##
##  A software suite designed for virus quasispecies analysis                  ##
##  See our website: <http://bioinfo.rjh.com.cn/labs/jhuang/tools/gap/>        ##
##                                                                             ##
##  Version 1.0                                                                ##
##                                                                             ##
##  Copyright (C) 2017 by Mingjie Wang, All rights reserved.                   ##
##  Contact:  huzai@sjtu.edu.cn                                                ##
##  Organization: Research Laboratory of Clinical Virology, Rui-jin Hospital,  ##
##  Shanghai Jiao Tong University, School of Medicine                          ##
##                                                                             ##
##  This file is a subprogram of GAP suite.                                    ##
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
##  License along with ViralFusionSeq; if not, see                             ##
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
printcol ("RawDataFiltraion","green");
print "\n";

##get workding directory
my $wk_dir = getcwd;
my $mainBin;
if ($RealBin =~ /(.*)\/bin/){
	$mainBin = $1;
}

####define command line arguments
my $help;
my $fq1;
my $fq2;
my $baseQual;
my $readLen;
my $outputDir;
my $pairEnded = 0;
my $maxN;
my $DateNow = `date +"%Y%m%d_%Hh%Mm%Ss"`;
chomp $DateNow;

GetOptions(
'1|fastq1|=s'      => \$fq1,
'2|fastq2|=s'      => \$fq2,
'q|baseQuality=s'  => \$baseQual,
'l|readLength|=s'  => \$readLen,
'o|outputDir|=s'   => \$outputDir,
'n|maxN|=s'        => \$maxN,
'h|help|'          => \$help,
);

##check command line arguments
if (defined $help){
	pod2usage(-verbose=>2,-exitval=>1);
}

if (defined $outputDir){
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
	$outputDir = File::Spec -> catfile($wk_dir,"qap_Results_for_RawDataFiltration_$DateNow");
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


if (defined $fq1){
	if (not -e $fq1){
		InfoError("Input fastq1 file $fq1 does NOT exist!",'red');
		pod2usage(-verbose=>0,-exitval=>1);
		exit;
	}
	$fq1 = abs_path($fq1);
}else{
	InfoError("Input fastq1 file must be provided!",'red');
	pod2usage(-verbose=>1,-exitval=>1);
	exit;
}

if (defined $fq2){
	if (not -e $fq2){
		InfoError("Input fastq2 file $fq2 does NOT exist!",'red');
		pod2usage(-verbose=>0,-exitval=>1);
		exit;
	}
	$fq2 = abs_path($fq2);
	$pairEnded = 1;
}


if (defined $baseQual){
	if(CheckPositiveInt($baseQual)){
		if($baseQual < 20){
			InfoWarn("You are using $baseQual for base quality filtration. It is too low, maybe a higher value would be better.","yellow");
		}else{
			#nothing
		}
	}else{
		InfoError("The base quanlity filtration cutoff value MUST be a positive integer.","red");
		exit;
	}
}else{
	InfoWarn("The base quanlity filtration cutoff value is not provided. The program will use \<30\> as default value.");
	$baseQual = 30;
}

if (defined $readLen){
	if(CheckPositiveInt($readLen)){
		if($readLen < 100){
			InfoWarn("You are using $readLen for read length filtration. It is too low, maybe a higher value would be better.","yellow");
		}else{
			#nothing
		}
	}else{
		InfoError("The read length filtration cutoff value MUST be a positive integer.","red");
		exit;
	}
}else{
	InfoWarn("The read length filtration cutoff value is not provided. The program will use \<250\> as default value.");
	$readLen = 250;
}

if (defined $maxN){
	if(CheckPositiveInt($maxN)){
		if ($maxN > $readLen){
			InfoError("The maximum number of base \'N\' allowed is to large. Please provide another one.");
			exit;
		}
	}else{
		InfoError("The number of base \'N\' should be a positive integer.");
		pod2usage(-verbose=>2,-exitval=>1);
		exit;
	}
}else{
	InfoWarn("The maximum number of base \'N\' is not provided. The program will use -n/--maxN 0 as default.");
	$maxN = 0;
}


##run cutadapt
#handle the naming of data files
my $fq1Name = basename($fq1);
my $fq2Name = basename($fq2) if $pairEnded;
my $tag = "fil";
my $fq1FilterFileName = addTagForRawData($fq1Name,$tag,1);
my $fq2FilterFileName = addTagForRawData($fq2Name,$tag,1) if $pairEnded;

#outputFile
my $fq1FilterFile = File::Spec -> catfile($outputDir, $fq1FilterFileName);
my $fq2FilterFile = File::Spec -> catfile($outputDir, $fq2FilterFileName);
my $cmd = "cutadapt --quality-cutoff $baseQual --max-n=$maxN --minimum-length $readLen -o $fq1FilterFile -p $fq2FilterFile $fq1 $fq2";
runcmd($cmd);

##run success
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
           



gap RawDataFiltration [options]

Use --help to see more information.

gap is still in development. If you have encounted any problem in usage, please feel no hesitation to cotact us.

=head1 DESCRIPTION

This script implements a function for raw sequencin data filtration. The script has B<several> mandatory options that MUST appear last. 

=head1 OPTIONS

=over 5

=item --fastq1,-1 F<FILE> [Required]

Path to next generation sequencing raw data. REQUIRED for both single-end or paired-end reads. Both compressed files and uncompressed files are allowed.

=item --fastq2,-2 F<FILE> [Optional]

Path to next generation sequencing raw data. REQUIRED for both paired-end reads. Both compressed files and uncompressed files are allowed.

=item --baseQuality,-q F<INTEGER> [Optional, 30]

Cut off Phread value for base quality filtration. If not provided, 30 would be used as default.

=item --readLength,-l F<INTEGER> [Optional, 250]

Cut off value for read length. If not provided, 250 would be used as default.

=item --maxN,-n F<INTEGER> [Optional,0]

The number of base 'N' allowed in sequencing reads. If not provided, 0 is used.

=item --outputDir,-o F<FILE> [Optional]

Path of the directory to storage filterd data files. If NOT provided, the program will generate a folder automatically.

=item --help,-h

Display this detailed help information.

=back

=head1 EXAMPLE

=over 5

gap RawDataFiltration -1 ngs_R1.fq.gz -2 ngs_R2.fq.gz -q 30 -l 250 -o ./filter

=back

=head1 AUTHOR

Mingjie Dr.Wang I<huzai@sjtu.edu.cn>

=head1 COPYRIGHT

Copyright (C) 2017, Mingjie Wang. All rights reserved.




