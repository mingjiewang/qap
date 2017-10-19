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
printcol ("SampleRename","green");
print "\n";


####define command line arguments
my $help;
my $sampleName;
my $sampleDir;
my $DateNow = `date +"%Y%m%d_%Hh%Mm%Ss"`;
chomp $DateNow;

GetOptions(
's|sampleName|=s'    => \$sampleName,
'd|sampleDir|=s'     => \$sampleDir,
'h|help|'            => \$help,
);

##check command line arguments
if (defined $help){
	pod2usage(-verbose=>2,-exitval=>1);
}

if(scalar(@ARGV) == 0){
	pod2usage(-verbose=>1,-exitval=>1);
}

if (defined $sampleName){
	if (not -e $sampleName){
 		InfoError("The file of sample names $sampleName does NOT exist.",'red');
 		exit;
	}
	$sampleName = abs_path($sampleName);
}else{
	InfoError("The file of sample names MUST be provided!",'red');
	pod2usage(-verbose=>2,-exitval=>1);
	exit;
}

if (defined $sampleDir){
	$sampleDir =~ s/\/$//;
	$sampleDir = abs_path($sampleDir) . "/";
	if (not -e $sampleDir){
 		InfoError("The directory for sample raw data $sampleDir does NOT exist.",'red');
 		exit;
	}
	$sampleDir = abs_path($sampleDir);
}else{
	InfoError("The directory for sample raw data MUST be provided!",'red');
	pod2usage(-verbose=>2,-exitval=>1);
	exit;
}


##read in name file
open N,$sampleName or die "Can NOT open sample Name file $sampleName:$!\n";
while (<N>){
	chomp;
	next if $_ !~ /^\w/;
	my $name1;
	my $name2;
	if (/(.*?)\s+(.*)/){
		$name1 = $1;
		$name2 = $2;
	}
	my $file1 = File::Spec -> catfile($sampleDir,$name1);
	my $file2 = File::Spec -> catfile($sampleDir,$name2);

	#renaming
	my $cmd = "mv $file1 $file2";
	Info("Renaming $name1 to $name2");
	system "$cmd";
}

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
           



qap SampleRename [options]

Use --help to see more information.

qap is still in development. If you have encounted any problem in usage, please feel no hesitation to cotact us.

=head1 DESCRIPTION

This script implements a function for renaming files in batch. The script has B<two> mandatory options that MUST appear last. 

First is the directory of the file contaning the name of all data files. Second is the directory storaging all data files. Both options are case B<sensitive> and are listed as follows:

=head1 OPTIONS

=over 5

=item --sampleName,-s F<FILE> [Required]

Path to the plain text file contaning all the file names to convert. The name file MUST be formatted accordingly: 1. No headline; 2. Two columns, the first contains the old names and the second contains the third name; 

3. Columns seperated by tabs; 4. Both the old names and new names should NOT contain space.

=item --sampleDir,-d F<FILE> [Required]

Path to the directory containg all the data files to rename.

=item --help,-h

Display this detailed help information.

=back

=head1 EXAMPLE

=over 5

qap rename -s sampleName.txt -d ./data/

=back

=head1 AUTHOR

Mingjie Dr.Wang I<huzai@sjtu.edu.cn>

=head1 COPYRIGHT

Copyright (C) 2017, Mingjie Wang. All rights reserved.


