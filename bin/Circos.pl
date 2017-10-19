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
use AppConfig;
use Getopt::Long;
use Pod::Usage;
use lib "$FindBin::Bin/../lib";
use Cwd qw/getcwd abs_path/;
use File::Basename;
use File::Copy;

####Use modules in this program####
use General;
use Mapper;
use Caller;

####---------------------------####
####The program begins here
####---------------------------####

##Show welcome
print "You are now running subprogram: ";
printcol ("Circos","green");
print "\n";

####Flash caches
$| = 1;

##get workding directory
my $wk_dir = getcwd;
my $mainBin;
if ($RealBin =~ /(.*)\/bin/){
	$mainBin = $1;
}

####define command line arguments
my $help;
my $outputDir;
my $confFile;

my $DateNow = `date +"%Y%m%d_%Hh%Mm%Ss"`;
chomp $DateNow;

GetOptions(
'o|outputDir|=s'     => \$outputDir,
'h|help|'            => \$help,
'c|conf|=s'          => \$confFile
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
	$outputDir = File::Spec -> catfile($wk_dir,"qap_Results_for_Circos_$DateNow");
	InfoWarn("The output directory is not provided!",'yellow');
	InfoWarn("Will mkdir \"$outputDir\" and use it as the output directory.",'yellow');
	
	if (!-e "$outputDir"){
		my $cmd = "mkdir -p $outputDir";
		system($cmd);
	}else{
		InfoError("Mkdir Failed! $outputDir already exists!","red");
		InfoError("Please specify another output directory using option -o/--outputDir\n");
		pod2usage(-verbose=>2,-exitval=>1);
		exit;
	}

}

my $confdir = File::Spec -> catfile($outputDir, 'conf');
makedir($confdir);

if(defined($confFile)){
	if($confFile =~ /null/i){
		my $makeconfScript = File::Spec -> catfile($mainBin,'bin','PerlScripts','makeCircosConf.pl');
		if(existFile($makeconfScript)){
			$confFile = File::Spec -> catfile($confdir, 'tmp.' . time() . ".conf");
			my $cmd = "perl $makeconfScript -o $confFile";
			system($cmd);
		}else{
			InfoError("Perl script $makeconfScript is missing. Please check.");
			exit(0);
		}
	}else{
		if(not existFile($confFile)){
			InfoError("Configuration file $confFile does NOT exist. Please check.");
			pod2usage(-verbose=>1,-exitval=>1);
			exit(0);
		}
	}
}else{
	InfoError("Configuration file MUST be specified.");
	pod2usage(-verbose=>1,-exitval=>1);
	exit;
}

##parse config file
sleep(1);
Info("Parsing input configuration file $confFile.");

my $config = AppConfig -> new();

#define arguments in config file
$config -> define("bedfile",{
	ARGCOUNT => AppConfig::ARGCOUNT_ONE,
	DEFAULT => "NA", #with default value
	ARGS => "=s",
});

$config -> define("trackData",{
	ARGCOUNT => AppConfig::ARGCOUNT_ONE,
	DEFAULT => "NA", #with default value
	ARGS => "=s",
});

$config -> define("dataType",{
	ARGCOUNT => AppConfig::ARGCOUNT_ONE,
	DEFAULT => "NA", #with default value
	ARGS => "=s",
});

$config -> define("plotType",{
	ARGCOUNT => AppConfig::ARGCOUNT_ONE,
	DEFAULT => "NA", #with default value
	ARGS => "=s",
});

#read config file
$config -> file($confFile);
my $bedfile = $config -> get('bedfile');
my $trackdata = $config -> get('trackData');
my $datatype = $config -> get('dataType');
my $plottype = $config -> get('plotType');

#check args from config file
if (not existFile($bedfile)){
	InfoError("The bed file $bedfile does NOT exist.");
	exit(0);
}

$trackdata =~ s/,$//;
my @trackdata = split ",",$trackdata;
my $numberOftrackdata = scalar(@trackdata);
$datatype =~ s/,$//;
my @datatype = split ",",$datatype;
my $numberOfdatatype = scalar(@datatype);
$plottype =~ s/,$//;
my @plottype = split ",",$plottype;
my $numberOfplottype = scalar(@plottype);
my @numTocheck = ($numberOfdatatype, $numberOfplottype, $numberOftrackdata);

if(not checkMultipleEqualNums(\@numTocheck)){
	InfoError("$numberOftrackdata track data, $numberOfdatatype data type and $numberOfplottype plot type detected from configuration file, which should be equal.");
	exit;
}

##start to write files used for CIRCOS program
#karyotype file
sleep(1);
Info("Writting karyotype files for Circos.");

my @col = qw/grey vdyellow blues-6-seq-1 bugn-3-seq bupu-4-seq-4 bupu-13-seq gnbu-5-seq-4 gnbu-13-seq-10 greens-6-seq-6 lred pred green pgreen lblue pblue purple ppurple orange lporange vvdpyellow vdgreen vdred greys-3-seq/;

open BED,$bedfile or die "Can NOT open $bedfile specified in configuration file:$!";
my @bed;
my %chrID;
my %chrStart;
my %chrEnd;
my %chrColor;

my $karyotypeFile = File::Spec -> catfile($confdir, "karyotype.txt");
open KA,">$karyotypeFile" or die "Can NOT output to file $karyotypeFile:$!";

my $i = 1;
while(<BED>){
	my @tmp = split /\s+/,$_;
	my ($chr,$start,$end) = @tmp;
	
	my $chrID = 'c' . $i;
	$chrID{$chr} = $chrID;
	$chrStart{$chr} = $start - 1;
	$chrEnd{$chr} = $end;
	$chrColor{$chr} = $col[($i - 1) % 23];

	print KA "chr - $chrID{$chr} $chr $chrStart{$chr} $chrEnd{$chr} $chrColor{$chr}\n";	
	$i++;
}
close KA;
close BED;

#start to write data files 
sleep(1);
Info("Preparing input track data for Circos.");

$i = 0;
my @trackdataNew;
for my $datafile (@trackdata){
	my $type = $datatype[$i];
	if($type eq 'link'){
		my $newdata = &convertToLinkData($datafile,\%chrID);
		push @trackdataNew,$newdata;
	}elsif($type eq 'region'){
		my $newdata = &convertToRegionData($datafile,\%chrID);
		push @trackdataNew,$newdata;
	}else{
		InfoError("The type of track data is illegal. Please check your configuration file.");
		exit;
	}
	
	$i++;
}

#start to write circos configuration file
sleep(1);
Info("Writing main configuration files for Circos.");

my $circosConf = File::Spec -> catfile($confdir,"circos.conf");
open CF,">$circosConf" or die "Can NOT output to configuration file $circosConf:$!";
my $ideogramConf = File::Spec -> catfile($confdir,"ideogram.conf");
open IF,">$ideogramConf" or die "Can NOT output to configuration file $ideogramConf:$!";
my $ticksConf = File::Spec -> catfile($confdir,"ticks.conf");
open TF,">$ticksConf" or die "Can NOT output to configuration file $ticksConf:$!";

print CF "karyotype = $karyotypeFile\n";
print CF "chromosomes_display_default = yes\n";
print CF "<<include $ideogramConf>>\n";
print CF "<<include $ticksConf>>\n\n";
print CF "<image>\n\n";
print CF "<<include etc/image.conf>>\n";
print CF "</image>\n";
##plots conf
print CF "<plots>\n";
for my $i (1..scalar(@plottype)){
	my $color;
	if(($i % 2) == 0){
		$color = 'blue';
	}else{
		$color = 'orange';
	}
	
	my $bgcolor;
	if(($i % 2) == 0){
		$bgcolor = 'vvlgrey';
	}else{
		$bgcolor = 'vlgrey';
	}
	
	my $ptype = $plottype[$i - 1];
	if($datatype[$i - 1] eq 'region'){
		print CF "<plot>\n";
		print CF "type = $ptype\n";
		print CF "file = $trackdataNew[$i - 1]\n\n";
		my $radius = (1 - 0.1 * $i) . "r";
		print CF "r1 = $radius + 100p\n";
		print CF "r0 = $radius\n\n";
		if($ptype eq 'histogram'){
			print CF "stroke_type = outline\n";
			print CF "thickness = 0\n";
			print CF "color = black\n";
			print CF "extend_bin = no\n";
			print CF "fill_color = $color\n";
			print CF "<backgrounds>\n";
			print CF "<background>\n";
			print CF "color = $bgcolor\n";
			print CF "</background>\n";
			print CF "</backgrounds>\n";
			print CF "</plot>\n\n";
		}elsif($ptype eq 'line'){
			print CF "orientation = out\n";
			print CF "thickness = 3\n";
			print CF "color = black\n";
			print CF "fill_color = vvlgrey\n";
			print CF "</plot>\n\n";
		}elsif($ptype eq 'heatmap'){
			print CF "color = blues-6-seq,vlgrey,ylgnbu-9-seq\n";
			print CF "</plot>\n\n";
		}else{
			InfoError("Plot type \'$ptype\' is illegal. Please check.");
		}
	}
	
}
print CF "</plots>\n";
##link conf
print CF "<links>\n" if isInArray('link',\@datatype);
for my $i (1..scalar(@plottype)){
	if($datatype[$i - 1] eq 'link'){
		print CF "flat = yes\n";
		print CF "stroke_color = dgrey\n";
		print CF "stoke_thickness = 1\n\n";
		print CF "<link>\n";
		print CF "file = $trackdataNew[$i - 1]\n";
		my $radius = (1 - (0.1 * $i) - 0.03) . "r";
		print CF "radius = $radius\n";
		print CF "bezier_radius = 0.01r\n";
		print CF "color = 55,204,38,0.3\n";
		if($plottype[$i - 1] eq 'linkline'){
			print CF "ribbon = no\n";
		}elsif($plottype[$i - 1] eq 'ribbon'){
			print CF "ribbon = yes\n";
		}else{
			#nothing
		}
		print CF "</link>\n\n";
	}
}
print CF "</links>\n" if isInArray('link',\@datatype);
##other conf
print CF "<<include etc/colors_fonts_patterns.conf>>\n";
print CF "<<include etc/housekeeping.conf>>\n";
print CF "data_out_of_range* = trim\n";

close CF;
##ideogram conf
sleep(1);
Info("Writing ideogram configuration files for Circos.");

print IF "<ideogram>\n\n";
print IF "<spacing>\n";
print IF "default = 0.005r\n";
print IF "</spacing>\n\n";
print IF "#Ideogram position, fill and outline\n";
print IF "radius = 0.90r\n";
print IF "thickness = 80p\n";
print IF "fill = yes\n";
print IF "stroke_color = dgrey\n";
print IF "stroke_thickness = 2p\n\n";
print IF "#Minimum definition for ideogram labels\n";
print IF "show_label = yes\n\n";
print IF "#see etc/fonts.conf for font configuration\n";
print IF "label_font = default\n";
print IF "label_radius = 1r - 60p\n";
print IF "label_size = 50\n";
print IF "label_parallel = yes\n";
print IF "\n</ideogram>\n";

close IF;
##ticks conf
sleep(1);
Info("Writing ticks configuration files for Circos.");

print TF "show_ticks = yes\n";
print TF "show_tick_labels = yes\n\n";
print TF "<ticks>\n";
print TF "radius = 1r\n";
print TF "color = black\n";
print TF "thickness = 2p\n\n";
print TF "multiplier = 1\n";
print TF "format = %d\n\n";
print TF "<tick>\n";
print TF "spacing = 20u\n";
print TF "size = 10p\n";
print TF "</tick>\n\n";
print TF "<tick>\n";
print TF "spacing = 100u\n";
print TF "size = 20p\n";
print TF "show_label = yes\n";
print TF "label_size = 20p\n";
print TF "label_offset = 10p\n";
print TF "format = %d\n";
print TF "</tick>\n\n";
print TF "<tick>\n";
print TF "spacing = 500u\n";
print TF "size = 20p\n";
print TF "show_label = yes\n";
print TF "label_size = 30p\n";
print TF "label_offset = 10p\n";
print TF "format = %d\n";
print TF "</tick>\n\n";
print TF "</ticks>\n";

close TF;

##start to run circos program
my $circos_excu = File::Spec -> catfile($mainBin,'bin','3rdPartyTools','circos','bin','circos');
if(not existFile($circos_excu)){
	InfoError("Circos main program $circos_excu is missing. Please check.");
	exit;
}

Info("Running Circos with $circosConf.");
my $cmd = "perl $circos_excu -conf $circosConf";
runcmd($cmd);

#mv file
my $outpng = File::Spec -> catfile($RealBin,"circos.png");
my $outsvg = File::Spec -> catfile($RealBin,"circos.svg");
my $newpng = File::Spec -> catfile($outputDir,"CircosPlot.png");
my $newsvg = File::Spec -> catfile($outputDir,"CircosPlot.svg");
move($outpng,$newpng);
move($outsvg,$newsvg);


##run success
Info("Program completed!",'green');


####sub program starts here
sub convertToLinkData {
	my $datafile = shift;
	my $chrID = shift;
	
	my %chrID = %$chrID;
	
	my $newdatafile = File::Spec -> catfile($confdir,basename($datafile));
	open ND,">$newdatafile" or die "Can NOT output to data file $newdatafile:$!";
	
	open D,"$datafile" or die "Can NOT open track data file $datafile:$!";
	while(<D>){
		my @tmp = split /\s+/,$_;
		
		my ($chr1,$start1,$end1,$chr2,$start2,$end2) = @tmp;
		
		print ND "$chrID{$chr1} $start1 $end1 $chrID{$chr2} $start2 $end2\n";
	}
	
	close ND;
	close D;
	
	return $newdatafile;
}

sub convertToRegionData {
	my $datafile = shift;
	my $chrID = shift;
	
	my %chrID = %$chrID;
	
	my $newdatafile = File::Spec -> catfile($confdir,basename($datafile));
	open ND,">$newdatafile" or die "Can NOT output to data file $newdatafile:$!";
	
	open D,"$datafile" or die "Can NOT open track data file $datafile:$!";
	while(<D>){
		my @tmp = split /\s+/,$_;
		
		my ($chr1,$start1,$end1,$value) = @tmp;
		
		print ND "$chrID{$chr1} $start1 $end1 $value\n";
	}
	
	close ND;
	close D;
	
	return $newdatafile;
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
           



qap Circos [options]

Use --help to see more information.

qap is still in development. If you have encounted any problem in usage, please feel no hesitation to cotact us.

=head1 DESCRIPTION

This script implements a function to draw Circos plot using formatted data. The script has B<several> mandatory options that MUST appear last. 

=head1 OPTIONS

=over 5

=item --conf,-c F<FILE> [Required]

Path to the configuration file for this program. If you already have one, please specify it by using --conf/-c <file path>. If you don't have one, you can generate a configuration file in a step-by-step manner by using --conf/-c null. See demo/circos.conf for demonstration. 

=item --outputDir,-o F<FILE> [Optional]

Path of the directory to storage result files. If NOT provided, the program will generate a folder automatically.

=item --help,-h

Display this detailed help information.

=back

=head1 EXAMPLE

=over 5

qap Circos -c test.conf -o ./Circos

=back

=head1 AUTHOR

Mingjie Dr.Wang I<huzai@sjtu.edu.cn>

=head1 COPYRIGHT

Copyright (C) 2017, Mingjie Wang. All rights reserved.
