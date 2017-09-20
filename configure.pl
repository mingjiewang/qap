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
use lib "$FindBin::Bin/lib";
use Cwd qw/getcwd abs_path/;
use File::Basename;

$| = 1;

####Use modules in this program####
use General;

####---------------------------####
####The program begins here
####---------------------------####

##check root status
#Info("Root user is required to check and install all the dependencies automatically. Make sure you are running this program as root user.");
#chomp(my $isRoot = `whoami`);
#if($isRoot eq 'root'){
#	#nothing
#}else{
#	InfoWarn("$isRoot, you are not the root user. Please login as root by using \"su root\" or contact the root user.");
#	Info("You can always check and install all the dependencies manually following online instructions (http://bioinfo.rjh.com.cn/labs/jhuang/tools/qap/)");
#	Info("Now exiting...");
#	exit;
#}

##Show welcome
Info("Running software configuration...");
print "\n";


##get workding directory
my $wk_dir = getcwd;
my $mainBin = $RealBin;
if ($RealBin =~ /(.*)\/bin/){
	$mainBin = $1;
}

####define command line arguments
my $help;
my $tmpDir;
my $clear;

my $DateNow = `date +"%Y%m%d_%Hh%Mm%Ss"`;
chomp $DateNow;

GetOptions(
'h|help|'             => \$help,
't|tmpDir|=s'         => \$tmpDir,
'c|clear|=s'          => \$clear
);



##check command line arguments
if (defined $help){
	pod2usage(-verbose=>2,-exitval=>1);
}

if (defined $tmpDir){
	$tmpDir = abs_path($tmpDir) . "/";
	if (not -e $tmpDir){
 		InfoWarn("The directory $tmpDir does NOT exist.",'yellow');
 		InfoWarn("Will mkdir $tmpDir and use it as the temporary log data directory.",'yellow');
		#pod2usage(-verbose=>0,-exitval=>1);
		#exit;
		if (!-e $tmpDir){
			mkdir $tmpDir or InfoDie("Cannot mkdir $tmpDir:$!","red");
		}else{
			InfoError("Mkdir Failed! Folder $tmpDir already exists!","red");
			InfoError("Please specify another temporary log data directory using option -t/--tmpDir");
			pod2usage(-verbose=>2,-exitval=>1);
			exit;
		}
	}
}else{
	$tmpDir = File::Spec -> catfile($wk_dir,"qap_tmpFiles_for_AutoInstall_$DateNow");
	$tmpDir .= "/";
	InfoWarn("The temporary log data directory is not provided!",'yellow');
	InfoWarn("Will mkdir \"$tmpDir\" and use it as the temporary log data directory.",'yellow');
	
	if (!-e "$tmpDir"){
		mkdir $tmpDir or InfoDie("Cannot mkdir $tmpDir:$!","red");
	}else{
		InfoError("Mkdir Failed! $tmpDir already exists!","red");
		InfoError("Please specify another temporary log data directory using option -t/--tmpDir");
		pod2usage(-verbose=>2,-exitval=>1);
		exit;
	}

}


##the core program starts here
my $cmd = '';
my $tmpOutFile = '';
my $autoInstallFile = File::Spec -> catfile($mainBin,"autoInstall");
my @softwareToInstall = ();
my $autoInstallCMD = '';

##start to write to autoinstall file
open AUTO,">$autoInstallFile" or die "Can not write to $autoInstallFile:$!";
print AUTO "#!/usr/bin/sh\n\n";
print AUTO "#####---------------------------------------------------#####\n";
my $DateNow2 = ` date +"%Y-%m-%d %H:%M:%S"`;
chomp $DateNow2;
print AUTO "##### Generated by qap at $DateNow2.\n";
print AUTO "##### Run this script to install denpendencies automatically.\n";
print AUTO "##### If you enconter any problems, please refer to our online instructions \n";
print AUTO "##### http://www.http://bioinfo.rjh.com.cn/labs/jhuang/tools/gap/.\n";
print AUTO "##### Wish you good luck! ^_^ \n";
print AUTO "#####---------------------------------------------------#####\n\n\n";


##check perl version
sleep(1);
Info("Checking Perl version");
#run check 
sleep(1);
$tmpOutFile = File::Spec -> catfile($tmpDir, "01.CheckPerlVersion.log");
open STDERR,">$tmpOutFile" or die "Cannot output to file $tmpOutFile:$!";
$cmd = "perl -v > $tmpOutFile";
system($cmd);
#check output
open T,$tmpOutFile or die "Cannot open $tmpOutFile:$!";
my @out1;
while (<T>){
	chomp;
	push @out1,$_;
}
close T;

my $out1 = join '',@out1;
if ($out1 =~ /This is perl (\d), version (\d+), subversion (\d+)/){
	my $mainversion = $1;
	my $version = $2;
	my $subversion = $3;
	
	if ($mainversion == 5 and $version >= 10){
		Info("Perl version ${mainversion}\.${version}\.${subversion}..........PASS",'green');
	}else{
		InfoError("Perl version MUST be > 5.10, current version ${mainversion}\.${version}\.${subversion}. Please update your Perl.");
		quit();
	}
}else{
	InfoError("Perl is not installed in your system or not in system PATH. Please check.");
	quit();
}



##check Perl threads
sleep(1);
Info("Checking threading status");

sleep(1);
$tmpOutFile = File::Spec -> catfile($tmpDir, "02.CheckPerlThreads.log");
open STDERR,">$tmpOutFile" or die "Cannot output to file $tmpOutFile:$!";
my $threads_usable = eval 'use threads; 1';
if ($threads_usable) {
	Info("Perl threading enabled..........PASS",'green');
	print STDERR "Perl threading enabled..........PASS";
} else {
	Info("No threading is possible. Please install perl module: threads or recompile perl with option -Dusethreads","red");
	print STDERR "No threading is possible. Please install perl module: threads or recompile perl with option -Dusethreads";
	quit();
}


##check Python version 
sleep(1);
Info("Checking Python version");
#run check 
sleep(1);
$tmpOutFile = File::Spec -> catfile($tmpDir, "03.CheckPythonVersion.log");
open STDERR,">$tmpOutFile" or die "Cannot output to file $tmpOutFile:$!";
$cmd = "python -V > $tmpOutFile";
system($cmd);
#check output
open T,$tmpOutFile or die "Cannot open $tmpOutFile:$!";
my @out2;
while (<T>){
	chomp;
	push @out2,$_;
}
close T;

my $out2 = join '',@out2;
if ($out2 =~ /Python (\d+)\.(\d+)\.(\d+)/){
	my $mainversion = $1;
	my $version = $2;
	my $subversion = $3;
	
	if($mainversion == 2){
		Info("Python version ${mainversion}\.${version}\.${subversion}..........PASS",'green');
	}else{
		InfoError("Python 2 is required. Please re-install python with version 2 (>=2.7.9)");
		quit();
	}
}else{
	InfoError("Python is not installed in your system or not in system PATH. Please check.");
	quit();
}


##check R version 
sleep(1);
Info("Checking R version");
#run check
sleep(1);
$tmpOutFile = File::Spec -> catfile($tmpDir, "04.CheckRVersion.log");
open STDERR,">$tmpOutFile" or die "Cannot output to file $tmpOutFile:$!";
$cmd = "R --version > $tmpOutFile";
system($cmd);
#check output
open T,$tmpOutFile or die "Cannot open $tmpOutFile:$!";
my @out3;
while (<T>){
	chomp;
	push @out3,$_;
}
close T;

my $out3 = join '',@out3;

if ($out3 =~ /R version (\d+)\.(\d+)\.(\d+)/){
	my $mainversion = $1;
	my $version = $2;
	my $subversion = $3;
	
	if ($mainversion == 3 and $version >= 3){
		Info("R version ${mainversion}\.${version}\.${subversion}.........PASS",'green');
	}else{
		InfoError("R >3.3.x is required. Please (re)install R software.");
		quit();
	}
}else{
	InfoError("R is not installed in your system or not in system PATH. Please check.");
	quit();
}


##check java version 
sleep(1);
Info("Checking Java version");
#run check
sleep(1);
$tmpOutFile = File::Spec -> catfile($tmpDir, "05.CheckJavaVersion.log");
open STDERR,">$tmpOutFile" or die "Cannot output to file $tmpOutFile:$!";
$cmd = "java -version > $tmpOutFile";
system($cmd);
#check output
open T,$tmpOutFile or die "Cannot open $tmpOutFile:$!";
my @out4;
while (<T>){
	chomp;
	push @out4,$_;
}
close T;

my $out4 = join '',@out4;

if ($out4 =~ /version \"(\d+)\.(\d+)\.(\d+)_\d+\"/){
	my $mainversion = $1;
	my $version = $2;
	my $subversion = $3;
	
	if ($version >= 6){
		Info("Java version ${mainversion}\.${version}\.${subversion}.........PASS",'green');
	}else{
		InfoError("Java Runtime Environment (JRE) >1.8.x is required. Please (re)install Java.");
		quit();
	}
}else{
	InfoError("Java is not installed in your system or not in system PATH. Please check.");
	quit();
}


##get software path
my $thirdPartyToolsPath = File::Spec -> catfile($mainBin, "bin", "3rdPartyTools");


##check python pip
my $pipPass = 0;
sleep(1);
Info("Checking Python pip");
#run check
sleep(1);
$tmpOutFile = File::Spec -> catfile($tmpDir, "06.CheckPip.log");
open STDERR,">$tmpOutFile" or die "Cannot output to file $tmpOutFile:$!";
$cmd = "pip -V > $tmpOutFile";
system($cmd);
#check output
open T,$tmpOutFile or die "Cannot open $tmpOutFile:$!";
my @out5;
while (<T>){
	chomp;
	push @out5,$_;
}
close T;

my $out5 = join '',@out5;

if ($out5 =~ /pip (\d+)\.(\d+)\.(\d+) from \//){
	my $mainversion = $1;
	my $version = $2;
	my $subversion = $3;
	
	if ($mainversion >= 7){
		Info("pip version ${mainversion}\.${version}\.${subversion}.........PASS",'green');
		$pipPass = 1;
	}else{
		InfoWarn("A newer version of pip is required. Installation CMD written to autoInstall file.");
		print STDERR "A newer version of pip is required. Will updating pip automatically.";
		
		$autoInstallCMD .= "##Upgrating python pip\n";
		$autoInstallCMD .= "pip install -U pip\n\n";
		
		push @softwareToInstall,"pip";
	}
}else{
	InfoError("pip is not installed in your system or not in system PATH. Please check.");
	print STDERR "pip is not installed in your system or not in system PATH.\n";
	
	$autoInstallCMD .= "##Installing python pip\n";
	##there are two ways to install pip, try to use both ways to make sure pip is installed.
	#method1--automatically install
	my $getpipScript = File::Spec -> catfile($mainBin, "bin", "PythonScripts", "get-pip.py");
	$autoInstallCMD .= "python $getpipScript\n\n";
	#method2--build from source
	my $pipnew = File::Spec -> catfile($thirdPartyToolsPath, "pip.tar.gz");
	my $pipold = File::Spec -> catfile($thirdPartyToolsPath, "pip");
	$autoInstallCMD .= "wget http://bioinfo.rjh.com.cn/labs/jhuang/tools/qap/download/pip.tar.gz -O $pipnew\n";
	$autoInstallCMD .= "cd $thirdPartyToolsPath\n";
	$autoInstallCMD .= "tar xvzf $pipnew\n";
	$autoInstallCMD .= "cd $pipold\n";
	$autoInstallCMD .= "python setup.py install\n";
	$autoInstallCMD .= "cd $mainBin\n\n";
	
	push @softwareToInstall,"pip";
	
}


##check bowtie2
sleep(1);
Info("Checking bowtie2");
#run check
sleep(1);
$tmpOutFile = File::Spec -> catfile($tmpDir, "07.CheckBowtie2.log");
open STDERR,">>$tmpOutFile" or die "Cannot output to file $tmpOutFile:$!";
my $bowtie2_excu = File::Spec -> catfile($mainBin, "bin", "3rdPartyTools", "bowtie2", "bowtie2");
my $bowtie2_build_excu = File::Spec -> catfile($mainBin, "bin", "3rdPartyTools", "bowtie2", "bowtie2-build");
$cmd = "$bowtie2_excu -h > $tmpOutFile";
system($cmd);
$cmd = "$bowtie2_build_excu -h >> $tmpOutFile";
system($cmd);

#check output
open T,$tmpOutFile or die "Cannot open $tmpOutFile:$!";
my @out6;
while (<T>){
	chomp;
	push @out6,$_;
}
close T;

my $out6 = join '',@out6;

if ($out6 =~ /Bowtie 2 version/){
	Info("Bowtie2.........PASS",'green');
}else{
	InfoError("Bowtie2 is not installed in your system or not in bin/ folder. Please check.");
	$autoInstallCMD .= "##Installing Bowtie2\n";
	my $bowtie2new = File::Spec -> catfile($mainBin, 'bin', '3rdPartyTools', 'bowtie2.tar.gz');
	$autoInstallCMD .= "wget http://bioinfo.rjh.com.cn/labs/jhuang/tools/qap/download/bowtie2.tar.gz -O $bowtie2new\n";
	my $bowtie2old = dirname($bowtie2_excu);
	$autoInstallCMD .= "rm -rf $bowtie2old\n";
	$autoInstallCMD .= "cd $thirdPartyToolsPath\n";
	$autoInstallCMD .= "tar xvzf $bowtie2new\n";
	$autoInstallCMD .= "cd $mainBin\n\n";
	
	push @softwareToInstall,'Bowtie2';
}


##check bwa
sleep(1);
Info("Checking BWA");
#run check
sleep(1);
$tmpOutFile = File::Spec -> catfile($tmpDir, "08.CheckBWA.log");
open STDERR,">$tmpOutFile" or die "Cannot output to file $tmpOutFile:$!";
my $bwa_excu = File::Spec -> catfile($mainBin, "bin", "3rdPartyTools", "bwa", "bwa");
$cmd = "$bwa_excu > $tmpOutFile";
system($cmd);

#check output
open T,$tmpOutFile or die "Cannot open $tmpOutFile:$!";
my @out7;
while (<T>){
	chomp;
	push @out7,$_;
}
close T;

my $out7 = join '',@out7;

if ($out7 =~ /Program: bwa \(alignment via Burrows-Wheeler transformation\)/){
	Info("BWA.........PASS",'green');
}else{
	InfoError("BWA is not installed in your system or not in bin/ folder. Please check.");
	$autoInstallCMD .= "##Installing BWA\n";
	my $bwanew = File::Spec -> catfile($mainBin, 'bin', '3rdPartyTools', 'bwa.tar.gz');
	$autoInstallCMD .= "wget http://bioinfo.rjh.com.cn/labs/jhuang/tools/qap/download/bwa.tar.gz -O $bwanew\n";
	my $bwaold = dirname($bwa_excu);
	$autoInstallCMD .= "rm -rf $bwaold\n";
	$autoInstallCMD .= "cd $thirdPartyToolsPath\n";
	$autoInstallCMD .= "tar xvzf $bwanew\n";
	$autoInstallCMD .= "cd $bwaold\n";
	$autoInstallCMD .= "make\n";
	$autoInstallCMD .= "cd $mainBin\n\n";
	
	push @softwareToInstall,'BWA';
}

##check samtools
sleep(1);
Info("Checking Samtools");
#run check
sleep(1);
$tmpOutFile = File::Spec -> catfile($tmpDir, "09.CheckSamtools.log");
open STDERR,">$tmpOutFile" or die "Cannot output to file $tmpOutFile:$!";
my $samtools_excu = File::Spec -> catfile($mainBin, "bin", "3rdPartyTools", "samtools", "samtools");
$cmd = "$samtools_excu > $tmpOutFile";
system($cmd);

#check output
open T,$tmpOutFile or die "Cannot open $tmpOutFile:$!";
my @out8;
while (<T>){
	chomp;
	push @out8,$_;
}
close T;

my $out8 = join '',@out8;

if ($out8 =~ /Program: samtools \(Tools for alignments in the SAM format\)/){
	Info("Samtools.........PASS",'green');
}else{
	InfoError("Samtools is not installed in your system or not in bin/ folder. Please check.");
	$autoInstallCMD .= "##Installing Samtools\n";
	my $samtoolsnew = File::Spec -> catfile($mainBin, 'bin', '3rdPartyTools', 'samtools.tar.gz');
	$autoInstallCMD .= "wget http://bioinfo.rjh.com.cn/labs/jhuang/tools/qap/download/samtools.tar.gz -O $samtoolsnew\n";
	my $samtoolsold = dirname($samtools_excu);
	$autoInstallCMD .= "rm -rf $samtoolsold\n";
	$autoInstallCMD .= "cd $thirdPartyToolsPath\n";
	$autoInstallCMD .= "tar xvzf $samtoolsnew\n";
	$autoInstallCMD .= "cd $samtoolsold\n";
	$autoInstallCMD .= "./configure\n";
	$autoInstallCMD .= "make\n";
	$autoInstallCMD .= "cd $mainBin\n\n";
	
	push @softwareToInstall,'Samtools';
}

##check Fastqc
sleep(1);
Info("Checking FastQC");
#run check
sleep(1);
$tmpOutFile = File::Spec -> catfile($tmpDir, "10.CheckFastQC.log");
open STDERR,">$tmpOutFile" or die "Cannot output to file $tmpOutFile:$!";
my $fastqc_excu = File::Spec -> catfile($mainBin, "bin", "3rdPartyTools", "fastqc", "fastqc");
$cmd = "$fastqc_excu -v > $tmpOutFile";
system($cmd);

#check output
open T,$tmpOutFile or die "Cannot open $tmpOutFile:$!";
my @out9;
while (<T>){
	chomp;
	push @out9,$_;
}
close T;

my $out9 = join '',@out9;

if ($out9 =~ /FastQC v[0-9\.]+/){
	Info("FastQC.........PASS",'green');
}else{
	InfoError("FastQC is not installed in your system or not in bin/ folder. Please check.");
	$autoInstallCMD .= "##Installing FastQC\n";
	my $fastqcnew = File::Spec -> catfile($mainBin, 'bin', '3rdPartyTools', 'fastqc.tar.gz');
	$autoInstallCMD .= "wget http://bioinfo.rjh.com.cn/labs/jhuang/tools/qap/download/fastqc.tar.gz -O $fastqcnew\n";
	my $fastqcold = dirname($fastqc_excu);
	$autoInstallCMD .= "rm -rf $fastqcold\n";
	$autoInstallCMD .= "cd $thirdPartyToolsPath\n";
	$autoInstallCMD .= "tar xvzf $fastqcnew\n";
	$autoInstallCMD .= "cd $fastqcold\n";
	$autoInstallCMD .= "chmod 755 fastqc\n";
	$autoInstallCMD .= "cd $mainBin\n\n";
	
	push @softwareToInstall,'FastQC';
}

##check cutadapt
sleep(1);
Info("Checking cutadapt");
#run check
sleep(1);
$tmpOutFile = File::Spec -> catfile($tmpDir, "11.CheckCutadapt.log");
open STDERR,">$tmpOutFile" or die "Cannot output to file $tmpOutFile:$!";
$cmd = "cutadapt -h > $tmpOutFile";
system($cmd);

#check output
open T,$tmpOutFile or die "Cannot open $tmpOutFile:$!";
my @out10;
while (<T>){
	chomp;
	push @out10,$_;
}
close T;

my $out10 = join '',@out10;

if ($out10 =~ /cutadapt version [0-9\.]+/){
	Info("cutadapt.........PASS",'green');
}else{
	InfoError("cutadapt is not installed in your system or not in system PATH. Please check.");
	if($pipPass){
		$autoInstallCMD .= "##Installing cutadapt\n";
		$autoInstallCMD .= "pip install --user --upgrade cutadapt\n\n";
	}else{
		InfoError("Automatic installation for cutadapt requires pip installation, however pip failed the check. Please fix pip then try again.")
	}
	
	push @softwareToInstall,'cutadapt';
}


##check FASTX-Toolkit
sleep(1);
Info("Checking FASTX-Toolkit");
#run check
sleep(1);
$tmpOutFile = File::Spec -> catfile($tmpDir, "12.CheckFastxToolkit.log");
open STDERR,">>$tmpOutFile" or die "Cannot output to file $tmpOutFile:$!";
my $fastx_trimmer_excu = File::Spec -> catfile($thirdPartyToolsPath, "fastx_toolkit", "fastx_trimmer");
my $fastx_splitter_excu = File::Spec -> catfile($thirdPartyToolsPath, "fastx_toolkit", "fastx_barcode_splitter.pl");
$cmd = "$fastx_trimmer_excu -h > $tmpOutFile";
system($cmd);
$cmd = "perl $fastx_splitter_excu -h >> $tmpOutFile";
system($cmd);

#check output
open T,$tmpOutFile or die "Cannot open $tmpOutFile:$!";
my @out11;
while (<T>){
	chomp;
	push @out11,$_;
}
close T;

my $out11 = join '',@out11;

if ($out11 =~ /usage\: fastx_trimmer.*Barcode Splitter, by Assaf Gordon/){
	Info("FASTX-Toolkit.........PASS",'green');
}else{
	InfoError("FASTX-Toolkit is not installed in your system or not in bin/ folder. Please check.");
	$autoInstallCMD .= "##Installing FASTX-Toolkit\n";
	my $fastxnew = File::Spec -> catfile($mainBin, 'bin', '3rdPartyTools', 'fastx_toolkit.tar.gz');
	$autoInstallCMD .= "wget http://bioinfo.rjh.com.cn/labs/jhuang/tools/qap/download/fastx_toolkit.tar.gz -O $fastxnew\n";
	my $fastxold = dirname($fastx_trimmer_excu);
	$autoInstallCMD .= "rm -rf $fastxold\n";
	$autoInstallCMD .= "cd $thirdPartyToolsPath\n";
	$autoInstallCMD .= "tar xvzf $fastxnew\n";
	$autoInstallCMD .= "cd $mainBin\n\n";
	
	push @softwareToInstall,'FASTX-Toolkit';
}


##check megacc
sleep(1);
Info("Checking MEGA");
#run check
sleep(1);
$tmpOutFile = File::Spec -> catfile($tmpDir, "13.CheckMEGA.log");
open STDERR,">>$tmpOutFile" or die "Cannot output to file $tmpOutFile:$!";
my $mega_excu = File::Spec -> catfile($thirdPartyToolsPath, "megacc", "megacc");
$cmd = "$mega_excu -h > $tmpOutFile";
system($cmd);

#check output
open T,$tmpOutFile or die "Cannot open $tmpOutFile:$!";
my @out12;
while (<T>){
	chomp;
	push @out12,$_;
}
close T;

my $out12 = join '',@out12;

if ($out12 =~ /-a --analysisOptions/){
	Info("MEGA.........PASS",'green');
}else{
	InfoError("MEGA is not installed in your system or not in bin/ folder. Please check.");
	$autoInstallCMD .= "##Installing MEGA\n";
	my $meganew = File::Spec -> catfile($mainBin, 'bin', '3rdPartyTools', 'megacc.tar.gz');
	$autoInstallCMD .= "wget http://bioinfo.rjh.com.cn/labs/jhuang/tools/qap/download/megacc.tar.gz -O $meganew\n";
	my $megaold = dirname($mega_excu);
	$autoInstallCMD .= "rm -rf $megaold\n";
	$autoInstallCMD .= "cd $thirdPartyToolsPath\n";
	$autoInstallCMD .= "tar xvzf $meganew\n";
	$autoInstallCMD .= "cd $mainBin\n\n";
	
	push @softwareToInstall,'MEGA';
}


##check clustalo
sleep(1);
Info("Checking Clustal Omega");
#run check
sleep(1);
$tmpOutFile = File::Spec -> catfile($tmpDir, "14.CheckClustalOmega.log");
open STDERR,">>$tmpOutFile" or die "Cannot output to file $tmpOutFile:$!";
my $clustalo_excu = File::Spec -> catfile($thirdPartyToolsPath, 'msa',"clustalo", "clustalo");
$cmd = "$clustalo_excu -h > $tmpOutFile";
system($cmd);

#check output
open T,$tmpOutFile or die "Cannot open $tmpOutFile:$!";
my @out13;
while (<T>){
	chomp;
	push @out13,$_;
}
close T;

my $out13 = join '',@out13;

if ($out13 =~ /-i, --in, --infile=/){
	Info("Clustal Omega.........PASS",'green');
}else{
	InfoError("Clustal Omega is not installed in your system or not in bin/ folder. Please check.");
	$autoInstallCMD .= "##Installing Clustal Omega\n";
	my $clustalonew = File::Spec -> catfile($mainBin, 'bin', '3rdPartyTools', 'msa','clustalo.tar.gz');
	$autoInstallCMD .= "wget http://bioinfo.rjh.com.cn/labs/jhuang/tools/qap/download/clustalo.tar.gz -O $clustalonew\n";
	my $clustaloold = dirname($clustalo_excu);
	$autoInstallCMD .= "rm -rf $clustaloold\n";
	$autoInstallCMD .= "cd $thirdPartyToolsPath/msa\n";
	$autoInstallCMD .= "tar xvzf $clustalonew\n";
	$autoInstallCMD .= "cd $mainBin\n\n";
	
	push @softwareToInstall,'Clustal Omega';
}


##check clustalx
sleep(1);
Info("Checking ClustalW");
#run check
sleep(1);
$tmpOutFile = File::Spec -> catfile($tmpDir, "15.CheckClustalW.log");
open STDERR,">>$tmpOutFile" or die "Cannot output to file $tmpOutFile:$!";
my $clustalw_excu = File::Spec -> catfile($thirdPartyToolsPath, 'msa',"clustalw", "clustalw");
$cmd = "$clustalw_excu -help > $tmpOutFile";
system($cmd);

#check output
open T,$tmpOutFile or die "Cannot open $tmpOutFile:$!";
my @out14;
while (<T>){
	chomp;
	push @out14,$_;
}
close T;

my $out14 = join '',@out14;

if ($out14 =~ /-INFILE=file.ext/){
	Info("ClustalW.........PASS",'green');
}else{
	InfoError("ClustalW is not installed in your system or not in bin/ folder. Please check.");
	$autoInstallCMD .= "##Installing ClustalX\n";
	my $clustalwnew = File::Spec -> catfile($mainBin, 'bin', '3rdPartyTools', 'msa','clustalw.tar.gz');
	$autoInstallCMD .= "wget http://bioinfo.rjh.com.cn/labs/jhuang/tools/qap/download/clustalw.tar.gz -O $clustalwnew\n";
	my $clustalwold = dirname($clustalw_excu);
	$autoInstallCMD .= "rm -rf $clustalwold\n";
	$autoInstallCMD .= "cd $thirdPartyToolsPath/msa\n";
	$autoInstallCMD .= "tar xvzf $clustalwnew\n";
	$autoInstallCMD .= "cd $mainBin\n\n";
	
	push @softwareToInstall,'ClustalW';
}


##check muscle
sleep(1);
Info("Checking MUSCLE");
#run check
sleep(1);
$tmpOutFile = File::Spec -> catfile($tmpDir, "16.CheckMUSCLE.log");
open STDERR,">>$tmpOutFile" or die "Cannot output to file $tmpOutFile:$!";
my $muscle_excu = File::Spec -> catfile($thirdPartyToolsPath, 'msa',"muscle", "muscle");
$cmd = "$muscle_excu > $tmpOutFile";
system($cmd);

#check output
open T,$tmpOutFile or die "Cannot open $tmpOutFile:$!";
my @out15;
while (<T>){
	chomp;
	push @out15,$_;
}
close T;

my $out15 = join '',@out15;

if ($out15 =~ /-in \<inputfile\>/){
	Info("MUSCLE.........PASS",'green');
}else{
	InfoError("MUSCLE is not installed in your system or not in bin/ folder. Please check.");
	$autoInstallCMD .= "##Installing MUSCLE\n";
	my $musclenew = File::Spec -> catfile($mainBin, 'bin', '3rdPartyTools', 'msa','muscle.tar.gz');
	$autoInstallCMD .= "wget http://bioinfo.rjh.com.cn/labs/jhuang/tools/qap/download/muscle.tar.gz -O $musclenew\n";
	my $muscleold = dirname($muscle_excu);
	$autoInstallCMD .= "rm -rf $muscleold\n";
	$autoInstallCMD .= "cd $thirdPartyToolsPath/msa\n";
	$autoInstallCMD .= "tar xvzf $musclenew\n";
	$autoInstallCMD .= "cd $mainBin\n\n";
	
	push @softwareToInstall,'MUSCLE';
}

##check gs
# gs must be in the system path

##check weblogo

##check perl modules
#Bio::Perl  cmd: cpan Bio::Perl ; perl -MCPAN -Mlocal::lib -e "CPAN::install(\"C/CJ/CJFIELDS/BioPerl-1.007001.tar.gz\")"
#Perl4::CoreLibs  to successful run perl scripts compile on perl4 

##check R packages
#bioconductor:Biostrings ; 
#base:gplots RColorBrewer ggplot2 seqinr optparse plyr xlsx stringr scatterplot3d devtools ggpubr ape
#github:ggbiplot
###library(devtools)
###install_github("vqv/ggbiplot")


##check Biopython
#pip install Biopython

##check gsl
#sudo apt-get upate -qq
#sudo apt-get install -y gsl-bin libgs10-dev
#sudo yum install -y gsl gsl-devel

##check shorah-master and other qsr programs

##check igv

##check circos modules 

##check makeblastdb,blastn

##get cmd file ready to run
print AUTO $autoInstallCMD;
$cmd = "chmod 755 $autoInstallFile";
system($cmd);


##sub program goes here
sub quit {
	InfoError("Aborting automatic check util this problem is solved. Please refer to $tmpDir for detailed error information");
	exit(0);
}



##run success
sleep(1);
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
           



gap MapReadsToRef [options]

Use --help to see more information.

gap is still in development. If you have encounted any problem in usage, please feel no hesitation to cotact us.

=head1 DESCRIPTION

This script implements a function to install all the dependencies (third-party softwares, modules, packages etc.) automatically.  

=head1 OPTIONS

=over 5

=item --temDir,-t F<FILE> [Optional]

Path to the directory to store the temporary log files.

=item --clear,-c S<STRING> [Optional]

Whether removal all the temporaty log files or not. Choose between 'T' (clear all data files) and 'F' (keep all data files). Default value is 'F'.

=item --help,-h

Display this detailed help information.

=back

=head1 EXAMPLE

=over 5

gap AutoInstall -t ~/tmp -c F

=back

=head1 AUTHOR

Mingjie Dr.Wang I<huzai@sjtu.edu.cn>

=head1 COPYRIGHT

Copyright (C) 2017, Mingjie Wang. All rights reserved.

