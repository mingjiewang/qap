#!/usr/bin/env perl

package Caller;

use diagnostics;
use strict;
use warnings;
use FindBin qw/$RealBin/;
use File::Spec;
use Getopt::Long;
use Pod::Usage;
use lib "$FindBin::Bin";
use Cwd qw/getcwd abs_path/;
use File::Basename;
use File::Copy;

####Use modules in this program####
use General;
use Mapper;

use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

our $VERSION     = 1.00;
our @ISA         = qw(Exporter);
our @EXPORT      = qw(mutationCallerPipeline prepareReference getMutationFromMSA);
our @EXPORT_OK   = qw(mutationCallerPipeline prepareReference getMutationFromMSA);
our %EXPORT_TAGS = ( DEFAULT => [qw(&mutationCallerPipeline &prepareReference &getMutationFromMSA)]);


##global variables
##get workding directory
my $wk_dir = getcwd;
my $mainBin;
if ($RealBin =~ /(.*)\/bin/){
	$mainBin = $1;
}

##subprogram goes here
sub preprocess {
	my $inputfile = shift;
	my $outputdir = shift;
	my $picard_excu = shift;
	my $samtools_excu = shift;
	my $threads = shift;
	
	my $outputfile = File::Spec -> catfile($outputdir, removeSamBamSuffix2(basename($inputfile),1) . ".rmdup.bam");
	Info("Start to mark duplicates for $inputfile.");
	
	#sort and index 
	#system("cp $inputfile $outputdir");
	copy($inputfile,$outputdir) or die "Can NOT copy $inputfile to $outputdir:$!";
	
	$inputfile = File::Spec -> catfile($outputdir,basename($inputfile));
	my $inputfilesort;
	if(isSamFile($inputfile)){
		sam2SortedAndIndexedBam($samtools_excu, $inputfile, 'pos', 1, $threads);
		$inputfilesort = $inputfile =~ s/\.sam$/.PosSorted.bam/r;
	}elsif(isBamFile($inputfile)){
		sortAndIndexBam($samtools_excu, $inputfile, 'pos', 1, $threads);
		$inputfilesort = $inputfile =~ s/\.bam$/.PosSorted.bam/r;
	}
	
	#fix read group
	my $inputfileWithRGfixed = $inputfilesort =~ s/\.bam$/.withRG.bam/r;
	fixReadGroup($inputfilesort,$inputfileWithRGfixed,$picard_excu);
	system("$samtools_excu index $inputfileWithRGfixed");
	
	#mark duplicates
	my $outmetrics = File::Spec -> catfile(dirname($outputfile),removeSamBamSuffix2(basename($inputfile),1) . ".metric");
	my $cmd = "java -Xmx10g -jar $picard_excu MarkDuplicates I=$inputfileWithRGfixed O=$outputfile CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=$outmetrics";
    runcmd($cmd);
    
    #remove tmp file
    system("rm -rf $inputfile");
    
    #return the path of output file
    return $inputfileWithRGfixed;
}

sub gatk {
	my $inputfile = shift;
	my $outputdir = shift;
	my $ref = shift;
	my $gatk_excu = shift;
	my $threads = shift;
	
	my $inputfileName = basename($inputfile);
	
	#outputdir
	my $gatkdir = File::Spec -> catfile($outputdir,"GATK");
	makedir($gatkdir);
	my $gatkUGdir = File::Spec -> catfile($gatkdir, "UnifiedGenotyper");
	makedir($gatkUGdir);
	my $gatkHCdir = File::Spec -> catfile($gatkdir, "HaplotypeCaller");
	makedir($gatkHCdir);
	
	#run gatk ug
	Info("Start to call variants for $inputfileName using GATK UnifiedGenotyper.");
	my $outputfile_UG = File::Spec -> catfile($gatkUGdir, removeSamBamSuffix2(basename($inputfile),1) . ".vcf");
	my $cmd = "java -Xmx10g -jar $gatk_excu -R $ref -dcov 10000 -nt $threads -T UnifiedGenotyper -glm BOTH --unsafe -I $inputfile -o $outputfile_UG";
	runcmd($cmd);
	
	#run gatk hc
	Info("Start to call variants for $inputfileName using GATK HaplotypeCaller.");
	my $outputfile_HC = File::Spec -> catfile($gatkHCdir, removeSamBamSuffix2(basename($inputfile),1) . ".vcf");
	$cmd = "java -Xmx10g -jar $gatk_excu -R $ref -T HaplotypeCaller --unsafe -I $inputfile -o $outputfile_HC";
	runcmd($cmd);
	
	#run vcf filter
	Info("Start to filter variants for $inputfileName called by GATK UnifiedGenotyper.");
	my $outputfile_UG_fil = File::Spec -> catfile($gatkUGdir, removeSamBamSuffix2(basename($inputfile),1) . ".fil.vcf");
	$cmd = "java -Xmx10g -jar $gatk_excu -R $ref -T VariantFiltration --variant $outputfile_UG --filterExpression \"MQ0 >= 4 \&\& \(\(MQ0 \/ \(1.0 * DP\)\) \> 0.1\)\" --filterName \"HARD_TO_VALIDATE\" --filterExpression \"DP \< 10 \" --filterName \"LowCoverage\" --filterExpression \"QD \< 1.5 \" --filterName \"LowQD\" --filterExpression \"QUAL \> 30.0 \&\& QUAL \< 50.0 \" --filterName \"LowQual\" --clusterWindowSize 10 -o $outputfile_UG_fil ";
	print("VCF Filter CMD: $cmd\n");
	system($cmd);
	
	Info("Start to filter variants for $inputfileName called by GATK HaplotypeCaller.");
	my $outputfile_HC_fil = File::Spec -> catfile($gatkHCdir, removeSamBamSuffix2(basename($inputfile),1) . ".fil.vcf");
	$cmd = "java -Xmx10g -jar $gatk_excu -R $ref -T VariantFiltration --variant $outputfile_HC --filterExpression \"DP \< 10 \" --filterName \"LowCoverage\" --filterExpression \"QD \< 1.5 \" --filterName \"LowQD\" --filterExpression \"QUAL \> 30.0 \&\& QUAL \< 50.0 \" --filterName \"LowQual\" --clusterWindowSize 10 -o $outputfile_HC_fil ";
	print("VCF Filter CMD: $cmd\n");
	system($cmd);

	return($outputfile_UG_fil,$outputfile_HC_fil);
}

sub lofreq {
	my $inputfile = shift;
	my $outputdir = shift;
	my $ref = shift;
	my $lofreq_excu = shift;
	my $threads = shift;
	
	my $inputfileName = basename($inputfile);
	
	Info("Start to call variants for $inputfileName using Lofreq.");

	#outputdir
	my $lofreqdir = File::Spec -> catfile($outputdir, "Lofreq");
	makedir($lofreqdir);
	
	#start to run lofreq
	my $outputvcf = File::Spec -> catfile($lofreqdir, removeSamBamSuffix2($inputfileName,1) . ".vcf");	
	my $cmd = "$lofreq_excu call-parallel --pp-threads $threads -f $ref --call-indels -o $outputvcf $inputfile";
	runcmd($cmd);	
	
	return($outputvcf);
}


sub mpileup {
	my $inputfile = shift;
	my $outputdir = shift;
	my $ref = shift;
	my $samtools_excu = shift;
	
	my $inputfileName = basename($inputfile);
	
	Info("Start to run mpileup for $inputfileName.");
	
	#outputdir
	my $mpileupdir = File::Spec -> catfile($outputdir, 'mpileup');
	makedir($mpileupdir);

	#run mpileup
	my $outputfile = File::Spec -> catfile($mpileupdir, removeSamBamSuffix2($inputfileName,1) . ".mpileup");
	my $cmd = "$samtools_excu mpileup -B -q 1 -f $ref $inputfile > $outputfile ";
	runcmd($cmd);
}

sub varscan {
	my $inputfile = shift;
	my $outputdir = shift;
	my $ref = shift;
	my $samtools_excu = shift;
	my $varscan_excu = shift;
	my $fmtOutputPerlScript = shift;
	
	
	my $inputfileName = basename($inputfile);
	
	Info("Start to call variants for $inputfileName using VarScan.");
	
	#outputdir
	my $varscandir = File::Spec -> catfile($outputdir,'VarScan');
	makedir($varscandir);
	
	#run mpileup
	mpileup($inputfile,$outputdir,$ref,$samtools_excu);
	my $mpileupdir = File::Spec -> catfile($outputdir, 'mpileup');
	my $mpileupfile = File::Spec -> catfile($mpileupdir, removeSamBamSuffix2($inputfileName,1) . ".mpileup");
	
	#run varscan
	my $outputvcf = File::Spec -> catfile($varscandir, removeSamBamSuffix2($inputfileName,1) . ".vcf");
	my $cmd = "java -Xmx10g -jar $varscan_excu mpileup2cns $mpileupfile --output-vcf 1 --variants > $outputvcf";
	runcmd($cmd);
	
	#format varscan output to GATK 
	my $outputvcf2 = $outputvcf =~ s/\.vcf$/.refmt.vcf/r;
	$cmd = "perl $fmtOutputPerlScript -i $outputvcf -o $outputvcf2";
	system($cmd);
	
	return($outputvcf2);
}

sub mutationCallerPipeline {
	my $inputfile = shift;
	my $outputdir = shift;
	my $ref = shift;
	my $program = shift;
	my $samtools_excu = shift;
	my $picard_excu = shift;
	my $gatk_excu = shift;
	my $lofreq_excu = shift;
	my $varscan_excu = shift;
	my $fmtOutputPerlScript = shift;
	my $getMutFreqPerlScript = shift;
	my $mergeMutFreqRScript = shift;
	my $threads = shift;
	
	my $inputfileName = basename($inputfile);
	
	Info("Start calling variants for $inputfileName using multiple callers.");
	
	#preprocessing
	my $rmdupfile = preprocess($inputfile,$outputdir,$picard_excu,$samtools_excu,$threads);
	
	#calling
	my @callingProgram = @$program;
	my @outputvcf;
	my @outputFreq;
	for my $p (@callingProgram){
		if(uc($p) eq 'GATK'){
			my ($gatkugvcf,$gatkhcvcf) = gatk($rmdupfile,$outputdir,$ref,$gatk_excu,$threads);
			
			push @outputvcf,$gatkugvcf;
			push @outputvcf,$gatkhcvcf;
			
			#extract freq
			Info("Extract mutation frequency for $gatkugvcf.");
			my $outputFreqName = basename($gatkugvcf) =~ s/\.vcf$/.gatkUG.mutFreq/r;
			my $freqdir = File::Spec -> catfile($outputdir, 'MutFreq');
			makedir($freqdir);
			my $outputFreq = File::Spec -> catfile($freqdir, $outputFreqName);
			my $cmd = "perl $getMutFreqPerlScript -i $gatkugvcf -o $outputFreq";
			system($cmd);
			push @outputFreq,$outputFreq;
			
			$outputFreqName = basename($gatkhcvcf) =~ s/\.vcf$/.gatkHC.mutFreq/r;
			$freqdir = File::Spec -> catfile($outputdir, 'MutFreq');
			makedir($freqdir);
			$outputFreq = File::Spec -> catfile($freqdir, $outputFreqName);
			$cmd = "perl $getMutFreqPerlScript -i $gatkhcvcf -o $outputFreq";
			system($cmd);
			push @outputFreq,$outputFreq;
			
		}elsif(uc($p) eq 'LOFREQ'){
			my $lofreqvcf = lofreq($rmdupfile,$outputdir,$ref,$lofreq_excu,$threads);
			
			push @outputvcf,$lofreqvcf;
			
			#extract freq
			Info("Extract mutation frequency for $lofreqvcf.");
			my $outputFreqName = basename($lofreqvcf) =~ s/\.vcf$/.lofreq.mutFreq/r;
			my $freqdir = File::Spec -> catfile($outputdir, 'MutFreq');
			makedir($freqdir);
			my $outputFreq = File::Spec -> catfile($freqdir, $outputFreqName);
			my $cmd = "perl $getMutFreqPerlScript -i $lofreqvcf -o $outputFreq";
			system($cmd);
			push @outputFreq,$outputFreq;
			
		}elsif(uc($p) eq 'VARSCAN'){
			my $varscanvcf = varscan($rmdupfile,$outputdir,$ref,$samtools_excu,$varscan_excu,$fmtOutputPerlScript);
			
			push @outputvcf,$varscanvcf;
			
			#extract freq
			Info("Extract mutation frequency for $varscanvcf.");
			my $outputFreqName = basename($varscanvcf) =~ s/\.vcf$/.varscan.mutFreq/r;
			my $freqdir = File::Spec -> catfile($outputdir, 'MutFreq');
			makedir($freqdir);
			my $outputFreq = File::Spec -> catfile($freqdir, $outputFreqName);
			my $cmd = "perl $getMutFreqPerlScript -i $varscanvcf -o $outputFreq";
			system($cmd);
			push @outputFreq,$outputFreq;
			
		}else{
			InfoError("Unkown calling program! Select among \'GATK\' or \'Lofreq\' or \'VarScan\'.");
		}
	}
	
	#generate vcf merege cmd
	my $mergevcfCMD;
	for my $v (@outputvcf){
		$mergevcfCMD .= " --variant $v ";
	}
	
	#merge the vcf
	Info("Start to merge vcf files for $inputfileName.");
	my $mergedvcf = File::Spec -> catfile($outputdir, removeSamBamSuffix2(basename($inputfile),1) . ".vcf");
	my $cmd = "java -Xmx10g -jar $gatk_excu -T CombineVariants -R $ref $mergevcfCMD -o $mergedvcf --unsafe --genotypemergeoption UNIQUIFY";
	runcmd($cmd);
	
	#calculate mutation frequency
	Info("Calculating mutation frequency for $inputfile.");
	my $outfreq = join ",",@outputFreq;
	my $mergefreq = File::Spec -> catfile($outputdir, "MutFreq", removeSamBamSuffix2(basename($inputfile),1) . ".mutFreq");
	my $outputxlsx = File::Spec -> catfile($outputdir, removeSamBamSuffix2(basename($inputfile),1) . ".xlsx");
	$cmd = "Rscript $mergeMutFreqRScript -i $outfreq -o $mergefreq -x $outputxlsx";
	runcmd($cmd);
	
	#move bam files
	my $bamdir = File::Spec -> catfile($outputdir, "Bam");
	makedir($bamdir);
	if (isBamFile($inputfile)){
		$cmd = "mv $outputdir/*.bam $outputdir/*.bai $outputdir/*.metric $bamdir";
	}else{
		$cmd = "mv $outputdir/*.bam $outputdir/*.bai $outputdir/*.metric $bamdir";
	}
	
	system($cmd);
	
}

sub prepareReference {
	my $ref = shift;
	my $picard_excu = shift;
	my $samtools_excu = shift;

	Info("Preparing reference files $ref");
	
	#create dict
	my $dictfile = removeFastaSuffix($ref) . ".dict";
	if(existFile($dictfile)){
		#nothing
	}else{
		my $cmd = "java -Xmx4g -jar $picard_excu CreateSequenceDictionary R=$ref O=$dictfile";
		runcmd($cmd);
	}
	
	#create index file
	my $indexfile = $ref . ".fai";
	if (existFile($indexfile)){
		#nothing
	}else{
		my $cmd = "$samtools_excu faidx $ref";
		runcmd($cmd);
	}
	
}

sub fixReadGroup {
	my $input = shift;
	my $output = shift;
	my $picard_excu = shift;
	
	Info("Fixing \@RG for $input.");
	
	my $cmd = "java -jar $picard_excu AddOrReplaceReadGroups INPUT=$input OUTPUT=$output SORT_ORDER=coordinate RGLB=VirusQuasispecies RGPL=illumina RGPU=ATCGAT RGSM=Sample RGCN=SomeLab";
	runcmd($cmd);
}

sub getMutationFromMSA{
	my $rscript = shift;
	my $refseqfile = shift;
	my $inputfile = shift;
	my $outputdir = shift;
	my $startpos = shift;
	my $cutoff = shift;
	
	#output xlsx file
	my $resdir = File::Spec -> catfile($outputdir, "result");
	makedir($resdir);
	my $outputfile = File::Spec -> catfile($resdir, removeSuffix(basename($inputfile)) . ".xlsx");
	
	#run cmd
	my $cmd = "Rscript $rscript -i $inputfile -o $outputfile -r $refseqfile -s $startpos -c $cutoff";
	runcmd($cmd);
}

1;

