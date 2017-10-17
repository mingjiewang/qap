#!perl

package NT;

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



####Use modules in this program####
use General;

use vars qw ($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);
our $VERSION         = 1.00;
our @ISA             = qw (Exporter);
our @EXPORT          = qw ($DEBUG_MODE translation getConsensusSeq getConsensusSeqWithInfo) ;
our @EXPORT_OK       = qw ($DEBUG_MODE translation getConsensusSeq getConsensusSeqWithInfo);
our %EXPORT_TAGS     = ( DEFAULT => [qw (&translation &getConsensusSeq &getConsensusSeqWithInfo) ] );

##Debug mode prints more information while running 
our $DEBUG_MODE = 1;

##subprograms begin here##
sub translation {
	my $dna = shift;  
	my $seqID = shift;
	
	my $protein = ' ';  
	my $codon;  
	 
	my $len = length $dna;
	
	for(my $i=0; $i<( $len -2);$i+=3)  
	{  
	    $codon=substr($dna,$i,3);  
	    $protein .= codon2aa($codon, $seqID);  
	}  
	
	$protein =~ s/\s//g;
	
	return $protein;   
} 

sub codon2aa {     
    my $codon = shift;
    my $seqID = shift;     
     
    $codon = uc $codon; #uc=uppercase;lc=lowercase    
    
    my(%genetic_code) = (     
    'TCA' => 'S',    # Serine     
    'TCC' => 'S',    # Serine     
    'TCG' => 'S',    # Serine     
    'TCT' => 'S',    # Serine     
    'TTC' => 'F',    # Phenylalanine     
    'TTT' => 'F',    # Phenylalanine     
    'TTA' => 'L',    # Leucine     
    'TTG' => 'L',    # Leucine     
    'TAC' => 'Y',    # Tyrosine      
    'TAT' => 'Y',    # Tyrosine     
    'TAA' => '*',    # Stop     
    'TAG' => '*',    # Stop     
    'TGC' => 'C',    # Cysteine     
    'TGT' => 'C',    # Cysteine     
    'TGA' => '*',    # Stop     
    'TGG' => 'W',    # Tryptophan     
    'CTA' => 'L',    # Leucine     
    'CTC' => 'L',    # Leucine     
    'CTG' => 'L',    # Leucine     
    'CTT' => 'L',    # Leucine     
    'CCA' => 'P',    # Proline     
    'CCC' => 'P',    # Proline     
    'CCG' => 'P',    # Proline     
    'CCT' => 'P',    # Proline     
    'CAC' => 'H',    # Histidine     
    'CAT' => 'H',    # Histidine     
    'CAA' => 'Q',    # Glutamine     
    'CAG' => 'Q',    # Glutamine     
    'CGA' => 'R',    # Arginine     
    'CGC' => 'R',    # Arginine     
    'CGG' => 'R',    # Arginine     
    'CGT' => 'R',    # Arginine     
    'ATA' => 'I',    # Isoleucine     
    'ATC' => 'I',    # Isoleucine     
    'ATT' => 'I',    # Isoleucine     
    'ATG' => 'M',    # Methionine     
    'ACA' => 'T',    # Threonine     
    'ACC' => 'T',    # Threonine     
    'ACG' => 'T',    # Threonine     
    'ACT' => 'T',    # Threonine     
    'AAC' => 'N',    # Asparagine     
    'AAT' => 'N',    # Asparagine     
    'AAA' => 'K',    # Lysine     
    'AAG' => 'K',    # Lysine     
    'AGC' => 'S',    # Serine     
    'AGT' => 'S',    # Serine     
    'AGA' => 'R',    # Arginine     
    'AGG' => 'R',    # Arginine     
    'GTA' => 'V',    # Valine     
    'GTC' => 'V',    # Valine     
    'GTG' => 'V',    # Valine     
    'GTT' => 'V',    # Valine     
    'GCA' => 'A',    # Alanine     
    'GCC' => 'A',    # Alanine     
    'GCG' => 'A',    # Alanine     
    'GCT' => 'A',    # Alanine         
    'GAC' => 'D',    # Aspartic Acid     
    'GAT' => 'D',    # Aspartic Acid     
    'GAA' => 'E',    # Glutamic Acid     
    'GAG' => 'E',    # Glutamic Acid     
    'GGA' => 'G',    # Glycine     
    'GGC' => 'G',    # Glycine     
    'GGG' => 'G',    # Glycine     
    'GGT' => 'G',    # Glycine     
    );     
     
    if(exists $genetic_code{$codon}){     
    	return $genetic_code{$codon};     
    }elsif($codon =~ /-/){
    	return '-'; #If there is deletion in the nucleotide sequence, translate the codon to '-' instead of translating the frameshift mutated sequences.
    }else{     
		InfoError("Sequence\[$seqID\]:Bad codon \"$codon\".");     
     	exit;     
    }     
} 

sub getConsensusSeq {
	my $inputfile = shift;
	my $rscript = shift;
	my $outputfile = shift;
	my $graphic = shift;
	my $gsProgram = shift;
	my $seqlogoProgram = shift;
	my $deletion = shift;
	
	#convert fasta file to two line format
	my $inputfile2Name = removeFastaSuffix(basename($inputfile)) . ".2line.fasta";
	my $inputfile2 = File::Spec -> catfile(dirname($outputfile),$inputfile2Name);
	formatFastaToTwoLineMode($inputfile,$inputfile2);
	
	#get fasta ready for r input
	#r input file
	my $rinputName = removeFastaSuffix(basename($inputfile)) . ".RInput";
	my $rinput = File::Spec -> catfile(dirname($outputfile), $rinputName);
	open R,">$rinput" or die "Can NOT output to $rinput:$!\n";
	
	open T,"$inputfile2" or die "Can NOT open $inputfile2:$!\n";
	while(my $line1 = <T>){
		chomp $line1;
		chomp (my $line2 = <T>);
		
		print R "$line2\n";
	}
	close T;
	close R;
	
	#get cs
	my $fileLabel = removeFastaSuffix(basename($inputfile));
	my $cmd = "Rscript $rscript -i $rinput -o $outputfile -l $fileLabel -d $deletion";
	system($cmd);
	
	#graphic weblogo
	if($graphic){
		# check folder 
		my $pdffolder = File::Spec -> catfile(dirname($outputfile), 'pdf');
		my $pngfolder = File::Spec -> catfile(dirname($outputfile), 'png');
		if(not -e $pdffolder){
			mkdir $pdffolder or die "Can NOT mkdir $pdffolder:$!";
		}
		if(not -e $pngfolder){
			mkdir $pngfolder or die "Can NOT mkdir $pngfolder:$!";
		}
		
		# pdf
		my $outgraphic1 = File::Spec -> catfile($pdffolder,removeFastaSuffix(basename($inputfile)) . ".pdf");
		my $cmd1 = "$seqlogoProgram -F pdf -f $inputfile2 -C 100 -h 4 -w 30 -c -e -n -p -Y -x Sequence -y Consensus > $outgraphic1";
		system($cmd1);
		
		# png
		my $outgraphic2 = File::Spec -> catfile($pngfolder,removeFastaSuffix(basename($inputfile)) . ".png");
		my $cmd2 = "$seqlogoProgram -F png -f $inputfile2 -C 100 -h 4 -w 30 -c -e -n -p -Y -x Sequence -y Consensus > $outgraphic2";
		system($cmd2);
		
	}
	
	#remove tmp file
	system("rm -rf $inputfile2 $rinput");
	
} 

sub getConsensusSeqWithInfo {
	my $inputfile = shift;
	my $rscript = shift;
	my $outputfile = shift;
	my $graphic = shift;
	my $gsProgram = shift;
	my $seqlogoProgram = shift;
	my $deletion = shift;
	
	Info("Calculating CS for $inputfile");
	
	#convert fasta file to two line format
	my $inputfile2Name = removeFastaSuffix(basename($inputfile)) . ".2line.fasta";
	my $inputfile2 = dirname($outputfile) . "/" . $inputfile2Name;
	formatFastaToTwoLineMode($inputfile,$inputfile2);
	
	#get fasta ready for r input
	#r input file
	my $rinputName = removeFastaSuffix(basename($inputfile)) . ".RInput";
	my $rinput = File::Spec -> catfile(dirname($outputfile), $rinputName);
	open R,">$rinput" or die "Can NOT output to $rinput:$!\n";
	
	open T,"$inputfile2" or die "Can NOT open $inputfile2:$!\n";
	while(my $line1 = <T>){
		chomp $line1;
		chomp (my $line2 = <T>);
		
		print R "$line2\n";
	}
	close T;
	close R;
	
	#get cs
	my $fileLabel = removeFastaSuffix(basename($inputfile));
	my $cmd = "Rscript $rscript -i $rinput -o $outputfile -l $fileLabel -d $deletion";
	runcmd($cmd);
	
	#graphic weblogo
	if($graphic){
		Info("Drawing graphics for $inputfile");
		
		# check folder 
		my $pdffolder = File::Spec -> catfile(dirname($outputfile), 'pdf');
		my $pngfolder = File::Spec -> catfile(dirname($outputfile), 'png');
		if(not -e $pdffolder){
			mkdir $pdffolder or die "Can NOT mkdir $pdffolder:$!";
		}
		if(not -e $pngfolder){
			mkdir $pngfolder or die "Can NOT mkdir $pngfolder:$!";
		}
		
		# pdf
		my $outgraphic1 = File::Spec -> catfile($pdffolder,removeFastaSuffix(basename($inputfile)) . ".pdf");
		my $cmd1 = "$seqlogoProgram -F pdf -f $inputfile2 -C 100 -h 4 -w 30 -c -e -n -p -Y -x Sequence -y Consensus > $outgraphic1";
		system($cmd1);
		
		# png
		my $outgraphic2 = File::Spec -> catfile($pngfolder,removeFastaSuffix(basename($inputfile)) . ".png");
		my $cmd2 = "$seqlogoProgram -F png -f $inputfile2 -C 100 -h 4 -w 30 -c -e -n -p -Y -x Sequence -y Consensus > $outgraphic2";
		system($cmd2);
		
	}
	
	#remove tmp file
	system("rm -rf $inputfile2 $rinput");
	
} 

1;

