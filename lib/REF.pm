#!perl

package REF;

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
our @EXPORT          = qw ($DEBUG_MODE copyRegion moveRegion subRegion) ;
our @EXPORT_OK       = qw ($DEBUG_MODE copyRegion moveRegion subRegion);
our %EXPORT_TAGS     = ( DEFAULT => [qw (&copyRegion &moveRegion &subRegion) ] );

##Debug mode prints more information while running 
our $DEBUG_MODE = 1;


##subprograms begin here##
sub copyRegion {
	##input fasta files must be in 2 line mode
	my $inputfile = shift;
	my $outputfile = shift;
	my $direction = shift; #'head2tail' or 'tail2head'
	my $length = shift;
	
	my @seqID;
	my @seq;
	open T,$inputfile or die "Can NOT open fasta file $inputfile:$!";
	while(my $line1 = <T>){
		chomp $line1;
		chomp (my $line2 = <T>);
		
		if ($line1 !~ /^>/){
			InfoError("The format of input file $inputfile should be \'fasta\'.");
			return 0;
		}
		
		$line1 =~ s/^\>//;
		push @seqID,$line1;
		push @seq,$line2;
	}
	close T;
	
	if (CheckPositiveInt($length) or $length == 0){
		#nothing
	}else{
		InfoError("The length of sequences to be copied should be a positive integer.");
		return 0;
	}
	
	my @seqnew;
	if($direction eq 'head2tail'){
		for my $seq (@seq){
			if($length >= length($seq)){
				InfoError("The length of copy region can not be greater than sequence lenght.");
				return 0;
			}
			my $seqnew = $seq . substr($seq,0,$length);
			push @seqnew,$seqnew;
		}
	}elsif($direction eq 'tail2head'){
		for my $seq (@seq){
			if($length >= length($seq)){
				InfoError("The length of copy region can not be greater than sequence lenght.");
				return 0;
			}
			my $seqnew = substr($seq,length($seq) - $length,$length) . $seq;
			push @seqnew,$seqnew;
		}
	}else{
		InfoError("Error direction when copying sequences. Please use \'head2tail\' or \'tail2head\'.");
		return 0;
	}
	
	open OUT,">$outputfile" or die "Can NOT output to file $outputfile:$!";
	for my $i (0..scalar(@seqID) - 1){
		print OUT ">$seqID[$i]\n";
		print OUT "$seqnew[$i]\n"; 
	}
	close OUT;
	
	my $cmd = "sed -i \'s\/\\r\/\/g\' $outputfile";
	system($cmd);
}

sub moveRegion {
	##input fasta files must be in 2 line mode
	my $inputfile = shift;
	my $outputfile = shift;
	my $direction = shift; #'head2tail' or 'tail2head'
	my $length = shift;
	
	my @seqID;
	my @seq;
	open T,$inputfile or die "Can NOT open fasta file $inputfile:$!";
	while(my $line1 = <T>){
		chomp $line1;
		chomp (my $line2 = <T>);
		
		if ($line1 !~ /^>/){
			InfoError("The format of input file $inputfile should be \'fasta\'.");
			return 0;
		}
		
		$line1 =~ s/^\>//;
		push @seqID,$line1;
		push @seq,$line2;
	}
	close T;
	
	if (CheckPositiveInt($length) or $length == 0){
		#nothing
	}else{
		InfoError("The length of sequences to be copied should be a positive integer.");
		return 0;
	}
	
	my @seqnew;
	if($direction eq 'head2tail'){
		for my $seq (@seq){
			if($length >= length($seq)){
				InfoError("The length of copy region can not be greater than sequence lenght.");
				return 0;
			}
			my $seqnew = substr($seq, $length, length($seq) - $length) . substr($seq,0,$length);
			push @seqnew,$seqnew;
		}
	}elsif($direction eq 'tail2head'){
		for my $seq (@seq){
			if($length >= length($seq)){
				InfoError("The length of copy region can not be greater than sequence lenght.");
				return 0;
			}
			my $seqnew = substr($seq,length($seq) - $length,$length) . substr($seq, 0, length($seq) - $length);
			push @seqnew,$seqnew;
		}
	}else{
		InfoError("Error direction when copying sequences. Please use \'head2tail\' or \'tail2head\'.");
		return 0;
	}
	
	open OUT,">$outputfile" or die "Can NOT output to file $outputfile:$!";
	for my $i (0..scalar(@seqID) - 1){
		print OUT ">$seqID[$i]\n";
		print OUT "$seqnew[$i]\n"; 
	}
	close OUT;
	
	my $cmd = "sed -i \'s\/\\r\/\/g\' $outputfile";
	system($cmd);
}

sub subRegion {
	my $inputfile = shift;
	my $outputfile = shift;
	my $start = shift;
	my $end = shift;
	
	my $seqlen = checkSeqLen($inputfile,dirname($outputfile));
	
	my @seqID;
	my @seq;
	open T,$inputfile or die "Can NOT open fasta file $inputfile:$!";
	while(my $line1 = <T>){
		chomp $line1;
		chomp (my $line2 = <T>);
		
		if ($line1 !~ /^>/){
			InfoError("The format of input file $inputfile should be \'fasta\'.");
			return 0;
		}
		
		$line1 =~ s/^\>//;
		push @seqID,$line1;
		push @seq,$line2;
	}
	close T;
	
	if (CheckPositiveInt($start) and CheckPositiveInt($end)){
		#nothing
	}else{
		InfoError("Start and end position should be a positive integer.");
		return 0;
	}
	
	if($start > $seqlen or $end > $seqlen){
		InfoError("Start or end position is greater than reference sequence length: ${seqlen}bp.");
		return 0;
	}
	
	if($end < $start){
		InfoError("Start position is greater than end position.");
		return 0;
	}
	
	my @seqnew;
	for my $seq (@seq){
		my $seqnew = substr($seq, $start - 1, $end - $start + 1);
		push @seqnew,$seqnew;
	}
	
	open OUT,">$outputfile" or die "Can NOT output to file $outputfile:$!";
	for my $i (0..scalar(@seqID) - 1){
		print OUT ">$seqID[$i]\n";
		print OUT "$seqnew[$i]\n"; 
	}
	close OUT;
	
	my $cmd = "sed -i \'s\/\\r\/\/g\' $outputfile";
	system($cmd);
}

1;

