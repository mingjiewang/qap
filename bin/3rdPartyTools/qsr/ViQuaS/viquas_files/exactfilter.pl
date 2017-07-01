#!/usr/bin/perl

use Bio::Seq;
use Bio::SeqIO;

use warnings;

require "getopts.pl";

use vars qw($opt_f);
&Getopts('f:');

local $now = time;

# Reading Reference Genome

my $referenceSeqIO = Bio::SeqIO->new(-file => "SIM/reference.fsa", '-format' => 'Fasta');
my $referenceSeq = ($referenceSeqIO->next_seq)->seq;

# Reading Samples
 
my $sampleSeqIO = Bio::SeqIO->new(-file => $opt_f, '-format' => 'Fasta');

while(my $nextSampleSeq = $sampleSeqIO->next_seq) 
{
	my $sampleSeq = $nextSampleSeq->seq;
	my $sampleID = $nextSampleSeq->display_id;
	isPerfectMatch($sampleID,$sampleSeq,$referenceSeq);
}

sub isPerfectMatch
{
	my($sampleID, $sample, $referenceGenome) = @_;

	my $sampleLength = length($sample);
	my $referenceLength = length($referenceGenome);
	my $numberIterations = $referenceLength - $sampleLength + 1;

	for (my $k = 0; $k < $numberIterations; $k++)
	{
		my $referenceSubString = substr($referenceGenome, $k, $sampleLength);

		if($referenceSubString eq $sample)
		{
			my $log = "commonreads.fsa";
			open (LOG, ">>$log") || die "Can't write to $log -- fatal\n";
			print LOG ">".$sampleID."\n";
			print LOG $sample."\n";

			my $log2 = "rdlen.txt";
			open (LOG, ">>$log2") || die "Can't write to $log -- fatal\n";
			print LOG length($sample)."\n";
			return;
		}	
	}

	my $log = "filtered.fsa";
	open (LOG, ">>$log") || die "Can't write to $log -- fatal\n";
	print LOG ">".$sampleID."\n";
	print LOG $sample."\n";

	return;
}

$now = time - $now;
#printf("\nTotal running time: %02d:%02d:%02d\n\n", int($now / 3600), int(($now % 3600) / 60), int($now % 60));

exit;
