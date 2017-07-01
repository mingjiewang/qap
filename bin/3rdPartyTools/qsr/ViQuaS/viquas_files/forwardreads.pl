#!/usr/bin/perl
use strict;
use Bio::Seq;
use Bio::SeqIO;

my($path, $names, $seqs ) = @ARGV;

my $opfile = $path."-readsnew.fa";

open names, "$names" or die $!;
open seqs, "$seqs" or die $!;

my @allnames;
my @allseqs;

while (<names>) {
    push (@allnames, $_);
}
while (<seqs>) {
    push (@allseqs, $_);
}
chomp(@allnames);
chomp(@allseqs);
for (my $i = 0; $i < @allnames; $i++)
	{
		my $newname = $allnames[$i];
		my $newseq = $allseqs[$i];
		
		my $newread = Bio::Seq->new(-id => $newname, -seq => $newseq);
		my $newout = Bio::SeqIO->new(-file => '>>'.$opfile, -format => 'Fasta');
		$newout -> write_seq($newread);
	}
exit;
