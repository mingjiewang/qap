#!/usr/bin/perl

use Bio::Seq;
use Bio::SeqIO;
use warnings;
require "getopts.pl";

use vars qw($opt_f);
&Getopts('f:');

my $file = $opt_f;
open my $FILEHANDLE, $file or die "Could not open $file: $!";
my $ContigFileName = "";
my $CFH;
while(my $line = <$FILEHANDLE>)  
{
	if(substr($line,0,1) eq '>')
	{
		@header=split('\|',substr($line,1,length($line)));
		$newcontigfile = $header[0].".txt";
		open $CFH, ">".$newcontigfile or die $!;	
		print $CFH $line;
	}
	else
	{
		print $CFH $line;
	}
}
exit;
