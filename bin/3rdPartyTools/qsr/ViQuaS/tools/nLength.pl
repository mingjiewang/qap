#!/usr/bin/perl
#Rene Warren 2002-2008

if ($#ARGV<1){
   die "$0 <fasta file needed> <N=?  (e.g. for the N50, write 50)>\n";
}  

if(! -e $ARGV[0]){
   die "$ARGV[0] doesn't exist -- fatal.\n";
}

if($ARGV[1] <1 || $ARGV[1] > 100){
   die "Please specify a valid number between 1-100 -- fatal.\n";
}

open (OUT, $ARGV[0]) || die "can't open $ARGV[0] for reading -- fatal.\n";
my ($seq,$prev) = ("","");
my $bins;

while (<OUT>){
   chomp;
   if (/\>(\S+)/){
      my $tig = $1; 
      if ($prev ne $tig && $prev ne ""){
         $bins->{$prev}=length($seq);
      }
      $seq='';
      $prev=$tig;
   }elsif(/^([ACGNTX]*)$/){
      $seq .= $_; 
   }    
}  
$bins->{$prev}=length($seq);
close OUT;

my $n = &calculateN($bins,$ARGV[1]);
print "N$ARGV[1] = $n bp ($ARGV[1]% of the bases in your assembly are in contigs of length $n bp or larger.)\n";

exit;

#------------------------------------------------------------------------
sub calculateN{
    # Ny length = the length L such that y% of all base-pairs are contained in contigs of this length or larger

    my ($bins,$portion) = @_; # reference to an array of contig sizes

    my ($n,$total_size, $contig_number, $tig_size); # total_size = total bases in assembly

    # get total size
    foreach my $tig (sort {$bins->{$a}<=>$bins->{$b}} keys %$bins){
        $total_size += $bins->{$tig};
    }

    my $size_subset = $total_size * ($portion / 100);

    my $running_count=0;

    # look for largest size contig containing size_subset% of bases

    foreach my $tig (sort {$bins->{$b}<=>$bins->{$a}} keys %$bins){
        $tig_size=$bins->{$tig};

        $running_count+=$tig_size;
        if ($running_count >= $size_subset){
            $n=$tig_size;
            last;
        }
    }

    return $n;
}#end sub

