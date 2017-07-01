#!/usr/bin/perl

use strict;
use Data::Dumper;
require "getopts.pl";

use vars qw($opt_f $opt_m $opt_o $opt_v $opt_r $opt_p $opt_d $opt_k $opt_a $opt_z $opt_e $opt_g $opt_s $opt_t $opt_b $opt_c $opt_x $opt_n $opt_h $opt_i $opt_w);
&Getopts('f:m:o:v:r:p:d:k:a:z:e:g:s:t:b:c:x:n:h:i:w:');

my ($base_overlap,$min_overlap,$verbose,$MIN_READ_LENGTH,$SEQ_SLIDE,$min_base_ratio,$paired,$min_links,$max_link_ratio,$contig_size_cutoff,$insert_stdev,$unpaired_file,$seed_file,$max_trim,$base_name,$tracked,$forcetrack,$max_count_trim,$min_tig_overlap,$npad_gaps,$ignorehead,$space_restriction,$min_depth_of_coverage)=(2,20,0,21,1,0.7,0,4,0.50,100,0.75,"no-g","no-s",0,"",0,0,10,20,0,0,0,0);

my $version = "[v3.8]";
my $per;
my $MAX = 0;
my $MAX_TOP = 1500; # this is the very maximum anchoring edge sequence that will be searched (designed for use with -s to prevent long searches)
my $TRACK_COUNT = 0;
my $illuminaLengthCutoff = 300; ### means all sequence reads greater than this are not illumina sequences

#-------------------------------------------------

if(! $opt_f || ! $opt_w){
   print "Usage: $0 $version\n";
   print "-f  File containing all the [paired (-p 1)] reads (required)\n";
   print "\t  With -p 1:\n";
   print "\t! Target insert size must be indicated at the end of the header line (e.g. :200 for a 200bp insert)\n";
   print "\t! Paired reads must be separated by \":\"\n";
   print "\t  >template_name:200\n\t  ACGACACTATGCATAAGCAGACGAGCAGCGACGCAGCACG:GCGCACGACGCAGCACAGCAGCAGACGAC\n";
   print "-w  Minimum depth of coverage allowed for contigs (e.g. -w 1 = process all reads [v3.7 behavior], required)\n";
   print "    *The assembly will stop when 50+ contigs with coverage < -w have been seen.*\n";
   print "-s  Fasta file containing sequences to use as seeds exclusively (specify only if different from read set, optional)\n";
   print "\t-i Apply read space restriction to seeds (TASR behavior) while -s option in use (-i 1 = yes, default = no , optional)\n"; 
   print "-m  Minimum number of overlapping bases with the seed/contig during overhang consensus build up (default -m $min_overlap)\n";
   print "-o  Minimum number of reads needed to call a base during an extension (default -o $base_overlap)\n";
   print "-r  Minimum base ratio used to accept a overhang consensus base (default -r $min_base_ratio)\n";
   print "-t  Trim up to -t base(s) on the contig end when all possibilities have been exhausted for an extension (default -t $max_trim, optional)\n";
   print "-c  Track base coverage and read position for each contig (default -c $tracked, optional)\n";
   print "-h  Ignore read name/header *will use less RAM if set to -h 1* (-h 1 = yes, default = no, optional)\n";
   print "-b  Base name for your output files (optional)\n";
   print "-z  Minimum contig size to track base coverage and read position (default -z $contig_size_cutoff, optional)\n";
   print "-p  Paired-end reads used? (-p 1 = yes, default = no, optional)\n";
   print "-v  Runs in verbose mode (-v 1 = yes, default = no, optional)\n";
   print "============ Options below only considered with -p 1 ============\n";
   print "-e  Error (%) allowed on mean distance   e.g. -e 0.75  == distance +/- 75% (default -e $insert_stdev, optional)\n";
   print "-k  Minimum number of links (read pairs) to compute scaffold (default -k $min_links, optional)\n";
   print "-a  Maximum link ratio between two best contig pairs *higher values lead to least accurate scaffolding* (default -a $max_link_ratio, optional)\n";
   print "-x  Minimum overlap required between contigs to merge adjacent contigs in a scaffold (default -x $min_tig_overlap, optional)\n";
   print "-n  N-pad gaps (-n 1 = yes, default = no $npad_gaps, optional)\n";
   die "-g  Fasta file containing unpaired sequence reads (optional)\n";
}

my $file = $opt_f;
$min_overlap = $opt_m if($opt_m);
$base_overlap = $opt_o if($opt_o);
$min_base_ratio = $opt_r if($opt_r);
$max_trim = $opt_t if($opt_t);
$verbose = $opt_v if($opt_v);
$paired = $opt_p if($opt_p);
$min_links = $opt_k if($opt_k);
$max_link_ratio = $opt_a if($opt_a);
$contig_size_cutoff = $opt_z if($opt_z);
$insert_stdev = $opt_e if($opt_e);
$unpaired_file = $opt_g if($opt_g);
$seed_file = $opt_s if($opt_s);
$base_name = $opt_b if($opt_b);
$tracked = $opt_c if($opt_c);
$min_tig_overlap = $opt_x if($opt_x);
$npad_gaps = $opt_n if($opt_n);
$ignorehead = $opt_h if($opt_h);
$min_depth_of_coverage = $opt_w if($opt_w); #### Lower bound on contig coverage

if($paired || $tracked){ $forcetrack = 1; }

my $display_unpaired_file = $1 if ($unpaired_file=~/([^\/]*)$/);
my $display_seed_file = $1 if ($seed_file=~/([^\/]*)$/);

#-------------------------------------------------

if(! -e $file){
   die "Invalid file: $file -- fatal\n";
}

if($opt_s && ! -e $opt_s){
   die "The file $opt_s you specified does not exist -- fatal\n";
}else{
   $space_restriction = $opt_i if($opt_i);
}
$base_name = "contigs";
### Naming output files
if ($base_name eq ""){

   $base_name = $file . ".ssake_m" . $min_overlap . "_o" . $base_overlap . "_r" . $min_base_ratio . "_t" . $max_trim . "_w" . $min_depth_of_coverage;

   if($paired){
      $base_name .= "_e" . $insert_stdev . "_k" . $min_links . "_a" . $max_link_ratio . "_z" . $contig_size_cutoff . "_x" . $min_tig_overlap  . "_g-" . $display_unpaired_file;
   }
   if($opt_s){
      $base_name .= "_s-" . $display_seed_file;
   }

   my $pid_num = getpgrp(0);
   $base_name .= "_pid" . $pid_num;
}

my $contig = $base_name .  ".fsa";
my $singlet = $base_name . ".singlets";
my $short = $base_name . ".short";
my $log = $base_name . ".log";
my $scaffold = $base_name . ".scaffolds" if ($paired);
my $mergedtigs = $base_name . ".mergedcontigs" if ($paired);
my $issues = $base_name . ".pairing_issues" if ($paired);
my $distribution = $base_name . ".pairing_distribution.csv" if ($paired);
my $covfile = $base_name . ".coverage.csv" if ($tracked);
my $rdpositionfile = $base_name . ".readposition" if ($tracked);
my $pileupfile = $base_name . ".pileup" if ($space_restriction);


open (LOG, ">$log") || die "Can't write to $log -- fatal\n";

if($min_overlap < 16 || $min_overlap > 100){
   my $outofbound_message = "-m must be a number between 16-100 ...Exiting.\n";
   print $outofbound_message;
   print LOG $outofbound_message;
   close LOG;
   exit;
}

if($base_overlap < 1){
   my $outofbound_message = "-o must be set to 1 or higher ...Exiting.\n";
   print $outofbound_message;
   print LOG $outofbound_message;
   close LOG;
   exit;
}

if($min_base_ratio <= 0.5 || $min_base_ratio > 1){
   my $outofbound_message = "-r must be a number between 0.51 and 1.00 ...Exiting.\n";
   print $outofbound_message;
   print LOG $outofbound_message;
   close LOG;
   exit;
}

#-------------------------------------------------

my $init_message = "\nRunning: $0 $version\n-f $file\n-s $seed_file\n\t-i $space_restriction\n-h $ignorehead\n-w $min_depth_of_coverage\n-m $min_overlap\n-o $base_overlap\n-r $min_base_ratio\n-t $max_trim\n";
if($tracked){$init_message .= "-c $tracked\nCoverage: $covfile\nRead position: $rdpositionfile\n";$init_message .= "Pileup: $pileupfile\n" if($space_restriction);}
if($forcetrack){$init_message .= "-z $contig_size_cutoff\n";}
if($paired){$init_message .= "-p $paired\n-e $insert_stdev\n-k $min_links\n-a $max_link_ratio\n-x $min_tig_overlap\nUnpaired reads (optional) -g $unpaired_file\nScaffolds: $scaffold\nMerged contigs: $mergedtigs\nPairing issues: $issues\nPairing distance distribution: $distribution\n";}
$init_message .= "\nContigs: $contig\nSinglets: $singlet\n\nExcluded reads: $short\nLog: $log\n";

print $init_message;
print LOG $init_message;

#-------------------------------------------------

my $date = `date`;
chomp($date);

my $reading_reads_message = "\n=>Reading sequences initiated $date\n";
print $reading_reads_message;
print LOG $reading_reads_message;

my $encoded = &encodeBases();

my ($seed,$seedsplit);

#-------------------------------------------------
### Allow user to specify a fasta file containing sequences to use as seeds, exclusively
if(-e $opt_s){
   my $use_seed_message = "Using seed sequence file $opt_s for this assembly.\nNote: ONLY sequences in $opt_s will be used as seeds (i.e. -f $opt_f and -g $opt_g will NOT be used as seeds, only used for extension)\n";
   print LOG $use_seed_message;
   print $use_seed_message if ($verbose);
   ($seed,$seedsplit) = &loadSeed($opt_s); 
}

my($set,$bin,$matepair);
($set,$bin,$matepair) = &readFasta($matepair,$set,$bin,$file,$short,$paired,$encoded,$seedsplit,$space_restriction,$ignorehead);
($set,$bin,$matepair) = &readFasta($matepair,$set,$bin,$unpaired_file,$short,0,$encoded,$seedsplit,$space_restriction,$ignorehead) if (-e $opt_g);

if(! $opt_s){
   $seed = $set;
   $TRACK_COUNT = 0;
}

my $seed_number_message = "Number of unique seed sequences: " . keys( %$seed ) . "\n";
printf $seed_number_message;
print LOG $seed_number_message;
#-------------------------------------------------

$date = `date`;
chomp($date);

my $ssake_start_message = "\n=>Sequence assembly initiated $date\n";
print $ssake_start_message;
print LOG $ssake_start_message;

#-------------------------------------------------
my ($sgl_count,$tig_count,$previous_index) = (1,1,0);

open (TIG, ">$contig") || die "Can't write to $contig -- fatal\n";
open (SIN, ">$singlet") || die "Can't write to $singlet -- fatal\n";
if ($paired){open (SC, ">$scaffold") || die "Can't write to $scaffold -- fatal\n";}
if($tracked){open (CF, ">$covfile") || die "Can't write to $covfile -- fatal\n";}
if($tracked){open (RP, ">$rdpositionfile") || die "Can't write to $rdpositionfile -- fatal\n";}
if($space_restriction){open (PU, ">$pileupfile") || die "Can't write to $pileupfile -- fatal\n";}

my ($tig_length,$track_all,$alternate);

eval{

my $status_bar = "+";
for(my $i=1;$i<=99;$i++){
   $per->{$i}++;
   my $ct = $i /10;
   if($ct == int($ct)){$status_bar .= $ct;}else{$status_bar .= "-";}
}

$status_bar .= "+ x 10 [% complete]";
print "$status_bar\n.";
my $keys_start = keys ( %$seed );

my $low_total = 0;

#--------------------------------------------
ASSEMBLY:
foreach my $seq (sort {$seed->{$b}{'count'}<=>$seed->{$a}{'count'}} keys %$seed){#cycle through the input [normalized] reads

   my $track;

   my ($pu_seed_name,$pu_seed_ori)=("","");
   if($space_restriction){
      $pu_seed_name = $seed->{$seq}{'seed_name'};
      $pu_seed_ori = $seed->{$seq}{'ori'};
   }

   if(defined $seed->{$seq}){#sequence read hasn't been used, is longer than 16 nt and the user-defined overlap minimum -m

      my $seed_name = "";
      if(defined $seed->{$seq}{'seed_name'}){$seed_name = "|seed:" . $seed->{$seq}{'seed_name'};}

      my $orig_mer = length($seq);

      if($forcetrack){
         $track->{$seq}{'start'} = 1;
         $track->{$seq}{'end'} = $orig_mer;
         $track->{$seq}{'cov'} = $seed->{$seq}{'count'};
         $track->{$seq}{'names'} = $seed->{$seq}{'names'};
      }

      #### Delete keys ref

      my $start_sequence = $seq;
      my $reads_needed = $seed->{$seq}{'count'};                       #tracks coverage
      my $total_bases = $orig_mer * $reads_needed;
      ($bin,$set,$seed) = deleteData($bin,$set,$seed,$seq,$encoded);   #remove k-mer from hash table and prefix tree
     
      print "\n\n>>> START SEED SEQUENCE :: $seq <<<\n\n" if ($verbose);

      ($seq, $set, $bin, $reads_needed, $total_bases, $track) = doExtension("3 prime", $orig_mer, $seq, $set, $bin, $reads_needed, $total_bases, $min_overlap, $base_overlap, $min_base_ratio, $verbose, $track, $forcetrack, $tig_count, $max_trim, $encoded, $matepair);
      ####end of 3' extension, beginning of 5' extension  (via 3' RC)
      my $seqrc = reverseComplement($seq);
      ($seqrc, $set, $bin, $reads_needed, $total_bases, $track) = doExtension("5 prime", $orig_mer, $seqrc, $set, $bin, $reads_needed, $total_bases, $min_overlap, $base_overlap, $min_base_ratio, $verbose, $track, $forcetrack, $tig_count, $max_trim, $encoded, $matepair);
      ####end of 5' extension
      my $leng = length($seqrc);
      my $reversetig = reverseComplement($seqrc);                   ### return to sequence, as provided 
      my $trimmed_length = length($start_sequence) - 2*($max_trim);

      if($leng > $trimmed_length){ ### commented out: && $start_sequence ne $seqrc && $start_sequence ne $reversetig
         $tig_length->{$tig_count} = $leng;
         my $cov =  $total_bases / $leng;
         $low_total++ if ($cov < $min_depth_of_coverage);
         last ASSEMBLY if ($low_total >= 50);             ### sudden termination based on expected depth ov coverage

         printf TIG ">contig%i|size%i|read%i|cov%.2f$seed_name\n%s\n", ($tig_count,$leng,$reads_needed,$cov,$reversetig);    #print contigs to file

         if($forcetrack && $leng >= $contig_size_cutoff){

            if($tracked){ ### only execute & report if user specifies -c
               printf CF ">contig%i|size%i|read%i|cov%.2f$seed_name\n", ($tig_count,$leng,$reads_needed,$cov);
               printf RP ">contig%i|size%i|read%i|cov%.2f$seed_name\n", ($tig_count,$leng,$reads_needed,$cov);

               my @tigarr=split(//,$reversetig);
               my $pileup;
               my $gpos=0;
               if($space_restriction){
                  foreach my $gbase(@tigarr){
                     $gpos++;
                     $pileup->{$gpos}{'tig'}=$gbase;
                  }
               }

               ### initialize all positions;
               my $hashcov;
               for (my $bpo=1;$bpo<=$leng;$bpo++){
                  $hashcov->{$bpo}=0;
               }

               foreach my $rd (sort {$track->{$a}{'start'}<=>$track->{$b}{'start'}} keys %$track){
                  my $rdlist = $track->{$rd}{'names'};
                  foreach my $rdoflist (keys %$rdlist){
                     printf RP "$rdoflist,$track->{$rd}{'start'},$track->{$rd}{'end'}\n" if($rdoflist ne "");
                 
                     my @covposition;
                     if($track->{$rd}{'start'} < $track->{$rd}{'end'}){### plus strand
                        my ($sss,$eee) = ($track->{$rd}{'start'},$track->{$rd}{'end'});
                        if($sss<0){$sss=1;}###<<
                        if($eee>$leng){$eee=$leng;}###<<
                        @covposition = ($sss .. $eee);

                        my @rdseq=split(//,$rd);
                        my @rdqua=split(//,$rdlist->{$rdoflist});
                        my $rdpos=0;

                        if($space_restriction){
                           #foreach my $vpos(@covposition){
                           for(my $vpos=$track->{$rd}{'start'};$vpos<=$track->{$rd}{'end'};$vpos++){
                              if($vpos>=1 && $vpos<=$gpos){
                                 #print "***$rdoflist***$seed->{$rd}{'seed_name'}***\n";
                                 if($rdoflist eq $pu_seed_name){
                                    my @seedarr = split(//,$pu_seed_ori);
                                    $pileup->{$vpos}{'interest'} = $seedarr[$rdpos] if($seedarr[$rdpos]=~/[acgt]/);
                                    #print "$pileup->{$vpos}{'interest'} $rdoflist eq $seed->{$rd}{'seed_name'}\n";
                                 }

                                 my $tmpq = $rdqua[$rdpos];
                                 $pileup->{$vpos}{'qua'} .= "$tmpq";
                                 #if($vpos == $track->{$rd}{'start'}){
                                 #   $pileup->{$vpos}{'seq'} .= "^";
                                 #}elsif($vpos == $track->{$rd}{'end'}){
                                 #   $pileup->{$vpos}{'seq'} .= "\$";
                                 #}else{
                                    if($pileup->{$vpos}{'tig'} eq $rdseq[$rdpos]){
                                       $pileup->{$vpos}{'seq'} .= ".";
                                    }else{
                                       $pileup->{$vpos}{'seq'} .= $rdseq[$rdpos];
                                    }
                                 #}
                              }
                              $rdpos++;
                           }
                        }
                     }else{### minus strand
                        my ($sss,$eee) = ($track->{$rd}{'end'},$track->{$rd}{'start'});
                        if($sss<0){$sss=1;}###<<
                        if($eee>$leng){$eee=$leng;}###<<
                        @covposition = ($sss .. $eee);

                        my @rdseq=split(//,reverseComplement($rd));
                        my @rdqua=split(//,reverse($rdlist->{$rdoflist}));
                        my $rdpos=0;
                        if($space_restriction){
                           #foreach my $vpos(@covposition){
                           for(my $vpos=$track->{$rd}{'end'};$vpos<=$track->{$rd}{'start'};$vpos++){
                              if($vpos>=1 && $vpos<=$gpos){
                                 my $tmpq = $rdqua[$rdpos];
                                 $pileup->{$vpos}{'qua'} .= "$tmpq";
                                 #if($vpos == $track->{$rd}{'start'}){
                                 #   $pileup->{$vpos}{'seq'} .= "^";
                                 #}elsif($vpos == $track->{$rd}{'end'}){
                                 #   $pileup->{$vpos}{'seq'} .= "\$";
                                 #}else{
                                    if($pileup->{$vpos}{'tig'} eq $rdseq[$rdpos]){
                                       $pileup->{$vpos}{'seq'} .= ",";
                                    }else{
                                       $pileup->{$vpos}{'seq'} .= lc($rdseq[$rdpos]);
                                    }
                                 #}
                              }
                              $rdpos++;
                           }
                        }#space restriction
                     }#plus/minus
                     foreach my $pss (@covposition){$hashcov->{$pss}++;}
                  }     
               }

               foreach my $pss (sort {$a<=>$b} keys %$hashcov){
                  $hashcov->{$pss}=$base_overlap if(! $hashcov->{$pss});
                  print CF "$hashcov->{$pss},";
               }
               print CF "\n";
              
               if($space_restriction){
                  foreach my $tigpos (sort {$a<=>$b} keys %$pileup){
                     my $depth = length($pileup->{$tigpos}{'seq'});
                     my $base = "";
                     if($pileup->{$tigpos}{'interest'} ne ""){
                        $base = $pileup->{$tigpos}{'interest'};
                     }else{
                        $base = $pileup->{$tigpos}{'tig'};
                     }
                     print PU "contig$tig_count\t$tigpos\t$base\t$depth\t$pileup->{$tigpos}{'seq'} $pileup->{$tigpos}{'qua'}\n";
                  }
                  print PU "\n";
               }
            }
            ($track_all,$alternate) = &trackReads($track,$track_all,$alternate,$tig_count);  ### all pairs from all contigs (track for scaffolding)
         }

         $tig_count++;
      }else{
         my $cov = $reads_needed;
         my $singlet_leng = length($start_sequence);
         printf SIN ">singlet%i|size%i|read%i|cov%.2f$seed_name\n%s\n", ($sgl_count,$singlet_leng,$reads_needed,$cov,$start_sequence);    #print singlets to file
         $sgl_count++;
      }
   }

   my $keys_left = keys( %$seed );
   my $index = (int((($keys_start-$keys_left)/$keys_start)*100));
   if(defined $per->{$index}){
      print "." x ($index - $previous_index);
      $|=1; ###clear buffer
      delete $per->{$index};
   }
   $previous_index = $index;

   last ASSEMBLY if (! $keys_left);
}
print ".";
};###end eval block

$date = `date`;
chomp($date);

if($@){
   my $message = $@;
   my $failure = "\nSomething went wrong running $0 $date\n$message\n";
   print $failure;
   print LOG $failure; 
}else{
   my $success = "\nContig assembly executed normally $date\n";
   print $success;
   print LOG $success;
}

close TIG;
close SIN;
close SHO;
if($tracked){
   close CF;
   close RP;
}
if($space_restriction){
   close PU;
}

#------------------------------------
$date = `date`;
chomp($date);

if($paired){

   my $sc_start_message = "\n=>Scaffolding initiated $date\n";
   print $sc_start_message;
   print LOG $sc_start_message;

   my $pair  = &pairContigs($matepair, $track_all, $tig_length, $issues, $distribution, $verbose);
   &buildScaffolds($pair, $tig_length, $contig_size_cutoff, $verbose);
 
   close SC;
   $date = `date`;
   chomp($date);

   my $sc_end_message = "\nScaffolding ended $date\n";
   print $sc_end_message;
   print LOG $sc_end_message;

   my $me_start_message = "\n=>Merging contigs initiated $date\n";
   print $me_start_message;
   print LOG $me_start_message;

   &forcefillGaps($scaffold, $contig, $mergedtigs, $verbose, $min_tig_overlap, $max_count_trim, $npad_gaps, $alternate, $matepair, $min_overlap, $base_overlap, $min_base_ratio, $forcetrack, $max_trim, $encoded);

   $date = `date`;
   chomp($date);
   my $me_end_message = "\nMerging contigs ended $date\n";
   print $me_end_message;
   print LOG $me_end_message;

}

close LOG;
exit;

#-----------------------
sub collectOverhang{

   my ($overhang,$newpass,$dangle,$set,$verbose) = @_;

   my @over = split(//,$dangle);
   my $ct_oh = 0;

   foreach my $bz(@over){
      $ct_oh++;                                                   ### tracks overhang position passed the seed
      $overhang->{$ct_oh}{$bz} += $set->{$newpass}{'count'};      ### reflects read coverage (often real duplicates)
      print "$ct_oh - $bz = $overhang->{$ct_oh}{$bz}\n" if($verbose);
   }

   return $overhang;
}

#------------------------------------
#Order and orient contigs into scaffolds
sub buildScaffolds{

   my ($pair, $tig_length, $contig_size_cutoff, $verbose) = @_;

   my $seen_it;
   my $sc_ct = 0;
 
   #print SC "Scaffold Number,Scaffold Size (only contig lengths considered),Scaffold Chain: e.g. _f127z7068k12a0.58m42_r3090z62k7r0.14m76_  means: contig127(+ strand=f), size 7068 (z) has 12 links (k), link ratio of 0.58 (a) and with a mean gap/overlap of 42nt (m)  with reverse (r) of contig3090 (size 62) on the right.\n";

   SEED:
   foreach my $tig (sort {$tig_length->{$b}<=>$tig_length->{$a}} keys %$tig_length){

      my $ftig = "f" . $tig;
      my $rtig = "r" . $tig;

      if(! defined $seen_it->{$tig}){##should prevent re-using a contig as seed if it's already been incorporated into a scaffold

         $sc_ct++;

         my $chainleft = "";
          
         my $ori_chainright = $ftig . "Z" . $tig_length->{$tig};
         my $chainright = $ori_chainright;
         my $total = $tig_length->{$tig};

         ($total, $chainright, $seen_it) = &computeLayout("R", $chainright, $ftig, $pair, $tig_length, $total, $seen_it, $tig);
         ($total, $chainleft, $seen_it) = &computeLayout("L", $chainleft, $rtig, $pair, $tig_length, $total, $seen_it, $tig);

         $seen_it->{$tig}++;

         delete $pair->{$ftig};
         delete $pair->{$rtig};
         delete $tig_length->{$tig};

         my $scaffold = $chainleft . $chainright;
         print SC "scaffold" . $sc_ct . ",$total,$scaffold\n" if($total >= $contig_size_cutoff);
      }
   }
}

#------------------------------------
# links contigs together into a chain - must satisfy user-defined criterions (-k -a)
sub computeLayout{

   my ($ext, $chain, $tig, $pair, $tig_length, $total, $seen_it, $orig_tig_number) = @_;

   my $orig_tig = $tig;
   my $extension = 1;

   EXTENSION:
   while($extension){

      my $tnum = $1 if($tig=~/[fr](\d+)/);
      my $tnumf = "f" . $tnum;
      my $tnumr = "r" . $tnum;

      if(! defined $seen_it->{$tnum}){

         $seen_it->{$tnum}++ if($tnumf ne $orig_tig);

         print "Attempt to extend $tig\n" if ($verbose);      
         my $list = $pair->{$tig};
         my ($match1,$link1,$gaps1,$match2,$link2,$gaps2,$cntloop)=("",0,0,"",0,0,0);

         LINK:
         foreach my $match (sort {$list->{$b}{'links'}<=>$list->{$a}{'links'}} keys %$list){

            if($cntloop){
               ($match2,$link2,$gaps2) = ($match,$list->{$match}{'links'},$list->{$match}{'gaps'});
               print "$tig links second best $match2 (links:$link2 total sz:$gaps2)\n" if ($verbose);
               last LINK;
            }else{
               ($match1,$link1,$gaps1) = ($match,$list->{$match}{'links'},$list->{$match}{'gaps'});
               print "$tig links best $match1 (links:$link1 total sz:$gaps1)\n" if ($verbose);
            }
            $cntloop++;
         }

         ###ratio
         my $ratio = 0.00;
         $ratio = $link2 / $link1 if ($link1);        ## relative ratio of the two most abundant contig pairs
         if ($ratio =~ /(\d+\.\d{2})/){$ratio = $1;}
         ###mean
         my $mean = 0;
         $mean = int($gaps1 / $link1) if ($link1);

         my $tempnum = $1 if($match1 =~ /[fr](\d+)/);

         #### Assessment
         if(defined $seen_it->{$tempnum} || $link1 < $min_links || $ratio > $max_link_ratio || $tempnum == $orig_tig_number){
            $extension = 0;
            print "defined seen_it->{ $tempnum } || $link1 < $min_links || $ratio > $max_link_ratio\n L1:$link1 L2:$link2  M1:$match1 M2:$match2 G1:$gaps1 G2:$gaps2 "  if ($verbose);

            last EXTENSION;
         }{### pass filter.. does this contig 
            print "$ext extension.  mean: $mean links:$link1 linkratio:$ratio\n" if ($verbose);

            if($ext eq "R"){
               $chain .= "k" . $link1 . "a" . $ratio . "m" . $mean . "_" . $match1 . "z" . $tig_length->{$tempnum};
            }else{
               my $temp_match = "";
               if($match1 =~ /^r(\d+)/){$temp_match = "f" . $1;}else{$temp_match = "r". $1;}            
               $chain = $temp_match . "z" . $tig_length->{$tempnum} . "k" . $link1 . "a" . $ratio . "m" . $mean . "_" . $chain;
            }   
            $total += $tig_length->{$tempnum};

            print "NEXT TIG TO LOOK AT= $match1\n" if ($verbose);
            $tig = $match1;
            $extension = 1; 
          
            print "Will flag $tnum as seen  (only if $tnumf != $orig_tig)." if ($verbose);
   
            if($tnumf ne $orig_tig){
               delete $pair->{$tnumf};
               delete $pair->{$tnumr};
               delete $tig_length->{$tnum};
            }else{
               delete $pair->{$tnumf};
            }
         }
      }else{
         print "NO MORE MATCH FOR $tig in hash: pair>>\n" if ($verbose);
         $extension = 0;
         last EXTENSION;
      }
   }### pair is defined
   return $total, $chain, $seen_it;
}

#------------------------------------
sub trackReads{

   my ($track, $track_all, $alternate, $tig_count) = @_;

   foreach my $rd (keys %$track){
      if(! defined $track_all->{$rd}){
         $track_all->{$rd}{'tig'}   = $tig_count;
         $track_all->{$rd}{'start'} = $track->{$rd}{'start'};
         $track_all->{$rd}{'end'}   = $track->{$rd}{'end'};

         $alternate->{$tig_count}{$track->{$rd}{'start'}}{$rd} = $track->{$rd}{'end'};
         delete $track->{$rd};
      }
   }
   return $track_all,$alternate;
}

#------------------------------------
sub getDistance{

   my ($insert_size, $length_i, $start_i, $start_j) = @_;

   # L  ------  --------- R
   # i    ->        <-    j
   #      ....  ......    insert_span
   #      ============    insert_size

   my $insert_span = ($length_i - $start_i) + $start_j;
   my $gap_or_overlap = $insert_size - $insert_span;

   return $gap_or_overlap;
}

#-----------------
#build contig pairs based on template information  -  must satisfy user-defined criterions (-d -e)
sub pairContigs{

   my ($matepair,$track,$tig_length,$issues,$distribution,$verbose) = @_;
   my ($ct_illogical, $ct_ok_contig, $ct_ok_pairs, $ct_problem_pairs, $ct_iz_issues, $ct_single, $ct_both)= (0,0,0,0,0,0,0);
   my $ct_illogical_hash;
   my $ct_ok_contig_hash;
   my $ct_ok_pairs_hash;
   my $ct_problem_pairs_hash;
   my $ct_iz_issues_hash;
   my $ct_single_hash;
   my $ct_both_hash;

   my ($pair,$err,$track_insert);

   print "Pairing contigs...\n" if ($verbose);

   open(PET, ">$issues") || die "Can't open $issues for writing -- fatal\n";

   foreach my $read_a (keys %$matepair){ 

      my $mateslist = $matepair->{$read_a};

      foreach my $read_b (keys %$mateslist){

         if(! $matepair->{$read_a}{$read_b}{'bt'} && ! $matepair->{$read_b}{$read_a}{'bt'}){

            ##2 lines below indicates this specific pair has been seen
            $matepair->{$read_a}{$read_b}{'bt'}=1;
            $matepair->{$read_b}{$read_a}{'bt'}=1;

            my $insert_size = $mateslist->{$read_b}{'is'};
            my $min_allowed = -1 * ($insert_stdev * $insert_size);
            my ($low_iz, $up_iz) = ($insert_size + $min_allowed, $insert_size - $min_allowed);

            print "Pair read1=$read_a read2=$read_b\n" if ($verbose);

            if(defined $track->{$read_a}{'tig'} && defined $track->{$read_b}{'tig'}){### both pairs assembled

               $ct_both++;
               $ct_both_hash->{$insert_size}++;

               my $tig_a = $track->{$read_a}{'tig'};
               my $tig_b = $track->{$read_b}{'tig'};

               my $ftig_a = "f" . $tig_a;
               my $ftig_b = "f" . $tig_b;

               my $rtig_a = "r" . $tig_a;
               my $rtig_b = "r" . $tig_b;

               my $A_length = $tig_length->{$tig_a};
               my $A_start = $track->{$read_a}{'start'};
               my $A_end = $track->{$read_a}{'end'};
 
               my $B_length = $tig_length->{$tig_b};
               my $B_start = $track->{$read_b}{'start'} ;
               my $B_end = $track->{$read_b}{'end'};

               if ($tig_a != $tig_b){####paired reads located on <> contigs

                  ####Determine most likely possibility
                  if ($track->{$read_a}{'start'} < $track->{$read_a}{'end'}){

                     if ($track->{$read_b}{'end'} < $track->{$read_b}{'start'}){####-> <- :::  A-> <-B  /  rB -> <- rA
                         my $d = &getDistance($insert_size, $A_length, $A_start, $B_start);
                         print "A-> <-B  WITH $tig_a -> <- $tig_b GAP $d A=$A_length ($A_start-$A_end) B=$B_length ($B_start-$B_end) Alen, Astart,Bstart\n" if($verbose);
                         if($d >= $min_allowed){
                            $pair->{$ftig_a}{$ftig_b}{'links'}++;
                            $pair->{$ftig_a}{$ftig_b}{'gaps'} += $d;                  
                            $pair->{$rtig_b}{$rtig_a}{'links'}++;
                            $pair->{$rtig_b}{$rtig_a}{'gaps'} += $d;
                            $ct_ok_pairs++;
                            $ct_ok_pairs_hash->{$insert_size}++;
                         }else{
                            my $err_pair = $ftig_a . "-". $ftig_b;
                            $err->{$err_pair}{'links'}++;
                            $err->{$err_pair}{'gaps'} += $d;
                            $ct_problem_pairs++;
                            $ct_problem_pairs_hash->{$insert_size}++;
                            print PET "Pairs unsatisfied in distance within a contig pair.  A-> <-B  WITH tig#$tig_a -> $d <- tig#$tig_b, A=$A_length nt (start:$A_start, end:$A_end) B=$B_length nt (start:$B_start, end:$B_end) CALCULATED DISTANCE APART: $d < $min_allowed\n";
                         }
                      }else{#### -> -> ::: A-> <-rB  / B-> <-rA 
                         my $rB_start = $B_length - $B_start;
                         my $d = &getDistance($insert_size, $A_length, $A_start, $rB_start);
                         print "A-> <-rB  WITH $tig_a -> <- r.$tig_b GAP $d A=$A_length ($A_start-$A_end) B=$B_length ($B_start-$B_end) Alen,Astart,rBstart\n" if($verbose);
                         if($d >= $min_allowed){
                            $pair->{$ftig_a}{$rtig_b}{'links'}++;
                            $pair->{$ftig_a}{$rtig_b}{'gaps'} += $d;
                            $pair->{$ftig_b}{$rtig_a}{'links'}++;
                            $pair->{$ftig_b}{$rtig_a}{'gaps'} += $d;
                            $ct_ok_pairs++;
                            $ct_ok_pairs_hash->{$insert_size}++;
                         }else{
                            my $err_pair = $ftig_a . "-". $rtig_b;
                            $err->{$err_pair}{'links'}++;
                            $err->{$err_pair}{'gaps'} += $d;
                            $ct_problem_pairs++;
                            $ct_problem_pairs_hash->{$insert_size}++;
                            print PET "Pairs unsatisfied in distance within a contig pair.  A-> <-rB  WITH tig#$tig_a -> $d <- tig#r.$tig_b, A=$A_length  nt (start:$A_start, end:$A_end) B=$B_length nt (start:$B_start, end:$B_end) CALCULATED DISTANCE APART: $d < $min_allowed\n";
                         }
                      }
                  }else{

                     if ($track->{$read_b}{'end'} > $track->{$read_b}{'start'}){####<-  -> ::: B-> <-A / rA -> <- rB
                        my $d = &getDistance($insert_size, $B_length, $B_start, $A_start);
                        print "B-> <-A  WITH $tig_b -> <- $tig_a GAP $d A=$A_length ($A_start-$A_end) B=$B_length ($B_start-$B_end) Blen,Bstart,Astart\n" if($verbose);
                        if($d >= $min_allowed){
                           $pair->{$ftig_b}{$ftig_a}{'links'}++;
                           $pair->{$ftig_b}{$ftig_a}{'gaps'} += $d;
                           $pair->{$rtig_a}{$rtig_b}{'links'}++;
                           $pair->{$rtig_a}{$rtig_b}{'gaps'} += $d;
                           $ct_ok_pairs++;
                           $ct_ok_pairs_hash->{$insert_size}++;
                        }else{
                           my $err_pair = $ftig_b . "-". $ftig_a;
                           $err->{$err_pair}{'links'}++;
                           $err->{$err_pair}{'gaps'} += $d;
                           $ct_problem_pairs++;
                           $ct_problem_pairs_hash->{$insert_size}++;
                           print PET "Pairs unsatisfied in distance within a contig pair.  B-> <-A  WITH tig#$tig_b -> $d <- tig#$tig_a, B=$B_length nt (start:$B_start, end:$B_end) A=$A_length nt (start:$A_start, end:$A_end) CALCULATED DISTANCE APART: $d < $min_allowed\n";
                        }
                     }else{                          ####<- <-  :::  rB-> <-A / rA-> <-B
                        my $rB_start = $B_length - $B_start;
                        my $d = &getDistance($insert_size, $B_length, $rB_start, $A_start);
                        print "rB-> <-A WITH r.$tig_b -> <- $tig_a GAP $d A=$A_length ($A_start-$A_end) B=$B_length ($B_start-$B_end) Blen,rBstart,Astart\n" if($verbose);
                        if($d >= $min_allowed){
                           $pair->{$rtig_b}{$ftig_a}{'links'}++;
                           $pair->{$rtig_b}{$ftig_a}{'gaps'} += $d;
                           $pair->{$rtig_a}{$ftig_b}{'links'}++;
                           $pair->{$rtig_a}{$ftig_b}{'gaps'} += $d;
                           $ct_ok_pairs++;
                           $ct_ok_pairs_hash->{$insert_size}++;
                        }else{
                           my $err_pair = $rtig_b . "-". $ftig_a;
                           $err->{$err_pair}{'links'}++;
                           $err->{$err_pair}{'gaps'} += $d;
                           $ct_problem_pairs++;
                           $ct_problem_pairs_hash->{$insert_size}++;
                           print PET "Pairs unsatisfied in distance within a contig pair.  rB-> <-A WITH tig#r.$tig_b -> $d <- tig#$tig_a, B=$B_length nt (start:$B_start, end:$B_end) A=$A_length nt (start:$A_start, end:$A_end) CALCULATED DISTANCE APART: $d < $min_allowed\n";
                        }
                     }
                  }
                  #print Dumper($pair);
               }else{###Clone, paired reads located on the same contig -- could be used to investigate misassemblies
           
                  print "Pair ($read_a and $read_b) located on same contig $tig_a ($A_length nt)\n" if ($verbose);
                  my $pet_size = 0;

                  if ($A_start > $B_start && ($B_start < $B_end) && ($A_start > $A_end)){    # B --> <-- A
                     $pet_size = $A_start - $B_start;
                     $track_insert->{$pet_size}++;
                     if($pet_size >= $low_iz && $pet_size <= $up_iz){
                        $ct_ok_contig++;
                        $ct_ok_contig_hash->{$insert_size}++;
                     }else{
                        print PET "Pairs unsatisfied in distance within a contig.  Pair ($read_a - $read_b) on contig $tig_a ($A_length nt) Astart:$A_start Aend:$A_end Bstart:$B_start Bend:$B_end CALCULATED DISTANCE APART: $pet_size\n";
                        $ct_iz_issues++;
                        $ct_iz_issues_hash->{$insert_size}++;
                     }
                  }elsif($B_start > $A_start && ($B_start > $B_end) && ($A_start < $A_end)){ # A --> <-- B
                     $pet_size = $B_start - $A_start;
                     $track_insert->{$pet_size}++;
                     if($pet_size >= $low_iz && $pet_size <= $up_iz){
                        $ct_ok_contig++;
                        $ct_ok_contig_hash->{$insert_size}++;
                     }else{
                        print PET "Pairs unsatisfied in distance within a contig.  Pair ($read_a - $read_b) on contig $tig_a ($A_length nt) Astart:$A_start Aend:$A_end Bstart:$B_start Bend:$B_end CALCULATED DISTANCE APART: $pet_size\n";
                        $ct_iz_issues++;
                        $ct_iz_issues_hash->{$insert_size}++;
                     }
                  }else{
                     $ct_illogical++;
                     $ct_illogical_hash->{$insert_size}++;
                     print PET "Pairs unsatisfied in pairing logic within a contig.  Pair ($read_a - $read_b) on contig $tig_a ($A_length nt) Astart:$A_start Aend:$A_end Bstart:$B_start Bend:$B_end\n";
                  }
               }
            }else{###both pairs assembled
               $ct_single++;
               $ct_single_hash->{$insert_size}++;
            }
         }#if unseen
      }#pairing read b
   }#read a

   ### summary of contig pair issues
   print PET "------------- Putative issues with contig pairing - Summary  ----------------\n";
   foreach my $err_pair (sort {$err->{$b}{'links'}<=>$err->{$a}{'links'}} keys %$err){
      my $mean_iz = 0;
      $mean_iz = $err->{$err_pair}{'gaps'} / $err->{$err_pair}{'links'} if ($err->{$err_pair}{'links'});
      print PET "Pair $err_pair has $err->{$err_pair}{'links'} links and mean distance = $mean_iz\n";
   }
   close PET;
 
   my $satisfied = $ct_ok_pairs + $ct_ok_contig;
   my $unsatisfied = $ct_problem_pairs + $ct_iz_issues + $ct_illogical;
   my $ct_both_reads = $ct_both * 2;

   print LOG "\n===========PAIRED-END READS STATS===========\n"; 
   print LOG "At least one sequence/pair missing from contigs >= $contig_size_cutoff bp (user-defined -z): $ct_single\n";
   print LOG "Assembled pairs: $ct_both ($ct_both_reads sequences)\n";
   print LOG "\tSatisfied in distance/logic within contigs (i.e. -> <-, distance on target: $ct_ok_contig\n";
   print LOG "\tUnsatisfied in distance within contigs (i.e. distance out-of-bounds): $ct_iz_issues\n";
   print LOG "\tUnsatisfied pairing logic within contigs (i.e. illogical pairing ->->, <-<- or <-->): $ct_illogical\n";
   print LOG "\t---\n";
   print LOG "\tSatisfied in distance/logic within a given contig pair (pre-scaffold): $ct_ok_pairs\n";
   print LOG "\tUnsatisfied in distance within a given contig pair (i.e. calculated distances out-of-bounds): $ct_problem_pairs\n";
   print LOG "\t---\n";
   print LOG "Total satisfied: $satisfied\tunsatisfied: $unsatisfied\n\nBreakdown by insert sizes:\n";

   foreach my $izopt(sort {$a<=>$b} keys %$ct_both_hash){
      print LOG "--------Reads with $izopt bp inserts--------\n";
      my $maopt = -1 * ($insert_stdev * $izopt);
      my ($low_izopt, $up_izopt) = ($izopt + $maopt, $izopt - $maopt);
      print LOG "MIN:$low_izopt MAX:$up_izopt as defined by $izopt * $insert_stdev\n";
      print LOG "At least one sequence/pair missing from contigs >= $contig_size_cutoff bp (user-defined -z): $ct_single_hash->{$izopt}\n";
      print LOG "Assembled pairs: $ct_both_hash->{$izopt}\n";
      print LOG "\tSatisfied in distance/logic within contigs (i.e. -> <-, distance on target: $ct_ok_contig_hash->{$izopt}\n";
      print LOG "\tUnsatisfied in distance within contigs (i.e. distance out-of-bounds): $ct_iz_issues_hash->{$izopt}\n";
      print LOG "\tUnsatisfied pairing logic within contigs (i.e. illogical pairing ->->, <-<- or <-->): $ct_illogical_hash->{$izopt}\n";
      print LOG "\t---\n";
      print LOG "\tSatisfied in distance/logic within a given contig pair (pre-scaffold): $ct_ok_pairs_hash->{$izopt}\n";
      print LOG "\tUnsatisfied in distance within a given contig pair (i.e. calculated distances out-of-bounds): $ct_problem_pairs_hash->{$izopt}\n";
   }
   print LOG "============================================\n";

   open (CSV, ">$distribution") || die "Can't open $distribution for writing -- fatal";

   foreach my $is (sort {$a<=>$b} keys %$track_insert){
      print CSV "$is,$track_insert->{$is}\n";
   }

   close CSV;
   return $pair;
}

#-----------------
# SSAKE contig extension
sub doExtension{

   my ($direction, $orig_mer, $seq, $set, $bin, $reads_needed, $total_bases, $min_overlap, $base_overlap, $min_base_ratio, $verbose, $track, $paired, $tig_count, $max_trim, $e, $matepair) = @_;

   my $extended = 1;
   my $trim_ct = 0;     #trim counter - keeps track of 3'-end trim

   if($orig_mer > $MAX){$orig_mer=$MAX;}  ### Deals with special cases where the seed sequences are different from the read set (and possibly very large) - goal here is not to increase sequence coverage of seed, but rather to extend it.

   TRIM:
   while($trim_ct <= $max_trim){
      while($extended){
         my $growing_tig_length = length($seq);
         my ($pos,$span) = (0,"");

         ### Added 19March08
         if($growing_tig_length >= $MAX){   # $seq is length of contig being extended -- if larger than largest read, make sure the largest read could align and all subsequent rds.
            $span = $MAX - $TRACK_COUNT;
         }else{
            $span = $growing_tig_length - $TRACK_COUNT;
         }

         my $overhang = {};
         my $overlapping_reads = {};
         my $long = 0;

         for (my $x=1;$x <= ($orig_mer * 2);$x++){
            ($overhang->{$x}{'A'},$overhang->{$x}{'C'},$overhang->{$x}{'G'},$overhang->{$x}{'T'}) = (0,0,0,0);
         }

         ### COLLECT SEQUENCES 
         while ($span >= $min_overlap){  # will slide the subseq, until the user-defined min overlap size

            $pos = $growing_tig_length - $span;
            print "MAX:$MAX, SPAN:$span, POS:$pos" if ($verbose);

            my $subseq = substr($seq, $pos, $span);              #make a sub-sequence of length l-(1..i) for searching
            my @s = $subseq =~ /\S{4}/g;
            my $subset = $bin->{$e->{$s[0]}}{$e->{$s[1]}}{$e->{$s[2]}}{$e->{$s[3]}}; #Will grab everything even the reverse complement ones

            print "####$direction SEARCH Position:$pos Span:$span - Subseq:$subseq Previous:$seq\n" if ($verbose);

            ### SEARCH -- this cycles through limited k-mer space
            foreach my $pass (sort {$subset->{$b} <=> $subset->{$a}} keys %$subset){

               if($pass =~ /^$subseq([ACGT]+)/){#### OVERHANG 
                  #can we align perfectly that subseq to another rd start?
                  my $dangle = $1;
                  print "\n", "=" x 80, "\n$direction'- FOUND sequence: $pass -> subset: $subseq -> overhang: $dangle\n", "=" x 80, "\n\n" if ($verbose);
                  #---------------------------------
                  my $psr;
                  my $pass_rc = reverseComplement($pass);

                  if(defined $set->{$pass}){
                     $psr->{$pass}{'start'} = $pos + 1;
                     $psr->{$pass}{'end'} = $pos + length($pass);
                  }

                  if(defined $set->{$pass_rc}){
                     $psr->{$pass_rc}{'start'} = $pos + length($pass_rc);
                     $psr->{$pass_rc}{'end'} = $pos + 1;
                  }

                  ###############################################
                  # CONSIDER CERTAIN READS FOR OVERLAP, PREFERABLY THOSE WITH LOGICAL MATES AND FWD READS IN TIG LARGE ENOUGH TO HAVE SUCH PAIRS
                  foreach my $newpass(keys %$psr){
                     if($paired && defined $matepair->{$newpass} && ($psr->{$newpass}{'end'} < $psr->{$newpass}{'start'}) && ($psr->{$newpass}{'start'} >= $set->{$newpass}{'grace'})){#paired, pairingRds, <---, outside grace

                        #print "$B_end <-- $B_start *** [upper limit is grace]***\n";
                        #     <--B
                        #========
                        my $mateshash = $matepair->{$newpass};

                        MATESEARCH:
                        foreach my $matchingread (keys %$mateshash){
                           if(defined $track->{$matchingread}){              #a mate has been found on this contig
                              my $insert_size = $mateshash->{$matchingread}{'is'};
                              my $min_allowed = -1 * ($insert_stdev * $insert_size);
                              my ($low_iz, $up_iz) = ($insert_size + $min_allowed, $insert_size - $min_allowed);
                              my $A_start = $track->{$matchingread}{'start'};
                              my $A_end = $track->{$matchingread}{'end'};
                              my $pet_size = $psr->{$newpass}{'start'} - $A_start;
                              print "\t$newpass ($psr->{$newpass}{'start'} - $psr->{$newpass}{'end'}) and $matchingread ($A_start - $A_end ) are on same tig#$tig_count  [$low_iz to $up_iz]\n" if($verbose);
                              if(($psr->{$newpass}{'start'} > $A_start) && ($A_start < $A_end) && ($pet_size >= $low_iz) && ($pet_size <= $up_iz)){ # A --> <-- B(candidate for extension)
                                 #print "\tTRACKING $psr->{$newpass}{'start'} > $A_start && ($psr->{$newpass}{'start'} > $psr->{$newpass}{'end'}) && ($A_start < $A_end)\n" if($verbose);
                                 $overhang = collectOverhang($overhang,$newpass,$dangle,$set,$verbose);
                                 $overlapping_reads->{$pass}++;
                                 $long=1 if(length($newpass) > $illuminaLengthCutoff);
                                 last MATESEARCH;  
                              }
                           }
                        }
                     }else{### not paired or paired (B-->): collect all overlaps
                        $overhang = collectOverhang($overhang,$newpass,$dangle,$set,$verbose);
                        $overlapping_reads->{$pass}++;
                        $long=1 if(length($newpass) > $illuminaLengthCutoff);
                     }
                  }###for $newpass
                  #----------------------------------
               }elsif($subseq =~ /$pass/){     #### EMBEDDED

                  my $complement_pass = reverseComplement($pass);

                  print "$pass found in $subseq ($set->{$pass}{'count'}) - deleting read: $pass and complement ($set->{$complement_pass}): $complement_pass\n\n" if ($verbose);

                  if(defined $set->{$pass}){
                     my $current_reads = $set->{$pass}{'count'};
                     my $current_bases = length($pass) * $current_reads;
                     $reads_needed += $current_reads;
                     $total_bases += $current_bases;
                     if($paired){
                        $track->{$pass}{'start'} = $pos + 1;
                        $track->{$pass}{'end'} = $pos + length($pass);
                        $track->{$pass}{'cov'} = $current_reads;
                        $track->{$pass}{'names'} = $set->{$pass}{'names'};
                     }
                     ($bin,$set,$seed) = deleteData($bin,$set,$seed,$pass,$e);
                     #print "EMBED .. $pass ($current_reads)\n";
                  }

                  if(defined $set->{$complement_pass}){
                     my $current_reads = $set->{$complement_pass}{'count'};
                     my $current_bases = length($complement_pass) * $current_reads;
                     $reads_needed += $current_reads;
                     $total_bases += $current_bases;
                     if($paired){ 
                        $track->{$complement_pass}{'end'} = $pos + 1;
                        $track->{$complement_pass}{'start'} = $pos + length($complement_pass);
                        $track->{$complement_pass}{'cov'} = $current_reads;
                        $track->{$complement_pass}{'names'} = $set->{$complement_pass}{'names'};
                     }
                     ($bin,$set,$seed) = deleteData($bin,$set,$seed,$complement_pass,$e);
                     #print "EMBED .. $complement_pass ($current_reads)\n";
                  }
               }  
            }
            $span--;
         }#while overlap >= user-defined -m minimum

         my $consensus = "";
         print "Finished Collecting Overlapping Reads - BUILDING CONSENSUS...\n" if ($verbose);
         print Dumper(%$overlapping_reads) if ($verbose);

         my $tmp_base_overlap = $base_overlap;
         if($long){
            $tmp_base_overlap = 2;
         }

         ### Build consensus
         CONSENSUS:
         foreach my $ohpos (sort {$a<=>$b} keys %$overhang){
            if($ohpos){

               my $coverage = $overhang->{$ohpos}{'A'}+$overhang->{$ohpos}{'C'}+$overhang->{$ohpos}{'G'}+$overhang->{$ohpos}{'T'};
               print "pos:$ohpos cov:$coverage A:$overhang->{$ohpos}{'A'} C:$overhang->{$ohpos}{'C'} G:$overhang->{$ohpos}{'G'} T:$overhang->{$ohpos}{'T'}\n" if($verbose);

               if ($coverage < $tmp_base_overlap){
                  print "COVERAGE BELOW THRESHOLD: $coverage < -o $tmp_base_overlap @ $ohpos :: will extend by: $consensus\n" if ($verbose);
                  last CONSENSUS;
               }
               my $baselist = $overhang->{$ohpos};

               my $ct_dna=0;
               my $previous_bz = "";

               BASE:
               foreach my $bz (sort {$baselist->{$b}<=>$baselist->{$a}} keys %$baselist){
                  #print "\t$ct_dna -> $bz..$baselist->{$previous_bz} > $baselist->{$bz}\n";
                  if($ct_dna){## the two most abundant bases at that position
                     #print "\t\t$ct_dna\n";
                     if($previous_bz ne "" && ($baselist->{$previous_bz} / $coverage) >= $min_base_ratio && $baselist->{$previous_bz} > $baselist->{$bz}){### a simple consensus btw top 2 
                        $consensus .= $previous_bz;                                         ### build consensus
                        print "Added base $previous_bz (cov = $baselist->{$previous_bz}) to $consensus **\n" if ($verbose);
                        last BASE;
                     }else{
                        print "ISSUES EXTENDING: best base = $previous_bz (cov=$baselist->{$previous_bz}) at $ohpos.  Second-Best: $bz (cov=$baselist->{$bz}) (ratio best=$baselist->{$previous_bz} / total=$coverage) >= $min_base_ratio (-r) -- will terminate with $consensus\n" if($verbose);
                        last CONSENSUS;
                     }
                  }
                  $previous_bz = $bz;                 
                  $ct_dna++;
               }
            }
         }
         $long = 0;
         ### deal with sequence reads making up the consensus/newly formed contig
         if($consensus ne ""){

            print "Will extend $seq\nwith: $consensus\n\n" if($verbose);
            my $temp_sequence = $seq . $consensus;
            my $position_buffer = 0;
            if($growing_tig_length > $MAX){
               $temp_sequence = substr($seq,$growing_tig_length-$MAX,$MAX) . $consensus;
               $position_buffer = $growing_tig_length-$MAX;
            }
            my $integral = 0;
            my $temp_sequence_portion = "-" x length($temp_sequence);

            foreach my $ro (keys %$overlapping_reads){
               my $or_pos = -99;
               while($temp_sequence =~ /$ro/g){     ### want the last position
                  $or_pos = pos($temp_sequence);
               } 
               if($or_pos > 0){

                  ###TRACK COVERAGE TO PREVENT FAULTY EXTENSIONS
                  my $linestring = '.' x length($ro);
                  substr ($temp_sequence_portion,$or_pos,length($ro),$linestring);

                  $or_pos += $position_buffer;
                  my $complement_ro = reverseComplement($ro);

                  print "$ro found in $seq ($set->{$ro}{'count'}) - deleting read: $ro and complement ($set->{$complement_ro}{'count'}): $complement_ro\n\n" if ($verbose); 
                  if(defined $set->{$ro}){          
                     my $current_reads = $set->{$ro}{'count'};  
                     #print "fwd SET:$current_reads BIN $subset->{$ro}\n";
                     my $current_bases = length($ro) * $current_reads;
                     $integral += $current_reads;
                     $reads_needed += $current_reads;
                     $total_bases += $current_bases;
                     if($paired){
                        $track->{$ro}{'start'} = $or_pos - length($ro) + 1;
                        $track->{$ro}{'end'} = $or_pos;
                        $track->{$ro}{'cov'} = $current_reads;
                        $track->{$ro}{'names'} = $set->{$ro}{'names'};
                        #my $lro=length($ro);
                        #print "OVER FWD\n$seq ($growing_tig_length)\n$temp_sequence\n$consensus\n$ro ($current_reads)\tend=$or_pos ($track->{$ro}{'start'}-$track->{$ro}{'end'})\n\n";
                     }
                     ($bin,$set,$seed) = deleteData($bin,$set,$seed,$ro,$e);
                  }

                  if(defined $set->{$complement_ro}){          
                     my $current_reads = $set->{$complement_ro}{'count'}; 
                     #print "rc SET:$current_reads BIN $subset_rc->{$complement_ro}\n";
                     my $current_bases = length($complement_ro) * $current_reads;
                     $integral += $current_reads;
                     $reads_needed += $current_reads;
                     $total_bases += $current_bases;
                     if($paired){
                        $track->{$complement_ro}{'end'} = $or_pos - length($ro) + 1;
                        $track->{$complement_ro}{'start'} = $or_pos;
                        $track->{$complement_ro}{'cov'} = $current_reads;
                        $track->{$complement_ro}{'names'} = $set->{$complement_ro}{'names'};
                        #my $lro=length($ro);
                        #print "OVER REV\n$seq ($growing_tig_length)\n$temp_sequence\n$consensus\n$complement_ro ($current_reads)\tstart=$or_pos($track->{$complement_ro}{'start'}-$track->{$complement_ro}{'end'})\n\n";
                     }
                     ($bin,$set,$seed) = deleteData($bin,$set,$seed,$complement_ro,$e);
                  }
               }
            }
            #if(! $integral){### no reads are found overlapping with the consensus might be indicative of low complexity regions -- Stop the extension
            if($integral < $tmp_base_overlap){
               print "No overlapping reads agree with the consensus or number of agreeing reads is lower than target coverage (tmp:$tmp_base_overlap  -o $base_overlap). Stopping extension" if ($verbose);
               $extended = 0;
            }else{
               ###Added R.Warren 5/5/2010
               if($temp_sequence_portion =~ /\.(\-*)?$/){### will not extend a contig with 3' consensus bases if no reads overlap them (mitigate assembly errors)
                  $consensus = substr($consensus,0,length($consensus)-length($1));
               }

               ###
               $seq .= $consensus;   ##### contig extension
               print "New Contig is: $seq\n" if ($verbose);
               $extended = 1;
            }
         }else{### no consensus built, will stop the extension
            $extended = 0;
         }

      }###while get the OK for extension

      $trim_ct++;
      if ($trim_ct <= $max_trim){
         last TRIM if (length($seq) <= $MIN_READ_LENGTH); #terminate assembly if trimming becomes too agressive
         $seq = substr($seq, 0, -1);
         $extended = 1;
         print "\n$direction EXTENSION ROUND $trim_ct COMPLETE UNTIL $max_trim nt TRIMMED OFF => TRIMMED SEQUENCE:$seq\n\n" if ($verbose);
      }
      
   }### while trimming within bounds
   #### Adjust the position if tracking paired reads in assembly
   if($paired){
      foreach my $rd (keys %$track){
         $track->{$rd}{'start'} = length($seq) - $track->{$rd}{'start'} + 1;
         $track->{$rd}{'end'} = length($seq) - $track->{$rd}{'end'} + 1;
      }
   }
   
   print "\n*** NOTHING ELSE TO BE DONE IN $direction - PERHAPS YOU COULD DECREASE THE MINIMUM OVERLAP -m (currently set to -m $min_overlap) ***\n\n" if ($verbose);

   return $seq, $set, $bin, $reads_needed, $total_bases, $track;
}


#-----------------------
sub deleteData {
   my ($bin,$set,$seed,$sequence,$e) = @_;
   
   my @o = $sequence =~ /\S{4}/g;
   my $comp_seq = reverseComplement($sequence);
   my @c = $comp_seq =~ /\S{4}/g;

   #remove k-mer from hash table and prefix tree
   delete $bin->{$e->{$o[0]}}{$e->{$o[1]}}{$e->{$o[2]}}{$e->{$o[3]}}{$sequence};
   delete $bin->{$e->{$c[0]}}{$e->{$c[1]}}{$e->{$c[2]}}{$e->{$c[3]}}{$comp_seq};
   delete $set->{$sequence};
   delete $seed->{$sequence};

   return $bin, $set, $seed;
}

#-----------------------
sub reverseComplement{
   $_ = shift;
   $_ = uc();
   tr/ATGC/TACG/;
   return (reverse());
}

#-----------------
sub readFasta{
   my ($matepair,$set,$bin,$file,$short,$paired,$encoded,$seedsplit,$space_restriction,$ignorehead) = @_;

   my $ctrd = 0;
   my $ctline = 0;
 
   my $head = "";
   my $insert_size = 0;

   open(IN,$file) || die "Can't open $file -- fatal\n";
   open (SHO, ">$short") || die "Can't write to $short -- fatal\n";

   print "Sequence reads loaded:\n";
   
   while(<IN>){
      chomp;
      $ctline++;
      if(/^([^\>]*)$/i){
         my $sdna = $1;

         if($paired){
            if($sdna =~/^([ACGT]*)\:([ACGT]*)$/i){### don't crash on Ns, but ignore see below
               my ($rd1,$rd2) = ($1,$2);
               my $head1 = $head . "1";
               my $head2 = $head . "2";

               $matepair->{$rd1}{$rd2}{'is'} = $insert_size;
               $matepair->{$rd2}{$rd1}{'is'} = $insert_size;
               $matepair->{$rd1}{$rd2}{'bt'} = 0;
               $matepair->{$rd2}{$rd1}{'bt'} = 0;

               $head1="" if($ignorehead);
               $head2="" if($ignorehead);
               ($set,$bin,$ctrd) = &loadSequence($set,$bin,$ctrd,$rd1,$encoded,$head1,$insert_size,$seedsplit,$space_restriction);
               ($set,$bin,$ctrd) = &loadSequence($set,$bin,$ctrd,$rd2,$encoded,$head2,$insert_size,$seedsplit,$space_restriction);
            }else{
               my $pairing_failure_message = "Input error at line #$ctline: The sequence \"$sdna\" is not in the right format for paired-end reads  -- Fatal\nMake sure your input is in the form (input sequences can be of variable lengths):\n\n>test\nGCTACGACTATGACATACAGT:GTAGATTGATCGCATGCACGCT\n\nWhere : separates paired reads.  Spaces, <<.>> or any characters other than A,C,G or T in your input file might have caused this error, including reads with Ns.\n";
               print $pairing_failure_message;
               print LOG $pairing_failure_message;
               close LOG;
               exit;
            }
         }else{
            $head="" if($ignorehead);

            ($set,$bin,$ctrd) = &loadSequence($set,$bin,$ctrd,$sdna,$encoded,$head,$insert_size,$seedsplit,$space_restriction) if ($sdna =~ /^([ACGT]*)$/i);
         }
      }elsif(/^\>(\S+)/){
         $head = $1;
         if($paired){
            ($head,$insert_size) = ($1,$2) if($head=~/(\S+)\:(\d+)$/);
         }

         if($head eq "" || ($paired && $insert_size == 0)){
            my $input_failure_message = "Input error at line #$ctline: $_ -- Either you forgot to name the mate pair template (head=$head), improperly formatted the insert size (iz=$insert_size) # or both.  Headers should look like: >templateA:200\nIf you are not using the -p 1 option, you can leave \":insert_size\" out.\n";
            print $input_failure_message;
            print LOG $input_failure_message;
            close LOG;
            exit;
         }
      }
   }
   close IN;
   close SHO;

   my $read_number_message = "\r$ctrd total sequences (" . keys( %$set ) . " unique)\n";
   printf $read_number_message;
   print LOG $read_number_message;

   return $set,$bin,$matepair;
}

#-----------------
### added 31Jan08 R.Warren
sub loadSeed{

   my $file = shift;
   my $seed;  
   my $seedsplit;

   open(IN,$file) || die "Can't open $file -- fatal\n";
   
   my ($subseq,$prev)=('','');

      while(<IN>){
      chomp;

      if (/\>(\S+)/){
         my $head=$1;
         my $subseq_length = length($subseq);
         $MAX=$subseq_length if ($subseq_length > $MAX);
         if($head ne $prev && $subseq ne '' && $subseq_length >= $MIN_READ_LENGTH && $subseq_length >= $min_overlap){
            my $ucsub = uc($subseq);
            $seed->{$ucsub}{'count'}++;
            $seed->{$ucsub}{'names'}{$prev}="";
            $seed->{$ucsub}{'seed_name'}=$prev;
            $seed->{$ucsub}{'ori'}=$subseq;
            for(my $pos==0;$pos<= ($subseq_length-15);$pos++){
               my $word=substr($ucsub,$pos,15);
               my $word_rc = reverseComplement($word);
               $seedsplit->{$word}=1;
               $seedsplit->{$word_rc}=1;
            }

            if($subseq=~/([NX])/i){print "WARNING: the fasta sequence >$prev in your seed file contains characters other than ACGT (i.e. $1) and may prevent proper contig extension.\n";}
         }
         $subseq='';
         $prev=$head;
      }elsif(/^([ACGTNX]*)$/i){
         $subseq .= $_;
      }
   }
   my $subseq_length = length($subseq);
   $MAX=$subseq_length if ($subseq_length > $MAX);
   if($subseq ne '' && $subseq_length >= $MIN_READ_LENGTH && $subseq_length >= $min_overlap){
      my $ucsub = uc($subseq);
      $seed->{$ucsub}{'count'}++;
      $seed->{$ucsub}{'names'}{$prev}="";
      $seed->{$ucsub}{'seed_name'}=$prev;
      $seed->{$ucsub}{'ori'}=$subseq;
      for(my $pos==0;$pos<= ($subseq_length-16);$pos++){
         my $word=substr($ucsub,$pos,16);
         my $word_rc = reverseComplement($word);
         $seedsplit->{$word}=1;
         $seedsplit->{$word_rc}=1;
      }

      if($subseq=~/([NX])/i){print "WARNING: the fasta sequence >$prev in your seed file contains characters other than ACGT (i.e. $1) and may prevent proper contig extension.\n";}
   }
 
   close IN;

   return $seed,$seedsplit;
}

#-----------------
sub loadSequence{

   my ($set,$bin,$ctrd,$seq,$e,$head,$insert_size,$seedsplit,$space_restriction) = @_;

   my $orig=uc($seq);
   my $orig_mer = length($orig);

   if ($orig ne '' && $orig_mer >= $MIN_READ_LENGTH && $orig_mer >= $min_overlap){

      my @f = $orig =~ /\S{4}/g;

      my $rc = reverseComplement($orig);
      my @r = $rc =~ /\S{4}/g;

      my $first_f = substr($orig,0,16);
      my $first_r = substr($rc,0,16); 

      ### added 31Jan08 R.Warren
      $MAX=$orig_mer if ($orig_mer > $MAX);
      my $min_allowed = -1 * ($insert_stdev * $insert_size);
      my ($low_iz, $up_iz) = ($insert_size + $min_allowed, $insert_size - $min_allowed);

      if(! $space_restriction || ($space_restriction && (defined $seedsplit->{$first_f} || defined $seedsplit->{$first_r}))){

         $set->{$orig}{'count'}++;
         $set->{$orig}{'names'}{$head}="";
         $set->{$orig}{'grace'} = $up_iz;
         $bin->{$e->{$f[0]}}{$e->{$f[1]}}{$e->{$f[2]}}{$e->{$f[3]}}{$orig}++;
         $bin->{$e->{$r[0]}}{$e->{$r[1]}}{$e->{$r[2]}}{$e->{$r[3]}}{$rc}++;

         $ctrd++;
         print "\r$ctrd";
         $|++;
      }
   }elsif($orig ne ''){
      if($orig_mer < $MIN_READ_LENGTH){
         print SHO "$seq\tInput sequence shorter than minimum read length allowed ($orig_mer < $MIN_READ_LENGTH nt)\n";
      }elsif($orig_mer < $min_overlap){
         print SHO "$seq\tInput sequence shorter than minimum overlap specified($orig_mer < -m $min_overlap)\n";
      }
   }
  
   $MAX = $MAX_TOP if ($MAX > $MAX_TOP);
   return $set,$bin,$ctrd;
}

#-----------------
sub encodeBases{

   my $encoded;

   my @pos1= ('A','C','G','T');
   my @pos2 = @pos1;
   my @pos3 = @pos1;
   my @pos4 = @pos1;

   my @chararr = ("À","Á","Â","Ã","Ä","Å","Ā","Ą","Ă","Æ","Ç","Ć","Č","Ĉ","Ċ","Ď","Đ","È","É","Ê","Ë","Ē","Ę","Ě","Ĕ","Ė","Ĝ","Ğ","Ġ","Ģ","Ĥ","Ħ","Ì","Í","Î","Ï","Ī","Ĩ","Ĭ","Į","İ","Ĳ","Ĵ","Ķ","Ł","Ľ","Ĺ","Ļ","Ŀ","Ñ","Ń","Ň","Ņ","Ŋ","Ò","Ó","Ô","Õ","Ö","Ø","Ō","Ő","Ŏ","Œ","Ŕ","Ř","Ŗ","Ś","Š","Ş","Ŝ","Ș","Ť","Ţ","Ŧ","Ț","Ù","Ú","Û","Ü","Ū","Ů","Ű","Ŭ","Ũ","Ų","Ŵ","Ý","Ŷ","Ÿ","Ź","Ž","Ż","à","á","â","ã","ä","å","ā","ą","ă","æ","ç","ć","č","ĉ","ċ","ď","đ","è","é","ê","ë","ē","ę","ě","ĕ","ė","ƒ","ĝ","ğ","ġ","ģ","ĥ","ħ","ì","í","î","ï","ī","ĩ","ĭ","į","ı","ĳ","ĵ","ķ","ĸ","ł","ľ","ĺ","ļ","ŀ","ñ","ń","ň","ņ","ŉ","ŋ","ò","ó","ô","õ","ö","ø","ō","ő","ŏ","œ","ŕ","ř","ŗ","ś","š","ş","ŝ","ș","ť","ţ","ŧ","ț","ù","ú","û","ü","ū","ů","ű","ŭ","ũ","ų","ŵ","ý","ÿ","ŷ","ž","ż","ź","a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","t","u","v","w","x","y","z","A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z","1","2","3","4","5","6","7","8","9","0","α","β","χ","δ","ε");
   #print "@chararr**\n";exit;
   my $el=0;
   foreach my $p1 (@pos1){
      foreach my $p2 (@pos2){
         foreach my $p3 (@pos3){
            foreach my $p4 (@pos4){
               my $quad = $p1.$p2.$p3.$p4;
               $encoded->{$quad}=$chararr[$el];
               #print "$quad .. $chararr[$el] .. $el :: $quad $encoded->{$quad}\n";
               $el++;
            } 
         }
      }
   }
   return $encoded;
}

#-----------------
sub forcefillGaps{

   my ($scaffold, $contig, $mergedtigs, $verbose, $min_tig_overlap, $max_count_trim, $npad_gaps, $alternate, $matepair, $min_overlap, $base_overlap, $min_base_ratio, $forcetrack, $max_trim, $encoded) = @_;

   ### STATIC
   my $chunk = 300;
   my $search_distance = 200;

   open(IN,$scaffold) || die "can't read $scaffold -- fatal\n";
   open(OUT,">$mergedtigs") || die "can't write to $mergedtigs -- fatal\n";

   my ($tot,$sct,$ct_merge) = (0,0,0);
   print "Loading mate pairs to force-fill gaps...\n";

   while(<IN>){### each line is a scaffold
      chomp;
      my $sc="";;
      my @a = split(/\,/);
      my @tig;

      if($a[2]=~/\_/){
         @tig = split(/\_/,$a[2]);
      }else{
         push @tig, $a[2];
      }

      $sct++;
      my ($ct,$mct) = (0,0);
      my ($prev,$word,$template) = ("NA","NA","NA");
      my ($seq,$prevseq,$headconcat,$prevEstimatedDistance) = ("","","","");
      my ($prev_miniset,$prev_minibin)=({},{});

      print "$_\n" if($verbose);

      foreach my $t (@tig){### each contig

        if($t=~/([fr])(\d+)z(\d+)(\S+)?/i){
            my $orient = $1;
            my $tnum=$2;
            my $head = $orient . $tnum;
            my $search = "tig" . $tnum;
            my $tlen = $3;
            my $other = $4;
            $tot += $tlen;

            my $estimatedDistance = "";
            $estimatedDistance = $1 if($other=~/m((\-)?\d+)/);
            print "\tSC $a[0] - TIG $ct.  pattern: $t search: $search totalTigSize: $tot Orientation: $orient Gap/Overlap estimated distance: $estimatedDistance\n" if($verbose);

            my $count_trim = 0;

            open(FA,$contig) || die "Can't read $contig -- fatal\n";
            READ:
            while(<FA>){
               chomp;
               if (/\>(\S+)/){
                  my $head=$1;
                  $seq =~ s/[BDEFHIJKLMOPQRSUVWXYZ]/N/g;
                  if ($prev=~/$search\|/i && $prev ne $head && $prev ne "NA"){
                     last READ;
                  }
                  $prev = $head;
                  $seq='';
               }elsif(/^(\S+)$/){
                  $seq.=uc($1);
               }
            }
            close FA;
            $seq = reverseComplement($seq) if($orient eq "f");  ###turn tig around for extension with previous reads
            my $track;
            my ($reads_needed,$total_bases,$tig_count)=(0,0,0);
            ($seq, $prev_miniset, $prev_minibin, $reads_needed, $total_bases, $track) = doExtension("GAP left", $MAX, $seq, $prev_miniset, $prev_minibin, $reads_needed, $total_bases, $min_overlap, $base_overlap, $min_base_ratio, $verbose, $track, $forcetrack, $tig_count, $max_trim, $encoded, $matepair) if($ct);###no need to extend first ($ct=0) contig on the left
            ($prev_miniset, $prev_minibin)=({},{});
            $seq = reverseComplement($seq);  ###re-change orientation for right extension (both fwd/rev tigs)

            #-------------------FORCE-FILL GAPS (Controlled contig extension) BEFORE MERGING
            #Inspect 2 contig edges for reads, collect partner and fill the gap with corresponding mates, while respecting the pairing logic:
            #even if mates have been used before (likely repeats)

            #existing(e) recruited(r)
            # e-->  <--r
            #   e-->  <--r
            #       r-->  <--e
            #  =====    =====   two contigs predicted to be linked
            # recruited (r) can be assembled elsewhere (repeats), but NOT within $search_distance

            my ($miniset,$minibin,$ctrd,$collect_candidates);
            my $head = "NA";
            my $next_tig = $ct+1;
            if(defined $tig[$next_tig]){### a linking contig exists on the right
               #collect reads from left contig
               my $leftCoordList = $alternate->{$tnum};
               if($orient eq "f"){      ### fwd contigs
                  my $instart = $tlen - $search_distance;
                  $instart = 0 if($instart<0);
                  my @window = ($instart .. $tlen); ### extend 5'->3' as-is
                  foreach my $coordstart(@window){
                     my $allrds = $leftCoordList->{$coordstart};
                     foreach my $partner (keys %$allrds){
                        if($coordstart < $allrds->{$partner}){ ### if start<end for a read means :: -->
                           my $partnerlist = $matepair->{$partner};###means mateseq should be <-- in the gap on the right
                           foreach my $mateseq(keys %$partnerlist){
                              print "* $partner partner mates with $mateseq on left fwd tig $tnum ($coordstart) [$instart-$tlen] . collecting $mateseq\n" if($verbose);
                              $collect_candidates->{$mateseq} = $partnerlist->{$mateseq};
                           }
                        }
                     }
                  }
               }else{                  ### rev contig
                  my @window = (0 .. $search_distance); ### right gap is on left of left reverse tig
                  foreach my $coordstart(@window){
                     my $allrds = $leftCoordList->{$coordstart};
                     foreach my $partner (keys %$allrds){
                        if($coordstart > $allrds->{$partner}){ ### reads are inversed:: <--
                           my $partnerlist = $matepair->{$partner};
                           foreach my $mateseq(keys %$partnerlist){
                              print "* $partner partner mates with $mateseq on left rev tig $tnum ($coordstart) [0-$search_distance]. collecting $mateseq\n" if($verbose);
                              $collect_candidates->{$mateseq} = $partnerlist->{$mateseq};
                           }
                        }
                     }
                  }
               }
               if($verbose){
                  my $numcollected = keys(%$collect_candidates);
                  print "$numcollected collected missing mates on gapleft\n";
               }

               #collect reads from right contig
               my ($rr_orient,$rr_tnum,$rr_tlen)=($1,$2,$3) if($tig[$next_tig] =~/([fr])(\d+)z(\d+)(\S+)?/i);### contig to the right
               my $rightCoordList = $alternate->{$rr_tnum};
               if($rr_orient eq "r"){
                  my $instart = $rr_tlen - $search_distance; ### reverse contig on the right, must collect missing mates on the right
                  $instart = 0 if($instart<0);
                  my @window = ($instart .. $rr_tlen);
                  foreach my $coordstart(@window){
                     my $allrds = $rightCoordList->{$coordstart};
                     READ:
                     foreach my $partner (keys %$allrds){
                        if(defined $collect_candidates->{$partner}){### assembled read in right tig was collected from left tig as missing mate..,ust delete
                           print "! previously collected $partner candidate is already assembled on right tig #$rr_tnum ... will delete.\n" if($verbose);
                           delete $collect_candidates->{$partner};
                           next READ;
                        }
                        if($coordstart < $allrds->{$partner}){###  read in this direction :: --->
                           my $partnerlist = $matepair->{$partner};
                           foreach my $mateseq(keys %$partnerlist){### implies missing mate in gap on the left (r tig makes it look left) is <---
                              $collect_candidates->{$mateseq} = $partnerlist->{$mateseq};
                              print "* About to collect $mateseq mating with assembled $partner on right rev tig $rr_tnum ($coordstart) [$instart-$rr_tlen] tig=$ct next=$next_tig string=$tig[$next_tig] orient=$rr_orient len=$rr_tlen\n" if($verbose);
                           }
                        }
                     }
                  }
               }else{###forward contig, gap is on the left
                  my @window = (0 .. $search_distance);
                  foreach my $coordstart(@window){
                     my $allrds = $rightCoordList->{$coordstart};
                     READ:
                     foreach my $partner (keys %$allrds){
                        if(defined $collect_candidates->{$partner}){
                           print "! previously collected $partner candidate is already assembled on right tig #$rr_tnum ... will delete.\n" if($verbose);
                           delete $collect_candidates->{$partner};
                           next READ;
                        }
                        if($coordstart > $allrds->{$partner}){#  read <--- assembled
                           my $partnerlist = $matepair->{$partner};
                           foreach my $mateseq(keys %$partnerlist){#looking for ---> mate
                              $collect_candidates->{$mateseq} = $partnerlist->{$mateseq};
                              print "* About to collect $mateseq mating with assembled $partner on right fwd tig $rr_tnum ($coordstart) [0-$search_distance] tig=$ct next=$next_tig string=$tig[$next_tig] orient=$rr_orient len=$rr_tlen\n" if($verbose);
                           }
                        }
                     }
                  }
               }
               if($verbose){
                  my $numcollected = keys(%$collect_candidates);
                  print "$numcollected collected revised after considering existing tig on right\n";
               }
               my $track;
               my ($reads_needed,$total_bases,$tig_count)=(0,0,0);
               foreach my $mateseq(keys %$collect_candidates){
                  ($miniset,$minibin,$ctrd) = &loadSequence($miniset,$minibin,$ctrd,$mateseq,$encoded,$head,$collect_candidates->{$mateseq}) if(length($mateseq) <= $illuminaLengthCutoff);
               }
               ($seq, $miniset, $minibin, $reads_needed, $total_bases, $track) = doExtension("GAP right", $MAX, $seq, $miniset, $minibin, $reads_needed, $total_bases, $min_overlap, $base_overlap, $min_base_ratio, $verbose, $track, $forcetrack, $tig_count, $max_trim, $encoded, $matepair);### extend first contig on the right
               $prev_miniset = $miniset;###used reads have been deleted
               $prev_minibin = $minibin;
               ($miniset,$minibin)=({},{});### flush custom prefix tree minibin and miniset

            }### another contig in the scaffold is true (means a gap to fill)
            print "\t$prev\n" if($verbose);

            #### CONTIG MERGE CODE ####
            if($word ne "NA"){
               #####
               if(length($seq)<=$chunk){
                  $template = $seq;
               }else{
                  $template = substr($seq,0,$chunk);
               }

               ##### word search
               my $dynamic_word = $word;

               SCAN:
               until($template =~ /$dynamic_word/){
                  
                  $dynamic_word = substr($dynamic_word,1,length($dynamic_word));###5' del 1bp at a time
                  if(length($dynamic_word) < $min_tig_overlap){
                     $count_trim++;
                     last SCAN if($count_trim >= $max_count_trim);
                     $dynamic_word = substr($word,0,length($word)-$count_trim); ###when overlap too short,reset then 3'del
                  }
               }

               if($seq =~ /^\S{0,$max_count_trim}$dynamic_word(.*)/){### will grab the left-most match which is ok
                  my $tail = $1;
                  my $all = "ERROR_";
                  while($prevseq =~ /^(.*)$dynamic_word/ig){
                     $all = $1;
                  }
                  print "$prevseq **** $all **** WORD:$word *** DWord:$dynamic_word *** COUNTTRIM:$count_trim\n" if ($verbose && $all=~/ERROR/);

                  $prevseq = $all . lc($dynamic_word) . $tail;
                  my $overlap = length($dynamic_word);
                  $ct_merge++;
                  print "$ct_merge. GROUNDS FOR MERGING ($overlap nt overlap) !!!\n" if($verbose);
                  $headconcat .= "+" . $prev;

               }else{### no overlaps

                  if(! $npad_gaps){### do not put Ns in gaps
                     print "No MERGE, will print previous sequence and memorize current.\n" if($verbose);
                     my $scsz = length($prevseq);

                     print OUT ">$a[0].$mct|size$scsz $headconcat\n$prevseq\n";
                     $prevseq = $seq;
                     $headconcat = $prev;
                     $mct++;
                  }else{           ### place Ns in gaps, n for predicted but undetected overlaps
                     ### ADDED RLW 5.MAR.2010
                     if($prevEstimatedDistance <= 0){
                        $prevseq .= "n" . $seq
                     }else{
                        $prevseq .= ("N" x $prevEstimatedDistance) . $seq;
                     }
                     $headconcat .= "+" . $prev;
                  }
               }
            }else{
               $prevseq = $seq;
               $headconcat = $prev;
               $mct++;
            }

            ##### For the next search
            if(length($seq)<=$chunk){
               $word = $seq;
            }else{
               $word = substr($seq,length($seq)-$chunk,$chunk); ### this will be the next word to search with
            }
            ###########################

            $prevEstimatedDistance = $estimatedDistance;

         }#tig regex
         $ct++;
      }#each tig
      my $scsz = length($prevseq);
      print OUT ">$a[0].$mct|size$scsz $headconcat\n$prevseq\n";
      $prevseq = '';
   }
   close IN;
   close OUT;
}

## We hope this code is useful to you -- Please send comments & suggestions to rwarren at bcgsc.ca