for my $f (glob "*.pl"){
	my $filename = $f;
	
	open (T,"$filename") || die "Can\'t open $filename:$!\n";
	
	my $outfile = $filename =~ s/\.pl$/.txt/r;
	
	my $id = $filename =~ s/\.pl$//r;
	
	open OUT,">./htmldoc/$outfile" or die "cannot output";
	
  	print OUT "<h4><b>USAGE</b></h4> <p style=\"font-size:1.3em;\">qap $id [options]</p><br/>";
    my $flag = 0;
	while (<T>) {
		chomp;
		$_ =~ s/F<//g;
		$_ =~ s/B<//g;
		$_ =~ s/>//g;
		$_ =~ s/<//g; 
		$_ =~ s/qap/qap/r;
		
		if (/qap is still in/){
			#print OUT "  $_\n";
			$flag = 1;
		}
		
		if (/=head1 DESC/){
			print OUT "  <h4><b>DESCRIPTION</b></h4>";
		}
		
		if($flag){
			if (/=head1 OPTIONS/){
				print OUT "  <h4><b>OPTIONS</b></h4>";
				next;
			}
			if($_ eq ''){
				next;
			}
			if(/=over/){
				next;
			}
			if(/\=\item --help/){
				last;
			}
			if(/=\item/){
				$_ =~ s/=\item/ /;
				print OUT "<b>$_</b><br/>";
				next;
			}	
			
			if(/Display this detailed help/){
				next;
			}
			if (/=back/){
				$flag = 0;
			}
			if(/=head1 DESCRIPTION/){
				next;
			}
			print OUT "  $_<br/>";
		}
		
	}
	
	close (T);
	
}
