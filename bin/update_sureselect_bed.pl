$oldbed = "/scratch/cqs/shengq1/references/sureselect/S0276129_Mouse_All_Exon_V1/S0276129_Regions.bed";
$newbed = "/scratch/cqs/shengq1/references/sureselect/S0276129_Mouse_All_Exon_V1/S0276129_Mouse_All_Exon_V1.bed";

open( ILL, ">$newbed" );

open( REF, $oldbed );
while (<REF>) {
  if ( $_ =~ m/^chrM/ ) {
    my $data = substr $_, 4;
    print ILL "MT$data";
  }
  elsif ( $_ =~ m/^chr/ ) {
    my $data = substr $_, 3;
    print ILL $data;
  }
  else {
    print ILL $_;
  }
}
close(REF);
close(ILL);
