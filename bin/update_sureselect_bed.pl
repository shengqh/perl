$oldbed     = "/scratch/cqs/shengq1/references/sureselect/S0276129_Mouse_All_Exon_V1/S0276129_Regions.bed";
$newbed     = "/scratch/cqs/shengq1/references/sureselect/S0276129_Mouse_All_Exon_V1/S0276129_mm9_All_Exon_V1_M.bed";
$newslimbed = "/scratch/cqs/shengq1/references/sureselect/S0276129_Mouse_All_Exon_V1/S0276129_mm9_All_Exon_V1_slim.bed";

open( ILL,  ">$newbed" );
open( SLIM, ">$newslimbed" );

open( REF, $oldbed );
while (<REF>) {
  if ( $_ =~ m/^chr/ ) {
    my $data = substr $_, 3;
    print ILL $data;

    if ( $_ =~ m/^chr[XYM]/ ) {
      next;
    }
    else {
      print SLIM $data;
    }
  }
  else {
    print ILL $_;
  }
}
close(REF);
close(ILL);
close(SLIM);
