my $curdepthFile = "H:/shengquanhu/projects/test/C2_1_TBX5_peaks.name.bed.depth.tmp";
my $depthFile    = "H:/shengquanhu/projects/test/C2_1_TBX5_peaks.name.bed.depth";

open( my $tmp,   $curdepthFile )     or die "Cannot open file $curdepthFile";
open( my $depth, "> " . $depthFile ) or die "Cannot open file $depthFile";
my $lastchr  = "";
my $lastpos  = 0;
my $lastfile = "";
my $header   = readline($tmp);
print $depth $header;

my @headers = split /\t/, $header;
my $zeroes = "\t0" x ( scalar(@headers) - 3 );

while ( my $line = <$tmp> ) {
  $line =~ s/\r|\n//g;
  my @parts = split /\t/, $line;

  if ( $lastfile ne $parts[ scalar(@parts) - 1 ] ) {
    $lastpos  = $parts[1];
    $lastfile = $parts[ scalar(@parts) - 1 ];
    print $depth $line, "\n";
    next;
  }

  my $gap = $parts[1] - $lastpos;
  if ( $gap > 2 ) {
    print $depth $parts[0], "\t", ( $lastpos + 1 ), $zeroes, "\t", $lastfile, "\n";
    print $depth $parts[0], "\t", ( $parts[1] - 1 ), $zeroes, "\t", $lastfile, "\n";
  }

  print $depth $line, "\n";
  $lastpos = $parts[1];
}
close $tmp;
close $depth;
