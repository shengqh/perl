use strict;
use warnings;

use CQS::HashUtils;

sub GetSampleGeneInfo {
  my ( $uncompressedfile, $keypattern, $ignorepattern ) = @_;

  open( my $readfile, $uncompressedfile )
    or die "Cannot open file $uncompressedfile : $!";

  my $line;
  my @samples;

  while ( !eof($readfile) ) {
    defined( $line = <$readfile> )
      or die "readline failed for $uncompressedfile: $!";

    chomp($line);

    if ( $line =~ /$keypattern/ ) {
      my @parts = split( '\t', $line );
      my $len = @parts;
      @samples = return_unique( @parts[ 1 .. $len - 1 ] );
      last;
    }
  }

  my @genes;
  my $genesymbol = "";
  while ( !eof($readfile) ) {
    defined( $line = <$readfile> )
      or die "readline failed for $uncompressedfile: $!";

    chomp($line);

    if ( defined($ignorepattern) && ( $line =~ /$ignorepattern/ ) ) {
      next;
    }

    if ( !( $line =~ /\S/ ) ) {
      last;
    }

    if ( ( length($line) > 0 ) && ( substr( $line, 0, 1 ) eq "!" ) ) {
      last;
    }

    if ( $line =~ /^(\S+)/ ) {
      push( @genes, $1 );
    }
  }

  my @uniquegenes = return_unique(@genes);
  my $genecount   = @uniquegenes;
  if ( $genecount < 5 ) {
    $genesymbol = join( ';', @uniquegenes );
  }
  else {
    my @names = @uniquegenes[ 0, 1, $genecount - 2, $genecount - 1 ];
    $genesymbol = join( ';', @names );
  }
  return ( 'sample' => \@samples, 'gene' => \@uniquegenes, 'symbol' => $genesymbol );
}

1;
