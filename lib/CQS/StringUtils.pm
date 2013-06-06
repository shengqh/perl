package CQS::StringUtils;

use strict;
use warnings;

require Exporter;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw( merge_string print_hash)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

sub merge_string {
  my ( $delimiter, @values ) = @_;
  my $first  = 1;
  my $result = "";
  foreach my $value (@values) {
    if ($first) {
      $result = $value;
      $first  = 0;
    }
    else {
      $result = $result . $delimiter . $value;
    }
  }
  return ($result);
}

sub print_hash {
  my $hash = shift;
  foreach my $k ( sort keys %{$hash} ) {
    print "$k => @{$hash->{$k}}\n";
  }
}

1;
