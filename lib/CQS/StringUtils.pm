package CQS::StringUtils;

use strict;
use warnings;

require Exporter;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw( merge_string)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

sub merge_string {
    my ($delimiter, @values) = @_;
    my $first = 1;
    my $result = "";
    foreach my $value (@values){
    	if($first){
    		$result = $value;
    		$first = 0;
    	}
    	else{
    		$result = $result . $delimiter . $value;
    	}
    }
    return ($result);
}

1;
