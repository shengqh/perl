package CQS::SystemUtils;

use strict;
use warnings;

require Exporter;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw( get_run_now)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';
sub get_run_now {
	my $runNow = 0;
	if ( $#ARGV >= 0 ) {
		my $isRunNow = $ARGV[0];
		return ( $isRunNow eq "y" );
	}
	else {
		return (false);
	}
}

1;
