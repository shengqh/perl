#!/usr/bin/perl
use strict;
use warnings;

use CQS::RNASeq;
use Config::Std;

if ( $#ARGV < 0 ) {
	print "Usage:\n";
	print "  tophat2_pbs config_file run_now[y/n]:\n";
	exit;
}

my $config_file = $ARGV[0];

if ( -e $config_file ) {
	my $run_now = 0;
	if ( $#ARGV > 0 ) {
		$run_now = $ARGV[1] eq 'y';
	}

	tophat2_by_pbs_config( $config_file, $run_now );
}
else {
	my $config = {};

	tophat2_create_config($config);

	foreach my $key ( sort keys %{$config} ) {
		print "$key\n";
		my $v = $config->{$key};
		foreach my $key2 ( keys %{$v} ) {
			print "\t$key2=$v->{$key2}\n";
		}
	}
	write_config $config, $config_file;
	print "Template config file has been saved to $config_file\n";
}

