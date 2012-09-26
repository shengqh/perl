#!/usr/bin/perl
use strict;
use warnings;

use CQS::RNASeq;
use CQS::SystemUtils;

my $runNow = get_run_now();

require "./TestConfig.pl";

tophat2_by_pbs( $main::config, $runNow );

1;
