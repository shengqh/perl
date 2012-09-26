#!/usr/bin/perl
use strict;
use warnings;

use CQS::RNASeq;
use CQS::SystemUtils;
my $runNow = get_run_now();

require "Config.pl";

my $section = "tophat2";

tophat2_by_pbs( $Task::config, $section, $runNow );

1;
