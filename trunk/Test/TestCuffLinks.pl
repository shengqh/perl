#!/usr/bin/perl
use strict;
use warnings;

use CQS::RNASeq;
use CQS::SystemUtils;
my $runNow = get_run_now();

require "Config.pl";

cufflinks_by_pbs( $Task::config, "cufflinks", $runNow );

cufflinks_by_pbs( $Task::config, "cufflinks2", $runNow );

1;
