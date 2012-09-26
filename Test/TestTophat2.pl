#!/usr/bin/perl
use strict;
use warnings;

use CQS::RNASeq;
use CQS::SystemUtils;
use CQS::FileUtils;
use XML::Simple;

my $runNow = get_run_now();

require "Config.pl";

my $section = "tophat2";
create_directory_or_die( $Task::config->{$section}{target_dir} );

my $configfile = $Task::config->{$section}{target_dir} . "/" . $Task::config->{general}{task_name} . ".xml";
XMLout( $Task::config, OutputFile => $configfile, RootName => "RNASeqPipeline" );

tophat2_by_pbs( $Task::config, $section, $runNow );

1;
