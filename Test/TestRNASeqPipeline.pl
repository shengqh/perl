#!/usr/bin/perl
use strict;
use warnings;

use CQS::QC;
use CQS::RNASeq;
use CQS::FileUtils;
use CQS::SystemUtils;

require "Config.pl";

tophat2_by_pbs( $Task::config, "tophat2", 0);

cufflinks_by_pbs( $Task::config, "cufflinks", 0 );

cufflinks_by_pbs( $Task::config, "cufflinks2", 0 );

cuffmerge_by_pbs( $Task::config, "cuffmerge", 0 );

cuffmerge_by_pbs( $Task::config, "cuffmerge2", 0 );

cuffdiff_by_pbs( $Task::config, "cuffdiff", 0 );

cuffdiff_by_pbs( $Task::config, "cuffdiff2", 0 );

cuffdiff_by_pbs( $Task::config, "cufflinks_cuffdiff", 0 );

1;
