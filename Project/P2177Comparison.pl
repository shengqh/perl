#!/usr/bin/perl
use strict;
use warnings;

use CQS::RNASeq;
use CQS::FileUtils;
use CQS::SystemUtils;

my @roots =
  ( "/scratch/cqs/shengq1/rnaseq/P2177_2/cufflinks_cuffdiff/result", "/scratch/cqs/shengq1/rnaseq/P2177_2/NG_cufflinks_cuffdiff/result", "/scratch/cqs/shengq1/rnaseq/P2177_2/NG_DEFAULT_cufflinks_cuffdiff/result" );

foreach my $root (@roots) {
	my $targetdir = create_directory_or_die( $root . "/comparison" );
	my $config    = {
		"RenameDiff" => {
			target_dir => $targetdir,
			root_dir   => $root
		},
	};

	copy_and_rename_cuffdiff_file( $config, "RenameDiff" );
}

1;

