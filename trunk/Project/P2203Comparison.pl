#!/usr/bin/perl
use strict;
use warnings;

use CQS::RNASeq;
use CQS::FileUtils;
use CQS::SystemUtils;

my $root;
if ( is_linux() ) {
	$root = "/scratch/cqs/shengq1/rnaseq/P2203/cuffdiff/result";
}
else {
	$root = "D:/projects/P2203/Cuffdiff";
}

my $targetdir = create_directory_or_die( $root . "/comparison" );
my $config    = {
	"RenameDiff" => {
		target_dir => $targetdir,
		root_dir   => $root
	},
	"MLN_vs_None_minus_R_MLN_vs_R" => {
		target_file => "${targetdir}/MLN_vs_None_minus_R_MLN_vs_R.gene_exp.diff.sig",
		operation   => "minus",
		files       => [ "${targetdir}/FLO-1_MLN_vs_FLO-1.gene_exp.diff.sig", "${targetdir}/FLO-1_R_MLN_vs_FLO-1_R.gene_exp.diff.sig" ]
	},
};

copy_and_rename_cuffdiff_file( $config, "RenameDiff" );

compare_cuffdiff( $config, "MLN_vs_None_minus_R_MLN_vs_R" );

1;

