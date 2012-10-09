#!/usr/bin/perl
use strict;
use warnings;

use File::Basename;
use File::Copy;
use CQS::RNASeq;
use CQS::FileUtils;
use CQS::SystemUtils;

my $root;
if ( is_linux() ) {
	$root = "/scratch/cqs/shengq1/rnaseq/P2277/cuffdiff/result";
}
else {
	$root = "D:/projects/P2277/Cuffdiff";
}

my $targetdir = create_directory_or_die( $root . "/comparison" );
my $config    = {
	"RenameDiff" => {
		target_dir => $targetdir,
		root_dir   => $root
	},
	"B_TAAS_LAP_BKM_minus_TAAS_LAP" => {
		target_file => "${targetdir}/B_TAAS_LAP_BKM_minus_TAAS_LAP.gene_exp.diff",
		operation   => "minus",
		files       => [ "${targetdir}/B_TAAS_LAP_BKM_vs_B_CON.gene_exp.diff", "${targetdir}/B_TAAS_LAP_vs_B_CON.gene_exp.diff" ]
	},
	"B_BKM_only" => {
		target_file => "${targetdir}/B_(TAAS_LAP_BKM_minus_TAAS_LAP)_intersect_(BKM_vs_CON).gene_exp.diff",
		operation   => "intersect",
		files       => [ "${targetdir}/B_TAAS_LAP_BKM_minus_TAAS_LAP.gene_exp.diff", "${targetdir}/B_BKM_vs_B_CON.gene_exp.diff" ]
	},
	"HCC_TAAS_LAP_BKM_minus_TAAS_LAP" => {
		target_file => "${targetdir}/HCC_TAAS_LAP_BKM_minus_TAAS_LAP.gene_exp.diff",
		operation   => "minus",
		files       => [ "${targetdir}/HCC_TAAS_LAP_BKM_vs_HCC_CON.gene_exp.diff", "${targetdir}/HCC_TAAS_LAP_vs_HCC_CON.gene_exp.diff" ]
	},
	"HCC_BKM_only" => {
		target_file => "${targetdir}/HCC_(TAAS_LAP_BKM_minus_TAAS_LAP)_intersect_(BKM_vs_CON).gene_exp.diff",
		operation   => "intersect",
		files       => [ "${targetdir}/HCC_TAAS_LAP_BKM_minus_TAAS_LAP.gene_exp.diff", "${targetdir}/HCC_BKM_vs_HCC_CON.gene_exp.diff" ]
	}
};

copy_and_rename_cuffdiff_file( $config, "RenameDiff" );

compare_cuffdiff( $config, "B_TAAS_LAP_BKM_minus_TAAS_LAP" );

compare_cuffdiff( $config, "B_BKM_only" );

compare_cuffdiff( $config, "HCC_TAAS_LAP_BKM_minus_TAAS_LAP" );

compare_cuffdiff( $config, "HCC_BKM_only" );

1;

