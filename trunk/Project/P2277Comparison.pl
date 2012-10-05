#!/usr/bin/perl
use strict;
use warnings;

use File::Basename;
use File::Copy;
use CQS::RNASeq;
use CQS::FileUtils;
use CQS::SystemUtils;
use List::Compare;

sub output_file {
	my ( $data1, $data2, $fileName, $header, @genes ) = @_;
	open OUT, ">$fileName" or die "Cannot create file $fileName";
	print OUT "$header\n";
	my @sortedgenes = sort @genes;
	for my $gene (@sortedgenes) {
		if ( defined $data1->{$gene} ) {
			print OUT "$data1->{$gene}\n";
		}
		if ( defined $data2->{$gene} ) {
			print OUT "$data2->{$gene}\n";
		}
	}
	close(OUT);
}

sub compare_cuffdiff {
	my ( $config, $section ) = @_;

	my $info = $config->{$section};

	my ( $file1, $file2 ) = @{ $info->{"files"} };

	my ( $data1, $header )  = read_cuffdiff_significant_genes($file1);
	my ( $data2, $header2 ) = read_cuffdiff_significant_genes($file2);

	my @genes1 = keys %{$data1};
	my @genes2 = keys %{$data2};

	my $lc = List::Compare->new( \@genes1, \@genes2 );

	my @resultgenes = ();
	if ( $info->{operation} eq "minus" ) {
		@resultgenes = $lc->get_Lonly();
	}
	elsif ( $info->{operation} eq "intersect" ) {
		@resultgenes = $lc->get_intersection();
	}
	else {
		die "Only minus or intersect is supported.";
	}

	my $resultFileName = $info->{"target_file"};

	output_file( $data1, $data2, $resultFileName, $header, @resultgenes );
}

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

