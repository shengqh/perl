#!/usr/bin/perl
use strict;
use warnings;

use File::Basename;
use CQS::RNASeq;
use CQS::FileUtils;
use CQS::SystemUtils;
use List::Compare;

sub output_file {
	my ( $data1, $data2, $fileName, $header, @genes ) = @_;
	open OUT, "<$fileName" or die "Cannot create file $fileName";
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

sub compare {
	my ( $targetdir, $dirmap ) = @_;

	my ( $dir1, $dir2 ) = keys %{$dirmap};

	my ( $data1, $header )  = read_cuffdiff_significant_genes("$dir1/gene_exp.diff");
	my ( $data2, $header2 ) = read_cuffdiff_significant_genes("$dir2/gene_exp.diff");

	my @genes1 = keys %{$data1};
	my @genes2 = keys %{$data2};

	my $lc = List::Compare->new( \@genes1, \@genes2 );

	my @common = $lc->get_intersection();
	my @lonly  = $lc->get_Lonly();
	my @ronly  = $lc->get_Ronly();

	my $filename1      = $dirmap->{$dir1};
	my $filename2      = $dirmap->{$dir2};
	my $commonFileName = "$targetdir/${filename1}_${filename2}.diff";
	my $lonlyFileName  = "$targetdir/${filename1}_only.diff";
	my $ronlyFileName  = "$targetdir/${filename2}_only.diff";

	output_file( $data1, $data2, $lonlyFileName,  $header, @lonly );
	output_file( $data1, $data2, $ronlyFileName,  $header, @ronly );
	output_file( $data1, $data2, $commonFileName, $header, @common );
}

my $root   = "/scratch/cqs/shengq1/rnaseq/P2277/cuffdiff/result";
my $dirmap = {
	"${root}/BC2" => "B_TAAS_LAP_BKM",
	"${root}/BC1" => "B_TAAS_LAP",
};

my $comparedir = create_directory_or_die( $root . "/comparison" );

compare( $comparedir, $dirmap );

1;

