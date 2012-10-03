#!/usr/bin/perl
use strict;
use warnings;

use CQS::RNASeq;
use CQS::FileUtils;

#my $cufflinksdir = "/scratch/cqs/shengq1/rnaseq/P2277/cufflinks/result";
my $cufflinksdir = "d:/tmp";

my @subdirs = list_directories($cufflinksdir);

@subdirs = sort @subdirs;

my $alldata = {};
my @genes   = ();
foreach my $subdir (@subdirs) {
	my $file = "${cufflinksdir}/${subdir}/genes.fpkm_tracking";
	my $data = read_cufflinks_fpkm($file);

	my @keys = keys %{$data};
	foreach my $gene (@keys) {
		if ( $gene =~ /^CUFF/ ) {
			delete( $data->{$gene} );
		}
	}
	$alldata->{$subdir} = $data;

	if ( scalar(@genes) == 0 ) {
		@genes = keys %{$data};
	}
	else {
		@genes = grep { exists $data->{$_} } @genes;
	}
}

my $resultfile = $cufflinksdir . "/fpkm.tsv";
open OUT, ">$resultfile" or die "Cannot create file $resultfile";

print OUT "gene";
foreach my $subdir (@subdirs) {
	print OUT "\t$subdir";
}
print OUT "\n";

@genes = sort @genes;

foreach my $gene (@genes) {
	print OUT $gene;
	foreach my $subdir (@subdirs) {
		print OUT "\t" . $alldata->{$subdir}->{$gene};
	}
	print OUT "\n";
}

close(OUT);

1;

