#!/usr/bin/perl
use strict;
use warnings;

use CQS::RNASeq;
use CQS::FileUtils;
use CQS::SystemUtils;

my $cufflinksdir;
if ( is_linux() ) {
	$cufflinksdir = "/scratch/cqs/shengq1/rnaseq/P2277/cufflinks/result";
}
else {
	$cufflinksdir = "d:/tmp";
}

my @subdirs = list_directories($cufflinksdir);

@subdirs = sort @subdirs;

my $alldata = {};
my @genes   = ();
foreach my $subdir (@subdirs) {
	my $file = "${cufflinksdir}/${subdir}/genes.fpkm_tracking";
	print "Reading $file ...\n";

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

sub save {
	my ( $resultfile, $subdirs_ref, $genes_ref, $alldata ) = @_;

	my @subdirs = @{$subdirs_ref};
	my @genes   = sort @{$genes_ref};

	print "Saving $resultfile ...\n";
	open OUT, ">$resultfile" or die "Cannot create file $resultfile";

	print OUT "gene";
	foreach my $subdir (@subdirs) {
		print OUT "\t$subdir";
	}
	print OUT "\n";

	foreach my $gene (@genes) {
		print OUT $gene;
		foreach my $subdir (@subdirs) {
			print OUT "\t" . $alldata->{$subdir}->{$gene};
		}
		print OUT "\n";
	}

	close(OUT);
}

sub filter {
	my ( $subdirs_ref, $genes_ref, $alldata, $mincount ) = @_;

	my @subdirs = @{$subdirs_ref};
	my @genes   = sort @{$genes_ref};

	my $allcount = scalar(@subdirs);

	my @result = ();
	foreach my $gene (@genes) {
		my $count = 0;
		foreach my $subdir (@subdirs) {
			if ( $alldata->{$subdir}->{$gene} > 0 ) {
				$count++;
			}
		}
		if ( $count >= $mincount ) {
			push( @result, $gene );
		}
	}

	return (@result);
}

my $resultfile = $cufflinksdir . "/fpkm.tsv";

save( $resultfile, \@subdirs, \@genes, $alldata );

my @counts = ( 2, 3, 4, 5 );
foreach my $i (@counts) {
	my @filteredgenes = filter( \@subdirs, \@genes, $alldata, $i );
	my $file = $cufflinksdir . "/fpkm_${i}.tsv";
	save( $resultfile, \@subdirs, \@filteredgenes, $alldata );
}

print "Done\n";

1;

