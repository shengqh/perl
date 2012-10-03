#!/usr/bin/perl
use strict;
use warnings;

use CQS::RNASeq;
use CQS::FileUtils;
use CQS::SystemUtils;

my $groups = {
	"B_CON"          => [ "P2277-12", "P2277-17", "P2277-22" ],
	"B_TAAS_LAP"     => ["P2277-15"],
	"B_TAAS_LAP_BKM" => ["P2277-18"],
	"B_LAP_BKM"  => [ "P2277-21", "P2277-24", ],
	"B_BKM"      => ["P2277-27"],
	"B_TAAS_BKM" => ["P2277-28"],
	"HCC_CON"          => [ "P2277-11", "P2277-19", "P2277-30", ],
	"HCC_LAP"          => ["P2277-09"],
	"HCC_TAAS_LAP"     => [ "P2277-10", "P2277-14", "P2277-20", ],
	"HCC_TAAS_LAP_BKM" => [ "P2277-13", "P2277-16", "P2277-29", ],
	"HCC_BKM"          => [ "P2277-23", "P2277-25", "P2277-26", ],
};

my $map = {};
foreach my $groupName ( sort keys %{$groups} ) {
	my @samples = @{ $groups->{$groupName} };
	foreach my $sample (@samples) {
		$map->{$sample} = $groupName . "_" . $sample;
	}
}

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
	if ( !defined $map->{$subdir} ) {
		next;
	}

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
	my ( $resultfile, $subdirs_ref, $genes_ref, $alldata, $map ) = @_;

	my @subdirs = @{$subdirs_ref};
	my @genes   = sort @{$genes_ref};

	print "Saving $resultfile ...\n";
	open OUT, ">$resultfile" or die "Cannot create file $resultfile";

	print OUT "gene";
	foreach my $subdir (@subdirs) {
		print OUT "\t" . $map->{$subdir};
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

save( $resultfile, \@subdirs, \@genes, $alldata, $map );

my @counts = ( 2, 3, 4, 5 );
foreach my $i (@counts) {
	my @filteredgenes = filter( \@subdirs, \@genes, $alldata, $i );
	my $file = $cufflinksdir . "/fpkm_${i}.tsv";
	save( $file, \@subdirs, \@filteredgenes, $alldata, $map );
}

print "Done\n";

1;

