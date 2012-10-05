#!/usr/bin/perl
use strict;
use warnings;

use CQS::RNASeq;
use CQS::FileUtils;
use CQS::SystemUtils;
use Statistics::Descriptive;

my $groups = {
	#"OTHER"          => [ "P2277-01", "P2277-02", "P2277-03", "P2277-04", "P2277-05", "P2277-06", "P2277-07", "P2277-08", "P2277-31" ],
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
		my $number = substr $sample, -2;
		$map->{$sample} = $groupName . "_" . $number;
	}
}

my $cufflinksdir;
if ( is_linux() ) {
	$cufflinksdir = "/scratch/cqs/shengq1/rnaseq/P2277/cufflinks/result";
}
else {
	$cufflinksdir = "D:/projects/P2277/cufflinks";
}

my @allsubdirs = list_directories($cufflinksdir);

my @subdirs = ();
foreach my $subdir (@allsubdirs) {
	if ( !defined $map->{$subdir} ) {
		next;
	}
	push( @subdirs, $subdir );
}

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
	my ( $resultfile, $subdirs_ref, $genes_ref, $alldata, $map ) = @_;

	my @subdirs   = @{$subdirs_ref};
	my @genes     = sort @{$genes_ref};
	my $genecount = scalar(@genes);

	print "Saving $resultfile ($genecount) ...\n";
	open OUT, ">$resultfile" or die "Cannot create file $resultfile";

	print OUT "gene";
	foreach my $subdir (@subdirs) {
		my $name = $map->{$subdir};
		print OUT "\t$name";
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
	my ( $subdirs_ref, $genes_ref, $alldata, $mincount, $top ) = @_;

	my @subdirs = @{$subdirs_ref};
	my @genes   = sort @{$genes_ref};

	my $allcount = scalar(@genes);

	my %sdmap = ();
	foreach my $gene (@genes) {
		my $count = 0;
		my $stat  = Statistics::Descriptive::Full->new();
		foreach my $subdir (@subdirs) {
			if ( $alldata->{$subdir}->{$gene} >= 1 ) {
				$count++;
			}
			$stat->add_data( $alldata->{$subdir}->{$gene} );
		}
		if ( $count >= $mincount ) {
			$sdmap{$gene} = $stat->standard_deviation();
		}
	}

	my @sortedkeys = ( sort { $sdmap{$b} <=> $sdmap{$a} } keys %sdmap );
	my $number = $top * scalar(@sortedkeys);

	my @result = ();
	for ( my $i = 0 ; $i < $number ; $i++ ) {
		push( @result, $sortedkeys[$i] );
	}

	return (@result);
}

my $resultfile = $cufflinksdir . "/fpkm.txt";

save( $resultfile, \@subdirs, \@genes, $alldata, $map );

my @tops = ( 0.05, 0.1, 0.2 );
foreach my $top (@tops) {
	my @filteredgenes = filter( \@subdirs, \@genes, $alldata, 1, $top );
	my $file = $cufflinksdir . "/fpkm_${top}.txt";
	save( $file, \@subdirs, \@filteredgenes, $alldata, $map );
}

print "Done\n";

1;

