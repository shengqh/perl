#!/usr/bin/perl
package CQS::RNASeq;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;

require Exporter;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(tophat2_by_pbs)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

use Cwd;

sub tophat2_by_pbs {
	my ( $genomeDb, $gtfFile, $gtfIndex, $tophat2param, $rootDir, $taskName, $refSampleNames, $refSampleFiles ) = @_;
	my @sampleNames = @{$refSampleNames};

	my @sampleFiles = @{$refSampleFiles};

	my $sampleNameCount = scalar(@sampleNames);
	my $sampleFileCount = scalar(@sampleFiles);

	my $isSingle = 1;
	if ( $sampleNameCount == $sampleFileCount ) {
		$isSingle = 1;
	}
	elsif ( $sampleNameCount * 2 == $sampleFileCount ) {
		$isSingle = 0;
	}
	else {
		die "Count of SampleName should be equals to count/half count of SampleFiles";
	}

	my $pathFile = '/home/shengq1/bin/path.txt';

	my ( $logDir, $pbsDir, $resultDir ) = init_dir($rootDir);

	my $tophatDir = $resultDir . "/tophat2";

	unless ( -e $tophatDir or mkdir($tophatDir) ) {
		die "Cannot create directory $tophatDir\n";
	}

	my ($pbsDesc) = get_pbs_desc();

	my $pbsFile = $pbsDir . "/${taskName}_tophat2.pbs";
	my $log     = $logDir . "/${taskName}_tophat2.log";

	open( OUT, ">$pbsFile" ) or die $!;
	print OUT $pbsDesc;
	print OUT "#PBS -o $log\n";
	print OUT "#PBS -j oe\n\n";
	print OUT "source $pathFile\n";
	print OUT "echo tophat2=`date` \n";

	my $gtfIndexFa = $gtfIndex . ".fa";
	for ( my $i = 0 ; $i < $sampleNameCount ; $i++ ) {
		my $sampleName = $sampleNames[$i];
		my $curDir     = $tophatDir . "/$sampleName";

		unless ( -e $curDir or mkdir($curDir) ) {
			die "Cannot create directory $curDir\n";
		}

		if ( ( $i == 0 ) && ( not -e $gtfIndexFa ) ) {
			print OUT "tophat2 $tophat2param -G $gtfFile --transcriptome-index=$gtfIndex -o $curDir $genomeDb ";
		}
		else {
			print OUT "tophat2 $tophat2param --transcriptome-index=$gtfIndex -o $curDir $genomeDb ";
		}
		
		if ($isSingle) {
			print OUT "$sampleFiles[$i]\n";
		}
		else {
			print OUT "$sampleFiles[$i*2] $sampleFiles[$i*2+1]\n";
		}

	}
	print OUT "echo finished=`date`\n";
	close OUT;

    print "$pbsFile\n";
	`qsub $pbsFile`;
}

1;
