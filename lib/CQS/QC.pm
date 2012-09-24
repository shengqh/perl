#!/usr/bin/perl
package CQS::QC;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::SystemUtils;
use CQS::FileUtils;

require Exporter;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(fastqc_by_pbs)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

use Cwd;

sub fastqc_by_pbs {
	my ( $rootDir, $refSampleNames, $refSampleFiles, $pbsParamHashRef, $runNow ) = @_;

	my @sampleNames = @{$refSampleNames};
	my @sampleFiles = @{$refSampleFiles};

	my $sampleNameCount = scalar(@sampleNames);

	my $pathFile = '/home/shengq1/bin/path.txt';

	my ( $logDir, $pbsDir, $resultDir ) = init_dir($rootDir);

	my ($pbsDesc) = get_pbs_desc($pbsParamHashRef);

	my $fastqcDir = create_directory_or_die( $resultDir . "/fastqc" );
	for ( my $index = 0 ; $index < $sampleNameCount ; $index++ ) {
		my $sampleName = $sampleNames[$index];

		my $pbsFile = $pbsDir . "/${sampleName}_fastqc.pbs";
		my $log     = $logDir . "/${sampleName}_fastqc.log";

		open( OUT, ">$pbsFile" ) or die $!;
		print OUT $pbsDesc;
		print OUT "#PBS -o $log\n";
		print OUT "#PBS -j oe\n\n";
		print OUT "source $pathFile\n";
		print OUT "echo fastqc=`date`\n";
		print OUT "fastqc -o $fastqcDir $sampleFiles[$index]\n";
		print OUT "echo finished=`date`\n";
		close OUT;

		if ($runNow) {
			`qsub $pbsFile`;
			print "$pbsFile submitted\n";
		}
		else {
			print "$pbsFile created\n";
		}
	}
}

1;
