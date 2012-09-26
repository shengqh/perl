#!/usr/bin/perl
package CQS::QC;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::Config;

require Exporter;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(fastqc_by_pbs)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

use Cwd;

sub fastqc_by_pbs {
	my ( $config, $section, $runNow ) = @_;

	my ( $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option ) = get_parameter( $config, $section );

	my %rawFiles = %{ $config->{ $config->{$section}{source} } };

	my $shfile = $pbsDir . "/submit.sh";
	open( SH, ">$shfile" ) or die "Cannot create $shfile";

	for my $groupName ( sort keys %rawFiles ) {
		my %sampleMap = %{ $rawFiles{$groupName} };
		for my $sampleName ( sort keys %sampleMap ) {
			my @sampleFiles = @{ $sampleMap{$sampleName} };

			my $pbsName = "${sampleName}_fastqc.pbs";
			my $pbsFile = "${pbsDir}/$pbsName";

			print SH "qsub ./$pbsName \n";

			my $log = "${logDir}/${sampleName}_fastqc.log";

			open( OUT, ">$pbsFile" ) or die $!;
			print OUT $pbsDesc;
			print OUT "#PBS -o $log\n";
			print OUT "#PBS -j oe\n\n";

			if ( -e $path_file ) {
				print OUT "source $path_file\n";
			}
			print OUT "echo fastqc=`date`\n";

			my $sampleCount = scalar(@sampleFiles);

			print OUT "fastqc $option -t $sampleCount -o $resultDir ";
			for my $sampleFile (@sampleFiles) {
				print OUT "$sampleFile ";
			}
			print OUT "\n";
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
	close(SH);
}

1;
