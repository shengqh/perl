#!/usr/bin/perl
package CQS::QC;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::ConfigUtils;

require Exporter;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(fastqc_by_pbs)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

use Cwd;

sub fastqc_by_pbs {
	my ( $config, $section ) = @_;

	my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option ) = get_parameter( $config, $section );

	my %rawFiles = %{ get_raw_files( $config, $section ) };

	my $shfile = $pbsDir . "/${task_name}.sh";
	open( SH, ">$shfile" ) or die "Cannot create $shfile";

	for my $sampleName ( sort keys %rawFiles ) {
		my @sampleFiles = @{ $rawFiles{$sampleName} };

		my $pbsName = "${sampleName}_fq.pbs";
		my $pbsFile = "${pbsDir}/$pbsName";

		print SH "qsub ./$pbsName \n";

		my $log = "${logDir}/${sampleName}_fq.log";

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

		print "$pbsFile created\n";
	}
	close(SH);
	print "!!!shell file $shfile created, you can run this shell file to submit all fastqc tasks.\n";
}

1;
