#!/usr/bin/perl
package CQS::CNV;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;

require Exporter;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(cnvnator conifer)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

use Cwd;

sub cnvnator {
	my ( $config, $section ) = @_;

	my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option ) = get_parameter( $config, $section );
	my $binsize = $config->{$section}{binsize} or die "define ${section}::binsize first";

	my %rawFiles = %{ get_raw_files( $config, $section ) };

	my $shfile = $pbsDir . "/${task_name}.sh";
	open( SH, ">$shfile" ) or die "Cannot create $shfile";
	print SH "type -P qsub &>/dev/null && export MYCMD=\"qsub\" || export MYCMD=\"bash\" \n";

	for my $sampleName ( sort keys %rawFiles ) {
		my @sampleFiles = @{ $rawFiles{$sampleName} };

		my $bamFile = $sampleFiles[0];

		my $pbsName = "${sampleName}_cnvnator.pbs";
		my $pbsFile = "${pbsDir}/$pbsName";

		print SH "\$MYCMD ./$pbsName \n";

		my $log = "${logDir}/${sampleName}_cnvnator.log";

		open( OUT, ">$pbsFile" ) or die $!;
		print OUT $pbsDesc;
		print OUT "#PBS -o $log\n";
		print OUT "#PBS -j oe\n\n";

		if ( -e $path_file ) {
			print OUT "source $path_file\n";
		}

		my $curDir   = create_directory_or_die( $resultDir . "/$sampleName" );
		my $rootFile = $sampleName . ".root";
		my $callFile = $sampleName . ".call";

		print OUT "cd $curDir\n\n";

		print OUT "if [ -s $callFile ]; then\n";
		print OUT "  echo job has already been done. if you want to do again, delete $callFile and submit job again.\n";
		print OUT "else\n";
		print OUT "  if [ ! -s $rootFile ]; then\n";
		print OUT "    echo extract=`date`\n";
		print OUT "    cnvnator -root $rootFile -tree $bamFile \n";
		print OUT "  fi\n";
		print OUT "  echo call=`date`\n";
		print OUT "  cnvnator -root $rootFile -his $binsize \n";
		print OUT "  cnvnator -root $rootFile -stat $binsize \n";
		print OUT "  cnvnator -root $rootFile -partition $binsize \n";
		print OUT "  cnvnator -root $rootFile -call $binsize > $callFile \n";
		print OUT "fi\n\n";

		print OUT "echo finished=`date`\n";
		close OUT;

		print "$pbsFile created\n";
	}
	close(SH);

	if ( is_linux() ) {
		chmod 0755, $shfile;
	}

	print "!!!shell file $shfile created, you can run this shell file to submit all cnvnator tasks.\n";

	#`qsub $pbsFile`;
}

sub conifer {
	my ( $config, $section ) = @_;

	my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option ) = get_parameter( $config, $section );
	my $conifer   = $config->{$section}{conifer} or die "define conifer program location first.\nconifer => \"location\"";
	my $probefile = $config->{$section}{probefile};
	my $probedef  = "";
	if ( defined $probefile ) {
		$probedef = "--probes $probefile";
	}

	my $bampattern = $config->{$section}{sorted_bam_replace_pattern};

	my %rawFiles = %{ get_raw_files( $config, $section ) };

	my $shfile = $pbsDir . "/${task_name}_rpkm.sh";
	open( SH, ">$shfile" ) or die "Cannot create $shfile";
	print SH "type -P qsub &>/dev/null && export MYCMD=\"qsub\" || export MYCMD=\"bash\" \n";

	create_directory_or_die( $resultDir . "/rpkm" );
	for my $sampleName ( sort keys %rawFiles ) {
		my @sampleFiles = @{ $rawFiles{$sampleName} };

		my $pbsName = "${sampleName}_rpkm.pbs";
		my $pbsFile = "${pbsDir}/$pbsName";

		print SH "\$MYCMD ./$pbsName \n";

		my $log = "${logDir}/${sampleName}_rpkm.log";

		open( OUT, ">$pbsFile" ) or die $!;
		print OUT $pbsDesc;
		print OUT "#PBS -o $log\n";
		print OUT "#PBS -j oe\n\n";

		if ( -e $path_file ) {
			print OUT "source $path_file\n";
		}
		print OUT "cd $resultDir\n\n";
		print OUT "echo rpkm=`date`\n";

		my $bamFile = $sampleFiles[0];
		if ( defined $bampattern ) {
			$bamFile =~ $bampattern;
		}

		print $bamFile . "\n";

		my $rpkm = "rpkm/" . $sampleName . ".rpkm";

		print OUT "if [ ! -s $rpkm ]; then\n";
		print OUT "  echo conifer=`date`\n";
		print OUT "  python $conifer rpkm $probedef --input $bamFile --output $rpkm \n";
		print OUT "fi\n";

		print OUT "echo finished=`date`\n";
		close OUT;

		print "$pbsFile created\n";
	}
	close(SH);

	if ( is_linux() ) {
		chmod 0755, $shfile;
	}

	my $pbsName   = "${task_name}_after_rpkm.pbs";
	my $pbsFile   = "${pbsDir}/$pbsName";
	my $hdf5File  = "${task_name}.hdf5";
	my $svalsFile = "${task_name}.svals";
	my $callFile  = "${task_name}.call";

	open( OUT, ">$pbsFile" ) or die $!;
	print OUT $pbsDesc;
	my $log = "${logDir}/${task_name}_after_rpkm.log";
	print OUT "#PBS -o $log\n";
	print OUT "#PBS -j oe\n\n";
	if ( -e $path_file ) {
		print OUT "source $path_file\n";
	}

	create_directory_or_die( $resultDir . "/call_images" );

	print OUT "cd $resultDir\n\n";

	print OUT "\n";
	print OUT "#2 analysis\n";
	print OUT "echo analyze=`date`\n";
	print OUT "python $conifer analyze $probedef --rpkm_dir rpkm/ --output $hdf5File --svd 6 --write_svals $svalsFile \n";
	print OUT "\n";
	print OUT "#3 call\n";
	print OUT "echo call=`date`\n";
	print OUT "python $conifer call --input $hdf5File --output $callFile \n";
	print OUT "\n";
	print OUT "#4 plot\n";
	print OUT "python $conifer plotcalls --input $hdf5File --calls $callFile --output call_images \n";
	close OUT;

	print "$pbsFile created\n";

	print "!!!shell file $shfile created, you can run this shell file to submit all conifer rpkm tasks.\n";

	#`qsub $pbsFile`;
}

1;
