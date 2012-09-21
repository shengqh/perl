#!/usr/bin/perl
package CQS::RNASeq;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::FileUtils;

require Exporter;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(tophat2_by_pbs_batch tophat2_by_pbs_individual cuffdiff_by_pbs)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

use Cwd;

sub tophat2_by_pbs_batch {
	my ( $genomeDb, $gtfFile, $gtfIndex, $tophat2param, $rootDir, $taskName, $refSampleNames, $refSampleFiles ) = @_;
	my @sampleNames = @{$refSampleNames};
	my @sampleFiles = @{$refSampleFiles};

	my $sampleNameCount = scalar(@sampleNames);

	my ($isSingle) = check_is_single( $sampleNameCount, @sampleFiles );

	my $pathFile = '/home/shengq1/bin/path.txt';

	my ( $logDir, $pbsDir, $resultDir ) = init_dir($rootDir);
	my $tophatDir = create_directory_or_die( $resultDir . "/tophat2" );
	my ($pbsDesc) = get_pbs_desc();

	my $pbsFile = $pbsDir . "/${taskName}_tophat2.pbs";
	my $log     = $logDir . "/${taskName}_tophat2.log";

	output_header( $pbsFile, $pbsDesc, $pathFile, $log );

	for ( my $index = 0 ; $index < $sampleNameCount ; $index++ ) {
		my $sampleName = $sampleNames[$index];
		output_tophat2_script( $genomeDb, $gtfFile, $gtfIndex, $tophat2param, $tophatDir, $sampleName, $index, $isSingle, @sampleFiles );
	}

	output_footer();

	print "$pbsFile\n";

	`qsub $pbsFile`;
}

sub tophat2_by_pbs_individual {
	my ( $genomeDb, $gtfFile, $gtfIndex, $tophat2param, $rootDir, $refSampleNames, $refSampleFiles, $runNow ) = @_;
	my @sampleNames = @{$refSampleNames};
	my @sampleFiles = @{$refSampleFiles};

	my $sampleNameCount = scalar(@sampleNames);

	my ($isSingle) = check_is_single( $sampleNameCount, @sampleFiles );

	my $pathFile = '/home/shengq1/bin/path.txt';

	my ( $logDir, $pbsDir, $resultDir ) = init_dir($rootDir);
	my $tophatDir = create_directory_or_die( $resultDir . "/tophat2" );
	my ($pbsDesc) = get_pbs_desc();

	for ( my $index = 0 ; $index < $sampleNameCount ; $index++ ) {
		my $sampleName = $sampleNames[$index];

		my $pbsFile = $pbsDir . "/${sampleName}_tophat2.pbs";
		my $log     = $logDir . "/${sampleName}_tophat2.log";

		output_header( $pbsFile, $pbsDesc, $pathFile, $log );
		output_tophat2_script( $genomeDb, $gtfFile, $gtfIndex, $tophat2param, $tophatDir, $sampleName, $index, $isSingle, @sampleFiles );
		output_footer();

		print "$pbsFile\n";

        if($runNow){
		  `qsub $pbsFile`;
        }
	}
}

sub cuffdiff_by_pbs {
	my ( $genomeFasta, $gtfFile, $cuffdiffparam, $rootDir, $taskName, $labels, @files ) = @_;

	my $pathFile = '/home/shengq1/bin/path.txt';

	my ( $logDir, $pbsDir, $resultDir ) = init_dir($rootDir);
	my $cuffdiffDir = create_directory_or_die( $resultDir . "/cuffdiff" );
	my ($pbsDesc) = get_pbs_desc();

	my $pbsFile = $pbsDir . "/${taskName}_cuffdiff.pbs";
	my $log     = $logDir . "/${taskName}_cuffdiff.log";

	output_header( $pbsFile, $pbsDesc, $pathFile, $log );
	print OUT "cuffdiff $cuffdiffparam -o $cuffdiffDir -L $labels ";

	if ( not( $genomeFasta eq "" ) ) {
		print OUT "-b $genomeFasta ";
	}

	print OUT " $gtfFile ";

	foreach my $file (@files) {
		print OUT "$file ";
	}
	print OUT "\n";

	output_footer();

	print "$pbsFile\n";

	#`qsub $pbsFile`;
}

sub output_tophat2_script {
	my ( $genomeDb, $gtfFile, $gtfIndex, $tophat2param, $tophatDir, $sampleName, $index, $isSingle, @sampleFiles ) = @_;
	my $curDir = $tophatDir . "/$sampleName";

	print OUT "echo tophat2=`date` \n";

	unless ( -e $curDir or mkdir($curDir) ) {
		die "Cannot create directory $curDir\n";
	}

	my $gtfIndexFile = $gtfIndex . ".rev.2.bt2";

	if ( -e $gtfFile ) {
		if ( ( $index == 0 ) && ( not -e $gtfIndexFile ) ) {
			print OUT "tophat2 $tophat2param -G $gtfFile --transcriptome-index=$gtfIndex -o $curDir $genomeDb ";
		}
		else {
			print OUT "tophat2 $tophat2param --transcriptome-index=$gtfIndex -o $curDir $genomeDb ";
		}
	}
	elsif ( -e $gtfIndexFile ) {
		print OUT "tophat2 $tophat2param --transcriptome-index=$gtfIndex -o $curDir $genomeDb ";
	}
	else {
		print OUT "tophat2 $tophat2param -o $curDir $genomeDb ";
	}

	if ($isSingle) {
		print OUT "$sampleFiles[$index]\n";
	}
	else {
		print OUT "$sampleFiles[$index*2] $sampleFiles[$index*2+1]\n";
	}
}

sub check_is_single() {
	my ( $sampleNameCount, @sampleFiles ) = @_;
	my $sampleFileCount = scalar(@sampleFiles);
	my $isSingle        = 1;
	if ( $sampleNameCount == $sampleFileCount ) {
		$isSingle = 1;
	}
	elsif ( $sampleNameCount * 2 == $sampleFileCount ) {
		$isSingle = 0;
	}
	else {
		die "Count of SampleName should be equals to count/half count of SampleFiles";
	}

	return ($isSingle);
}

sub output_header {
	my ( $pbsFile, $pbsDesc, $pathFile, $log ) = @_;
	open( OUT, ">$pbsFile" ) or die $!;
	print OUT $pbsDesc;
	print OUT "#PBS -o $log\n";
	print OUT "#PBS -j oe\n\n";
	print OUT "source $pathFile\n";
}

sub output_footer() {
	print OUT "echo finished=`date`\n";
	close OUT;
}

1;
