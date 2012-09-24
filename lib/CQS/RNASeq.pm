#!/usr/bin/perl
package CQS::RNASeq;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::FileUtils;

require Exporter;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(tophat2_by_pbs_batch tophat2_by_pbs_individual cufflinks_by_pbs cuffdiff_by_pbs)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

use Cwd;

sub tophat2_parse_and_check_parameters {
	my ($refParamHash) = @_;

	my %paramHash = %{$refParamHash};

	my $rootDir      = $paramHash{"root_dir"}      or die "define root_dir first";
	my $genomeDb     = $paramHash{"genome_db"}     or die "define genome_db first";
	my $tophat2param = $paramHash{"tophat2_param"} or die "define tophat2_param first";
	my $gtfFile  = $paramHash{"gtf_file"};     #optional parameter
	my $gtfIndex = $paramHash{"gtf_index"};    #optional parameter
	my $pathFile = $paramHash{"path_file"};    #optional parameter
	
	print "root_dir = $rootDir\n";
    print "genome_db = $genomeDb\n";
    print "tophat2_param = $tophat2param\n";
    print "gtf_file = $gtfFile\n";
    print "gtf_index = $gtfIndex\n";
    print "path_file = $pathFile\n";

	if ( defined $gtfFile ) {
		if ( !-e $gtfFile ) {
			die "gtf_file $gtfFile defined but not exists!";
		}

		if ( !defined $gtfIndex ) {
			die "gtf_file was defined but gtf_index was not defined, you should defined gtf_index to cache the parsing result.";
		}
	}
	if ( ( defined $pathFile ) && ( !-s $pathFile ) ) {
		die "path_file $pathFile defined but not exists!";
	}

	return ( $rootDir, $genomeDb, $tophat2param, $gtfFile, $gtfIndex, $pathFile );
}

sub tophat2_by_pbs_batch {
	my ( $refParamHash, $refSampleNames, $refSampleFiles, $refPbsParamHash, $runNow ) = @_;

	my ( $rootDir, $genomeDb, $tophat2param, $gtfFile, $gtfIndex, $pathFile ) = tophat2_parse_and_check_parameters($refParamHash);

	my $taskName = $refParamHash->{"task_name"} or die "task_name is not defined.";

	my @sampleNames = @{$refSampleNames};
	my @sampleFiles = @{$refSampleFiles};

	my $sampleNameCount = scalar(@sampleNames);

	my ($isSingle) = check_is_single( $sampleNameCount, @sampleFiles );

	my ( $logDir, $pbsDir, $resultDir ) = init_dir($rootDir);
	my $tophatDir = create_directory_or_die( $resultDir . "/tophat2" );
	my ($pbsDesc) = get_pbs_desc($refPbsParamHash);

	my $pbsFile = $pbsDir . "/${taskName}_tophat2.pbs";
	my $log     = $logDir . "/${taskName}_tophat2.log";

	output_header( $pbsFile, $pbsDesc, $pathFile, $log );

	for ( my $index = 0 ; $index < $sampleNameCount ; $index++ ) {
		my $sampleName = $sampleNames[$index];
		output_tophat2_script( $genomeDb, $gtfFile, $gtfIndex, $tophat2param, $tophatDir, $sampleName, $index, $isSingle, @sampleFiles );
	}

	output_footer();

	if ($runNow) {
		`qsub $pbsFile`;
		print "$pbsFile submitted\n";
	}
	else {
		print "$pbsFile created\n";
	}
}

sub tophat2_by_pbs_individual {
	my ( $refParamHash, $refSampleNames, $refSampleFiles, $refPbsParamHash, $runNow ) = @_;

	my ( $rootDir, $genomeDb, $tophat2param, $gtfFile, $gtfIndex, $pathFile ) = tophat2_parse_and_check_parameters($refParamHash);

	my @sampleNames     = @{$refSampleNames};
	my @sampleFiles     = @{$refSampleFiles};
	my $sampleNameCount = scalar(@sampleNames);
	my ($isSingle) = check_is_single( $sampleNameCount, @sampleFiles );

	my ( $logDir, $pbsDir, $resultDir ) = init_dir($rootDir);
	my $tophatDir = create_directory_or_die( $resultDir . "/tophat2" );
	my ($pbsDesc) = get_pbs_desc($refPbsParamHash);

	for ( my $index = 0 ; $index < $sampleNameCount ; $index++ ) {
		my $sampleName = $sampleNames[$index];

		my $pbsFile = $pbsDir . "/${sampleName}_tophat2.pbs";
		my $log     = $logDir . "/${sampleName}_tophat2.log";

		output_header( $pbsFile, $pbsDesc, $pathFile, $log );
		output_tophat2_script( $genomeDb, $gtfFile, $gtfIndex, $tophat2param, $tophatDir, $sampleName, 0, $isSingle, @sampleFiles );
		output_footer();

		if ($runNow) {
			`qsub $pbsFile`;
			print "$pbsFile submitted\n";
		}
		else {
			print "$pbsFile created\n";
		}
	}
}

sub cufflinks_by_pbs {
	my ( $cufflinksparam, $rootDir, $refSampleNames, $refSampleFiles, $refPbsParamHash, $runNow ) = @_;

	my @sampleNames     = @{$refSampleNames};
	my @sampleFiles     = @{$refSampleFiles};
	my $sampleNameCount = scalar(@sampleNames);

	my $pathFile = '/home/shengq1/bin/path.txt';

	my ( $logDir, $pbsDir, $resultDir ) = init_dir($rootDir);
	my $cufflinkDir = create_directory_or_die( $resultDir . "/cufflinks" );
	my ($pbsDesc) = get_pbs_desc($refPbsParamHash);

	for ( my $index = 0 ; $index < $sampleNameCount ; $index++ ) {
		my $sampleName = $sampleNames[$index];
		my $sampleFile = $sampleFiles[$index];

		my $pbsFile = $pbsDir . "/${sampleName}_cufflinks.pbs";
		my $log     = $logDir . "/${sampleName}_cufflinks.log";

		output_header( $pbsFile, $pbsDesc, $pathFile, $log );

		my $curDir = create_directory_or_die( $cufflinkDir . "/$sampleName" );

		print OUT "echo cufflinks=`date` \n";
		print OUT "cufflinks $cufflinksparam -o $curDir $sampleFile \n";

		output_footer();

		if ($runNow) {
			`qsub $pbsFile`;
			print "$pbsFile submitted\n";
		}
		else {
			print "$pbsFile created\n";
		}
	}
}

sub cuffdiff_by_pbs {
	my ( $genomeFasta, $gtfFile, $cuffdiffparam, $rootDir, $taskName, $labels, $refFiles, $refPbsParamHash, $runNow ) = @_;

	my @files = @{$refFiles};

	my $pathFile = '/home/shengq1/bin/path.txt';

	my ( $logDir, $pbsDir, $resultDir ) = init_dir($rootDir);
	my $cuffdiffDir = create_directory_or_die( $resultDir . "/cuffdiff" );
	my ($pbsDesc) = get_pbs_desc($refPbsParamHash);

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

	if ($runNow) {
		`qsub $pbsFile`;
		print "$pbsFile submitted\n";
	}
	else {
		print "$pbsFile created\n";
	}
}

sub output_tophat2_script {
	my ( $genomeDb, $gtfFile, $gtfIndex, $tophat2param, $tophatDir, $sampleName, $index, $isSingle, @sampleFiles ) = @_;

    print "genome_db = $genomeDb\n";
    print "tophat2_param = $tophat2param\n";
    print "gtf_file = $gtfFile\n";
    print "gtf_index = $gtfIndex\n";

	my $curDir = create_directory_or_die( $tophatDir . "/$sampleName" );

	print OUT "echo tophat2=`date` \n";

	my $hasGtfFile = ( defined $gtfFile ) && ( -e $gtfFile );
	my $hasIndexFile = defined($gtfIndex) && (-e ($gtfIndex . ".rev.2.bt2"));

	print "hasGtfFile = $hasGtfFile\n";
	print "hasIndexFile = $hasIndexFile\n";

	my $tophat2file = $curDir . "/accepted_hits.bam";

	print OUT "if [ -s $tophat2file ];\n";
	print OUT "then\n";
	print OUT "  echo job has already been done. if you want to do again, delete accepted_hits.bam and submit job again.\n";
	print OUT "else\n";
	if ($hasGtfFile) {
		if ( ( $index == 0 ) && ( !$hasIndexFile ) ) {
			print OUT "  tophat2 $tophat2param -G $gtfFile --transcriptome-index=$gtfIndex -o $curDir $genomeDb ";
		}
		else {
			print OUT "  tophat2 $tophat2param --transcriptome-index=$gtfIndex -o $curDir $genomeDb ";
		}
	}
	elsif ($hasIndexFile) {
		print OUT "  tophat2 $tophat2param --transcriptome-index=$gtfIndex -o $curDir $genomeDb ";
	}
	else {
		print OUT "  tophat2 $tophat2param -o $curDir $genomeDb ";
	}

	if ($isSingle) {
		print OUT "$sampleFiles[$index]\n";
	}
	else {
		print OUT "$sampleFiles[$index*2] $sampleFiles[$index*2+1]\n";
	}
	print OUT "fi\n";
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
	if ( defined $pathFile ) {
		print OUT "source $pathFile\n";
	}
}

sub output_footer() {
	print OUT "echo finished=`date`\n";
	close OUT;
}

1;
