#!/usr/bin/perl
package CQS::RNASeq;

use strict;
use warnings;
use File::Basename;
use Config::Std;
use CQS::PBS;
use CQS::FileUtils;
use CQS::StringUtils;

require Exporter;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = (
	'all' => [
		qw(tophat2_by_pbs_batch tophat2_by_pbs_individual tophat2_by_pbs_config tophat2_create_config cufflinks_by_pbs cuffdiff_by_pbs)
	]
);

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

use Cwd;

sub tophat2_parse_and_check_parameters {
	my ($refParamHash) = @_;
	my %paramHash = %{$refParamHash};

	my $root_dir  = $paramHash{"root_dir"}  or die "define root_dir first";
	my $genome_db = $paramHash{"genome_db"} or die "define genome_db first";
	my $tophat2_param = $paramHash{"tophat2_param"}
	  or die "define tophat2_param first";
	my $gtf_file  = $paramHash{"gtf_file"};     #optional parameter
	my $gtfIndex  = $paramHash{"gtf_index"};    #optional parameter
	my $path_file = $paramHash{"path_file"};    #optional parameter

	#print "root_dir = $root_dir\n";
	#print "genome_db = $genome_db\n";
	#print "tophat2_param = $tophat2_param\n";
	#print "gtf_file = $gtf_file\n";
	#print "gtf_index = $gtfIndex\n";
	#print "path_file = $path_file\n";

	if ( defined $gtf_file ) {
		if ( !-e $gtf_file ) {
			die "gtf_file $gtf_file defined but not exists!";
		}

		if ( !defined $gtfIndex ) {
			die
"gtf_file was defined but gtf_index was not defined, you should defined gtf_index to cache the parsing result.";
		}
	}
	if ( defined $path_file ) {
		if ( !-e $path_file ) {
			die "path_file $path_file defined but not exists!";
		}
	}

	return (
		$root_dir, $genome_db, $tophat2_param,
		$gtf_file, $gtfIndex,  $path_file
	);
}

sub tophat2_by_pbs_batch {
	my (
		$refParamHash,    $refSampleNames, $refSampleFiles,
		$refPbsParamHash, $runNow
	) = @_;

	my (
		$root_dir, $genome_db, $tophat2_param,
		$gtf_file, $gtfIndex,  $path_file
	) = tophat2_parse_and_check_parameters($refParamHash);

	my $taskName = $refParamHash->{"task_name"}
	  or die "task_name is not defined.";

	my @sampleNames = @{$refSampleNames};
	my @sampleFiles = @{$refSampleFiles};

	my $sampleNameCount = scalar(@sampleNames);

	my ($isSingle) = check_is_single( $sampleNameCount, @sampleFiles );

	my ( $logDir, $pbsDir, $resultDir ) = init_dir($root_dir);
	my $tophatDir = create_directory_or_die( $resultDir . "/tophat2" );
	my ($pbsDesc) = get_pbs_desc($refPbsParamHash);

	my $pbsFile = $pbsDir . "/${taskName}_tophat2.pbs";
	my $log     = $logDir . "/${taskName}_tophat2.log";

	output_header( $pbsFile, $pbsDesc, $path_file, $log );

	for ( my $index = 0 ; $index < $sampleNameCount ; $index++ ) {
		my $sampleName = $sampleNames[$index];
		output_tophat2_script(
			$genome_db,     $gtf_file,  $gtfIndex,
			$tophat2_param, $tophatDir, $sampleName,
			$index,         $isSingle,  @sampleFiles
		);
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

sub check_config_bool_value {
	my ( $config, $key, $value ) = @_;
	if ( !defined $config->{$key} ) {
		return (0);
	}

	return ( $config->{$key} eq $value );
}

sub tophat2_check_config {
	my ($config) = @_;

	my $root_dir = $config->{general}{root_dir}
	  or die "define general::root_dir first";
	my $genome_db = $config->{general}{bowtie2_index}
	  or die "define general::bowtie2_index first";
	my $gtf_file = $config->{general}{transcript_gtf};    #optional parameter
	my $gtf_index =
	  $config->{general}{transcript_gtf_index};           #optional parameter
	my $path_file     = $config->{general}{path_file};    #optional parameter
	my $tophat2_param = $config->{option}{tophat2}
	  or die "define option::tophat2 first";

	if ( defined $gtf_file ) {
		if ( !-e $gtf_file ) {
			die "transcript_gtf $gtf_file defined but not exists!";
		}

		if ( !defined $gtf_index ) {
			die
"transcript_gtf was defined but transcript_gtf_index was not defined, you should defined transcript_gtf_index to cache the parsing result.";
		}
	}
	if ( defined $path_file ) {
		if ( !$path_file eq "" ) {
			if ( !-e $path_file ) {
				die "path_file $path_file defined but not exists!";
			}
		}
	}

	return (
		$root_dir, $genome_db, $tophat2_param,
		$gtf_file, $gtf_index, $path_file
	);
}

sub tophat2_create_config {
	my ($config) = @_;
	if ( !$config->{general}{root_dir} ) {
		$config->{general}{root_dir} = "root_dir";
	}
	if ( !$config->{general}{bowtie2_index} ) {
		$config->{general}{bowtie2_index} = "bowtie2_index";
	}
	if ( !$config->{general}{transcript_gtf} ) {
		$config->{general}{transcript_gtf} = "transcript_gtf";
	}
	if ( !$config->{general}{transcript_gtf_index} ) {
		$config->{general}{transcript_gtf_index} = "transcript_gtf_index";
	}
	if ( !$config->{general}{path_file} ) {
		$config->{general}{path_file} = "path_file";
	}
	if ( !$config->{option}{tophat2} ) {
		$config->{option}{tophat2} = "tophat2_option";
	}
	return ($config);
}

sub tophat2_by_pbs_config {
	my ( $paramFile, $runNow ) = @_;

	read_config $paramFile => my %config;

	my (
		$root_dir, $genome_db, $tophat2_param,
		$gtf_file, $gtfIndex,  $path_file
	) = tophat2_parse_and_check_parameters( \%config );
}

sub tophat2_by_pbs_individual {
	my (
		$refParamHash,    $refSampleNames, $refSampleFiles,
		$refPbsParamHash, $runNow
	) = @_;

	my (
		$root_dir, $genome_db, $tophat2_param,
		$gtf_file, $gtfIndex,  $path_file
	) = tophat2_parse_and_check_parameters($refParamHash);

	my @sampleNames     = @{$refSampleNames};
	my @sampleFiles     = @{$refSampleFiles};
	my $sampleNameCount = scalar(@sampleNames);

	if ( $sampleNameCount > 1 ) {
		if ( -e $gtf_file ) {
			if ( defined $gtfIndex ) {
				if ( !-e ( $gtfIndex . ".rev.2.bt2" ) ) {
					die
"Gtf file defined but index file has not been built, you should only run one job to build index first!";
				}
			}
		}
	}

	my ($isSingle) = check_is_single( $sampleNameCount, @sampleFiles );

	my ( $logDir, $pbsDir, $resultDir ) = init_dir($root_dir);
	my $tophatDir = create_directory_or_die( $resultDir . "/tophat2" );
	my ($pbsDesc) = get_pbs_desc($refPbsParamHash);

	for ( my $index = 0 ; $index < $sampleNameCount ; $index++ ) {
		my $sampleName = $sampleNames[$index];

		my $pbsFile = $pbsDir . "/${sampleName}_tophat2.pbs";
		my $log     = $logDir . "/${sampleName}_tophat2.log";

		output_header( $pbsFile, $pbsDesc, $path_file, $log );
		output_tophat2_script(
			$genome_db,     $gtf_file,  $gtfIndex,
			$tophat2_param, $tophatDir, $sampleName,
			0,              $isSingle,  @sampleFiles
		);
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
	my ( $config, $runNow ) = @_;

	my $root_dir = $config->{general}{root_dir}
	  or die "define general::root_dir first";
	my $cufflinksparam = $config->{cufflinks}{option}
	  or die "define cufflinks::option first";

	my $path_file = $config->{general}{path_file};
	my $refPbs = $config->{pbs} or die "define pbs parameters first";

	my ( $logDir, $pbsDir, $resultDir ) = init_dir($root_dir);
	my $cufflinkDir = create_directory_or_die( $resultDir . "/cufflinks" );
	my ($pbsDesc) = get_pbs_desc($refPbs);

	for ( keys %{ $config->{tophat2_result} } ) {
		my $sampleName = $_;
		my $sampleFile = $config->{tophat2_result}{$_};
		my $pbsFile    = $pbsDir . "/${sampleName}_cufflinks.pbs";
		my $log        = $logDir . "/${sampleName}_cufflinks.log";

		output_header( $pbsFile, $pbsDesc, $path_file, $log );

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
	my ( $config, $runNow ) = @_;

	my $root_dir = $config->{general}{root_dir}
	  or die "define general::root_dir first";
	my $genome_db = $config->{general}{bowtie2_index}
	  or die "define general::bowtie2_index first";
	my $genome_fasta = $genome_db . ".fa";

	my $task_name = $config->{general}{task_name}
	  or die "define general::task_name first";
	my $cuffdiffparam = $config->{cuffdiff}{option}
	  or die "define cuffdiff::option first";

	my $path_file = $config->{general}{path_file};

	my $refPbs = $config->{pbs} or die "define pbs parameters first";

	my $gtf_file  = $config->{general}{transcript_gtf};
	my $gtf_index = $config->{general}{transcript_gtf_index};

	if ( ( !defined($gtf_file) ) && ( !defined($gtf_index) ) ) {
		die
"define general::transcript_gtf or general::transcript_gtf_index first";
	}

	if ( defined $gtf_index ) {
		my $gff_file = $gtf_index . ".gff";
		if ( -e $gff_file ) {
			$gtf_file = $gff_file;
		}
	}

	my @labels = ();
	my @files  = ();
	for ( keys %{ $config->{files} } ) {
		push( @labels, $_ );
		my @gfiles = @{ $config->{files}{$_} };
		push( @files, merge_string( ",", @gfiles ) );
	}

	my $labels = merge_string( ",", @labels );

	my ( $logDir, $pbsDir, $resultDir ) = init_dir($root_dir);
	my $cuffdiffDir = create_directory_or_die( $resultDir . "/cuffdiff" );
	my ($pbsDesc) = get_pbs_desc($refPbs);

	my $pbsFile = $pbsDir . "/${task_name}_cuffdiff.pbs";
	my $log     = $logDir . "/${task_name}_cuffdiff.log";
	output_header( $pbsFile, $pbsDesc, $path_file, $log );

	print OUT "cuffdiff $cuffdiffparam -o $cuffdiffDir -L $labels ";

	print OUT "-b $genome_fasta ";

	print OUT " $gtf_file ";

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

sub cufflinks_by_pbs_old {
	my (
		$cufflinksparam, $root_dir,        $refSampleNames,
		$refSampleFiles, $refPbsParamHash, $runNow
	) = @_;

	my @sampleNames     = @{$refSampleNames};
	my @sampleFiles     = @{$refSampleFiles};
	my $sampleNameCount = scalar(@sampleNames);

	my $path_file = '/home/shengq1/bin/path.txt';

	my ( $logDir, $pbsDir, $resultDir ) = init_dir($root_dir);
	my $cufflinkDir = create_directory_or_die( $resultDir . "/cufflinks" );
	my ($pbsDesc) = get_pbs_desc($refPbsParamHash);

	for ( my $index = 0 ; $index < $sampleNameCount ; $index++ ) {
		my $sampleName = $sampleNames[$index];
		my $sampleFile = $sampleFiles[$index];

		my $pbsFile = $pbsDir . "/${sampleName}_cufflinks.pbs";
		my $log     = $logDir . "/${sampleName}_cufflinks.log";

		output_header( $pbsFile, $pbsDesc, $path_file, $log );

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

sub output_tophat2_script {
	my (
		$genome_db,     $gtf_file,  $gtfIndex,
		$tophat2_param, $tophatDir, $sampleName,
		$index,         $isSingle,  @sampleFiles
	) = @_;

	my $curDir = create_directory_or_die( $tophatDir . "/$sampleName" );

	print OUT "echo tophat2=`date` \n";

	my $hasgtf_file = -e $gtf_file;

	my $hasIndexFile = 0;
	if ( defined $gtfIndex ) {
		if ( -e ( $gtfIndex . ".rev.2.bt2" ) ) {
			$hasIndexFile = 1;
		}
	}

	my $tophat2file = $curDir . "/accepted_hits.bam";

	print OUT "if [ -s $tophat2file ];\n";
	print OUT "then\n";
	print OUT
"  echo job has already been done. if you want to do again, delete accepted_hits.bam and submit job again.\n";
	print OUT "else\n";
	if ($hasgtf_file) {
		if ( ( $index == 0 ) && ( !$hasIndexFile ) ) {
			print OUT
"  tophat2 $tophat2_param -G $gtf_file --transcriptome-index=$gtfIndex -o $curDir $genome_db ";
		}
		else {
			print OUT
"  tophat2 $tophat2_param --transcriptome-index=$gtfIndex -o $curDir $genome_db ";
		}
	}
	elsif ($hasIndexFile) {
		print OUT
"  tophat2 $tophat2_param --transcriptome-index=$gtfIndex -o $curDir $genome_db ";
	}
	else {
		print OUT "  tophat2 $tophat2_param -o $curDir $genome_db ";
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
		die
"Count of SampleName should be equals to count/half count of SampleFiles";
	}

	return ($isSingle);
}

sub output_header {
	my ( $pbsFile, $pbsDesc, $path_file, $log ) = @_;
	open( OUT, ">$pbsFile" ) or die $!;
	print OUT $pbsDesc;
	print OUT "#PBS -o $log\n";
	print OUT "#PBS -j oe\n\n";
	if ( defined $path_file ) {
		print OUT "source $path_file\n";
	}
}

sub output_footer() {
	print OUT "echo finished=`date`\n";
	close OUT;
}

1;
