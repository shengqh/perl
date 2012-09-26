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

our %EXPORT_TAGS = ( 'all' => [qw(tophat2_by_pbs get_tophat2_result cufflinks_by_pbs cuffmerge_by_pbs cuffdiff_by_pbs)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

use Cwd;

sub file_exists {
	my $file   = shift;
	my $result = 0;
	if ( defined($file) ) {
		$result = -e $file;
	}
	return ($result);
}

sub transcript_gtf_index_exists {
	my $transcript_gtf_index = shift;
	my $result               = 0;
	if ( defined($transcript_gtf_index) ) {
		my $file = $transcript_gtf_index . ".rev.1.bt2";
		$result = -e $file;
	}
	return ($result);
}

sub output_tophat2 {
	my ( $bowtie2_index, $transcript_gtf, $transcript_gtf_index, $tophat2_param, $tophatDir, $sampleName, $index, @sampleFiles ) = @_;

	my $curDir = create_directory_or_die( $tophatDir . "/$sampleName" );

	print OUT "echo tophat2=`date` \n";

	my $hasgtf_file = file_exists($transcript_gtf);

	my $hasIndexFile = transcript_gtf_index_exists($transcript_gtf_index);

	my $tophat2file = $curDir . "/accepted_hits.bam";

	print OUT "if [ -s $tophat2file ];\n";
	print OUT "then\n";
	print OUT "  echo job has already been done. if you want to do again, delete accepted_hits.bam and submit job again.\n";
	print OUT "else\n";
	if ($hasgtf_file) {
		if ( ( $index == 0 ) && ( !$hasIndexFile ) ) {
			print OUT "  tophat2 $tophat2_param -G $transcript_gtf --transcriptome-index=$transcript_gtf_index -o $curDir $bowtie2_index ";
		}
		else {
			print OUT "  tophat2 $tophat2_param --transcriptome-index=$transcript_gtf_index -o $curDir $bowtie2_index ";
		}
	}
	elsif ($hasIndexFile) {
		print OUT "  tophat2 $tophat2_param --transcriptome-index=$transcript_gtf_index -o $curDir $bowtie2_index ";
	}
	else {
		print OUT "  tophat2 $tophat2_param -o $curDir $bowtie2_index ";
	}

	for my $sampleFile (@sampleFiles) {
		print OUT "$sampleFile ";
	}
	print OUT "\nfi\n\n";
}

#get parameter which indicates a file. If required, not defined or not exists, die. If defined but not exists, die.
#returned file either undefined or exists.
sub get_param_file {
	my ( $file, $name, $required ) = @_;

	my $result = $file;

	if ($required) {
		if ( !defined $file ) {
			die "$name was not defined!";
		}

		if ( !-e $file ) {
			die "$name $file defined but not exists!";
		}
	}
	else {
		if ( defined($file) ) {
			if ( $file eq "" ) {
				undef($result);
			}
			elsif ( !-e $file ) {
				die "$name $file defined but not exists!";
			}
		}
	}
	return ($result);
}

sub get_parameter {
	my ( $config, $section ) = @_;

	my $path_file = get_param_file( $config->{general}{path_file}, "path_file", 0 );
	my $refPbs     = $config->{$section}{pbs}        or die "define ${section}::pbs parameters first";
	my $target_dir = $config->{$section}{target_dir} or die "define ${section}::target_dir parameters first";
	my ( $logDir, $pbsDir, $resultDir ) = init_dir($target_dir);
	my ($pbsDesc) = get_pbs_desc($refPbs);
	my $option = $config->{$section}{option} or die "define ${section}::option first";

	return ( $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option );
}

sub tophat2_by_pbs {
	my ( $config, $section, $runNow ) = @_;

	my ( $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option ) = get_parameter( $config, $section );

	my $bowtie2_index = $config->{general}{bowtie2_index} or die "define general::bowtie2_index first";

	my $batchmode = $config->{$section}{batchmode};
	my $task_name = $config->{general}{task_name};
	if ( !defined($batchmode) ) {
		$batchmode = 0;
	}
	if ($batchmode) {
		if ( !defined $task_name ) {
			die "define general::task_name for batchmode first!";
		}
	}

	my $sampleNameCount = 0;
	my %fqFiles         = %{ $config->{ $config->{$section}{source} } };
	while ( my ( $groupName, $sampleMap ) = each(%fqFiles) ) {
		$sampleNameCount = $sampleNameCount + scalar( keys %{$sampleMap} );
	}

	my $transcript_gtf = get_param_file( $config->{general}{transcript_gtf}, "transcript_gtf", 0 );
	my $transcript_gtf_index = $config->{general}{transcript_gtf_index};

	if ( -e $transcript_gtf ) {
		if ( !defined $transcript_gtf_index ) {
			die "transcript_gtf was defined but transcript_gtf_index was not defined, you should defined transcript_gtf_index to cache the parsing result.";
		}

		if ( ( !$batchmode ) && ( $sampleNameCount > 1 ) ) {
			if ( !-e ( $transcript_gtf_index . ".rev.1.bt2" ) ) {
				die "transcript_gtf was defined but transcript_gtf_index has not been built, you should run only one job to build index first!";
			}
		}
	}
	elsif ( defined $transcript_gtf_index ) {
		if ( !transcript_gtf_index_exists($transcript_gtf_index) ) {
			die "transcript_gtf_index $transcript_gtf_index defined but not exists!";
		}
	}

	if ($batchmode) {
		my $pbsFile = $pbsDir . "/${task_name}_tophat2.pbs";
		my $log     = $logDir . "/${task_name}_tophat2.log";

		output_header( $pbsFile, $pbsDesc, $path_file, $log );

		my $index = 0;
		for my $groupName ( sort keys %fqFiles ) {
			my %sampleMap = %{ $fqFiles{$groupName} };
			for my $sampleName ( sort keys %sampleMap ) {
				my @sampleFiles = @{ $sampleMap{$sampleName} };
				output_tophat2( $bowtie2_index, $transcript_gtf, $transcript_gtf_index, $option, $resultDir, $sampleName, $index, @sampleFiles );
				$index++;
			}
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
	else {
		my $shfile = $pbsDir . "/submit.sh";
		open( SH, ">$shfile" ) or die "Cannot create $shfile";

		for my $groupName ( sort keys %fqFiles ) {
			my %sampleMap = %{ $fqFiles{$groupName} };
			for my $sampleName ( sort keys %sampleMap ) {
				my @sampleFiles = @{ $sampleMap{$sampleName} };

				my $pbsFile = $pbsDir . "/${sampleName}_tophat2.pbs";
				my $log     = $logDir . "/${sampleName}_tophat2.log";

				output_header( $pbsFile, $pbsDesc, $path_file, $log );
				output_tophat2( $bowtie2_index, $transcript_gtf, $transcript_gtf_index, $option, $resultDir, $sampleName, 0, @sampleFiles );
				output_footer();

				print SH "qsub $pbsFile \n";
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
}

#get expected tophat2 result based on tophat2 definition
sub get_tophat2_bam {
	my ( $config, $section ) = @_;
	my $tophat_dir = $config->{$section}{target_dir} or die "${section}::target_dir not defined.";
	my ( $logDir, $pbsDir, $resultDir ) = init_dir( $tophat_dir, 0 );
	my %fqFiles = %{ $config->{ $config->{$section}{source} } };
	my $result  = {};
	for my $groupName ( sort keys %fqFiles ) {
		my %sampleMap = %{ $fqFiles{$groupName} };
		for my $sampleName ( sort keys %sampleMap ) {
			$result->{$groupName}{$sampleName} = "${resultDir}/${sampleName}/accepted_hits.bam";
		}
	}
	return ($result);
}

sub get_tophat2_map {
	my ( $config, $section ) = @_;

	my $result;

	if ( defined $config->{$section}{sourcefiles} ) {
		$result = $config->{$section}{sourcefiles};
	}
	elsif ( defined $config->{$section}{source} ) {
		$result = get_tophat2_bam( $config, $config->{$section}{source} );
	}
	else {
		die "either {$section}::sourcefiles or {$section}::source should be defined.";
	}

	return ($result);
}

sub cufflinks_by_pbs {
	my ( $config, $section, $runNow ) = @_;

	my ( $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option ) = get_parameter( $config, $section );

	my $tophat2map = get_tophat2_map( $config, $section );
	for my $groupName ( sort keys %{$tophat2map} ) {
		my %sampleMap = %{ $tophat2map->{$groupName} };
		for my $sampleName ( sort keys %sampleMap ) {
			my $tophat2File = $sampleMap{$sampleName};
			my $pbsFile     = $pbsDir . "/${sampleName}_cufflinks.pbs";
			my $log         = $logDir . "/${sampleName}_cufflinks.log";

			output_header( $pbsFile, $pbsDesc, $path_file, $log );

			my $curDir = create_directory_or_die( $resultDir . "/$sampleName" );

			print OUT "echo cufflinks=`date` \n";
			print OUT "cufflinks $option -o $curDir $tophat2File \n";

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
}

sub get_cufflinks_gtf {
	my ( $config, $section ) = @_;

	#get cufflinks root directory
	my $cufflinks_dir = $config->{$section}{target_dir} or die "${section}::target_dir not defined.";

	#get cufflinks result directory
	my ( $logDir, $pbsDir, $resultDir ) = init_dir( $cufflinks_dir, 0 );

	my $tophat2map = get_tophat2_map( $config, $section );

	my @result = ();

	#get expected gtf file list
	for my $groupName ( sort keys %{$tophat2map} ) {
		my %sampleMap = %{ $tophat2map->{$groupName} };
		for my $sampleName ( sort keys %sampleMap ) {
			push( @result, "${resultDir}/${sampleName}/transcripts.gtf" );
		}
	}

	#return expected gtf file list
	return ( \@result );
}

sub cuffmerge_by_pbs {
	my ( $config, $section, $runNow ) = @_;

	my ( $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option ) = get_parameter( $config, $section );

	my $task_name = $config->{general}{task_name} or die "define general::task_name first";
	my $transcript_gtf = get_param_file( $config->{general}{transcript_gtf}, "transcript_gtf", 0 );
	my $bowtie2_index = $config->{general}{bowtie2_index} or die "define general::bowtie2_index first";
	my $bowtie2_fasta = get_param_file( $bowtie2_index . ".fa", "bowtie2_fasta", 1 );

	my $assembliesfile = get_param_file( $config->{$section}{assemblies_file}, "${section}::assemblies_file", 0 );
	if ( !defined $assembliesfile ) {
		my $cufflinks_gtf = get_cufflinks_gtf( $config, $config->{$section}{source} );
		$assembliesfile = $target_dir . "/assemblies.txt";
		open( OUT, ">$assembliesfile" ) or die $!;
		for my $gtf ( @{$cufflinks_gtf} ) {
			print OUT "${gtf}\n";
		}
		close OUT;
	}

	my $pbsFile = $pbsDir . "/${task_name}_cuffmerge.pbs";
	my $log     = $logDir . "/${task_name}_cuffmerge.log";

	output_header( $pbsFile, $pbsDesc, $path_file, $log );

	print OUT "echo cuffmerge=`date` \n";

	my $gtfparam = "";
	if ($transcript_gtf) {
		$gtfparam = "-g $transcript_gtf";
	}

	print OUT "cuffmerge $option $gtfparam -s $bowtie2_fasta -o $resultDir $assembliesfile \n";

	output_footer();

	if ($runNow) {
		`qsub $pbsFile`;
		print "$pbsFile submitted\n";
	}
	else {
		print "$pbsFile created\n";
	}
}

sub cuffdiff_by_pbs {
	my ( $config, $section, $runNow ) = @_;

	my ( $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option ) = get_parameter( $config, $section );

	my $task_name     = $config->{general}{task_name}     or die "define general::task_name first";
	my $bowtie2_index = $config->{general}{bowtie2_index} or die "define general::bowtie2_index first";
	my $bowtie2_fasta = get_param_file( $bowtie2_index . ".fa", "bowtie2_fasta", 1 );

	my $transcript_gtf = get_param_file( $config->{general}{transcript_gtf}, "transcript_gtf", 1 );

	my $tophat2map = get_tophat2_map( $config, $section );

	my @labels = ();
	my @files  = ();
	for my $groupName ( sort keys %{$tophat2map} ) {
		push( @labels, $groupName );

		my @gfiles    = ();
		my %sampleMap = %{ $tophat2map->{$groupName} };
		for my $sampleName ( sort keys %sampleMap ) {
			my $tophat2File = $sampleMap{$sampleName};
			push( @gfiles, $tophat2File );
		}
		push( @files, merge_string( ",", @gfiles ) );
	}

	my $labels = merge_string( ",", @labels );

	my $pbsFile = $pbsDir . "/${task_name}_cuffdiff.pbs";
	my $log     = $logDir . "/${task_name}_cuffdiff.log";
	output_header( $pbsFile, $pbsDesc, $path_file, $log );
	print OUT "cuffdiff $option -o $resultDir -L $labels -b $bowtie2_fasta $transcript_gtf ";

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
	my ( $genome_db, $gtf_file, $gtfIndex, $tophat2_param, $tophatDir, $sampleName, $index, $isSingle, @sampleFiles ) = @_;

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
	print OUT "  echo job has already been done. if you want to do again, delete accepted_hits.bam and submit job again.\n";
	print OUT "else\n";
	if ($hasgtf_file) {
		if ( ( $index == 0 ) && ( !$hasIndexFile ) ) {
			print OUT "  tophat2 $tophat2_param -G $gtf_file --transcriptome-index=$gtfIndex -o $curDir $genome_db ";
		}
		else {
			print OUT "  tophat2 $tophat2_param --transcriptome-index=$gtfIndex -o $curDir $genome_db ";
		}
	}
	elsif ($hasIndexFile) {
		print OUT "  tophat2 $tophat2_param --transcriptome-index=$gtfIndex -o $curDir $genome_db ";
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
		die "Count of SampleName should be equals to count/half count of SampleFiles";
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
