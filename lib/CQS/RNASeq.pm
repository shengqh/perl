#!/usr/bin/perl
package CQS::RNASeq;

use strict;
use warnings;
use File::Basename;
use File::Copy;
use List::Compare;
use CQS::PBS;
use CQS::FileUtils;
use CQS::StringUtils;
use CQS::ConfigUtils;
use CQS::SystemUtils;

require Exporter;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = (
	'all' => [
		qw(tophat2_by_pbs get_tophat2_result cufflinks_by_pbs cuffmerge_by_pbs cuffdiff_by_pbs read_cufflinks_fpkm read_cuffdiff_significant_genes copy_and_rename_cuffdiff_file compare_cuffdiff miso_by_pbs)
	]
);

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

	my $has_gtf_file = file_exists($transcript_gtf);

	my $has_index_file = transcript_gtf_index_exists($transcript_gtf_index);

	my $tophat2file = $curDir . "/accepted_hits.bam";

	print OUT "if [ -s $tophat2file ];\n";
	print OUT "then\n";
	print OUT "  echo job has already been done. if you want to do again, delete accepted_hits.bam and submit job again.\n";
	print OUT "else\n";
	if ($has_gtf_file) {
		print OUT "  tophat2 $tophat2_param -G $transcript_gtf --transcriptome-index=$transcript_gtf_index -o $curDir $bowtie2_index ";
	}
	elsif ($has_index_file) {
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

sub tophat2_by_pbs {
	my ( $config, $section ) = @_;

	my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option ) = get_parameter( $config, $section );

	my $bowtie2_index = $config->{general}{bowtie2_index} or die "define general::bowtie2_index first";

	my $batchmode = $config->{$section}{batchmode};

	if ( !defined($batchmode) ) {
		$batchmode = 0;
	}

	my %fqFiles = %{ get_raw_files( $config, $section ) };
	my $sampleNameCount = 0;
	while ( my ( $sampleName, $sampleFiles ) = each(%fqFiles) ) {
		$sampleNameCount = $sampleNameCount + scalar( @{$sampleFiles} );
	}

	my $transcript_gtf = get_param_file( $config->{$section}{transcript_gtf}, "transcript_gtf", 0 );
	my $transcript_gtf_index;
	if ( defined $transcript_gtf ) {
		$transcript_gtf_index = $config->{$section}{transcript_gtf_index};
	}
	elsif ( defined $config->{$section}{transcript_gtf_ref} ) {
		$transcript_gtf = get_param_file( $config->{general}{ $config->{$section}{transcript_gtf_ref} }, "general::$config->{$section}{transcript_gtf_ref}", 1 );
		if ( defined $config->{$section}{transcript_gtf_index_ref} ) {
			$transcript_gtf_index = $config->{general}{ $config->{$section}{transcript_gtf_index_ref} };
		}
	}

	if ( ( defined $transcript_gtf ) && ( -e $transcript_gtf ) ) {
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
		my $pbsFile = $pbsDir . "/${task_name}_th2.pbs";
		my $log     = $logDir . "/${task_name}_th2.log";

		output_header( $pbsFile, $pbsDesc, $path_file, $log );

		my $index = 0;
		for my $sampleName ( sort keys %fqFiles ) {
			my @sampleFiles = @{ $fqFiles{$sampleName} };
			output_tophat2( $bowtie2_index, $transcript_gtf, $transcript_gtf_index, $option, $resultDir, $sampleName, $index, @sampleFiles );
			$index++;
		}

		output_footer();

		print "$pbsFile created\n";
	}
	else {
		my $shfile = $pbsDir . "/${task_name}.submit";
		open( SH, ">$shfile" ) or die "Cannot create $shfile";

		for my $sampleName ( sort keys %fqFiles ) {
			my @sampleFiles = @{ $fqFiles{$sampleName} };

			my $pbsName = "${sampleName}_th2.pbs";
			my $pbsFile = $pbsDir . "/$pbsName";
			my $log     = $logDir . "/${sampleName}_th2.log";

			output_header( $pbsFile, $pbsDesc, $path_file, $log );
			output_tophat2( $bowtie2_index, $transcript_gtf, $transcript_gtf_index, $option, $resultDir, $sampleName, 0, @sampleFiles );
			output_footer();

			print SH "qsub ./$pbsName \n";
			print SH "echo $pbsName was submitted. \n\n";
			print "$pbsFile created\n";
		}

		print SH "exit 0\n";
		close(SH);

		if ( is_linux() ) {
			chmod 0755, $shfile;
		}
		print "!!!shell file $shfile created, you can run this shell file to submit all tophat2 tasks.\n";
	}
}

#get expected tophat2 result based on tophat2 definition
sub get_tophat2_map {
	my ( $config, $section ) = @_;

	if ( defined $config->{$section}{source} ) {
		return $config->{$section}{source};
	}

	if ( defined $config->{$section}{source_ref} ) {
		my $result = {};

		my $tophatsection = $config->{$section}{source_ref};
		my $tophat_dir = $config->{$tophatsection}{target_dir} or die "${tophatsection}::target_dir not defined.";
		my ( $logDir, $pbsDir, $resultDir ) = init_dir( $tophat_dir, 0 );
		my %fqFiles = %{ get_raw_files( $config, $tophatsection ) };
		for my $sampleName ( keys %fqFiles ) {
			$result->{$sampleName} = "${resultDir}/${sampleName}/accepted_hits.bam";
		}

		return $result;
	}

	die "define ${section}::source or ${section}::source_ref first!";
}

sub cufflinks_by_pbs {
	my ( $config, $section ) = @_;

	my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option ) = get_parameter( $config, $section );

	my $transcript_gtf = get_param_file( $config->{$section}{transcript_gtf}, "transcript_gtf", 0 );
	my $gtf = "";
	if ( defined $transcript_gtf ) {
		$gtf = "-g $transcript_gtf";
	}

	my %tophat2map = %{ get_tophat2_map( $config, $section ) };

	my $shfile = $pbsDir . "/${task_name}.submit";
	open( SH, ">$shfile" ) or die "Cannot create $shfile";

	for my $sampleName ( sort keys %tophat2map ) {
		my $tophat2File = $tophat2map{$sampleName};

		my $pbsName = "${sampleName}_clinks.pbs";
		my $pbsFile = $pbsDir . "/$pbsName";

		print SH "if [ -s ${resultDir}/${sampleName}/transcripts.gtf ];\n";
		print SH "then\n";
		print SH "  echo job has already been done. if you want to do again, delete ${sampleName}/transcripts.gtf and submit job again.\n";
		print SH "else\n";
		print SH "  if [ ! -s $tophat2File ];\n";
		print SH "  then";
		print SH "    echo tophat2 of ${sampleName} has not finished, ignore current job. \n";
		print SH "  else\n";
		print SH "    qsub ./$pbsName \n";
		print SH "    echo $pbsName was submitted. \n";
		print SH "  fi\n";
		print SH "fi\n";

		my $log = $logDir . "/${sampleName}_clinks.log";

		output_header( $pbsFile, $pbsDesc, $path_file, $log );

		my $curDir = create_directory_or_die( $resultDir . "/$sampleName" );

		print OUT "echo cufflinks=`date` \n";

		print OUT "cufflinks $option $gtf -o $curDir $tophat2File \n";

		output_footer();

		print "$pbsFile created\n";
	}
	print SH "exit 0\n";
	close(SH);

	if ( is_linux() ) {
		chmod 0755, $shfile;
	}
	print "!!!shell file $shfile created, you can run this shell file to submit all cufflinks tasks.\n";
}

sub get_cufflinks_result_gtf {
	my ( $config, $section ) = @_;

	#get cufflinks root directory
	my $cufflinks_dir = $config->{$section}{target_dir} or die "${section}::target_dir not defined.";

	#get cufflinks result directory
	my ( $logDir, $pbsDir, $resultDir ) = init_dir( $cufflinks_dir, 0 );

	my %tophat2map = %{ get_tophat2_map( $config, $section ) };

	my @result = ();

	#get expected gtf file list
	for my $sampleName ( sort keys %tophat2map ) {
		push( @result, "${resultDir}/${sampleName}/transcripts.gtf" );
	}

	#return expected gtf file list
	return ( \@result );
}

sub get_assemblies_file {
	my ( $config, $section, $target_dir ) = @_;

	my $result = get_param_file( $config->{$section}{source}, "${section}::source", 0 );

	if ( defined $result ) {
		return $result;
	}

	my $cufflinkssection = $config->{$section}{source_ref};
	if ( defined $cufflinkssection ) {
		my $cufflinks_gtf = get_cufflinks_result_gtf( $config, $cufflinkssection );
		$result = $target_dir . "/assemblies.txt";
		open( OUT, ">$result" ) or die $!;
		for my $gtf ( @{$cufflinks_gtf} ) {
			print OUT "${gtf}\n";
		}
		close OUT;

		return $result;
	}

	die "define ${section}::source or ${section}::source_ref first!";
}

sub cuffmerge_by_pbs {
	my ( $config, $section ) = @_;

	my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option ) = get_parameter( $config, $section );

	my $transcript_gtf = get_param_file( $config->{general}{transcript_gtf}, "transcript_gtf", 0 );
	my $bowtie2_index = $config->{general}{bowtie2_index} or die "define general::bowtie2_index first";
	my $bowtie2_fasta = get_param_file( $bowtie2_index . ".fa", "bowtie2_fasta", 1 );

	my $assembliesfile = get_assemblies_file( $config, $section, $resultDir );

	my $pbsFile = $pbsDir . "/${task_name}_cmerge.pbs";
	my $log     = $logDir . "/${task_name}_cmerge.log";

	output_header( $pbsFile, $pbsDesc, $path_file, $log );

	print OUT "echo cuffmerge=`date` \n";

	my $gtfparam = "";
	if ($transcript_gtf) {
		$gtfparam = "-g $transcript_gtf";
	}

	print OUT "cuffmerge $option $gtfparam -s $bowtie2_fasta -o $resultDir $assembliesfile \n";

	output_footer();

	print "$pbsFile created\n";

	open( FILE, $assembliesfile ) or die("Unable to open file $assembliesfile");
	my @data = <FILE>;
	close(FILE);

	my $shfile = $pbsDir . "/${task_name}.submit";
	open( SH, ">$shfile" ) or die "Cannot create $shfile";
	for my $gtf (@data) {
		chomp($gtf);
		print SH "if [ ! -s $gtf ];\n";
		print SH "then\n";
		print SH "  echo $gtf is not exists, cannot submit the job.\n";
		print SH "  exit;\n";
		print SH "fi;\n\n";
	}
	print SH "qsub $pbsFile\n";
	print SH "exit 0\n";
	close(SH);

	if ( is_linux() ) {
		chmod 0755, $shfile;
	}
	print "!!!shell file $shfile created, you can run this shell file to submit cuffmerge task.\n";
}

sub get_cuffdiff_gtf {
	my ( $config, $section ) = @_;

	my $transcript_gtf = get_param_file( $config->{$section}{transcript_gtf}, "transcript_gtf", 0 );
	if ( defined $transcript_gtf ) {
		return ($transcript_gtf);
	}

	my $cuffmergesection = $config->{$section}{transcript_gtf_ref};
	if ( defined $cuffmergesection ) {
		my $cuffmerge_target_dir = $config->{$cuffmergesection}{target_dir};
		my ( $logDir, $pbsDir, $resultDir ) = init_dir($cuffmerge_target_dir);
		return ("${resultDir}/merged.gtf");
	}

	$transcript_gtf = get_param_file( $config->{general}{transcript_gtf}, "transcript_gtf", 0 );
	if ( defined $transcript_gtf ) {
		return ($transcript_gtf);
	}

	die "define ${section}::transcript_gtf or ${section}::transcript_gtf_ref or general::transcript_gtf first!";
}

sub get_cuffdiff_groups {
	my ( $config, $section ) = @_;

	my $result = $config->{$section}{groups};

	if ( defined $result ) {
		return $result;
	}

	if ( !defined $config->{$section}{groups_ref} ) {
		die "define ${section}::groups or ${section}::groups_ref first!";
	}

	$result = $config->{ $config->{$section}{groups_ref} };
	if ( !defined $result ) {
		die "$config->{$section}{groups_ref} was not defined!";
	}

	return $result;
}

sub get_cuffdiff_pairs {
	my ( $config, $section ) = @_;

	my $result = $config->{$section}{pairs};

	if ( defined $result ) {
		return $result;
	}

	if ( !defined $config->{$section}{pairs_ref} ) {
		die "define ${section}::pairs or ${section}::pairs_ref first!";
	}

	$result = $config->{ $config->{$section}{pairs_ref} };
	if ( !defined $result ) {
		die "$config->{$section}{pairs_ref} was not defined!";
	}

	return $result;
}

sub cuffdiff_by_pbs {
	my ( $config, $section ) = @_;

	my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option ) = get_parameter( $config, $section );

	my $bowtie2_index = $config->{general}{bowtie2_index} or die "define general::bowtie2_index first";
	my $bowtie2_fasta = get_param_file( $bowtie2_index . ".fa", "bowtie2_fasta", 1 );

	my $transcript_gtf = get_cuffdiff_gtf( $config, $section );

	my $tophat2map = get_tophat2_map( $config, $section );

	my $groups = get_cuffdiff_groups( $config, $section );

	my $pairs = get_cuffdiff_pairs( $config, $section );

	my @labels = ();
	my @files  = ();

	my %tpgroups = ();
	for my $groupName ( sort keys %{$groups} ) {
		my @samples = @{ $groups->{$groupName} };
		my @gfiles  = ();
		foreach my $sampleName (@samples) {
			my $tophat2File = $tophat2map->{$sampleName};
			push( @gfiles, $tophat2File );
		}
		$tpgroups{$groupName} = \@gfiles;

		#print " $groupName => $tpgroups{$groupName} \n";
	}

	my $shfile = $pbsDir . "/${task_name}.submit";
	open( SH, ">$shfile" ) or die "Cannot create $shfile";

	print SH "if [ ! -s $transcript_gtf ];\n";
	print SH "then\n";
	print SH "  echo $transcript_gtf is not exists! all job will be ignored.\n";
	print SH "  exit 1\n";
	print SH "fi\n\n";

	for my $pairName ( sort keys %{$pairs} ) {
		my @groupNames = @{ $pairs->{$pairName} };

		my $pbsName = "${pairName}_cdiff.pbs";

		my $pbsFile = $pbsDir . "/$pbsName";
		my $log     = $logDir . "/${pairName}_cdiff.log";

		my $curDir = create_directory_or_die( $resultDir . "/$pairName" );

		my $labels = merge_string( ",", @groupNames );

		output_header( $pbsFile, $pbsDesc, $path_file, $log );
		print OUT "cuffdiff $option -o $curDir -L $labels -b $bowtie2_fasta $transcript_gtf ";

		my @conditions = ();
		foreach my $groupName (@groupNames) {
			my @bamfiles = @{ $tpgroups{$groupName} };
			my $bams = merge_string( ",", @bamfiles );
			print OUT "$bams ";

			foreach my $bam (@bamfiles) {
				push( @conditions, "[ -s $bam ]" );
			}
		}
		print OUT "\n";

		output_footer();

		print "$pbsFile created. \n";

		my $condition = merge_string( " && ", @conditions );
		print SH "if [ -s ${curDir}/gene_exp.diff ];\n";
		print SH "then\n";
		print SH "  echo job has already been done. if you want to do again, delete ${curDir}/gene_exp.diff and submit job again.\n";
		print SH "else\n";
		print SH "  if $condition;\n";
		print SH "  then\n";
		print SH "    qsub ./$pbsName \n";
		print SH "    echo $pbsName was submitted. \n";
		print SH "  else\n";
		print SH "    echo some required file not exists! $pbsName will be ignored.\n";
		print SH "  fi\n";
		print SH "fi\n\n";
	}

	print SH "exit 0\n";
	close(SH);

	if ( is_linux() ) {
		chmod 0755, $shfile;
	}
	print "!!!shell file $shfile created, you can run this shell file to submit cuffdiff task.\n";
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

sub read_cufflinks_fpkm {
	my $genefile = shift;
	open GF, "<$genefile" or die "Cannot open file $genefile";
	my $header = <GF>;
	my @headers = split( /\t/, $header );
	my ($geneindex) = grep { $headers[$_] eq "gene_id" } 0 .. $#headers;
	my ($fpkmindex) = grep { $headers[$_] eq "FPKM" } 0 .. $#headers;

	my $result = {};
	while ( my $line = <GF> ) {
		my @values = split( /\t/, $line );
		$result->{ $values[$geneindex] } = $values[$fpkmindex];
	}
	close(GF);

	return $result;
}

sub read_cuffdiff_significant_genes {
	my $file = shift;

	my $result = {};
	open IN, "<$file" or die "Cannot open file $file";
	my $header = <IN>;
	chomp($header);
	while ( my $line = <IN> ) {
		chomp($line);
		my @part = split( /\t/, $line );
		if ( $part[13] eq "yes" ) {
			$result->{ $part[2] } = $line;
		}
	}
	close IN;

	return ( $result, $header );
}

#use CQS::RNASeq;
#my $root = "/home/xxx/cuffdiff/result";
#my $targetdir = create_directory_or_die( $root . "/comparison" );
#my $config    = {
#    "RenameDiff" => {
#        target_dir => $targetdir,
#        root_dir   => $root
#    }
#}
#copy_and_rename_cuffdiff_file($config, "RenameDiff");
sub copy_and_rename_cuffdiff_file {
	my ( $config, $section ) = @_;
	my $targetdir = $config->{$section}{"target_dir"} or die "define ${section}::target_dir first";
	my $dir       = $config->{$section}{"root_dir"}   or die "define ${section}::root_dir first";
	if ( !-d $dir ) {
		die "directory $dir is not exists";
	}

	my @subdirs = list_directories($dir);
	if ( 0 == scalar(@subdirs) ) {
		die "$dir has no sub CuffDiff directories";
	}

	my @filenames = ( "gene_exp.diff", "splicing.diff" );
	for my $subdir (@subdirs) {
		foreach my $filename (@filenames) {
			my $file = "${dir}/${subdir}/${filename}";
			print "Dealing " . $file . "\n";
			if ( -s $file ) {
				open IN, "<$file" or die "Cannot open file $file";
				my $line = <IN>;
				$line = <IN>;
				close(IN);

				my @parts = split( /\t/, $line );
				my $targetname = "${targetdir}/" . $parts[4] . "_vs_" . $parts[5] . "." . ${filename};

				copy( $file, $targetname ) or die "copy failed : $!";

				my $target_sign_name = $targetname . ".sig";
				my $cmd              = "cat $targetname | awk '\$14==\"yes\" || \$14==\"significant\"' > $target_sign_name";

				print $cmd . "\n";
				`$cmd`;
			}
		}
	}
}

sub output_compare_cuffdiff_file {
	my ( $data1, $data2, $fileName, $header, @genes ) = @_;
	open OUT, ">$fileName" or die "Cannot create file $fileName";
	print OUT "$header\n";
	my @sortedgenes = sort @genes;
	for my $gene (@sortedgenes) {
		if ( defined $data1->{$gene} ) {
			print OUT "$data1->{$gene}\n";
		}
		if ( defined $data2->{$gene} ) {
			print OUT "$data2->{$gene}\n";
		}
	}
	close(OUT);
}

sub compare_cuffdiff {
	my ( $config, $section ) = @_;

	my $info = $config->{$section};

	my ( $file1, $file2 ) = @{ $info->{"files"} };

	my ( $data1, $header )  = read_cuffdiff_significant_genes($file1);
	my ( $data2, $header2 ) = read_cuffdiff_significant_genes($file2);

	my @genes1 = keys %{$data1};
	my @genes2 = keys %{$data2};

	my $lc = List::Compare->new( \@genes1, \@genes2 );

	my @resultgenes = ();
	if ( $info->{operation} eq "minus" ) {
		@resultgenes = $lc->get_Lonly();
	}
	elsif ( $info->{operation} eq "intersect" ) {
		@resultgenes = $lc->get_intersection();
	}
	else {
		die "Only minus or intersect is supported.";
	}

	my $resultFileName = $info->{"target_file"};

	output_compare_cuffdiff_file( $data1, $data2, $resultFileName, $header, @resultgenes );
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

sub miso_by_pbs {
	my ( $config, $section ) = @_;

	my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option ) = get_parameter( $config, $section );

	my $indexed = $config->{$section}{indexed_gff};

	my %tophat2map = %{ get_tophat2_map( $config, $section ) };

	my $shfile = $pbsDir . "/${task_name}.submit";
	open( SH, ">$shfile" ) or die "Cannot create $shfile";

	for my $sampleName ( sort keys %tophat2map ) {
		my $tophat2File      = $tophat2map{$sampleName};
		my $tophat2indexFile = $tophat2File . ".bai";

		my $pbsName = "${sampleName}_clinks.pbs";
		my $pbsFile = $pbsDir . "/$pbsName";

		print SH "  if [ ! -s $tophat2File ];\n";
		print SH "  then";
		print SH "    echo tophat2 of ${sampleName} has not finished, ignore current job. \n";
		print SH "  else\n";
		print SH "    qsub ./$pbsName \n";
		print SH "    echo $pbsName was submitted. \n";
		print SH "  fi\n";

		my $log = $logDir . "/${sampleName}_miso.log";

		output_header( $pbsFile, $pbsDesc, $path_file, $log );

		my $curDir = create_directory_or_die( $resultDir . "/$sampleName" );

		print OUT "echo MISO=`date` \n";

		print OUT "samtools index $tophat2File \n";

		print OUT "run_events_analysis.py --compute-genes-psi $indexed $tophat2File --output-dir $curDir --read-len 35 \n";

		output_footer();

		print "$pbsFile created\n";
	}
	print SH "exit 0\n";
	close(SH);

	if ( is_linux() ) {
		chmod 0755, $shfile;
	}
	print "!!!shell file $shfile created, you can run this shell file to submit all miso tasks.\n";
}

1;
