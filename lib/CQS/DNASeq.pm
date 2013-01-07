#!/usr/bin/perl
package CQS::DNASeq;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;

require Exporter;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(bwa_by_pbs_single bwa_by_pbs_double samtools_index)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

use Cwd;

sub bwa_by_pbs_single {
	my ( $config, $section ) = @_;

	my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option ) = get_parameter( $config, $section );

	my $faFile = get_param_file( $config->{$section}{fasta_file}, "fasta_file", 1 );

	die "define ${section}::option_samse first" if ( !defined $config->{$section}{option_samse} );
	my $option_samse = $config->{$section}{option_samse};

	my %rawFiles = %{ get_raw_files( $config, $section ) };

	my $shfile = $pbsDir . "/${task_name}.sh";
	open( SH, ">$shfile" ) or die "Cannot create $shfile";
	print SH "type -P qsub &>/dev/null && export MYCMD=\"qsub\" || export MYCMD=\"bash\" \n";

	for my $sampleName ( sort keys %rawFiles ) {
		my @sampleFiles = @{ $rawFiles{$sampleName} };

		my $sampleFile1 = $sampleFiles[0];

		my ( $sampleName1, $directories1, $suffix1 ) = fileparse($sampleFile1);
		my $saiFile1      = $sampleName1 . ".sai";
		my $samFile       = $sampleName . ".sam";
		my $bamFile       = $sampleName . ".bam";
		my $sortedBamFile = $sampleName . "_sort";

		my $pbsName = "${sampleName}_bwa.pbs";
		my $pbsFile = "${pbsDir}/$pbsName";

		print SH "\$MYCMD ./$pbsName \n";

		my $log = "${logDir}/${sampleName}_bwa.log";

		open( OUT, ">$pbsFile" ) or die $!;
		print OUT $pbsDesc;
		print OUT "#PBS -o $log\n";
		print OUT "#PBS -j oe\n\n";

		if ( -e $path_file ) {
			print OUT "source $path_file\n";
		}
		print OUT "echo bwa=`date`\n";

		my $curDir = create_directory_or_die( $resultDir . "/$sampleName" );

		#my $tag="'\@RG\tID:$sample\tLB:$sample\tSM:$sample\tPL:ILLUMINA'";
		print OUT "cd $curDir\n\n";

		print OUT "if [ -s ${sortedBamFile}.bam ]; then\n";
		print OUT "  echo job has already been done. if you want to do again, delete ${sortedBamFile}.bam and submit job again.\n";
		print OUT "else\n";
		print OUT "  if [ ! -s $bamFile ]; then\n";
		print OUT "    if [ ! -s $samFile ]; then\n";
		print OUT "      if [ ! -s $saiFile1 ]; then\n";
		print OUT "        echo sai1=`date` \n";
		print OUT "        bwa aln $option $faFile $sampleFile1 >$saiFile1 \n";
		print OUT "      fi\n";
		print OUT "      echo aln=`date` \n";
		print OUT "      bwa samse $option_samse $faFile $saiFile1 $sampleFile1 > $samFile \n";
		print OUT "    fi\n";
		print OUT "    echo sam2bam=`date`\n";
		print OUT "    samtools view -b -S $samFile -o $bamFile\n";
		print OUT "  fi\n";
		print OUT "  echo sortbam=`date`\n";
		print OUT "  samtools sort $bamFile $sortedBamFile\n";
		print OUT "  echo bamstat=`date`\n";
		print OUT "  samtools flagstat ${sortedBamFile}.bam > ${sortedBamFile}.bam.stat\n";
		print OUT "fi\n\n";

		print OUT "echo finished=`date`\n";
		close OUT;

		print "$pbsFile created\n";
	}
	close(SH);

	if ( is_linux() ) {
		chmod 0755, $shfile;
	}

	print "!!!shell file $shfile created, you can run this shell file to submit all bwa tasks.\n";

	#`qsub $pbsFile`;
}

sub bwa_by_pbs_double {
	my ( $config, $section ) = @_;

	my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option ) = get_parameter( $config, $section );

	my $faFile = get_param_file( $config->{$section}{fasta_file}, "fasta_file", 1 );
	my $inserts = $config->{$section}{estimate_insert};

	die "define ${section}::option_sampe first" if ( !defined $config->{$section}{option_sampe} );
	my $option_sampe = $config->{$section}{option_sampe};

	my %rawFiles = %{ get_raw_files( $config, $section ) };

	my $shfile = $pbsDir . "/${task_name}.sh";
	open( SH, ">$shfile" ) or die "Cannot create $shfile";
	print SH "type -P qsub &>/dev/null && export MYCMD=\"qsub\" || export MYCMD=\"bash\" \n";

	for my $sampleName ( sort keys %rawFiles ) {
		my @sampleFiles = @{ $rawFiles{$sampleName} };

		my $sampleFile1 = $sampleFiles[0];
		my $sampleFile2 = $sampleFiles[1];

		my ( $sampleName1, $directories1, $suffix1 ) = fileparse($sampleFile1);
		my $saiFile1 = $sampleName1 . ".sai";
		my ( $sampleName2, $directories2, $suffix2 ) = fileparse($sampleFile2);
		my $saiFile2      = $sampleName2 . ".sai";
		my $samFile       = $sampleName . ".sam";
		my $bamFile       = $sampleName . ".bam";
		my $sortedBamFile = $sampleName . "_sort";

		my $pbsName = "${sampleName}_bwa.pbs";
		my $pbsFile = "${pbsDir}/$pbsName";

		print SH "\$MYCMD ./$pbsName \n";

		my $log = "${logDir}/${sampleName}_bwa.log";

		open( OUT, ">$pbsFile" ) or die $!;
		print OUT $pbsDesc;
		print OUT "#PBS -o $log\n";
		print OUT "#PBS -j oe\n\n";

		if ( -e $path_file ) {
			print OUT "source $path_file\n";
		}
		print OUT "echo bwa=`date`\n";

		my $curDir = create_directory_or_die( $resultDir . "/$sampleName" );

		#my $tag="'\@RG\tID:$sample\tLB:$sample\tSM:$sample\tPL:ILLUMINA'";
		print OUT "cd $curDir\n\n";

		print OUT "if [ -s ${sortedBamFile}.bam ]; then\n";
		print OUT "  echo job has already been done. if you want to do again, delete ${sortedBamFile}.bam and submit job again.\n";
		print OUT "else\n";
		print OUT "  if [ ! -s $bamFile ]; then\n";
		print OUT "    if [ ! -s $samFile ]; then\n";
		print OUT "      if [ ! -s $saiFile1 ]; then\n";
		print OUT "        echo sai1=`date` \n";
		print OUT "        bwa aln $option $faFile $sampleFile1 >$saiFile1 \n";
		print OUT "      fi\n";
		print OUT "      if [ ! -s $saiFile2 ]; then\n";
		print OUT "        echo sai2=`date` \n";
		print OUT "        bwa aln $option $faFile $sampleFile2 >$saiFile2 \n";
		print OUT "      fi\n";
		print OUT "      echo aln=`date` \n";
		print OUT "      bwa sampe $option_sampe $faFile $saiFile1 $saiFile2 $sampleFile1 $sampleFile2 > $samFile \n";
		print OUT "    fi\n";
		print OUT "    echo sam2bam=`date`\n";
		print OUT "    samtools view -b -S $samFile -o $bamFile\n";
		print OUT "  fi\n";
		print OUT "  echo sortbam=`date`\n";
		print OUT "  samtools sort $bamFile $sortedBamFile\n";
		print OUT "  echo bamstat=`date`\n";
		print OUT "  samtools flagstat ${sortedBamFile}.bam > ${sortedBamFile}.bam.stat\n";

		if ($inserts) {
			print OUT "  echo insertsize=`date`\n";
			print OUT "  samtools view ${sortedBamFile}.bam | awk 'and (\$2, 0x0002) && and (\$2, 0x0040)' | cut -f 9 | sed 's/^-//' > ${sortedBamFile}.len \n";
			print OUT
"  sort -n ${sortedBamFile}.len | awk ' { x[NR]=\$1; s+=\$1; } END {mean=s/NR; for (i in x){ss+=(x[i]-mean)^2}; sd=sqrt(ss/NR); if(NR %2) {median=x[(NR+1)/2];}else{median=(x[(NR/2)]+x[(NR/2)+1])/2.0;} print \"mean=\"mean \"; stdev=\"sd \"; median=\"median }' > ${sortedBamFile}.inserts \n";
		}
		print OUT "fi\n\n";

		print OUT "echo finished=`date`\n";
		close OUT;

		print "$pbsFile created\n";
	}
	close(SH);

	if ( is_linux() ) {
		chmod 0755, $shfile;
	}

	print "!!!shell file $shfile created, you can run this shell file to submit all bwa tasks.\n";

	#`qsub $pbsFile`;
}

sub samtools_index {
	my ( $config, $section ) = @_;

	my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option ) = get_parameter( $config, $section );

	my %rawFiles = %{ get_raw_files( $config, $section ) };

	my $shfile = $pbsDir . "/${task_name}.sh";
	open( SH, ">$shfile" ) or die "Cannot create $shfile";
	print SH "type -P qsub &>/dev/null && export MYCMD=\"qsub\" || export MYCMD=\"bash\" \n";

	for my $sampleName ( sort keys %rawFiles ) {
		my @sampleFiles = @{ $rawFiles{$sampleName} };

		my $pbsName = "${sampleName}_index.pbs";
		my $pbsFile = "${pbsDir}/$pbsName";

		print SH "\$MYCMD ./$pbsName \n";

		my $log = "${logDir}/${sampleName}_index.log";

		open( OUT, ">$pbsFile" ) or die $!;
		print OUT $pbsDesc;
		print OUT "#PBS -o $log\n";
		print OUT "#PBS -j oe\n\n";

		if ( -e $path_file ) {
			print OUT "source $path_file\n";
		}
		print OUT "echo index=`date`\n";

		my $bamFile = $sampleFiles[0];
		my $bamIndexFile = $bamFile . ".bai";
		print OUT "if [ ! -s $bamIndexFile ]; then\n";
		print OUT "  echo samtools_index=`date`\n";
		print OUT "  samtools index $bamFile \n";
		print OUT "fi\n";
		print OUT "echo finished=`date`\n";
		close OUT;

		print "$pbsFile created\n";
	}
	close(SH);

	if ( is_linux() ) {
		chmod 0755, $shfile;
	}

	print "!!!shell file $shfile created, you can run this shell file to submit all samtools index tasks.\n";

	#`qsub $pbsFile`;
}

1;
