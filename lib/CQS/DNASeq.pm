#!/usr/bin/perl
package CQS::DNASeq;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;

require Exporter;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(bwa_by_pbs_single bwa_by_pbs_double)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

use Cwd;

sub bwa_by_pbs_single {
	my ( $faFile, $sampleFile, $rootDir ) = @_;

	my ( $sampleName, $directories, $suffix ) = fileparse($sampleFile);

	$sampleName =~ s/\.(\w+)$//;

	my $pathFile = '/data/cqs/bin/path.txt';

	my ( $logDir, $pbsDir, $resultDir ) = init_dir($rootDir);

	my ($pbsDesc) = get_pbs_desc();

	my $saiFile       = $sampleName . ".sai";
	my $samFile       = $sampleName . ".sam";
	my $bamFile       = $sampleName . ".bam";
	my $sortedBamFile = $sampleName . "_sort";
	my $log           = $logDir . "/" . $sampleName . ".log";

	if ( -e $resultDir . "/" . $sortedBamFile ) {
		next;
	}

	#my $tag="'\@RG\tID:$sample\tLB:$sample\tSM:$sample\tPL:ILLUMINA'";
	my $pbsFile = $pbsDir . "/${sampleName}.pbs";
	print "$pbsDir\n";
	open( OUT, ">$pbsFile" ) or die $!;
	print OUT $pbsDesc;
	print OUT "#PBS -o $log\n";
	print OUT "#PBS -j oe\n\n";
	print OUT "source $pathFile\n";
	print OUT "cd $resultDir\n\n";

	if ( !( -e $resultDir . "/" . $bamFile ) ) {
		if ( !( -e $resultDir . "/" . $samFile ) ) {
			if ( !( -e $resultDir . "/" . $sampleFile ) ) {
				print OUT "echo sai=`date` \n";
				print OUT "bwa aln -q 15 $faFile $sampleFile >$saiFile \n";
			}
			print OUT "echo aln=`date` \n";
			print OUT "bwa sampe -n 3 $faFile $saiFile $sampleFile > $samFile \n";
		}
		print OUT "echo sam2bam=`date`\n";
		print OUT "samtools view -b -S $samFile -o $bamFile\n";
	}
	print OUT "echo sortbam=`date`\n";
	print OUT "samtools sort $bamFile $sortedBamFile\n";
	print OUT "echo bamstat=`date`\n";
	print OUT "samtools flagstat $sortedBamFile.bam > $sortedBamFile.bam.stat\n";
	print OUT "echo finished=`date`\n";
	close OUT;

	#`qsub $pbsFile`;
}

sub bwa_by_pbs_double {
	my ( $config, $section ) = @_;

	my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option ) = get_parameter( $config, $section );

	my $faFile = get_param_file( $config->{$section}{fasta_file}, "fasta_file", 1 );
	my $inserts = $config->{$section}{estimate_insert};

	my %rawFiles = %{ get_raw_files( $config, $section ) };

	my $shfile = $pbsDir . "/${task_name}.sh";
	open( SH, ">$shfile" ) or die "Cannot create $shfile";

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

		print SH "qsub ./$pbsName \n";

		my $log = "${logDir}/${sampleName}_bwa.log";

		open( OUT, ">$pbsFile" ) or die $!;
		print OUT $pbsDesc;
		print OUT "#PBS -o $log\n";
		print OUT "#PBS -j oe\n\n";

		if ( -e $path_file ) {
			print OUT "source $path_file\n";
		}
		print OUT "echo bwa=`date`\n";

		#my $tag="'\@RG\tID:$sample\tLB:$sample\tSM:$sample\tPL:ILLUMINA'";
		print OUT "cd $resultDir\n\n";

		print OUT "if [ -s ${sortedBamFile}.bam ]; then\n";
		print OUT "  echo job has already been done. if you want to do again, delete ${sortedBamFile}.bam and submit job again.\n";
		print OUT "else\n";
		print OUT "  if [ ! -s $bamFile ]; then\n";
		print OUT "    if [ ! -s $samFile ]; then\n";
		print OUT "      if [ ! -s $saiFile1 ]; then\n";
		print OUT "        echo sai1=`date` \n";
		print OUT "        bwa aln -q 15 $faFile $sampleFile1 >$saiFile1 \n";
		print OUT "      fi\n";
		print OUT "      if [ ! -s $saiFile2 ]; then\n";
		print OUT "        echo sai2=`date` \n";
		print OUT "        bwa aln -q 15 $faFile $sampleFile2 >$saiFile2 \n";
		print OUT "      fi\n";
		print OUT "      echo aln=`date` \n";
		print OUT "      bwa sampe -n 3 $faFile $saiFile1 $saiFile2 $sampleFile1 $sampleFile2 > $samFile \n";
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
			print OUT "  sort -n ${sortedBamFile}.len | awk ' { x[NR]=\$1; s+=$1; } END {mean=s/NR; for (i in x){ss+=(x[i]-mean)^2}; sd=sqrt(ss/NR); if(NR %2) {median=x[(NR+1)/2];}else{median=(x[(NR/2)]+x[(NR/2)+1])/2.0;} print \"mean=\"mean \"; stdev=\"sd; \"; median=\"median }' \n";
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

1;
