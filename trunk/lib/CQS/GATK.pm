#!/usr/bin/perl
package CQS::GATK;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;

require Exporter;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(RealignerTargetCreator IndelRealigner)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

use Cwd;

sub RealignerTargetCreator {
  my ( $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option ) = get_parameter( $config, $section );

  my $faFile = get_param_file( $config->{$section}{fasta_file}, "fasta_file", 1 );
  my $vcfFile = get_param_file( $config->{$section}{vcf_file}, "vcf_file", 1 );
  my $gatk_jar = get_param_file( $config->{$section}{gatk_jar}, "gatk_jar", 1 );
    my $markDuplicates_jar = get_param_file( $config->{$section}{markDuplicates_jar}, "markDuplicates_jar", 1 );

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $shfile = $pbsDir . "/${task_name}.sh";
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  print SH "type -P qsub &>/dev/null && export MYCMD=\"qsub\" || export MYCMD=\"bash\" \n";

  for my $sampleName ( sort keys %rawFiles ) {
    my @sampleFiles = @{ $rawFiles{$sampleName} };

    my $sampleFile1 = $sampleFiles[0];
    my $intervalFile = $sampleFile1 . ".intervals";
    my $realignedFile = change_extension($sampleFile1,".realigned.bam");
    my $redupFile = change_extension($realignedFile, ".redup.bam");

    my $pbsName = "${sampleName}_refine.pbs";
    my $pbsFile = "${pbsDir}/$pbsName";

    print SH "\$MYCMD ./$pbsName \n";

    my $log = "${logDir}/${sampleName}_refine.log";

    open( OUT, ">$pbsFile" ) or die $!;
    print OUT $pbsDesc;
    print OUT "#PBS -o $log\n";
    print OUT "#PBS -j oe\n\n";

    if ( -e $path_file ) {
      print OUT "source $path_file\n";
    }
    print OUT "echo bwa=`date`\n";

    my $curDir = create_directory_or_die( $resultDir . "/$sampleName" );

    print OUT "cd $curDir\n";
    print OUT "mkdir tmpdir\n\n";

    print OUT "echo RealignerTargetCreator=`date` \n";
    print OUT "java $option -jar $gatk_jar -I $sampleFile1 -R $faFile -T RealignerTargetCreator -o $intervalFile --known $vcfFile \n\n";
    
    print OUT "echo IndelRealigner=`date` \n";
    print OUT "java $option -Djava.io.tmpdir=tmpdir -jar $gatk_jar -I $sampleFile1 -R $faFile -T IndelRealigner -targetIntervals $intervalFile -o $realignedFile -known $vcfFile --consensusDeterminationModel KNOWNS_ONLY -LOD 0.4 \n\n";
    
    print OUT "echo RemoveDuplicate=`date` \n";
    print OUT "java $option -jar $markDuplicates_jar I=$realignedFile O=$redupFile M=duplicateMetricsFile VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true REMOVE_DUPLICATES=true \n\n";

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


sub IndelRealigner {
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

1;
