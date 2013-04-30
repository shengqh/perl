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

our %EXPORT_TAGS = ( 'all' => [qw(refine_bam_file)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

use Cwd;

sub refine_bam_file {
	my ( $config, $section ) = @_;

	my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option ) = get_parameter( $config, $section );

	my $faFile             = get_param_file( $config->{$section}{fasta_file},         "fasta_file",         1 );
	my $vcfFile            = get_param_file( $config->{$section}{vcf_file},           "vcf_file",           1 );
	my $gatk_jar           = get_param_file( $config->{$section}{gatk_jar},           "gatk_jar",           1 );
	my $markDuplicates_jar = get_param_file( $config->{$section}{markDuplicates_jar}, "markDuplicates_jar", 1 );
	my $thread_count       = $config->{$section}{thread_count};
	if ( !defined($thread_count) ) {
		$thread_count = 1;
	}

	my %rawFiles = %{ get_raw_files( $config, $section ) };

	my $shfile = $pbsDir . "/${task_name}.sh";
	open( SH, ">$shfile" ) or die "Cannot create $shfile";
	print SH "type -P qsub &>/dev/null && export MYCMD=\"qsub\" || export MYCMD=\"bash\" \n";

	for my $sampleName ( sort keys %rawFiles ) {
		my @sampleFiles = @{ $rawFiles{$sampleName} };

		my $sampleFile1   = $sampleFiles[0];
		my $redupFile     = change_extension( $sampleFile1, ".redup.bam" );
		my $intervalFile  = $redupFile . ".intervals";
		my $realignedFile = change_extension( $redupFile, ".realigned.bam" );
		my $grpFile       = $realignedFile . ".grp";
		my $recalFile     = change_extension( $sampleFile1, ".recal.bam" );

		my $pbsName = "${sampleName}_refine.pbs";
		my $pbsFile = "${pbsDir}/$pbsName";

		print SH "\$MYCMD ./$pbsName \n";

		my $log    = "${logDir}/${sampleName}_refine.log";
		my $curDir = create_directory_or_die( $resultDir . "/$sampleName" );

		open( OUT, ">$pbsFile" ) or die $!;

		print OUT "
$pbsDesc
#PBS -o $log
#PBS -j oe
";

		if ( -e $path_file ) {
			print OUT "source $path_file\n";
		}

		print OUT "
echo bwa=`date`
cd $curDir

if [ ! -s tmpdir ]; then
  mkdir tmpdir
fi

if [ ! -e $redupFile ]; then
  echo RemoveDuplicate=`date` 
  java $option -jar $markDuplicates_jar I=$sampleFile1 O=$redupFile M=${redupFile}.matrix VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true REMOVE_DUPLICATES=true 
fi

if [ -e $redupFile && ! -e $intervalFile ]; then
  echo RealignerTargetCreator=`date` 
  java $option -jar $gatk_jar -I $redupFile -R $faFile -T RealignerTargetCreator -o $intervalFile --known $vcfFile -nt $thread_count
fi

if [ -e $intervalFile && ! -e $realignedFile ]; then
  echo IndelRealigner=`date` 
  java $option -Djava.io.tmpdir=tmpdir -jar $gatk_jar -I $sampleFile1 -R $faFile -T IndelRealigner -targetIntervals $intervalFile -o $realignedFile -known $vcfFile --consensusDeterminationModel KNOWNS_ONLY -LOD 0.4 
fi

if [ -e $realignedFile && ! -e $grpFile ]; then
  echo BaseRecalibrator=`date` 
  java $option -jar $gatk_jar -R $faFile -I $realignedFile -T BaseRecalibrator -knownSites $vcfFile -o $grpFile
fi

if [ -e $grpFile && ! -e $recalFile ]; then
  echo TableRecalibration=`date`
  java $option -jar $gatk_jar -R $faFile -I $realignedFile -T PrintReads -o $recalFile -BQSR $grpFile
fi

echo finished=`date`
";

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
