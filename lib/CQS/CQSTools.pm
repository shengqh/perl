#!/usr/bin/perl
package CQS::CQSTools;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;

require Exporter;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(mirna_count)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

use Cwd;

#get expected alignment result based on alignment definition
sub get_bam_map {
  my ( $config, $section ) = @_;

  my ( $result, $issource ) = get_raw_files2( $config, $section );
  if ($issource) {
    retun $result;
  }

  my $alignsection = $config->{$section}{source_ref};
  my $align_dir = $config->{$alignsection}{target_dir} or die "${$alignsection}::target_dir not defined.";
  my ( $logDir, $pbsDir, $resultDir ) = init_dir( $align_dir, 0 );
  my %fqFiles = %{$result};

  my $alignresult = {};
  for my $sampleName ( keys %fqFiles ) {
    my $sam = "${resultDir}/${sampleName}/${sampleName}_sorted.bam";
    $alignresult->{$sampleName} = $sam;
  }
  return $alignresult;
}

sub mirna_count {
  my ( $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $cqsFile = get_param_file( $config->{$section}{cqs_tools}, "cqs_tools", 1 );
  my $gffFile = get_param_file( $config->{$section}{gff_file},  "gff_file",  1 );

  my %rawFiles = %{ get_bam_map( $config, $section ) };
  my $fasta_format = $config->{$section}{fasta_format};
  if ( !defined $fasta_format ) {
    $fasta_format = 0;
  }

  if ($fasta_format) {
    $option = $option . " --fasta";
  }

  my $shfile = $pbsDir . "/${task_name}_count.sh";
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  if ($sh_direct) {
    print SH "export MYCMD=\"bash\" \n";
  }
  else {
    print SH "type -P qsub &>/dev/null && export MYCMD=\"qsub\" || export MYCMD=\"bash\" \n";
  }

  for my $sampleName ( sort keys %rawFiles ) {
    my $bamFile = $rawFiles{$sampleName};

    my $curDir   = dirname($bamFile);
    my $fileName = basename($bamFile);

    my $countFile = $fileName . ".count";

    my $pbsName = "${sampleName}_count.pbs";
    my $pbsFile = "${pbsDir}/$pbsName";

    print SH "\$MYCMD ./$pbsName \n";

    my $log = "${logDir}/${sampleName}_count.log";

    open( OUT, ">$pbsFile" ) or die $!;
    print OUT "$pbsDesc
#PBS -o $log
#PBS -j oe

$path_file

cd $curDir

if [ -s $countFile ]; then
  echo job has already been done. if you want to do again, delete $countFile and submit job again.
  exit 0
fi

mono-sgen $cqsFile mirna_count $option -i $fileName -g $gffFile

echo finished=`date`

exit 1 
";

    close OUT;

    print "$pbsFile created \n";
  }
  close(SH);

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all bwa tasks.\n";

  #`qsub $pbsFile`;
}

1;
