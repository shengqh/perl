#!/usr/bin/perl
package CQS::IdenticalQueryBuilder;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::Task;

our @ISA = qw(CQS::Task);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name} = "IdenticalQueryBuilder";
  bless $self, $class;
  return $self;
}

sub generateScript {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $cqstools  = $config->{$section}{cqstools}  or die "define ${section}::cqstools first";
  my $extension = $config->{$section}{extension} or die "define ${section}::extension first";

  my $merge_result = $config->{$section}{merge_result};
  if ( !defined $merge_result ) {
    $merge_result = 0;
  }
  my $minlen = $config->{$section}{minlen};
  if ( defined $minlen ) {
    $minlen = "-l $minlen";
  }
  else {
    $minlen = "";
  }

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $shfile = $pbsDir . "/${task_name}_IQB.sh";
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  if ($sh_direct) {
    print SH "export MYCMD=\"bash\" \n";
  }
  else {
    print SH "type -P qsub &>/dev/null && export MYCMD=\"qsub\" || export MYCMD=\"bash\" \n";
  }

  for my $sampleName ( sort keys %rawFiles ) {
    my @sampleFiles = @{ $rawFiles{$sampleName} };
    my $finalFile   = $sampleName . $extension;

    my $pbsName = "${sampleName}_IQB.pbs";
    my $pbsFile = "${pbsDir}/$pbsName";

    print SH "\$MYCMD ./$pbsName \n";

    my $log = "${logDir}/${sampleName}_IQB.log";

    open( OUT, ">$pbsFile" ) or die $!;
    print OUT "$pbsDesc
#PBS -o $log
#PBS -j oe

$path_file

cd $resultDir

";
    if ( scalar(@sampleFiles) == 1 ) {
      print OUT "mono $cqstools fastq_identical -i $sampleFiles[0] $minlen -o $finalFile \n";
    }
    else {
      my $outputFiles = "";
      for my $sampleFile (@sampleFiles) {
        my $fileName = basename($sampleFile);
        my $outputFile = change_extension( $fileName, $extension );
        $outputFiles = $outputFiles . " " . $outputFile;
        print OUT "mono $cqstools fastq_identical -i $sampleFiles[0] $minlen -o $outputFile \n";
      }

      if ($merge_result) {
        print OUT "
cat $outputFiles > $finalFile
 
rm $outputFiles
";
      }
    }
    print OUT "
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

sub getExpectResult {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $extension = $config->{$section}{extension} or die "define ${section}::extension first";
  my $merge_result = $config->{$section}{merge_result};
  if ( !defined $merge_result ) {
    $merge_result = 0;
  }

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $result = {};
  for my $sampleName ( keys %rawFiles ) {
    my @sampleFiles = @{ $rawFiles{$sampleName} };
    my $finalFile   = $sampleName . $extension;
    my @resultFiles = ();

    if ( scalar(@sampleFiles) == 1 || $merge_result ) {
      push( @resultFiles, $resultDir . "/" . $finalFile );
    }
    else {
      for my $sampleFile (@sampleFiles) {
        my $fileName = basename($sampleFile);
        my $outputFile = change_extension( $fileName, $extension );
        push( @resultFiles, $resultDir . "/" . $outputFile );
      }
    }

    $result->{$sampleName} = \@resultFiles;
  }
  return $result;
}

1;
