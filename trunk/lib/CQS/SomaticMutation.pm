#!/usr/bin/perl
package CQS::Samtools;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::ConfigUtils;

require Exporter;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(call_wsmdetector)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

use Cwd;

sub call_wsmdetector {
  my ( $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option ) = get_parameter( $config, $section );

  my $wsmfile = get_param_file( $config->{$section}{execute_file}, "execute_file", 1 );
  my $fafile = get_param_file( $config->{$section}{reference_sequence}, "reference_sequence", 1 );
  my $minimumEventCount = $config->{$section}{minimum_event_count};
  my $minimumMappingQuality = $config->{$section}{minimum_mapping_quality};

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $shfile = $pbsDir . "/${task_name}.submit";
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  print SH "type -P qsub &>/dev/null && export MYCMD=\"qsub\" || export MYCMD=\"bash\" \n";

  for my $sampleName ( sort keys %rawFiles ) {
    my @sampleFiles = @{ $rawFiles{$sampleName} };

    my $pbsName = "${sampleName}_wsmdetector.pbs";
    my $pbsFile = "${pbsDir}/$pbsName";

    print SH "\$MYCMD ./$pbsName \n";

    my $log = "${logDir}/${sampleName}_wsmdetector.log";

    open( OUT, ">$pbsFile" ) or die $!;
    print OUT $pbsDesc;
    print OUT "#PBS -o $log\n";
    print OUT "#PBS -j oe\n\n";

    if ( defined $path_file ) {
      if ( -e $path_file ) {
        print OUT "source $path_file \n\n";
      }
    }
    
    print OUT "cd $resultDir \n\n";
    print OUT "echo wsmdetector=`date` \n\n";

    my $sampleCount = scalar(@sampleFiles);

    for my $sampleFile (@sampleFiles) {
      my $bamindex = $sampleFile . ".bai";

      print OUT "if [ ! -s $bamindex ]; \n";
      print OUT "then \n";
      print OUT "  samtools index $sampleFile \n";
      print OUT "fi \n\n";
    }
    
    print OUT "mono $wsmfile -g $fafile";
    if ( defined $minimumEventCount ) {
      if ( $minimumEventCount > 0 ) {
        print OUT " -c $minimumEventCount"
      }
    }
    
    if ( defined $minimumMappingQuality ) {
      if ( $minimumMappingQuality > 0 ) {
        print OUT " -q $minimumMappingQuality"
      }
    }
    
    my $first = 1;
    for my $sampleFile (@sampleFiles) {
      if($first){
        print OUT " -b $sampleFile";
        $first = 0;
      }
      else{
        print OUT ",$sampleFile";
      }
    }

    print OUT " -o $resultDir > ${sampleName}.summary \n\n";
    print OUT "echo finished=`date` \n";
    close OUT;

    print "$pbsFile created \n";
  }

  close(SH);

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all samtools mpileup tasks.\n";
}

1;
