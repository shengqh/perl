#!/usr/bin/perl
package CQS::SomaticMutation;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::ConfigUtils;

require Exporter;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(rsmc)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

use Cwd;

sub rsmc {
  my ( $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option,$sh_direct ) = get_parameter( $config, $section );

  my $rsmcfile = get_param_file( $config->{$section}{execute_file}, "execute_file", 1 );
  my $source_type = $config->{$section}{source_type} or die "source_type is not defined in $section";

  my $annovarBuildver = $config->{$section}{annovar_buildver} or die "annovar_buildver is not defined in $section";
  $option = $option . " --annovar --annovar_buildver $annovarBuildver ";

  my $rnaediting_db = $config->{$section}{rnaediting_db};
  if(defined $rnaediting_db) {
    $option = $option . " --rnaediting --rnaediting_db $rnaediting_db ";
  }

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $mpileupfile      = "";
  my $fafile           = "";
  my $mpileupParameter = "";
  my $isbam            = lc($source_type) eq "bam";
  if ($isbam) {
    $fafile = get_param_file( $config->{$section}{mpileup_sequence}, "mpileup_sequence", 1 );
    $mpileupParameter = $config->{$section}{mpileup_option};
    if ( defined $mpileupParameter ) {
      if ( $mpileupParameter eq "" ) {
        undef($$mpileupParameter);
      }
    }
  }
  else {
    $mpileupfile = get_param_file( $config->{$section}{mpileup_file}, "mpileup_file", 1 );
  }

  my $shfile = $pbsDir . "/${task_name}.submit";
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  if ($sh_direct) {
    print SH "export MYCMD=\"bash\" \n";
  }
  else {
    print SH "type -P qsub &>/dev/null && export MYCMD=\"qsub\" || export MYCMD=\"bash\" \n";
  }

  for my $sampleName ( sort keys %rawFiles ) {
    my @sampleFiles = @{ $rawFiles{$sampleName} };

    my $pbsName = "rsmc_${sampleName}.pbs";
    my $pbsFile = "${pbsDir}/$pbsName";

    print SH "\$MYCMD ./$pbsName \n";

    my $log = "${logDir}/rsmc_${sampleName}.log";

    open( OUT, ">$pbsFile" ) or die $!;
    print OUT "$pbsDesc
#PBS -o $log
#PBS -j oe

$path_file 

echo rsmc=`date` 
";

    my $sampleCount = scalar(@sampleFiles);
    my $curDir      = create_directory_or_die( $resultDir . "/$sampleName" );

    if ($isbam) {
      for my $sampleFile (@sampleFiles) {
        my $bamindex = $sampleFile . ".bai";

        print OUT "if [ ! -s $bamindex ]; \n";
        print OUT "then \n";
        print OUT "  samtools index $sampleFile \n";
        print OUT "fi \n\n";
      }

      if ( defined $mpileupParameter ) {
        print OUT "samtools mpileup -f $fafile $mpileupParameter";
        for my $sampleFile (@sampleFiles) {
          print OUT " $sampleFile";
        }
        print OUT " | mono $rsmcfile all -t console $option";
      }
      else {
        print OUT "mono $rsmcfile all -t bam -f $fafile $option";

        my $first = 1;
        for my $sampleFile (@sampleFiles) {
          if ($first) {
            print OUT " -b $sampleFile";
            $first = 0;
          }
          else {
            print OUT ",$sampleFile";
          }
        }
      }
    }
    else {
      print OUT "mono $rsmcfile all -t mpileup -m $mpileupfile $option";
    }

    print OUT " -o $curDir \n\n";
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

sub muTect {
  my ( $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option,$sh_direct ) = get_parameter( $config, $section );

  my $rsmcfile = get_param_file( $config->{$section}{execute_file}, "execute_file", 1 );
  my $source_type = $config->{$section}{source_type} or die "source_type is not defined in $section";

  my $annovarBuildver = $config->{$section}{annovar_buildver} or die "annovar_buildver is not defined in $section";
  $option = $option . " --annovar --annovar_buildver $annovarBuildver ";

  my $rnaediting_db = $config->{$section}{rnaediting_db};
  if(defined $rnaediting_db) {
    $option = $option . " --rnaediting --rnaediting_db $rnaediting_db ";
  }

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $mpileupfile      = "";
  my $fafile           = "";
  my $mpileupParameter = "";
  my $isbam            = lc($source_type) eq "bam";
  if ($isbam) {
    $fafile = get_param_file( $config->{$section}{mpileup_sequence}, "mpileup_sequence", 1 );
    $mpileupParameter = $config->{$section}{mpileup_option};
    if ( defined $mpileupParameter ) {
      if ( $mpileupParameter eq "" ) {
        undef($$mpileupParameter);
      }
    }
  }
  else {
    $mpileupfile = get_param_file( $config->{$section}{mpileup_file}, "mpileup_file", 1 );
  }

  my $shfile = $pbsDir . "/${task_name}.submit";
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  if ($sh_direct) {
    print SH "export MYCMD=\"bash\" \n";
  }
  else {
    print SH "type -P qsub &>/dev/null && export MYCMD=\"qsub\" || export MYCMD=\"bash\" \n";
  }

  for my $sampleName ( sort keys %rawFiles ) {
    my @sampleFiles = @{ $rawFiles{$sampleName} };

    my $pbsName = "rsmc_${sampleName}.pbs";
    my $pbsFile = "${pbsDir}/$pbsName";

    print SH "\$MYCMD ./$pbsName \n";

    my $log = "${logDir}/rsmc_${sampleName}.log";

    open( OUT, ">$pbsFile" ) or die $!;
    print OUT "$pbsDesc
#PBS -o $log
#PBS -j oe

$path_file 

echo rsmc=`date` 
";

    my $sampleCount = scalar(@sampleFiles);
    my $curDir      = create_directory_or_die( $resultDir . "/$sampleName" );

    if ($isbam) {
      for my $sampleFile (@sampleFiles) {
        my $bamindex = $sampleFile . ".bai";

        print OUT "if [ ! -s $bamindex ]; \n";
        print OUT "then \n";
        print OUT "  samtools index $sampleFile \n";
        print OUT "fi \n\n";
      }

      if ( defined $mpileupParameter ) {
        print OUT "samtools mpileup -f $fafile $mpileupParameter";
        for my $sampleFile (@sampleFiles) {
          print OUT " $sampleFile";
        }
        print OUT " | mono $rsmcfile all -t console $option";
      }
      else {
        print OUT "mono $rsmcfile all -t bam -f $fafile $option";

        my $first = 1;
        for my $sampleFile (@sampleFiles) {
          if ($first) {
            print OUT " -b $sampleFile";
            $first = 0;
          }
          else {
            print OUT ",$sampleFile";
          }
        }
      }
    }
    else {
      print OUT "mono $rsmcfile all -t mpileup -m $mpileupfile $option";
    }

    print OUT " -o $curDir \n\n";
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
