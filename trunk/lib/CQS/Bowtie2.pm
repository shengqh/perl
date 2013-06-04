#!/usr/bin/perl
package CQS::Bowtie2;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::Task;
use CQS::NGSCommon;

our @ISA = qw(CQS::Task);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name} = "Bowtie2";
  bless $self, $class;
  return $self;
}

sub generateScript {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $faFile = get_param_file( $config->{$section}{fasta_file}, "fasta_file", 1 );
  my $bowtie2_index = $config->{$section}{bowtie2_index} or die "define ${section}::bowtie2_index first";
  my $samonly = $config->{$section}{samonly};
  if ( !defined $samonly ) {
    $samonly = 0;
  }

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $shfile = $pbsDir . "/${task_name}.sh";
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  if ($sh_direct) {
    print SH "export MYCMD=\"bash\" \n";
  }
  else {
    print SH "type -P qsub &>/dev/null && export MYCMD=\"qsub\" || export MYCMD=\"bash\" \n";
  }

  for my $sampleName ( sort keys %rawFiles ) {
    my @sampleFiles = @{ $rawFiles{$sampleName} };
    my $samFile     = $sampleName . ".sam";
    my $bamFile     = $sampleName . ".bam";

    my $indent = "";
    my $tag    = "--sam-RG ID:$sampleName --sam-RG LB:$sampleName --sam-RG SM:$sampleName --sam-RG PL:ILLUMINA";

    my $fastqs = join( ',', @sampleFiles );
    my $bowtie2_aln_command = "bowtie2 $option -x $bowtie2_index -U $fastqs $tag";

    my ( $bamSortedFile, $bamSortedPrefix ) = get_sorted_bam( $bamFile, $indent );

    my $sam2bam_command = get_sam2bam_command( $samFile, $bamFile, $indent );
    my $sort_index_command = get_sort_index_command( $bamFile, $bamSortedPrefix, $indent );
    my $stat_command = get_stat_command( $bamSortedFile, $indent );

    my $pbsName = "${sampleName}_bowtie2.pbs";
    my $pbsFile = "${pbsDir}/$pbsName";
    my $curDir  = create_directory_or_die( $resultDir . "/$sampleName" );
    my $log     = "${logDir}/${sampleName}_bowtie2.log";

    print SH "\$MYCMD ./$pbsName \n";

    open( OUT, ">$pbsFile" ) or die $!;
    print OUT "$pbsDesc
#PBS -o $log
#PBS -j oe

$path_file

cd $curDir
";
    if ($samonly) {
      print OUT "
if [ -s $samFile ]; then
  echo job has already been done. if you want to do again, delete $samFile and submit job again.
  exit 0
fi

$bowtie2_aln_command -S $samFile
";
    }
    else {
      print OUT "
if [ -s $bamSortedFile ]; then
  echo job has already been done. if you want to do again, delete $bamSortedFile and submit job again.
  exit 0
fi

$bowtie2_aln_command | samtools view -S -b - | samtools sort - $bamSortedPrefix

$stat_command
";
    }

    print OUT "
echo finished=`date`

exit 1;
";

    close OUT;

    print "$pbsFile created\n";
  }
  close(SH);

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all bwa tasks.\n";
}

sub getExpectResult {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $samonly = $config->{$section}{samonly};
  if ( !defined $samonly ) {
    $samonly = 0;
  }

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $result = {};
  for my $sampleName ( sort keys %rawFiles ) {
    my $samFile = $sampleName . ".sam";
    my $curDir  = $resultDir . "/$sampleName";

    my $finalFile;
    if ($samonly) {
      $finalFile = $samFile;
    }
    else {
      my $bamFile = $sampleName . ".bam";
      my ( $bamSortedFile, $bamSortedPrefix ) = get_sorted_bam($bamFile);
      $finalFile = $bamSortedFile;
    }
    my @resultFiles = ();
    push( @resultFiles, $curDir . "/" . $finalFile );

    $result->{$sampleName} = \@resultFiles;
  }

  return $result;
}

1;
