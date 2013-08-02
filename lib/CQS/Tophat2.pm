#!/usr/bin/perl
package CQS::Tophat2;

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
  $self->{_name} = "Tophat2";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  $option = $option . " --keep-fasta-order";

  my $bowtie2_index = $config->{$section}{bowtie2_index} or die "define ${section}::bowtie2_index first";
  my %fqFiles = %{ get_raw_files( $config, $section ) };

  my $transcript_gtf = get_param_file( $config->{$section}{transcript_gtf}, "${section}::transcript_gtf", 0 );
  my $transcript_gtf_index;
  if ( defined $transcript_gtf ) {
    $transcript_gtf_index = $config->{$section}{transcript_gtf_index};
  }
  elsif ( defined $config->{general}{transcript_gtf} ) {
    $transcript_gtf = get_param_file( $config->{general}{transcript_gtf}, "general::transcript_gtf", 1 );
    $transcript_gtf_index = $config->{general}{transcript_gtf_index};
  }

  my $has_gtf_file   = file_exists($transcript_gtf);
  my $has_index_file = transcript_gtf_index_exists($transcript_gtf_index);

  if ( $has_gtf_file && !defined $transcript_gtf_index ) {
    die "transcript_gtf was defined but transcript_gtf_index was not defined, you should defined transcript_gtf_index to cache the parsing result.";
  }

  if ( defined $transcript_gtf_index && !$has_index_file ) {
    die "transcript_gtf_index $transcript_gtf_index defined but not exists!";
  }

  my $shfile = $pbsDir . "/${task_name}.submit";
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  print SH get_run_command($sh_direct);

  for my $sampleName ( sort keys %fqFiles ) {
    my @sampleFiles = @{ $fqFiles{$sampleName} };
    my $samples = join( " ", @sampleFiles );

    my $pbsName = "${sampleName}_th2.pbs";
    my $pbsFile = $pbsDir . "/$pbsName";
    my $log     = $logDir . "/${sampleName}_th2.log";

    my $curDir      = create_directory_or_die( $resultDir . "/$sampleName" );
    my $rgline      = "--rg-id $sampleName --rg-sample $sampleName --rg-library $sampleName";
    my $tophat2file = "accepted_hits.bam";

    my $gtfstr = "";
    if ($has_gtf_file) {
      $gtfstr = "-G $transcript_gtf --transcriptome-index=$transcript_gtf_index";
    }
    elsif ($has_index_file) {
      $gtfstr = "--transcriptome-index=$transcript_gtf_index";
    }

    open( OUT, ">$pbsFile" ) or die $!;

    print OUT "$pbsDesc
#PBS -o $log
#PBS -j oe

$path_file

cd $curDir 

if [ -s $tophat2file ]; then
  echo job has already been done. if you want to do again, delete ${curDir}/${tophat2file} and submit job again.
  exit 1;
fi

tophat2 $option $rgline $gtfstr -o . $bowtie2_index $samples

samtools index $tophat2file

1;
";
    close(OUT);

    print SH "\$MYCMD ./$pbsName \n";
    print "$pbsFile created\n";
  }
  print SH "exit 0\n";
  close(SH);

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all tasks.\n";
}

sub result {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $result = {};
  for my $sampleName ( keys %rawFiles ) {
    my $finalFile   = $resultDir . "/" . $sampleName . "/accepted_hits.bam";
    my @resultFiles = ();
    push( @resultFiles, $finalFile );
    $result->{$sampleName} = \@resultFiles;
  }
  return $result;
}

1;
