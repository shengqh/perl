#!/usr/bin/perl
package CQS::VarScan2;

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
  $self->{_name} = "VarScan2";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $varscan2_jar = get_param_file( $config->{$section}{VarScan2_jar},  "VarScan2_jar",  1 );
  my $faFile     = get_param_file( $config->{$section}{fasta_file},  "fasta_file",  1 );

  my $min_coverage    = $config->{$section}{min_coverage}    or die "min_coverage is not defined in $section!";
  my $somatic_p_value = $config->{$section}{somatic_p_value} or die "somatic_p_value is not defined in $section!";

  my $annovarParameter = $config->{$section}{annovar_param} or die "annovar_param is not defined in $section";
  $option = $option . " " . $annovarParameter;

  my $annovarDB = $config->{$section}{annovar_db} or die "annovar_db is not defined in $section";

  my $rawFiles = get_raw_files( $config, $section );
  my $groups = get_raw_files( $config, $section, "groups" );
  my %group_sample_map = ();
  for my $groupName ( sort keys %{$groups} ) {
    my @samples = @{ $groups->{$groupName} };
    my @gfiles  = ();
    my $index   = 0;
    foreach my $sampleName (@samples) {
      my @bamFiles = @{ $rawFiles->{$sampleName} };
      push( @gfiles, $bamFiles[0] );
    }
    $group_sample_map{$groupName} = \@gfiles;
  }

  my $shfile = $pbsDir . "/${task_name}_mt.submit";
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  if ($sh_direct) {
    print SH "export MYCMD=\"bash\" \n";
  }
  else {
    print SH "type -P qsub &>/dev/null && export MYCMD=\"qsub\" || export MYCMD=\"bash\" \n";
  }

  $pbsDesc =~ /\=(\d+)gb/;
  my $gb = $1;

  for my $groupName ( sort keys %group_sample_map ) {
    my @sampleFiles = @{ $group_sample_map{$groupName} };
    my $sampleCount = scalar(@sampleFiles);

    if ( $sampleCount != 2 ) {
      die "SampleFile should be normal,tumor paired.";
    }

    my $curDir      = create_directory_or_die( $resultDir . "/$groupName" );

    my $normal = $sampleFiles[0];
    my $tumor  = $sampleFiles[1];
    
    for my $sampleFile (@sampleFiles) {
      my $bamindex = $sampleFile . ".bai";
      if ( !-s $bamindex ) {
        die "bam file must be indexed : $sampleFile";
      }
    }

    my $normalfile = basename($normal);
    my $tumorfile  = basename($tumor);

    my $normal_mpileup = "${normalfile}.mpileup";
    my $tumor_mpileup  = "${tumorfile}.mpileup";

    my $snpvcf = "${groupName}.snp.vcf";
    my $passinput = "${groupName}.somatic.pass.avinput";
    my $annovar   = "${groupName}.somatic.pass.annovar";

    my $pbsName = "${groupName}_vs2.pbs";
    my $pbsFile = "${pbsDir}/$pbsName";

    print SH "\$MYCMD ./$pbsName \n";

    my $log = "${logDir}/${groupName}_vs2.log";

    open( OUT, ">$pbsFile" ) or die $!;
    print OUT "$pbsDesc
#PBS -o $log
#PBS -j oe

$path_file 

echo varscan2=`date` 

cd $curDir

if [ ! -s $normal_mpileup ]; then
  echo NORMAL_MPILEUP=`date`
  samtools mpileup -q 20 -f $faFile $normal > $normal_mpileup
fi

if [ ! -s $tumor_mpileup ]; then
  echo TUMOR_MPILEUP=`date`
  samtools mpileup -q 20 -f $faFile $tumor > $tumor_mpileup
fi

if [ ! -s $snpvcf ]; then
  java -Xmx${gb}g -jar $varscan2_jar somatic $option $normal_mpileup $tumor_mpileup $groupName --output-vcf --somatic-p-value $somatic_p_value --min-coverage $min_coverage --strand-filter
fi

java -Xmx${gb}g -jar $varscan2_jar processSomatic $snpvcf --p-value $somatic_p_value

convert2annovar.pl -format vcf4 ${snpvcf}.Somatic.hc -includeinfo > $passinput

summarize_annovar.pl $annovarParameter --outfile $annovar $passinput $annovarDB

echo finished=`date` \n";
    close OUT;

    print "$pbsFile created \n";
  }

  close(SH);

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all VarScan2 tasks.\n";
}

sub result {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $groups = get_raw_files( $config, $section, "groups" );

  my $result = {};
  for my $groupName ( keys %{$groups} ) {
    my @resultFiles = ();
    my $curDir      = $resultDir . "/$groupName";
    my $snpvcf = "${groupName}.snp.vcf";
    my $annovar   = "${groupName}.somatic.pass.annovar";
    push( @resultFiles, "$curDir/${annovar}.genome_summary.csv" );
    push( @resultFiles, "$curDir/${snpvcf}.Somatic.hc" );
    
    $result->{$groupName} = \@resultFiles;
  }
  return $result;
}

1;
