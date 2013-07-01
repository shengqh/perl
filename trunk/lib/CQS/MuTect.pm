#!/usr/bin/perl
package CQS::MuTect;

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
  $self->{_name} = "MuTect";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $muTect_jar = get_param_file( $config->{$section}{muTect_jar},  "muTect_jar",  1 );
  my $faFile     = get_param_file( $config->{$section}{fasta_file},  "fasta_file",  1 );
  my $cosmicfile = get_param_file( $config->{$section}{cosmic_file}, "cosmic_file", 1 );
  my $dbsnpfile  = get_param_file( $config->{$section}{dbsnp_file},  "dbsnp_file",  1 );

  my $annovarParameter = $config->{$section}{annovar_param} or die "annovar_param is not defined in $section";
  $option = $option . " " . $annovarParameter;

  my $java_option = $config->{$section}{java_option};
  if ( !defined $java_option ) {
    $java_option = "";
  }

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

  for my $groupName ( sort keys %group_sample_map ) {
    my @sampleFiles = @{ $group_sample_map{$groupName} };
    my $sampleCount = scalar(@sampleFiles);
    my $curDir      = create_directory_or_die( $resultDir . "/$groupName" );

    if ( $sampleCount != 2 ) {
      die "SampleFile should be normal,tumor paired.";
    }

    my $normal = $sampleFiles[0];
    my $tumor  = $sampleFiles[1];

    my $out       = "${groupName}.somatic.out";
    my $vcf       = "${groupName}.somatic.vcf";
    my $passvcf   = "${groupName}.somatic.pass.vcf";
    my $passinput = "${groupName}.somatic.pass.avinput";
    my $annovar   = "${groupName}.somatic.pass.annovar";
    my $result    = "${groupName}.somatic.pass.annovar.genome_summary.csv";

    my $pbsName = "${groupName}_mt.pbs";
    my $pbsFile = "${pbsDir}/$pbsName";
    print SH "\$MYCMD ./$pbsName \n";

    my $log = "${logDir}/${groupName}_mt.log";

    open( OUT, ">$pbsFile" ) or die $!;
    print OUT "$pbsDesc
#PBS -o $log
#PBS -j oe

$path_file 

echo muTect=`date` 

cd $curDir

if [ !-s ${normal}.bai ]; then
  samtools index ${normal}
fi

if [ !-s ${tumor}.bai ]; then
  samtools index ${tumor}
fi

if [ ! -s $vcf ]; then
  java $java_option -jar $muTect_jar $option --analysis_type MuTect --reference_sequence $faFile --cosmic $cosmicfile --dbsnp $dbsnpfile --input_file:normal $normal --input_file:tumor $tumor -o $out --coverage_file ${groupName}.coverage.txt --vcf $vcf
fi 

if [[ -s $vcf && ! -s $passvcf ]]; then
  grep -v REJECT $vcf > $passvcf
fi

if [[ -s $passvcf && ! -s $result ]]; then
  convert2annovar.pl -format vcf4 $passvcf -includeinfo > $passinput

  summarize_annovar.pl $annovarParameter --outfile $annovar $passinput $annovarDB
fi

echo finished=`date` \n";
    close OUT;

    print "$pbsFile created \n";
  }

  close(SH);

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all MuTect tasks.\n";
}

sub result {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $groups = get_raw_files( $config, $section, "groups" );

  my $result = {};
  for my $groupName ( keys %{$groups} ) {
    my @resultFiles = ();
    my $curDir      = $resultDir . "/$groupName";
    push( @resultFiles, "$curDir/${groupName}.somatic.pass.annovar.genome_summary.csv" );
    push( @resultFiles, "$curDir/${groupName}.somatic.pass.vcf" );
    $result->{$groupName} = \@resultFiles;
  }
  return $result;
}


1;
