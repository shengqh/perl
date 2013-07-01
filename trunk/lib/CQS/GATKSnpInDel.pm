#!/usr/bin/perl
package CQS::GATKSnpInDel;

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
  $self->{_name} = "GATKSnpInDel";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $faFile   = get_param_file( $config->{$section}{fasta_file}, "fasta_file", 1 );
  my $gatk_jar = get_param_file( $config->{$section}{gatk_jar},   "gatk_jar",   1 );

  my @vcfFiles = @{ $config->{$section}{vcf_files} };
  my $knownvcf = "";
  foreach my $vcf (@vcfFiles) {
    if ( $knownvcf eq "" ) {
      $knownvcf = "-D $vcf";
    }
    else {
      $knownvcf = $knownvcf . " -comp $vcf";
    }
  }

  my $java_option = $config->{$section}{java_option};
  if ( !defined $java_option ) {
    $java_option = "";
  }

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $shfile = $pbsDir . "/${task_name}_snp.sh";
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  print SH "type -P qsub &>/dev/null && export MYCMD=\"qsub\" || export MYCMD=\"bash\" \n";

  for my $sampleName ( sort keys %rawFiles ) {
    my $curDir       = create_directory_or_die( $resultDir . "/$sampleName" );
    my $listfilename = "${sampleName}.list";
    my $listfile     = $curDir . "/$listfilename";
    open( LIST, ">$listfile" ) or die "Cannot create $listfile";
    my @sampleFiles = @{ $rawFiles{$sampleName} };
    foreach my $sampleFile (@sampleFiles) {
      print LIST $sampleFile . "\n";
    }
    close(LIST);

    my $snpOut  = $sampleName . "_snp.vcf";
    my $snpStat = $sampleName . "_snp.stat";

    my $pbsName = "${sampleName}_snp.pbs";

    print SH "\$MYCMD ./$pbsName \n";

    my $log = "${logDir}/${sampleName}_snp.log";

    my $pbsFile = "${pbsDir}/$pbsName";
    open( OUT, ">$pbsFile" ) or die $!;
    print OUT "
$pbsDesc
#PBS -o $log
#PBS -j oe

$path_file

cd $curDir

echo SNP=`date` 
java -jar $java_option $gatk_jar $option -T UnifiedGenotyper -R $faFile -I $listfilename $knownvcf --out $snpOut -metrics $snpStat -glm SNP

echo finished=`date`
";
    close OUT;
    print "$pbsFile created\n";

    my $indelOut  = $sampleName . "_indel.vcf";
    my $indelStat = $sampleName . "_indel.stat";

    $pbsName = "${sampleName}_id.pbs";

    print SH "\$MYCMD ./$pbsName \n";

    $log = "${logDir}/${sampleName}_id.log";

    $pbsFile = "${pbsDir}/$pbsName";
    open( OUT, ">$pbsFile" ) or die $!;
    print OUT "
$pbsDesc
#PBS -o $log
#PBS -j oe

$path_file

cd $curDir

echo InDel=`date` 
java -jar $java_option $gatk_jar -T UnifiedGenotyper -R $faFile -I $listfilename $knownvcf --out $indelOut -metrics $indelStat -glm INDEL $option

echo finished=`date`
";
    close OUT;

    print "$pbsFile created\n";
  }
  close(SH);

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all GATK SnpInDel tasks.\n";
}

sub result {
  my ( $self, $config, $section ) = @_;
  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );
  my $result = {};

  my %rawFiles = %{ get_raw_files( $config, $section ) };
  for my $sampleName ( sort keys %rawFiles ) {
    my $curDir      = $resultDir . "/$sampleName";
    my $snpOut      = $sampleName . "_snp.vcf";
    my $indelOut    = $sampleName . "_indel.vcf";
    my @resultFiles = ();
    push( @resultFiles, "${curDir}/${snpOut}" );
    push( @resultFiles, "${curDir}/${indelOut}" );
    $result->{$sampleName} = \@resultFiles;
  }
  return $result;
}

1;
