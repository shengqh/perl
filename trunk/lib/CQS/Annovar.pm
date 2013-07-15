#!/usr/bin/perl
package CQS::Annovar;

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
  $self->{_name} = "Annovar";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $buildver = $config->{$section}{buildver} or die "buildver is not defined in $section";
  $option = "-buildver $buildver $option";
    
  my $annovarDB = $config->{$section}{annovar_db} or die "annovar_db is not defined in $section";
  my $isvcf = $config->{$section}{isvcf};
  if ( !defined $isvcf ) {
    $isvcf = 0;
  }

  my $rawFiles = get_raw_files( $config, $section );

  my $shfile = $pbsDir . "/${task_name}_ann.submit";
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  if ($sh_direct) {
    print SH "export MYCMD=\"bash\" \n";
  }
  else {
    print SH "type -P qsub &>/dev/null && export MYCMD=\"qsub\" || export MYCMD=\"bash\" \n";
  }

  for my $sampleName ( sort keys %{$rawFiles} ) {
    my @sampleFiles = @{ $rawFiles->{$sampleName} };

    my $pbsName = "${sampleName}_ann.pbs";

    my $pbsFile = $pbsDir . "/$pbsName";
    my $log     = $logDir . "/${sampleName}_ann.log";

    my $curDir = create_directory_or_die( $resultDir . "/$sampleName" );

    open( OUT, ">$pbsFile" ) or die $!;
    print OUT "$pbsDesc
#PBS -o $log
#PBS -j oe

$path_file

cd $curDir
";

    for my $sampleFile (@sampleFiles) {
      my ($filename, $dir) = fileparse($sampleFile);
      
      if($dir eq $curDir){
        $sampleFile = $filename;
      }
      
      my $passinput = change_extension( $filename, ".avinput" );
      my $annovar   = change_extension( $filename, ".annovar" );
      my $result    = "${annovar}.${buildver}_multianno.csv";
      my $vcf       = "";
      if ($isvcf) {
        $vcf = "convert2annovar.pl -format vcf4 $sampleFile -includeinfo > $passinput
        ";
      }
      print OUT "
if [ ! -s $result ]; then
  $vcf
  
  table_annovar.pl $passinput $annovarDB $option --csvout --outfile $annovar 
fi
";
    }
    print OUT "
echo finished=`date`

exit 1
";
    close(OUT);

    print "$pbsFile created. \n";

    print SH "\$MYCMD ./$pbsName \n";
  }
  print SH "exit 1\n";
  close(SH);

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }
  print "!!!shell file $shfile created, you can run this shell file to submit Annovar tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $buildver = $config->{$section}{buildver} or die "buildver is not defined in $section";
  my $rawFiles = get_raw_files( $config, $section );

  my $result = {};
  for my $sampleName ( sort keys %{$rawFiles} ) {
    my @sampleFiles = @{ $rawFiles->{$sampleName} };
    my $curDir      = $resultDir . "/$sampleName";
    my @resultFiles = ();
    for my $sampleFile (@sampleFiles) {
      my $annovar = change_extension( $sampleFile, ".annovar" );
      my $result    = "${annovar}.${buildver}_multianno.csv";
      push( @resultFiles, $curDir . "/$result" );
    }
    $result->{$sampleName} = filter_array( \@resultFiles, $pattern );
  }
  return $result;
}

1;
