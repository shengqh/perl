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

our %EXPORT_TAGS = ( 'all' => [qw(rsmc muTect varscan2)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

use Cwd;

sub rsmc {
  my ( $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $rsmcfile = get_param_file( $config->{$section}{execute_file}, "execute_file", 1 );
  my $source_type = $config->{$section}{source_type} or die "source_type is not defined in $section";

  my $annovarBuildver = $config->{$section}{annovar_buildver} or die "annovar_buildver is not defined in $section";
  $option = $option . " --annovar --annovar_buildver $annovarBuildver ";

  my $rnaediting_db = $config->{$section}{rnaediting_db};
  if ( defined $rnaediting_db ) {
    $option = $option . " --rnaediting --rnaediting_db $rnaediting_db ";
  }

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $mpileupfile      = "";
  my $fafile           = "";
  my $mpileupParameter = "";
  my $isbam            = lc($source_type) eq "bam";
  if ($isbam) {
    $fafile = get_param_file( $config->{$section}{fasta_file}, "fasta_file (for mpileup)", 1 );
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
    my $sampleCount = scalar(@sampleFiles);
    if ( $sampleCount != 2 ) {
      die "SampleFile should be normal,tumor paired.";
    }

    my $curDir = create_directory_or_die( $resultDir . "/$sampleName" );

    my $normal = $sampleFiles[0];
    my $tumor  = $sampleFiles[1];

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

    if ($isbam) {
      for my $sampleFile (@sampleFiles) {
        my $bamindex = $sampleFile . ".bai";
        print OUT "if [ ! -s $bamindex ]; then
  samtools index $sampleFile 
fi
";
      }

      if ( defined $mpileupParameter ) {
        print OUT "samtools mpileup -f $fafile $mpileupParameter $normal $tumor | mono $rsmcfile all -t console $option";
      }
      else {
        print OUT "mono $rsmcfile all -t bam -f $fafile $option -b $normal,$tumor";
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

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $executefile = get_param_file( $config->{$section}{execute_file}, "execute_file (muTect jar file)", 1 );
  my $faFile      = get_param_file( $config->{$section}{fasta_file},   "fasta_file",                     1 );
  my $cosmicfile  = get_param_file( $config->{$section}{cosmic_file},  "cosmic_file",                    1 );
  my $dbsnpfile   = get_param_file( $config->{$section}{dbsnp_file},   "dbsnp_file",                     1 );

  my $annovarParameter = $config->{$section}{annovar_param} or die "annovar_param is not defined in $section";
  $option = $option . " " . $annovarParameter;

  my $annovarDB = $config->{$section}{annovar_db} or die "annovar_db is not defined in $section";

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $shfile = $pbsDir . "/${task_name}.submit";
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  if ($sh_direct) {
    print SH "export MYCMD=\"bash\" \n";
  }
  else {
    print SH "type -P qsub &>/dev/null && export MYCMD=\"qsub\" || export MYCMD=\"bash\" \n";
  }

  $pbsDesc =~ /\=(\d+)gb/;
  my $gb = $1;

  for my $sampleName ( sort keys %rawFiles ) {
    my @sampleFiles = @{ $rawFiles{$sampleName} };
    my $sampleCount = scalar(@sampleFiles);
    my $curDir      = create_directory_or_die( $resultDir . "/$sampleName" );

    if ( $sampleCount != 2 ) {
      die "SampleFile should be normal,tumor paired.";
    }

    my $normal = $sampleFiles[0];
    my $tumor  = $sampleFiles[1];

    my $checkindex = "";
    for my $sampleFile (@sampleFiles) {
      my $bamindex = $sampleFile . ".bai";

    }

    my $out       = "${sampleName}.somatic.out";
    my $vcf       = "${sampleName}.somatic.vcf";
    my $passvcf   = "${sampleName}.somatic.pass.vcf";
    my $passinput = "${sampleName}.somatic.pass.avinput";
    my $annovar   = "${sampleName}.somatic.pass.annovar";

    my $pbsName = "muTect_${sampleName}.pbs";
    my $pbsFile = "${pbsDir}/$pbsName";

    print SH "\$MYCMD ./$pbsName \n";

    my $log = "${logDir}/muTect_${sampleName}.log";

    open( OUT, ">$pbsFile" ) or die $!;
    print OUT "$pbsDesc
#PBS -o $log
#PBS -j oe

$path_file 

echo muTect=`date` 

if [ ! -s ${normal}.bai ]; then
  samtools index $normal
fi

if [ ! -s ${tumor}.bai ]; then
  samtools index $tumor
fi

cd $curDir
    
java -Xmx${gb}g -jar $executefile --analysis_type MuTect --reference_sequence $faFile --cosmic $cosmicfile --dbsnp $dbsnpfile --input_file:normal $normal --input_file:tumor $tumor -o $out --coverage_file ${sampleName}.coverage.txt --vcf $vcf 

grep -v REJECT $vcf > $passvcf

convert2annovar.pl -format vcf4 $passvcf -includeinfo > $passinput

summarize_annovar.pl $annovarParameter --outfile $annovar $passinput $annovarDB

echo finished=`date` \n";
    close OUT;

    print "$pbsFile created \n";
  }

  close(SH);

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all tasks.\n";
}

sub varscan2 {
  my ( $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $executefile = get_param_file( $config->{$section}{execute_file}, "execute_file (varscan2 jar file)", 1 );
  my $faFile = get_param_file( $config->{$section}{fasta_file}, "fasta_file", 1 );

  my $min_coverage    = $config->{$section}{min_coverage}    or die "min_coverage is not defined in $section!";
  my $somatic_p_value = $config->{$section}{somatic_p_value} or die "somatic_p_value is not defined in $section!";

  my $annovarParameter = $config->{$section}{annovar_param} or die "annovar_param is not defined in $section";

  my $annovarDB = $config->{$section}{annovar_db} or die "annovar_db is not defined in $section";

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $shfile = $pbsDir . "/${task_name}.submit";
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  if ($sh_direct) {
    print SH "export MYCMD=\"bash\" \n";
  }
  else {
    print SH "type -P qsub &>/dev/null && export MYCMD=\"qsub\" || export MYCMD=\"bash\" \n";
  }

  $pbsDesc =~ /\=(\d+)gb/;
  my $gb = $1;

  for my $sampleName ( sort keys %rawFiles ) {
    my @sampleFiles = @{ $rawFiles{$sampleName} };
    my $sampleCount = scalar(@sampleFiles);
    my $curDir      = create_directory_or_die( $resultDir . "/$sampleName" );

    if ( $sampleCount != 2 ) {
      die "SampleFile should be normal,tumor paired.";
    }

    my $normal = $sampleFiles[0];
    my $tumor  = $sampleFiles[1];

    my $normalfile = basename($normal);
    my $tumorfile  = basename($tumor);

    my $normal_mpileup = "${normalfile}.mpileup";
    my $tumor_mpileup  = "${tumorfile}.mpileup";

    my $out       = "${sampleName}.somatic.out";
    my $vcf       = "${sampleName}.somatic.vcf";
    my $passvcf   = "${sampleName}.somatic.pass.vcf";
    my $passinput = "${sampleName}.somatic.pass.avinput";
    my $annovar   = "${sampleName}.somatic.pass.annovar";

    my $pbsName = "varscan2_${sampleName}.pbs";
    my $pbsFile = "${pbsDir}/$pbsName";

    print SH "\$MYCMD ./$pbsName \n";

    my $log = "${logDir}/varscan2_${sampleName}.log";

    open( OUT, ">$pbsFile" ) or die $!;
    print OUT "$pbsDesc
#PBS -o $log
#PBS -j oe

$path_file 

echo varscan2=`date` 

if [ ! -s ${normal}.bai ]; then
  samtools index $normal
fi

if [ ! -s ${tumor}.bai ]; then
  samtools index $tumor
fi

cd $curDir

if [ ! -s $normal_mpileup ]; then
  echo NORMAL_MPILEUP=`date`
  samtools mpileup -q 20 -f $faFile $normal > $normal_mpileup
fi

if [ ! -s $tumor_mpileup ]; then
  echo TUMOR_MPILEUP=`date`
  samtools mpileup -q 20 -f $faFile $tumor > $tumor_mpileup
fi

java -Xmx${gb}g -jar $executefile somatic $option $normal_mpileup $tumor_mpileup $sampleName --output-vcf --somatic-p-value $somatic_p_value --min-coverage $min_coverage --strand-filter

java -Xmx${gb}g -jar $executefile processSomatic ${sampleName}.snp.vcf --p-value $somatic_p_value

convert2annovar.pl -format vcf4 ${sampleName}.snp.vcf.Somatic.hc -includeinfo > ${sampleName}.somatic.avinput

summarize_annovar.pl $annovarParameter --outfile $annovar ${sampleName}.somatic.avinput $annovarDB

echo finished=`date` \n";
    close OUT;

    print "$pbsFile created \n";
  }

  close(SH);

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all tasks.\n";
}

1;
