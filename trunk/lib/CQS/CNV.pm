#!/usr/bin/perl
package CQS::CNV;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::DNASeq;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;

require Exporter;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(cnvnator conifer cnmops)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

use Cwd;

sub cnvnator {
  my ( $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option ) = get_parameter( $config, $section );
  my $binsize        = $config->{$section}{binsize}        or die "define ${section}::binsize first";
  my $chromosome_dir = $config->{$section}{chromosome_dir} or die "define ${section}::chromosome_dir first";

  my $genome    = $config->{$section}{genome};
  my $genomestr = "";

  if ( defined $genome ) {
    $genomestr = "-genome " . $genome;
  }

  my $isbamsorted = $config->{$section}{isbamsorted};
  if ( !defined($isbamsorted) ) {
    $isbamsorted = 0;
  }

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $shfile = $pbsDir . "/${task_name}.sh";
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  print SH "type -P qsub &>/dev/null && export MYCMD=\"qsub\" || export MYCMD=\"bash\" \n";

  for my $sampleName ( sort keys %rawFiles ) {
    my @sampleFiles = @{ $rawFiles{$sampleName} };

    my $bamFile = $sampleFiles[0];

    if ( !$isbamsorted ) {
      ( $bamFile, my $bamSorted ) = get_sorted_bam($bamFile);
    }
    my $pbsName = "cnvnator_${sampleName}.pbs";
    my $pbsFile = "${pbsDir}/$pbsName";

    print SH "\$MYCMD ./$pbsName \n";

    my $log = "${logDir}/cnvnator_${sampleName}.log";

    open( OUT, ">$pbsFile" ) or die $!;
    print OUT $pbsDesc;
    print OUT "#PBS -o $log\n";
    print OUT "#PBS -j oe\n\n";

    if ( -e $path_file ) {
      print OUT "source $path_file\n";
    }

    my $curDir   = create_directory_or_die( $resultDir . "/$sampleName" );
    my $rootFile = $sampleName . ".root";
    my $callFile = $sampleName . ".call";

    print OUT "cd $curDir\n\n";

    print OUT "if [ -s $callFile ]; then\n";
    print OUT "  echo job has already been done. if you want to do again, delete $callFile and submit job again.\n";
    print OUT "else\n";
    print OUT "  if [ ! -s $rootFile ]; then\n";
    print OUT "    echo \"EXTRACTING READ MAPPING FROM BAM/SAM FILES =\" `date`\n";
    print OUT "    cnvnator -root $rootFile $genomestr -unique -tree $bamFile \n";
    print OUT "  fi\n";
    print OUT "  echo \"GENERATING HISTOGRAM =\" `date`\n";
    print OUT "  cnvnator -root $rootFile -d $chromosome_dir -his $binsize \n";
    print OUT "  echo \"CALCULATING STATISTICS =\" `date`\n";
    print OUT "  cnvnator -root $rootFile -stat $binsize \n";
    print OUT "  echo \"RD SIGNAL PARTITIONING =\" `date`\n";
    print OUT "  cnvnator -root $rootFile -partition $binsize \n";
    print OUT "  echo \"CNV CALLING =\" `date`\n";
    print OUT "  cnvnator -root $rootFile -call $binsize > $callFile \n";
    print OUT "fi\n\n";

    print OUT "echo finished=`date`\n";
    close OUT;

    print "$pbsFile created\n";
  }
  close(SH);

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all cnvnator tasks.\n";

  #`qsub $pbsFile`;
}

sub conifer {
  my ( $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option ) = get_parameter( $config, $section );
  my $conifer   = $config->{$section}{conifer} or die "define conifer program location first.\nconifer => \"location\"";
  my $probefile = $config->{$section}{probefile};
  my $probedef  = "";
  if ( defined $probefile ) {
    $probedef = "--probes $probefile";
  }

  my $isbamsorted = $config->{$section}{isbamsorted};
  if ( !defined($isbamsorted) ) {
    $isbamsorted = 0;
  }

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $shfile = $pbsDir . "/${task_name}_rpkm.sh";
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  print SH "type -P qsub &>/dev/null && export MYCMD=\"qsub\" || export MYCMD=\"bash\" \n";

  create_directory_or_die( $resultDir . "/rpkm" );
  for my $sampleName ( sort keys %rawFiles ) {
    my @sampleFiles = @{ $rawFiles{$sampleName} };

    my $pbsName = "${sampleName}_rpkm.pbs";
    my $pbsFile = "${pbsDir}/$pbsName";

    print SH "\$MYCMD ./$pbsName \n";

    my $log = "${logDir}/${sampleName}_rpkm.log";

    open( OUT, ">$pbsFile" ) or die $!;
    print OUT $pbsDesc;
    print OUT "#PBS -o $log\n";
    print OUT "#PBS -j oe\n\n";

    if ( -e $path_file ) {
      print OUT "source $path_file\n";
    }
    print OUT "cd $resultDir\n\n";
    print OUT "echo rpkm=`date`\n";

    my $bamFile = $sampleFiles[0];

    if ( !$isbamsorted ) {
      ( $bamFile, my $bamSorted ) = get_sorted_bam($bamFile);

      #print $bamFile . "\n";
    }
    my $rpkm = "rpkm/" . $sampleName . ".rpkm";

    print OUT "if [ ! -s $rpkm ]; then\n";
    print OUT "  echo conifer=`date`\n";
    print OUT "  python $conifer rpkm $probedef --input $bamFile --output $rpkm \n";
    print OUT "fi\n";

    print OUT "echo finished=`date`\n";
    close OUT;

    print "$pbsFile created\n";
  }
  close(SH);

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  my $pbsName   = "${task_name}_after_rpkm.pbs";
  my $pbsFile   = "${pbsDir}/$pbsName";
  my $hdf5File  = "${task_name}.hdf5";
  my $svalsFile = "${task_name}.svals";
  my $callFile  = "${task_name}.call";

  open( OUT, ">$pbsFile" ) or die $!;
  print OUT $pbsDesc;
  my $log = "${logDir}/${task_name}_after_rpkm.log";
  print OUT "#PBS -o $log\n";
  print OUT "#PBS -j oe\n\n";
  if ( -e $path_file ) {
    print OUT "source $path_file\n";
  }

  create_directory_or_die( $resultDir . "/call_images" );

  print OUT "cd $resultDir\n\n";

  print OUT "\n";
  print OUT "#2 analysis\n";
  print OUT "echo analyze=`date`\n";
  print OUT "python $conifer analyze $probedef --rpkm_dir rpkm/ --output $hdf5File --svd 6 --write_svals $svalsFile \n";
  print OUT "\n";
  print OUT "#3 call\n";
  print OUT "echo call=`date`\n";
  print OUT "python $conifer call --input $hdf5File --output $callFile \n";
  print OUT "\n";
  print OUT "#4 plot\n";
  print OUT "python $conifer plotcalls --input $hdf5File --calls $callFile --output call_images \n";
  close OUT;

  print "$pbsFile created\n";

  print "!!!shell file $shfile created, you can run this shell file to submit all conifer rpkm tasks.\n";

  #`qsub $pbsFile`;
}

sub cnmops {
  my ( $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option ) = get_parameter( $config, $section );
  my $probefile = $config->{$section}{probefile};

  my $isbamsorted = $config->{$section}{isbamsorted};
  if ( !defined($isbamsorted) ) {
    $isbamsorted = 0;
  }

  my $pairmode = $config->{$section}{pairmode};
  if ( !defined($pairmode) ) {
    $pairmode = "unpaired";
  }

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $rfile = $pbsDir . "/cnmops_${task_name}.r";
  open( R, ">$rfile" ) or die "Cannot create $rfile";
  print R "library(cn.mops) \n";
  print R "setwd(\"$resultDir\") \n";
  print R "SampleNames <- c( \n";
  my $isfirst = 1;
  for my $sampleName ( sort keys %rawFiles ) {
    if ($isfirst) {
      print R "\"$sampleName\"\n";
      $isfirst = 0;
    }
    else {
      print R ",\"$sampleName\"\n";
    }
  }
  print R ") \n";
  print R "BAMFiles <- c( \n";

  $isfirst = 1;
  for my $sampleName ( sort keys %rawFiles ) {
    my @sampleFiles = @{ $rawFiles{$sampleName} };
    my $bamFile     = $sampleFiles[0];

    if ( !$isbamsorted ) {
      ( $bamFile, my $bamSorted ) = get_sorted_bam($bamFile);
    }

    if ($isfirst) {
      print R "\"$bamFile\"\n";
      $isfirst = 0;
    }
    else {
      print R ",\"$bamFile\"\n";
    }
  }
  print R ") \n";

  if ( defined $probefile ) {
    print R "segments <- read.table(\"$probefile\", sep=\"\\t\", as.is=TRUE, header=T) \n";
    print R "gr <- GRanges(segments[,1], IRanges(segments[,2],segments[,3]), gene=segments[,4]) \n";
    print R "x <- getSegmentReadCountsFromBAM(BAMFiles, GR=gr, sampleNames=SampleNames, mode=\"$pairmode\") \n";
    print R "save(x, file=\"${task_name}_x_getSegmentReadCountsFromBAM.Rdata\") \n";
    print R "resCNMOPS <- exomecn.mops(x) \n";
    print R "save(resCNMOPS, file=\"${task_name}_resCNMOPS_exomecn.mops.Rdata\") \n";
  }
  else {
    print R "x <- getReadCountsFromBAM(BAMFiles, sampleNames=SampleNames, mode=\"$pairmode\") \n";
    print R "save(x, file=\"${task_name}_x_getReadCountsFromBAM.Rdata\") \n";
    print R "resCNMOPS <- cn.mops(x) \n";
    print R "save(resCNMOPS, file=\"${task_name}_resCNMOPS_cn.mops.Rdata\") \n";
  }

  close R;

  my $pbsFile = "${pbsDir}/cnmops_${task_name}.pbs";
  my $log     = "${logDir}/cnmops_${task_name}.log";

  open( OUT, ">$pbsFile" ) or die $!;
  print OUT $pbsDesc;
  print OUT "#PBS -o $log\n";
  print OUT "#PBS -j oe\n\n";

  if ( -e $path_file ) {
    print OUT "source $path_file\n";
  }

  print OUT "cd $pbsDir\n\n";
  print OUT "echo cnmops=`date`\n";
  print OUT "R --vanilla < $rfile \n";
  print OUT "echo finished=`date`\n";
  close OUT;

  print "$pbsFile created\n";
}

1;
