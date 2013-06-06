#!/usr/bin/perl
package CQS::RNASeq;

use strict;
use warnings;
use File::Copy;
use File::Basename;
use List::Compare;
use CQS::PBS;
use CQS::FileUtils;
use CQS::StringUtils;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::NGSCommon;

require Exporter;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = (
  'all' => [
    qw(tophat2 call_tophat2 tophat2_by_pbs get_tophat2_result call_RNASeQC cufflinks_by_pbs cuffmerge_by_pbs cuffdiff_by_pbs read_cufflinks_fpkm read_cuffdiff_significant_genes copy_and_rename_cuffdiff_file compare_cuffdiff miso_by_pbs novoalign shrimp2)
  ]
);

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

use Cwd;

sub get_sorted_bam_prefix {
  my ($oldbam) = @_;
  my ( $filename, $dirs, $suffix ) = fileparse( $oldbam, qr/\.[^.]*/ );
  return ( $filename . "_sorted" );
}

sub tophat2_by_pbs {
  my ( $config, $section ) = @_;
  my $obj = instantiate("Tophat2");
  $obj->perform( $config, $section );
}

sub tophat2 {
  tophat2_by_pbs(@_);
}

sub call_tophat2 {
  tophat2_by_pbs(@_);
}

#get expected tophat2 result based on tophat2 definition
sub get_tophat2_map {
  my ( $config, $section ) = @_;

  my ( $result, $issource ) = get_raw_files2( $config, $section );
  if ($issource) {
    retun $result;
  }

  my $tophatsection = $config->{$section}{source_ref};
  my $tophat_dir = $config->{$tophatsection}{target_dir} or die "${tophatsection}::target_dir not defined.";
  my ( $logDir, $pbsDir, $resultDir ) = init_dir( $tophat_dir, 0 );
  my %fqFiles = %{$result};

  my $tpresult = {};
  for my $sampleName ( keys %fqFiles ) {
    my $bam = "${resultDir}/${sampleName}/accepted_hits.bam";
    $tpresult->{$sampleName} = $bam;
  }
  return $tpresult;
}

sub call_RNASeQC {
  my ( $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option ) = get_parameter( $config, $section );

  my $transcript_gtf = get_param_file( $config->{$section}{transcript_gtf}, "transcript_gtf", 1 );
  my $genome_fasta   = get_param_file( $config->{$section}{genome_fasta},   "genome_fasta",   1 );
  my $rnaseqc_jar    = get_param_file( $config->{$section}{rnaseqc_jar},    "rnaseqc_jar",    1 );

  my %tophat2map = %{ get_tophat2_map( $config, $section ) };

  my $pbsFile = "${pbsDir}/${task_name}_RNASeQC.pbs";
  my $log     = $logDir . "/${task_name}_RNASeQC.log";
  output_header( $pbsFile, $pbsDesc, $path_file, $log );

  my $sampleFile = $resultDir . "/sample.list";
  open( SF, ">$sampleFile" ) or die "Cannot create $sampleFile";
  print SF "Sample ID\tBam File\tNotes\n";

  for my $sampleName ( sort keys %tophat2map ) {
    my $tophat2File   = $tophat2map{$sampleName};
    my $sortedBamFile = get_sorted_bam($tophat2File);

    print OUT "if [ ! -e $sortedBamFile ];\n";
    print OUT "then\n";

    my ( $filename, $dirs, $suffix ) = fileparse( $tophat2File, qr/\.[^.]*/ );
    my $sortedBamPrefix = get_sorted_bam_prefix($tophat2File);

    print OUT "  cd $dirs \n";
    print OUT "  samtools sort ${filename}${suffix} $sortedBamPrefix \n";
    print OUT "  samtools index ${sortedBamPrefix}.bam \n";
    print OUT "fi\n\n";

    print SF "${sampleName}\t${sortedBamFile}\t${sampleName}\n";
  }

  close SF;

  print OUT "cd $resultDir \n";
  print OUT "java -jar $rnaseqc_jar -s $sampleFile -t $transcript_gtf -r $genome_fasta -o . \n";
  output_footer();

  print "$pbsFile created\n";
}

sub cufflinks_by_pbs {
  my ( $config, $section ) = @_;
  my $obj = instantiate("Cufflinks");
  $obj->perform( $config, $section );
}

sub get_cufflinks_result_gtf {
  my ( $config, $section ) = @_;

  #get cufflinks root directory
  my $cufflinks_dir = $config->{$section}{target_dir} or die "${section}::target_dir not defined.";

  #get cufflinks result directory
  my ( $logDir, $pbsDir, $resultDir ) = init_dir( $cufflinks_dir, 0 );

  my %tophat2map = %{ get_tophat2_map( $config, $section ) };

  my @result = ();

  #get expected gtf file list
  for my $sampleName ( sort keys %tophat2map ) {
    push( @result, "${resultDir}/${sampleName}/transcripts.gtf" );
  }

  #return expected gtf file list
  return ( \@result );
}

sub get_assemblies_file {
  my ( $config, $section, $target_dir ) = @_;

  my $result = get_param_file( $config->{$section}{source}, "${section}::source", 0 );

  if ( defined $result ) {
    return $result;
  }

  my $cufflinkssection = $config->{$section}{source_ref};
  if ( defined $cufflinkssection ) {
    my $cufflinks_gtf = get_cufflinks_result_gtf( $config, $cufflinkssection );
    $result = $target_dir . "/assemblies.txt";
    open( OUT, ">$result" ) or die $!;
    for my $gtf ( @{$cufflinks_gtf} ) {
      print OUT "${gtf}\n";
    }
    close OUT;

    return $result;
  }

  die "define ${section}::source or ${section}::source_ref first!";
}

sub cuffmerge_by_pbs {
  my ( $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option ) = get_parameter( $config, $section );

  my $transcript_gtf = get_param_file( $config->{general}{transcript_gtf}, "transcript_gtf", 0 );
  my $bowtie2_index = $config->{general}{bowtie2_index} or die "define general::bowtie2_index first";
  my $bowtie2_fasta = get_param_file( $bowtie2_index . ".fa", "bowtie2_fasta", 1 );

  my $assembliesfile = get_assemblies_file( $config, $section, $resultDir );

  my $pbsFile = $pbsDir . "/${task_name}_cmerge.pbs";
  my $log     = $logDir . "/${task_name}_cmerge.log";

  output_header( $pbsFile, $pbsDesc, $path_file, $log );

  print OUT "echo cuffmerge=`date` \n";

  my $gtfparam = "";
  if ($transcript_gtf) {
    $gtfparam = "-g $transcript_gtf";
  }

  print OUT "cuffmerge $option $gtfparam -s $bowtie2_fasta -o $resultDir $assembliesfile \n";

  output_footer();

  print "$pbsFile created\n";

  open( FILE, $assembliesfile ) or die("Unable to open file $assembliesfile");
  my @data = <FILE>;
  close(FILE);

  my $shfile = $pbsDir . "/${task_name}.submit";
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  print SH "type -P qsub &>/dev/null && export MYCMD=\"qsub\" || export MYCMD=\"bash\" \n";
  for my $gtf (@data) {
    chomp($gtf);
    print SH "if [ ! -s $gtf ];\n";
    print SH "then\n";
    print SH "  echo $gtf is not exists, cannot submit the job.\n";
    print SH "  exit;\n";
    print SH "fi;\n\n";
  }
  print SH "\$MYCMD $pbsFile\n";
  print SH "exit 0\n";
  close(SH);

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }
  print "!!!shell file $shfile created, you can run this shell file to submit cuffmerge task.\n";
}

sub get_cuffdiff_gtf {
  my ( $config, $section ) = @_;

  my $transcript_gtf = get_param_file( $config->{$section}{transcript_gtf}, "transcript_gtf", 0 );
  if ( defined $transcript_gtf ) {
    return ($transcript_gtf);
  }

  my $cuffmergesection = $config->{$section}{transcript_gtf_ref};
  if ( defined $cuffmergesection ) {
    my $cuffmerge_target_dir = $config->{$cuffmergesection}{target_dir};
    my ( $logDir, $pbsDir, $resultDir ) = init_dir($cuffmerge_target_dir);
    return ("${resultDir}/merged.gtf");
  }

  $transcript_gtf = get_param_file( $config->{general}{transcript_gtf}, "transcript_gtf", 0 );
  if ( defined $transcript_gtf ) {
    return ($transcript_gtf);
  }

  die "define ${section}::transcript_gtf or ${section}::transcript_gtf_ref or general::transcript_gtf first!";
}

sub cuffdiff_by_pbs {
  my ( $config, $section ) = @_;
  my $obj = instantiate("Cuffdiff");
  $obj->perform( $config, $section );
}

sub output_tophat2_script {
  my ( $genome_db, $gtf_file, $gtfIndex, $tophat2_param, $tophatDir, $sampleName, $index, $isSingle, @sampleFiles ) = @_;

  my $curDir = create_directory_or_die( $tophatDir . "/$sampleName" );

  print OUT "echo tophat2=`date` \n";

  my $hasgtf_file = -e $gtf_file;

  my $hasIndexFile = 0;
  if ( defined $gtfIndex ) {
    if ( -e ( $gtfIndex . ".rev.2.bt2" ) ) {
      $hasIndexFile = 1;
    }
  }

  my $tophat2file = $curDir . "/accepted_hits.bam";

  print OUT "if [ -s $tophat2file ];\n";
  print OUT "then\n";
  print OUT "  echo job has already been done. if you want to do again, delete accepted_hits.bam and submit job again.\n";
  print OUT "else\n";
  if ($hasgtf_file) {
    if ( ( $index == 0 ) && ( !$hasIndexFile ) ) {
      print OUT "  tophat2 $tophat2_param -G $gtf_file --transcriptome-index=$gtfIndex -o $curDir $genome_db ";
    }
    else {
      print OUT "  tophat2 $tophat2_param --transcriptome-index=$gtfIndex -o $curDir $genome_db ";
    }
  }
  elsif ($hasIndexFile) {
    print OUT "  tophat2 $tophat2_param --transcriptome-index=$gtfIndex -o $curDir $genome_db ";
  }
  else {
    print OUT "  tophat2 $tophat2_param -o $curDir $genome_db ";
  }

  if ($isSingle) {
    print OUT "$sampleFiles[$index]\n";
  }
  else {
    print OUT "$sampleFiles[$index*2] $sampleFiles[$index*2+1]\n";
  }
  print OUT "fi\n";
}

sub check_is_single() {
  my ( $sampleNameCount, @sampleFiles ) = @_;
  my $sampleFileCount = scalar(@sampleFiles);
  my $isSingle        = 1;
  if ( $sampleNameCount == $sampleFileCount ) {
    $isSingle = 1;
  }
  elsif ( $sampleNameCount * 2 == $sampleFileCount ) {
    $isSingle = 0;
  }
  else {
    die "Count of SampleName should be equals to count/half count of SampleFiles";
  }

  return ($isSingle);
}

sub read_cufflinks_fpkm {
  my $genefile = shift;
  open GF, "<$genefile" or die "Cannot open file $genefile";
  my $header = <GF>;
  my @headers = split( /\t/, $header );
  my ($geneindex) = grep { $headers[$_] eq "gene_id" } 0 .. $#headers;
  my ($fpkmindex) = grep { $headers[$_] eq "FPKM" } 0 .. $#headers;

  my $result = {};
  while ( my $line = <GF> ) {
    my @values = split( /\t/, $line );
    $result->{ $values[$geneindex] } = $values[$fpkmindex];
  }
  close(GF);

  return $result;
}

sub read_cuffdiff_significant_genes {
  my $file = shift;

  my $result = {};
  open IN, "<$file" or die "Cannot open file $file";
  my $header = <IN>;
  chomp($header);
  while ( my $line = <IN> ) {
    chomp($line);
    my @part = split( /\t/, $line );
    if ( $part[13] eq "yes" ) {
      $result->{ $part[2] } = $line;
    }
  }
  close IN;

  return ( $result, $header );
}

#use CQS::RNASeq;
#my $root = "/home/xxx/cuffdiff/result";
#my $targetdir = create_directory_or_die( $root . "/comparison" );
#my $config    = {
#    "RenameDiff" => {
#        target_dir => $targetdir,
#        root_dir   => $root
#    }
#}
#copy_and_rename_cuffdiff_file($config, "RenameDiff");
sub copy_and_rename_cuffdiff_file {
  my ( $config, $section ) = @_;
  my $dir = $config->{$section}{"root_dir"} or die "define ${section}::root_dir first";
  if ( !-d $dir ) {
    die "directory $dir is not exists";
  }

  my $targetdir = $config->{$section}{"target_dir"} or die "define ${section}::target_dir first";
  if ( !-d $targetdir ) {
    create_directory_or_die($targetdir);
  }

  my $gene_only = $config->{$section}{"gene_only"};
  if ( !defined($gene_only) ) {
    $gene_only = 0;
  }

  my @subdirs = list_directories($dir);
  if ( 0 == scalar(@subdirs) ) {
    die "$dir has no sub CuffDiff directories";
  }

  my @filenames;
  if ($gene_only) {
    @filenames = ("gene_exp.diff");
  }
  else {
    @filenames = ( "gene_exp.diff", "splicing.diff" );
  }

  for my $subdir (@subdirs) {
    foreach my $filename (@filenames) {
      my $file = "${dir}/${subdir}/${filename}";
      if ( -s $file ) {
        print "Dealing " . $file . "\n";
        open IN, "<$file" or die "Cannot open file $file";
        my $line = <IN>;
        $line = <IN>;
        close(IN);

        if ( defined($line) ) {
          my @parts      = split( /\t/, $line );
          my $partcount  = scalar(@parts);
          my $targetname = "${targetdir}/${subdir}.${filename}";

          copy( $file, $targetname ) or die "copy failed : $!";

          my $target_sign_name = $targetname . ".sig";
          my $cmd;
          if ($gene_only) {
            $cmd = "cat $targetname | awk '(\$3 != \"-\") && (\$$partcount==\"yes\" || \$$partcount==\"significant\")' > $target_sign_name";
          }
          else {
            $cmd = "cat $targetname | awk '\$$partcount==\"yes\" || \$$partcount==\"significant\"' > $target_sign_name";
          }

          print $cmd . "\n";
          `$cmd`;
        }
      }
    }
  }
}

sub output_compare_cuffdiff_file {
  my ( $data1, $data2, $fileName, $header, @genes ) = @_;
  open OUT, ">$fileName" or die "Cannot create file $fileName";
  print OUT "$header\n";
  my @sortedgenes = sort @genes;
  for my $gene (@sortedgenes) {
    if ( defined $data1->{$gene} ) {
      print OUT "$data1->{$gene}\n";
    }
    if ( defined $data2->{$gene} ) {
      print OUT "$data2->{$gene}\n";
    }
  }
  close(OUT);
}

sub compare_cuffdiff {
  my ( $config, $section ) = @_;

  my $info = $config->{$section};

  my ( $file1, $file2 ) = @{ $info->{"files"} };

  my ( $data1, $header )  = read_cuffdiff_significant_genes($file1);
  my ( $data2, $header2 ) = read_cuffdiff_significant_genes($file2);

  my @genes1 = keys %{$data1};
  my @genes2 = keys %{$data2};

  my $lc = List::Compare->new( \@genes1, \@genes2 );

  my @resultgenes = ();
  if ( $info->{operation} eq "minus" ) {
    @resultgenes = $lc->get_Lonly();
  }
  elsif ( $info->{operation} eq "intersect" ) {
    @resultgenes = $lc->get_intersection();
  }
  else {
    die "Only minus or intersect is supported.";
  }

  my $resultFileName = $info->{"target_file"};

  output_compare_cuffdiff_file( $data1, $data2, $resultFileName, $header, @resultgenes );
}

sub output_header {
  my ( $pbsFile, $pbsDesc, $path_file, $log ) = @_;
  open( OUT, ">$pbsFile" ) or die $!;
  print OUT "$pbsDesc
PBS -o $log
#PBS -j oe

$path_file
";
}

sub output_footer() {
  print OUT "echo finished=`date`\n";
  close OUT;
}

sub miso_by_pbs {
  my ( $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option ) = get_parameter( $config, $section );

  my $gff3file = get_param_file( $config->{$section}{gff3_file}, "gff3_file", 1 );
  my $gff3index = $gff3file . "indexed";

  my %tophat2map = %{ get_tophat2_map( $config, $section ) };

  my $shfile = $pbsDir . "/${task_name}.submit";
  open( SH, ">$shfile" ) or die "Cannot create $shfile";

  print SH "if [ ! -d $gff3index ];\n";
  print SH "then\n";
  print SH "  index_gff.py --index $gff3file $gff3index \n";
  print SH "fi\n\n";

  for my $sampleName ( sort keys %tophat2map ) {
    my $tophat2File      = $tophat2map{$sampleName};
    my $tophat2indexFile = $tophat2File . ".bai";

    my $pbsName = "${sampleName}_miso.pbs";
    my $pbsFile = $pbsDir . "/$pbsName";

    print SH "if [ ! -s $tophat2File ];\n";
    print SH "then";
    print SH "  echo tophat2 of ${sampleName} has not finished, ignore current job. \n";
    print SH "else\n";
    print SH "  qsub ./$pbsName \n";
    print SH "  echo $pbsName was submitted. \n";
    print SH "fi\n";

    my $log = $logDir . "/${sampleName}_miso.log";

    output_header( $pbsFile, $pbsDesc, $path_file, $log );

    my $curDir = create_directory_or_die( $resultDir . "/$sampleName" );

    print OUT "echo MISO=`date` \n";

    print OUT "if [ ! -e $tophat2indexFile ];\n";
    print OUT "then\n";
    print OUT "  samtools index $tophat2File \n";
    print OUT "fi\n";

    print OUT "run_events_analysis.py --compute-genes-psi $gff3index $tophat2File --output-dir $curDir --read-len 35 \n";

    output_footer();

    print "$pbsFile created\n";
  }
  print SH "\nexit 0\n";
  close(SH);

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }
  print "!!!shell file $shfile created, you can run this shell file to submit all miso tasks.\n";
}

sub novoalign {
  my ( $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option ) = get_parameter( $config, $section );

  my $novoindex = get_param_file( $config->{$section}{novoindex}, "novoindex", 1 );

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $shfile = $pbsDir . "/${task_name}.sh";
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  print SH "type -P qsub &>/dev/null && export MYCMD=\"qsub\" || export MYCMD=\"bash\" \n";

  for my $sampleName ( sort keys %rawFiles ) {
    my @sampleFiles = @{ $rawFiles{$sampleName} };

    my $sampleFile = $sampleFiles[0];

    my $samFile         = $sampleName . ".sam";
    my $bamFile         = $sampleName . ".bam";
    my $sortedBamPrefix = $sampleName . "_sorted";
    my $sortedBamFile   = $sortedBamPrefix . ".bam";

    my $pbsName = "${sampleName}_nalign.pbs";
    my $pbsFile = "${pbsDir}/$pbsName";

    print SH "\$MYCMD ./$pbsName \n";

    my $log    = "${logDir}/${sampleName}_nalign.log";
    my $curDir = create_directory_or_die( $resultDir . "/$sampleName" );

    open( OUT, ">$pbsFile" ) or die $!;
    print OUT "$pbsDesc
#PBS -o $log
#PBS -j oe

$path_file

echo novoalign=`date`

cd $curDir

novoalign -d $novoindex -f $sampleFile $option -o SAM > $samFile

samtools view -b -S $samFile -o $bamFile 

samtools sort $bamFile $sortedBamPrefix 

samtools index $sortedBamFile 

samtools flagstat $sortedBamFile > ${sortedBamFile}.stat 

echo finished=`date` 
";

    close OUT;

    print "$pbsFile created \n";
  }
  close(SH);

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all bwa tasks.\n";

  #`qsub $pbsFile`;
}

sub shrimp2 {
  my ( $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $genome_index = $config->{$section}{genome_index} or die "define ${section}::genome_index first";
  die "genome index ${genome_index}.genome not exist" if ( !-e "${genome_index}.genome" );
  my $is_mirna   = $config->{$section}{is_mirna}   or die "define ${section}::is_mirna first";
  my $output_bam = $config->{$section}{output_bam} or die "define ${section}::output_bam first";
  my $mirna = "-M mirna" if $is_mirna or "";

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
    my $sampleFile = $rawFiles{$sampleName};

    my $shrimpFile      = $sampleName . ".shrimp";
    my $samFile         = $sampleName . ".sam";
    my $bamFile         = $sampleName . ".bam";
    my $sortedBamPrefix = $sampleName . "_sorted";
    my $sortedBamFile   = $sortedBamPrefix . ".bam";

    my $pbsName = "${sampleName}_shrimp2.pbs";
    my $pbsFile = "${pbsDir}/$pbsName";

    print SH "\$MYCMD ./$pbsName \n";

    my $log    = "${logDir}/${sampleName}_shrimp2.log";
    my $curDir = create_directory_or_die( $resultDir . "/$sampleName" );

    open( OUT, ">$pbsFile" ) or die $!;

    print OUT "$pbsDesc
#PBS -o $log
#PBS -j oe

$path_file

echo shrimp2=`date`

cd $curDir
";

    if ($output_bam) {
      print OUT "gmapper -L $genome_index $sampleFile $mirna $option --extra-sam-fields >$samFile

if [ -s $samFile ]; then
  samtools view -b -S $samFile -o $bamFile 
  samtools sort $bamFile $sortedBamPrefix 
  samtools index $sortedBamFile 
  samtools flagstat $sortedBamFile > ${sortedBamFile}.stat 
fi

echo finished=`date` 
";
    }
    else {
      print OUT "gmapper -L $genome_index $sampleFile $mirna $option --pretty >$shrimpFile
      
echo finished=`date` 
";
    }

    close OUT;

    print "$pbsFile created \n";
  }
  close(SH);

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all shrimp2 tasks.\n";

  #`qsub $pbsFile`;
}

1;
