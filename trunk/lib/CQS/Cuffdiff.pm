#!/usr/bin/perl
package CQS::Cuffdiff;

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
  $self->{_name} = "Cuffdiff";
  bless $self, $class;
  return $self;
}

sub perform {
  my ($self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option ) = get_parameter( $config, $section );

  my $bowtie2_index = $config->{$section}{bowtie2_index} or die "define ${section}::bowtie2_index first";
  my $bowtie2_fasta = get_param_file( $bowtie2_index . ".fa", "bowtie2_fasta", 1 );

  my $transcript_gtf = get_cuffdiff_gtf( $config, $section );

  my $tophat2map = get_tophat2_map( $config, $section );

  my $groups = get_cuffdiff_groups( $config, $section );

  my $pairs = get_cuffdiff_pairs( $config, $section );

  my @labels = ();
  my @files  = ();

  my %tpgroups = ();
  for my $groupName ( sort keys %{$groups} ) {
    my @samples = @{ $groups->{$groupName} };
    my @gfiles  = ();
    foreach my $sampleName (@samples) {
      my $tophat2File = $tophat2map->{$sampleName};
      push( @gfiles, $tophat2File );
    }
    $tpgroups{$groupName} = \@gfiles;

    #print " $groupName => $tpgroups{$groupName} \n";
  }

  my $shfile = $pbsDir . "/${task_name}.submit";
  open( SH, ">$shfile" ) or die "Cannot create $shfile";

  print SH "
if [ ! -s $transcript_gtf ]; then
  echo $transcript_gtf is not exists! all job will be ignored.
  exit 1
fi

type -P qsub &>/dev/null && export MYCMD=\"qsub\" || export MYCMD=\"bash\" 
";

  for my $pairName ( sort keys %{$pairs} ) {
    my @groupNames = @{ $pairs->{$pairName} };

    my $pbsName = "${pairName}_cdiff.pbs";

    my $pbsFile = $pbsDir . "/$pbsName";
    my $log     = $logDir . "/${pairName}_cdiff.log";

    my $curDir = create_directory_or_die( $resultDir . "/$pairName" );

    my $labels = merge_string( ",", @groupNames );

    output_header( $pbsFile, $pbsDesc, $path_file, $log );
    print OUT "cuffdiff $option -o $curDir -L $labels -b $bowtie2_fasta $transcript_gtf ";

    my @conditions = ();
    foreach my $groupName (@groupNames) {
      my @bamfiles = @{ $tpgroups{$groupName} };
      my $bams = merge_string( ",", @bamfiles );
      print OUT "$bams ";

      foreach my $bam (@bamfiles) {
        push( @conditions, "[ -s $bam ]" );
      }
    }
    print OUT "\n";

    output_footer();

    print "$pbsFile created. \n";

    my $condition = merge_string( " && ", @conditions );
    print SH "
if [ -s ${curDir}/gene_exp.diff ];then
  echo job has already been done. if you want to do again, delete ${curDir}/gene_exp.diff and submit job again.
else
  if $condition;then
    \$MYCMD ./$pbsName 
    echo $pbsName was submitted. 
  else
    echo some required file not exists! $pbsName will be ignored.
  fi
fi

";
  }

  print SH "exit 0\n";
  close(SH);

  my $sigfile = $pbsDir . "/${task_name}_sig.pl";
  open( SH, ">$sigfile" ) or die "Cannot create $sigfile";

  print SH "#!/usr/bin/perl
use strict;
use warnings;

use CQS::RNASeq;

my \$config = {
  rename_diff => {
    target_dir => \"${target_dir}/result/comparison\",
    root_dir   => \"${target_dir}/result\",
    gene_only  => 0
  },
};

copy_and_rename_cuffdiff_file(\$config, \"rename_diff\");

1;

";
  close(SH);

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }
  print "!!!shell file $shfile created, you can run this shell file to submit cuffdiff task.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $fasta_format = $config->{$section}{fasta_format};
  if ( !defined $fasta_format ) {
    $fasta_format = 0;
  }

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $result = {};
  for my $sampleName ( keys %rawFiles ) {
    my $curDir = $resultDir . "/$sampleName";

    my @bamFiles  = @{ $rawFiles{$sampleName} };
    my $bamFile   = $bamFiles[0];
    my $fileName  = basename($bamFile);
    my $countFile = $curDir . "/" . $fileName . ".count";

    my @resultFiles = ();
    if ( !defined $pattern || $countFile =~ m/$pattern/ ) {
      push( @resultFiles, $countFile );
    }

    my $unmapped;
    if ($fasta_format) {
      $unmapped = change_extension( $countFile, ".unmapped.fasta" );
    }
    else {
      $unmapped = change_extension( $countFile, ".unmapped.fastq" );
    }
    if ( !defined $pattern || $unmapped =~ m/$pattern/ ) {
      push( @resultFiles, $unmapped );
    }

    $result->{$sampleName} = \@resultFiles;
  }
  return $result;
}

1;
