#!/usr/bin/perl
package CQS::PerformSmallRNA;

use strict;
use warnings;
use Pipeline::SmallRNA;

require Exporter;
our @ISA = qw(Exporter);

our %EXPORT_TAGS = (
  'all' => [
    qw(hg19_genome hg20_genome mm10_genome rn5_genome cel235_genome getSmallRNADefinition performSmallRNA_hg19 performSmallRNATask_hg19 performSmallRNA_hg20 performSmallRNATask_hg20 performSmallRNA_mm10 performSmallRNATask_mm10 performSmallRNA_rn5 performSmallRNATask_rn5 performSmallRNA_cel235 performSmallRNATask_cel235)
  ]
);

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

sub hg19_genome {
  return {

    #genome database
    mirbase_count_option  => "-p hsa",
    coordinate            => "/scratch/cqs/shengq1/references/smallrna/hg19_miRBase20_ucsc-tRNA_ensembl75.bed",
    coordinate_fasta      => "/scratch/cqs/shengq1/references/smallrna/hg19_miRBase20_ucsc-tRNA_ensembl75.bed.fa",
    bowtie1_index         => "/scratch/cqs/shengq1/references/hg19_16569_MT/bowtie_index_1.1.2/hg19_16569_MT",
    bowtie1_miRBase_index => "/data/cqs/shengq1/reference/miRBase20/bowtie_index_1.1.1/mature.dna",
    gsnap_index_directory => "/scratch/cqs/shengq1/references/hg19_16569_MT/gsnap_index_k14_2015-06-23/",
    gsnap_index_name      => "hg19_16569_MT",
    star_index_directory => "/scratch/cqs/shengq1/references/hg19_16569_MT/STAR_index_v37.75_2.4.2a_sjdb49"
  };
}

sub hg20_genome {
  return {

    #genome database
    mirbase_count_option  => "-p hsa",
    coordinate            => "/scratch/cqs/shengq1/references/smallrna/hg20_miRBase21_ucsc-tRNA_ensembl78.bed",
    coordinate_fasta      => "/scratch/cqs/shengq1/references/smallrna/hg20_miRBase21_ucsc-tRNA_ensembl78.bed.fa",
    bowtie1_index         => "/scratch/cqs/shengq1/references/hg20/bowtie_index_1.1.0/hg20",
    bowtie1_miRBase_index => "/data/cqs/shengq1/reference/miRBase21/bowtie_index_1.1.1/mature.dna",
    gsnap_index_directory => "/scratch/cqs/shengq1/references/hg20/gsnap_index_k14_2015-06-23/",
    gsnap_index_name      => "hg20",
    star_index_directory => "/scratch/cqs/shengq1/references/hg20/STAR_index_v38.81_2.4.2a_sjdb49"
  };
}

sub mm10_genome {
  return {

    #genome database
    mirbase_count_option  => "-p mmu",
    coordinate            => "/scratch/cqs/shengq1/references/smallrna/mm10_miRBase21_ucsc-tRNA_ensembl78.bed",
    coordinate_fasta      => "/scratch/cqs/shengq1/references/smallrna/mm10_miRBase21_ucsc-tRNA_ensembl78.bed.fa",
    bowtie1_index         => "/scratch/cqs/shengq1/references/mm10/bowtie_index_1.1.1/mm10",
    bowtie1_miRBase_index => "/data/cqs/shengq1/reference/miRBase21/bowtie_index_1.1.1/mature.dna",
    gsnap_index_directory => "/scratch/cqs/shengq1/references/mm10/gsnap_index_k14_2015-06-23/",
    gsnap_index_name      => "mm10",
    star_index_directory => "/scratch/cqs/shengq1/references/mm10/STAR_index_v38.81_2.4.2a_sjdb49"
  };
}

sub rn5_genome {
  return {

    #genome database
    mirbase_count_option  => "-p rno",
    coordinate            => "/scratch/cqs/shengq1/references/smallrna/rn5_miRBase21_ucsc-tRNA_ensembl78.bed",
    coordinate_fasta      => "/scratch/cqs/shengq1/references/smallrna/rn5_miRBase21_ucsc-tRNA_ensembl78.bed.fa",
    bowtie1_index         => "/scratch/cqs/shengq1/references/rn5/bowtie_index_1.1.1/rn5",
    bowtie1_miRBase_index => "/data/cqs/shengq1/reference/miRBase21/bowtie_index_1.1.1/mature.dna",
    gsnap_index_directory => "/scratch/cqs/shengq1/references/rn5/gsnap_index_k14_2015-06-23/",
    gsnap_index_name      => "rn5",
    star_index_directory => "/scratch/cqs/shengq1/references/rn5/STAR_index_v79_2.4.2a_sjdb49"
  };
}

sub cel235_genome {
  return {

    #genome database
    mirbase_count_option  => "-p cel",
    coordinate            => "/scratch/cqs/shengq1/references/smallrna/cel235_miRBase21_ensembl78.bed",
    coordinate_fasta      => "/scratch/cqs/shengq1/references/smallrna/cel235_miRBase21_ensembl78.bed.fa",
    bowtie1_index         => "/scratch/cqs/zhangp2/reference/wormbase/bowtie_index_1.1.0/Caenorhabditis_elegans.WBcel235.dna.toplevel",
    bowtie1_miRBase_index => "/data/cqs/shengq1/reference/miRBase21/bowtie_index_1.1.1/mature.dna",
    gsnap_index_directory => "/scratch/cqs/shengq1/references/cel235/gsnap_index_k14_2015-06-23/",
    gsnap_index_name      => "cel235",
    star_index_directory => "/scratch/cqs/shengq1/references/cel235/STAR_index_v78_2.4.2a_sjdb49"
  };
}

sub getSmallRNADefinition {
  my ( $userdef, $genome ) = @_;

  my $cluster = "slurm";

  if ( defined $userdef->{cluster} ) {
    $cluster = $userdef->{cluster};
  }

  my $min_read_length = 16;
  if ( defined $userdef->{min_read_length} ) {
    $min_read_length = $userdef->{min_read_length};
  }

  my $smallrnacount_option = "-s";
  if ( defined $userdef->{smallrnacount_option} ) {
    $smallrnacount_option = $userdef->{smallrnacount_option};
  }

  my $def = {

    #General options
    task_name  => $userdef->{task_name},
    email      => $userdef->{email},
    target_dir => $userdef->{target_dir},
    max_thread => $userdef->{max_thread},
    cluster    => $cluster,

    #Data
    files => $userdef->{files},

    #Default software parameter (don' t change it except you really know it ) 
    bowtie1_option_1mm => "-a -m 100 --best --strata -v 1 -p 8",
    bowtie1_option_pm    => "-a -m 100 --best --strata -v 0 -p 8",
    smallrnacount_option => $smallrnacount_option,
    min_read_length      => $min_read_length,

    #Software and miRBase database options
    samtools => "/scratch/cqs/shengq1/local/bin/samtools",
    cqstools => "/home/shengq1/cqstools/CQS.Tools.exe",

    #genome database
    mirbase_count_option => $genome->{mirbase_count_option},

    fastq_remove_N => $userdef->{fastq_remove_N},
    adapter        => $userdef->{adapter},
    run_cutadapt => $userdef->{run_cutadapt},

    coordinate            => $genome->{coordinate},
    coordinate_fasta      => $genome->{coordinate_fasta},
    bowtie1_index         => $genome->{bowtie1_index},
    bowtie1_miRBase_index => $genome->{bowtie1_miRBase_index},
    
    gsnap_index_directory => $genome->{gsnap_index_directory},
    gsnap_index_name      => $genome->{gsnap_index_name},
    star_index_directory => $genome->{star_index_directory}
  };

  return $def;
}

sub performSmallRNA_hg19 {
  my ( $userdef, $perform ) = @_;
  my $def = getSmallRNADefinition( $userdef, hg19_genome() );

  my $config = performSmallRNA( $def, $perform );
  return $config;
}

sub performSmallRNATask_hg19 {
  my ( $userdef, $task ) = @_;
  my $def = getSmallRNADefinition( $userdef, hg19_genome() );

  performSmallRNATask( $def, $task );
}

sub performSmallRNA_hg20 {
  my ( $userdef, $perform ) = @_;
  my $def = getSmallRNADefinition( $userdef, hg20_genome() );

  my $config = performSmallRNA( $def, $perform );
  return $config;
}

sub performSmallRNATask_hg20 {
  my ( $userdef, $task ) = @_;
  my $def = getSmallRNADefinition( $userdef, hg20_genome() );

  performSmallRNATask( $def, $task );
}

sub performSmallRNA_mm10 {
  my ( $userdef, $perform ) = @_;
  my $def = getSmallRNADefinition( $userdef, mm10_genome() );

  my $config = performSmallRNA( $def, $perform );
  return $config;
}

sub performSmallRNATask_mm10 {
  my ( $userdef, $task ) = @_;
  my $def = getSmallRNADefinition( $userdef, mm10_genome() );

  performSmallRNATask( $def, $task );
}

sub performSmallRNA_rn5 {
  my ( $userdef, $perform ) = @_;
  my $def = getSmallRNADefinition( $userdef, rn5_genome() );

  my $config = performSmallRNA( $def, $perform );
  return $config;
}

sub performSmallRNATask_rn5 {
  my ( $userdef, $task ) = @_;
  my $def = getSmallRNADefinition( $userdef, rn5_genome() );

  performSmallRNATask( $def, $task );
}

sub performSmallRNA_cel235 {
  my ( $userdef, $perform ) = @_;
  my $def = getSmallRNADefinition( $userdef, cel235_genome() );

  my $config = performSmallRNA( $def, $perform );
  return $config;
}

sub performSmallRNATask_cel235 {
  my ( $userdef, $task ) = @_;
  my $def = getSmallRNADefinition( $userdef, cel235_genome() );

  performSmallRNATask( $def, $task );
}

1;
