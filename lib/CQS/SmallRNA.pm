#!/usr/bin/perl
package CQS::SmallRNA;

use strict;
use warnings;
use Pipeline::SmallRNA;

require Exporter;
our @ISA = qw(Exporter);

our %EXPORT_TAGS =
  ( 'all' =>
    [qw(performSmallRNA_hg19 performSmallRNATask_hg19 performSmallRNA_hg20 performSmallRNATask_hg20 performSmallRNA_mm10 performSmallRNATask_mm10 performSmallRNA_rn5 performSmallRNATask_rn5 performSmallRNA_cel235 performSmallRNATask_cel235)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

my $hg19_genome = {

  #genome database
  mirbase_count_option  => "-p hsa",
  coordinate            => "/scratch/cqs/shengq1/references/smallrna/hg19_mirBase20_ucsc-tRNA_ensembl75.bed",
  coordinate_fasta      => "/scratch/cqs/shengq1/references/smallrna/hg19_mirBase20_ucsc-tRNA_ensembl75.bed.fa",
  bowtie1_index         => "/data/cqs/guoy1/reference/hg19/bowtie_index_hg19_rCRS_1.0.0/hg19_rCRS",
  bowtie1_miRBase_index => "/data/cqs/shengq1/reference/miRBase20/bowtie_index_1.1.1/mature.dna",
};

my $hg20_genome = {

  #genome database
  mirbase_count_option  => "-p hsa",
  coordinate            => "/scratch/cqs/shengq1/references/smallrna/hg20_mirBase21_ucsc-tRNA_ensembl78.bed",
  coordinate_fasta      => "/scratch/cqs/shengq1/references/smallrna/hg20_mirBase21_ucsc-tRNA_ensembl78.bed.fa",
  bowtie1_index         => "/scratch/cqs/shengq1/references/hg20/bowtie_index_1.1.0/hg20",
  bowtie1_miRBase_index => "/data/cqs/shengq1/reference/miRBase21/bowtie_index_1.1.1/mature.dna",
};

my $mm10_genome = {

  #genome database
  mirbase_count_option  => "-p mmu",
  coordinate            => "/scratch/cqs/shengq1/references/smallrna/mm10_miRBase21_ucsc-tRNA_ensembl78.bed",
  coordinate_fasta      => "/scratch/cqs/shengq1/references/smallrna/mm10_miRBase21_ucsc-tRNA_ensembl78.bed.fa",
  bowtie1_index         => "/data/cqs/shengq1/reference/mm10/bowtie_index/mm10",
  bowtie1_miRBase_index => "/data/cqs/shengq1/reference/miRBase21/bowtie_index_1.1.1/mature.dna",
};

my $rn5_genome = {

  #genome database
  mirbase_count_option  => "-p rno",
  coordinate            => "/scratch/cqs/shengq1/references/smallrna/rn5_miRBase21_ucsc-tRNA_ensembl78.bed",
  coordinate_fasta      => "/scratch/cqs/shengq1/references/smallrna/rn5_miRBase21_ucsc-tRNA_ensembl78.bed.fa",
  bowtie1_index         => "/data/cqs/shengq1/reference/rn5/bowtie_index_1.1.0/rn5",
  bowtie1_miRBase_index => "/data/cqs/shengq1/reference/miRBase21/bowtie_index_1.1.1/mature.dna",
};

my $cel235_genome = {

  #genome database
  mirbase_count_option  => "-p cel",
  coordinate            => "/scratch/cqs/shengq1/references/smallrna/cel235_miRBase21_ensembl78.bed",
  coordinate_fasta      => "/scratch/cqs/shengq1/references/smallrna/cel235_miRBase21_ensembl78.bed.fa",
  bowtie1_index         => "/scratch/cqs/zhangp2/reference/wormbase/bowtie_index_1.1.0/Caenorhabditis_elegans.WBcel235.dna.toplevel",
  bowtie1_miRBase_index => "/data/cqs/shengq1/reference/miRBase21/bowtie_index_1.1.1/mature.dna",
};

sub getDefinition {
  my ( $userdef, $genome ) = @_;

  my $cluster = "slurm";

  if ( defined $userdef->{cluster} ) {
    $cluster = $userdef->{cluster};
  }

  my $min_read_length = 16;
  if ( defined $userdef->{min_read_length} ) {
    $min_read_length = $userdef->{min_read_length};
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

    #Default software parameter (don't change it except you really know it)
    bowtie1_option_1mm         => "-a -m 100 --best --strata -v 1 -p 8",
    bowtie1_option_pm          => "-a -m 100 --best --strata -v 0 -p 8",
    mirnacount_option          => "-s",                                         #ignore score
    smallrnacount_option       => "-s --min_overlap 0.5 --length --sequence",
    mirna_overlap_count_option => "-s --min_overlap 0.5 --gtf_key miRNA",
    min_read_length            => $min_read_length,

    #Software and miRBase database options
    samtools => "/scratch/cqs/shengq1/local/bin/samtools",
    cqstools => "/home/shengq1/cqstools/CQS.Tools.exe",

    #genome database
    mirbase_count_option => $genome->{mirbase_count_option},

    fastq_remove_N => $userdef->{fastq_remove_N},
    adapter        => $userdef->{adapter},

    coordinate            => $genome->{coordinate},
    coordinate_fasta      => $genome->{coordinate_fasta},
    bowtie1_index         => $genome->{bowtie1_index},
    bowtie1_miRBase_index => $genome->{bowtie1_miRBase_index}
  };

  return $def;
}

sub performSmallRNA_hg19 {
  my ($userdef) = @_;
  my $def = getDefinition( $userdef, $hg19_genome );

  performSmallRNA($def);
}

sub performSmallRNATask_hg19 {
  my ( $userdef, $task ) = @_;
  my $def = getDefinition( $userdef, $hg19_genome );

  performSmallRNATask( $def, $task );
}

sub performSmallRNA_hg20 {
  my ($userdef) = @_;
  my $def = getDefinition( $userdef, $hg20_genome );

  performSmallRNA($def);
}

sub performSmallRNATask_hg20 {
  my ( $userdef, $task ) = @_;
  my $def = getDefinition( $userdef, $hg20_genome );

  performSmallRNATask( $def, $task );
}

sub performSmallRNA_mm10 {
  my ($userdef) = @_;
  my $def = getDefinition( $userdef, $mm10_genome );

  performSmallRNA($def);
}

sub performSmallRNATask_mm10 {
  my ( $userdef, $task ) = @_;
  my $def = getDefinition( $userdef, $mm10_genome );

  performSmallRNATask( $def, $task );
}

sub performSmallRNA_rn5 {
  my ($userdef) = @_;
  my $def = getDefinition( $userdef, $rn5_genome );
  performSmallRNA($def);
}

sub performSmallRNATask_rn5 {
  my ( $userdef, $task ) = @_;
  my $def = getDefinition( $userdef, $rn5_genome );

  performSmallRNATask( $def, $task );
}

sub performSmallRNA_cel235 {
  my ($userdef) = @_;
  my $def = getDefinition( $userdef, $cel235_genome );
  performSmallRNA($def);
}

sub performSmallRNATask_cel235 {
  my ( $userdef, $task ) = @_;
  my $def = getDefinition( $userdef, $cel235_genome );

  performSmallRNATask( $def, $task );
}

1;
