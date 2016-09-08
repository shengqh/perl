#!/usr/bin/perl
use strict;
use warnings;

use CQS::PerformSmallRNA;
use CQS::FileUtils;
use Pipeline::SmallRNAUtils;
use Pipeline::SmallRNA;

my $genome = mm10_genome();

my $root = create_directory_or_die("/scratch/cqs/shengq1/vickers/20160708_smallRNA_3018-KCV-31-32_mouse_test_contig/");

my $userdef = {
  email                     => "quanhu.sheng\@vanderbilt.edu",
  max_thread                => 8,
  cqstools                  => "/home/shengq1/cqstools/CQS.Tools.exe",
  sequencetask_run_time     => 12,
  table_vis_group_text_size => 14,

  #Default software parameter (don't change it except you really know it)
  fastq_remove_N        => 0,
  remove_sequences      => "'CCACGTTCCCGTGG;ACAGTCCGACGATC'",
  search_unmapped_reads => 1,
  blast_unmapped_reads  => 1,
  search_not_identical  => 0,

  #General options
  task_name  => "withcontig",
  target_dir => $root . "with_contig",

  #Data
  files => {
    "HDL_09_WT" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-32-i21_GTTTCG_L002_R1_001.fastq.gz"],
    "HDL_10_WT" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-32-i22_CGTACG_L002_R1_001.fastq.gz"],
    "HDL_11_WT" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-32-i23_GAGTGG_L002_R1_001.fastq.gz"],
    "HDL_12_WT" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-32-i24_GGTAGC_L002_R1_001.fastq.gz"],
  }
};

my $defWithContig = getSmallRNADefinition( $userdef, $genome );
performSmallRNA($defWithContig);

$userdef->{target_dir}   = $root . "without_contig";
$userdef->{task_name}    = "without_contig";
$genome->{bowtie1_index} = "/scratch/cqs/shengq1/references/mm10_sorted_M_slim/bowtie_index_1.1.2/mm10";

my $defWithoutContig = getSmallRNADefinition( $userdef, $genome );
performSmallRNA($defWithoutContig);

1;
