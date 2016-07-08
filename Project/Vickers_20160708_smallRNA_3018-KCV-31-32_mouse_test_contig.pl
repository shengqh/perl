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

  #General options
  task_name  => "3018-31-32-withcontig",
  email      => "quanhu.sheng\@vanderbilt.edu",
  target_dir => $root . "with_contig",
  max_thread => 8,

  #Default software parameter (don't change it except you really know it)
  fastq_remove_N => 0,

  #Data
  files => {
    "HDL_09_WT"       => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-32-i21_GTTTCG_L002_R1_001.fastq.gz"],
    "HDL_10_WT"       => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-32-i22_CGTACG_L002_R1_001.fastq.gz"],
    "HDL_11_WT"       => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-32-i23_GAGTGG_L002_R1_001.fastq.gz"],
    "HDL_12_WT"       => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-32-i24_GGTAGC_L002_R1_001.fastq.gz"],
  }
};

my $defWithContig = getSmallRNADefinition( $userdef, $genome );
performSmallRNA( $defWithContig);

$userdef->{target_dir} =  $root . "without_contig";
$genome->{bowtie1_index} = "/scratch/cqs/shengq1/references/mm10_sorted_M_slim/bowtie_index_1.1.2/mm10";

my $defWithoutContig = getSmallRNADefinition( $userdef, $genome );
performSmallRNA( $defWithoutContig);

1;
