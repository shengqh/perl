#!/usr/bin/perl
use strict;
use warnings;

use CQS::PerformSmallRNA;

my $def = {

  #General options
  task_name  => "3018-KCV-10",
  email      => "quanhu.sheng\@vanderbilt.edu",
  target_dir => "/scratch/cqs/shengq1/vickers/20150224_smallRNA_3018-KCV-10_human",
  max_thread => 8,

  #Default software parameter (don't change it except you really know it)
  fastq_remove_N => 1,

  #Data
  files => {
    "3018-KCV-10-37" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-10_Human_Cow/3018-KCV-10-37_CGGAAT_L001_R1_001.fastq.gz"],
    "3018-KCV-10-38" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-10_Human_Cow/3018-KCV-10-38_CTAGCT_L001_R1_001.fastq.gz"],
    "3018-KCV-10-39" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-10_Human_Cow/3018-KCV-10-39_CTATAC_L001_R1_001.fastq.gz"],
    "3018-KCV-10-40" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-10_Human_Cow/3018-KCV-10-40_CTCAGA_L001_R1_001.fastq.gz"],
    "3018-KCV-10-41" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-10_Human_Cow/3018-KCV-10-41_GACGAC_L001_R1_001.fastq.gz"],
    "3018-KCV-10-42" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-10_Human_Cow/3018-KCV-10-42_TAATCG_L001_R1_001.fastq.gz"],
    "3018-KCV-10-43" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-10_Human_Cow/3018-KCV-10-43_TACAGC_L001_R1_001.fastq.gz"],
    "3018-KCV-10-44" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-10_Human_Cow/3018-KCV-10-44_TATAAT_L001_R1_001.fastq.gz"],
    "3018-KCV-10-45" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-10_Human_Cow/3018-KCV-10-45_TCATTC_L001_R1_001.fastq.gz"],
    "3018-KCV-10-46" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-10_Human_Cow/3018-KCV-10-46_TCCCGA_L001_R1_001.fastq.gz"],
    "3018-KCV-10-47" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-10_Human_Cow/3018-KCV-10-47_TCGAAG_L001_R1_001.fastq.gz"],
    "3018-KCV-10-48" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-10_Human_Cow/3018-KCV-10-48_TCGGCA_L001_R1_001.fastq.gz"],
    }
};

performSmallRNA_hg19($def);

1;

