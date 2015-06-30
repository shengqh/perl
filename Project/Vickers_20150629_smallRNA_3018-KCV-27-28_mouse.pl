#!/usr/bin/perl
use strict;
use warnings;

use CQS::PerformSmallRNA;

my $def = {

  #General options
  task_name  => "3018-27-28",
  email      => "quanhu.sheng\@vanderbilt.edu",
  target_dir => "/scratch/cqs/shengq1/vickers/20150629_smallRNA_3018-KCV-27-28_mouse",
  max_thread => 8,

  #Default software parameter (don't change it except you really know it)
  fastq_remove_N => 0,

  #Data
  files => {
    "i01" => ["/gpfs21/scratch/vantage_repo/Vickers/3018/3018-KCV-27/3018-KCV-27-i1_ATCACG_L003_R1_001.fastq.gz"],
    "i02" => ["/gpfs21/scratch/vantage_repo/Vickers/3018/3018-KCV-27/3018-KCV-27-i2_CGATGT_L003_R1_001.fastq.gz"],
    "i03" => ["/gpfs21/scratch/vantage_repo/Vickers/3018/3018-KCV-27/3018-KCV-27-i3_TTAGGC_L003_R1_001.fastq.gz"],
    "i04" => ["/gpfs21/scratch/vantage_repo/Vickers/3018/3018-KCV-27/3018-KCV-27-i4_TGACCA_L003_R1_001.fastq.gz"],
    "i05" => ["/gpfs21/scratch/vantage_repo/Vickers/3018/3018-KCV-28/3018-KCV-28-i5_ACAGTG_L004_R1_001.fastq.gz"],
    "i06" => ["/gpfs21/scratch/vantage_repo/Vickers/3018/3018-KCV-28/3018-KCV-28-i6_GCCAAT_L004_R1_001.fastq.gz"],
    "i07" => ["/gpfs21/scratch/vantage_repo/Vickers/3018/3018-KCV-28/3018-KCV-28-i7_CAGATC_L004_R1_001.fastq.gz"],
    "i08" => ["/gpfs21/scratch/vantage_repo/Vickers/3018/3018-KCV-28/3018-KCV-28-i8_ACTTGA_L004_R1_001.fastq.gz"],
    "i09" => ["/gpfs21/scratch/vantage_repo/Vickers/3018/3018-KCV-27/3018-KCV-27-i9_GATCAG_L003_R1_001.fastq.gz"],
    "i10" => ["/gpfs21/scratch/vantage_repo/Vickers/3018/3018-KCV-27/3018-KCV-27-i10_TAGCTT_L003_R1_001.fastq.gz"],
    "i11" => ["/gpfs21/scratch/vantage_repo/Vickers/3018/3018-KCV-27/3018-KCV-27-i11_GGCTAC_L003_R1_001.fastq.gz"],
    "i12" => ["/gpfs21/scratch/vantage_repo/Vickers/3018/3018-KCV-27/3018-KCV-27-i12_CTTGTA_L003_R1_001.fastq.gz"],
    "i13" => ["/gpfs21/scratch/vantage_repo/Vickers/3018/3018-KCV-28/3018-KCV-28-i13_AGTCAA_L004_R1_001.fastq.gz"],
    "i14" => ["/gpfs21/scratch/vantage_repo/Vickers/3018/3018-KCV-28/3018-KCV-28-i14_AGTTCC_L004_R1_001.fastq.gz"],
    "i15" => ["/gpfs21/scratch/vantage_repo/Vickers/3018/3018-KCV-28/3018-KCV-28-i15_ATGTCA_L004_R1_001.fastq.gz"],
    "i16" => ["/gpfs21/scratch/vantage_repo/Vickers/3018/3018-KCV-28/3018-KCV-28-i16_CCGTCC_L004_R1_001.fastq.gz"],
    "i17" => ["/gpfs21/scratch/vantage_repo/Vickers/3018/3018-KCV-27/3018-KCV-27-i17_GTAGAG_L003_R1_001.fastq.gz"],
    "i18" => ["/gpfs21/scratch/vantage_repo/Vickers/3018/3018-KCV-27/3018-KCV-27-i18_GTCCGC_L003_R1_001.fastq.gz"],
    "i19" => ["/gpfs21/scratch/vantage_repo/Vickers/3018/3018-KCV-27/3018-KCV-27-i19_GTGAAA_L003_R1_001.fastq.gz"],
    "i20" => ["/gpfs21/scratch/vantage_repo/Vickers/3018/3018-KCV-27/3018-KCV-27-i20_GTGGCC_L003_R1_001.fastq.gz"],
    "i21" => ["/gpfs21/scratch/vantage_repo/Vickers/3018/3018-KCV-28/3018-KCV-28-i21_GTTTCG_L004_R1_001.fastq.gz"],
    "i22" => ["/gpfs21/scratch/vantage_repo/Vickers/3018/3018-KCV-28/3018-KCV-28-i22_CGTACG_L004_R1_001.fastq.gz"],
    "i23" => ["/gpfs21/scratch/vantage_repo/Vickers/3018/3018-KCV-28/3018-KCV-28-i23_GAGTGG_L004_R1_001.fastq.gz"],
    "i24" => ["/gpfs21/scratch/vantage_repo/Vickers/3018/3018-KCV-28/3018-KCV-28-i24_GGTAGC_L004_R1_001.fastq.gz"],
  }
};

performSmallRNA_mm10($def);

1;
