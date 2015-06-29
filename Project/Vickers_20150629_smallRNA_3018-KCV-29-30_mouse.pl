#!/usr/bin/perl
use strict;
use warnings;

use CQS::PerformSmallRNA;

my $def = {

  #General options
  task_name  => "3018-29-30",
  email      => "quanhu.sheng\@vanderbilt.edu",
  target_dir => "/scratch/cqs/shengq1/vickers/20150629_smallRNA_3018-KCV-29-30_mouse",
  max_thread => 8,

  #Default software parameter (don't change it except you really know it)
  fastq_remove_N => 0,

  #Data
  files => {
    "i01" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-30/3018-KCV-30-i1_ATCACG_L006_R1_001.fastq.gz"],
    "i10" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-30/3018-KCV-30-i10_TAGCTT_L006_R1_001.fastq.gz"],
    "i11" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-30/3018-KCV-30-i11_GGCTAC_L006_R1_001.fastq.gz"],
    "i12" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-30/3018-KCV-30-i12_CTTGTA_L006_R1_001.fastq.gz"],
    "i02" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-30/3018-KCV-30-i2_CGATGT_L006_R1_001.fastq.gz"],
    "i03" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-30/3018-KCV-30-i3_TTAGGC_L006_R1_001.fastq.gz"],
    "i04" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-30/3018-KCV-30-i4_TGACCA_L006_R1_001.fastq.gz"],
    "i05" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-30/3018-KCV-30-i5_ACAGTG_L006_R1_001.fastq.gz"],
    "i06" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-30/3018-KCV-30-i6_GCCAAT_L006_R1_001.fastq.gz"],
    "i07" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-30/3018-KCV-30-i7_CAGATC_L006_R1_001.fastq.gz"],
    "i08" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-30/3018-KCV-30-i8_ACTTGA_L006_R1_001.fastq.gz"],
    "i09" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-30/3018-KCV-30-i9_GATCAG_L006_R1_001.fastq.gz"],
    "i13" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-29/3018-KCV-29-i13_AGTCAA_L005_R1_001.fastq.gz"],
    "i14" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-29/3018-KCV-29-i14_AGTTCC_L005_R1_001.fastq.gz"],
    "i15" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-29/3018-KCV-29-i15_ATGTCA_L005_R1_001.fastq.gz"],
    "i16" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-29/3018-KCV-29-i16_CCGTCC_L005_R1_001.fastq.gz"],
    "i17" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-29/3018-KCV-29-i17_GTAGAG_L005_R1_001.fastq.gz"],
    "i18" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-29/3018-KCV-29-i18_GTCCGC_L005_R1_001.fastq.gz"],
    "i19" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-29/3018-KCV-29-i19_GTGAAA_L005_R1_001.fastq.gz"],
    "i20" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-29/3018-KCV-29-i20_GTGGCC_L005_R1_001.fastq.gz"],
    "i21" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-29/3018-KCV-29-i21_GTTTCG_L005_R1_001.fastq.gz"],
    "i22" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-29/3018-KCV-29-i22_CGTACG_L005_R1_001.fastq.gz"],
    "i23" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-29/3018-KCV-29-i23_GAGTGG_L005_R1_001.fastq.gz"],
    "i24" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-29/3018-KCV-29-i24_GGTAGC_L005_R1_001.fastq.gz"],
    "i45" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-29/3018-KCV-29-i45_TCATTC_L005_R1_001.fastq.gz"],
    "i46" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-29/3018-KCV-29-i46_TCCCGA_L005_R1_001.fastq.gz"],
    "i47" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-30/3018-KCV-30-i47_TCGAAG_L006_R1_001.fastq.gz"],
    "i48" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-30/3018-KCV-30-i48_TCGGCA_L006_R1_001.fastq.gz"],
  }
};

performSmallRNA_mm10($def);

1;
