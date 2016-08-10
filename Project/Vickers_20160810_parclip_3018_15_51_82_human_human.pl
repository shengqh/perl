#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;
use Pipeline::SmallRNAUtils;
use Pipeline::ParclipSmallRNA;
use CQS::PerformSmallRNA;
use Data::Dumper;

use Hash::Merge qw( merge );

my $userdef = merge(
  {

    #General options
    task_name  => "3018_15_51_82",
    email      => "quanhu.sheng\@vanderbilt.edu",
    target_dir => "/scratch/cqs/shengq1/vickers/20160810_parclip_3018_15_51_82_human_human/",
    max_thread => 8,
    cluster    => "slurm",

    #Default software parameter (don't change it except you really know it)
    fastq_remove_N      => 0,
    remove_sequences    => "'CCACGTTCCCGTGG;ACAGTCCGACGATC'",
    fastq_remove_random => 0,                                       #for trueseq kit
    adapter             => "TGGAATTCTCGGGTGCCAAGG",
    cqstools            => "/home/shengq1/cqstools/cqstools.exe",

    #Data
    files => {
      "3018-KCV-15-15"  => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-15_parclip/3018-KCV-15-15_ATGTCA_L006_R1_001.fastq.gz"],
      "3018-KCV-15-36"  => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-15_parclip/3018-KCV-15-36_CCAACA_L006_R1_001.fastq.gz"],
      "3018-KCV-15-37"  => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-15_parclip/3018-KCV-15-37_CGGAAT_L006_R1_001.fastq.gz"],
      "3018-KCV-15-46"  => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-15_parclip/3018-KCV-15-46_TCCCGA_L006_R1_001.fastq.gz"],
      "3018-KCV-15-47"  => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-15_parclip/3018-KCV-15-47_TCGAAG_L006_R1_001.fastq.gz"],
      "S6_T"            => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-51/3018/3018-KCV-51-i11_S1_R1_001.fastq.gz"],
      "S6_B"            => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-51/3018/3018-KCV-51-i12_S2_R1_001.fastq.gz"],
      "S6_M"            => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-51/3018/3018-KCV-51-i13_S3_R1_001.fastq.gz"],
      "S7_T"            => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-51/3018/3018-KCV-51-i21_S4_R1_001.fastq.gz"],
      "S7_B"            => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-51/3018/3018-KCV-51-i22_S5_R1_001.fastq.gz"],
      "S7_M"            => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-51/3018/3018-KCV-51-i23_S6_R1_001.fastq.gz"],
      "S7_T_no_4SU"     => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-51/3018/3018-KCV-51-i24_S7_R1_001.fastq.gz"],
      "S7_B_no_4SU"     => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-51/3018/3018-KCV-51-i25_S8_R1_001.fastq.gz"],
      "S7_M_no_4SU"     => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-51/3018/3018-KCV-51-i26_S9_R1_001.fastq.gz"],
      "S8_T"            => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-51/3018/3018-KCV-51-i31_S10_R1_001.fastq.gz"],
      "S8_B"            => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-51/3018/3018-KCV-51-i32_S11_R1_001.fastq.gz"],
      "S8_M"            => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-51/3018/3018-KCV-51-i33_S12_R1_001.fastq.gz"],
      "S8_T_no_4SU"     => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-51/3018/3018-KCV-51-i34_S13_R1_001.fastq.gz"],
      "S8_B_no_4SU"     => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-51/3018/3018-KCV-51-i35_S14_R1_001.fastq.gz"],
      "S8_M_no_4SU"     => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-51/3018/3018-KCV-51-i36_S15_R1_001.fastq.gz"],
      "T1_T_T"          => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-82/3018/3018-KCV-82-i1_S1_R1_001.fastq.gz"],
      "B1_T_B"          => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-82/3018/3018-KCV-82-i2_S2_R1_001.fastq.gz"],
      "M1_T_M"          => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-82/3018/3018-KCV-82-i3_S3_R1_001.fastq.gz"],
      "T2_T_T"          => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-82/3018/3018-KCV-82-i4_S4_R1_001.fastq.gz"],
      "B2_T_B"          => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-82/3018/3018-KCV-82-i5_S5_R1_001.fastq.gz"],
      "M2_T_M"          => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-82/3018/3018-KCV-82-i6_S6_R1_001.fastq.gz"],
      "T3_T_T"          => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-82/3018/3018-KCV-82-i7_S7_R1_001.fastq.gz"],
      "B3_T_B"          => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-82/3018/3018-KCV-82-i8_S8_R1_001.fastq.gz"],
      "M3_T_M"          => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-82/3018/3018-KCV-82-i9_S9_R1_001.fastq.gz"],
      "KW8_EC_Podocyte" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-82/3018/3018-KCV-82-i38_S10_R1_001.fastq.gz"],
    },
  },
  hg19_3utr()
);

my $def = getSmallRNADefinition( $userdef, hg19_genome() );

#print Dumper($def);

performParclipSmallRNA($def);

1;
