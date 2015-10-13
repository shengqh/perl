#!/usr/bin/perl
use strict;
use warnings;

use CQS::PerformSmallRNA;
use CQS::ClassFactory;

my $def = {

  #General options
  task_name  => "3018-KCV-35-38",
  email      => "quanhu.sheng\@vanderbilt.edu",
  target_dir => "/scratch/cqs/shengq1/vickers/20150709_smallRNA_3018-KCV-35-38_human",
  max_thread => 8,

  #Default software parameter (don't change it except you really know it)
  fastq_remove_N => 0,

  #Data
  files => {
    "Lean_381230"  => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-35-i1_ATCACG_L005_R1_001.fastq.gz"],
    "Lean_62790"   => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-35-i2_CGATGT_L005_R1_001.fastq.gz"],
    "Lean_502220"  => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-35-i3_TTAGGC_L005_R1_001.fastq.gz"],
    "Lean_511060"  => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-35-i4_TGACCA_L005_R1_001.fastq.gz"],
    "T2D_97D"      => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-35-i5_ACAGTG_L005_R1_001.fastq.gz"],
    "T2D_74D"      => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-35-i6_GCCAAT_L005_R1_001.fastq.gz"],
    "T2D_119D"     => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-35-i7_CAGATC_L005_R1_001.fastq.gz"],
    "T2D_85D"      => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-35-i8_ACTTGA_L005_R1_001.fastq.gz"],
    "T2D_81D"      => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-35-i9_GATCAG_L005_R1_001.fastq.gz"],
    "Obese_504158" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-35-i21_GTTTCG_L005_R1_001.fastq.gz"],
    "Obese_275760" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-35-i22_CGTACG_L005_R1_001.fastq.gz"],
    "Obese_366342" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-35-i23_GAGTGG_L005_R1_001.fastq.gz"],
    "Obese_206240" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-35-i24_GGTAGC_L005_R1_001.fastq.gz"],
    "Obese_504307" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-35-i25_ACTGAT_L005_R1_001.fastq.gz"],
    "T2D_40D"      => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-36-i1_ATCACG_L006_R1_001.fastq.gz"],
    "T2D_78D"      => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-36-i2_CGATGT_L006_R1_001.fastq.gz"],
    "T2D_91D"      => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-36-i3_TTAGGC_L006_R1_001.fastq.gz"],
    "T2D_99D"      => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-36-i4_TGACCA_L006_R1_001.fastq.gz"],
    "Lean_18100"   => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-36-i5_ACAGTG_L006_R1_001.fastq.gz"],
    "Lean_327163"  => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-36-i6_GCCAAT_L006_R1_001.fastq.gz"],
    "Lean_368895"  => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-36-i7_CAGATC_L006_R1_001.fastq.gz"],
    "Lean_504930"  => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-36-i8_ACTTGA_L006_R1_001.fastq.gz"],
    "Lean_92162"   => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-36-i9_GATCAG_L006_R1_001.fastq.gz"],
    "Obese_77740"  => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-36-i26_ATGAGC_L006_R1_001.fastq.gz"],
    "Obese_510244" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-36-i27_ATTCCT_L006_R1_001.fastq.gz"],
    "Obese_502285" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-36-i28_CAAAAG_L006_R1_001.fastq.gz"],
    "Obese_511462" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-36-i29_CAACTA_L006_R1_001.fastq.gz"],
    "Obese_509217" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-36-i30_CACCGG_L006_R1_001.fastq.gz"],
    "Lean_511502"  => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-37-i10_TAGCTT_L007_R1_001.fastq.gz"],
    "Lean_505327"  => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-37-i11_GGCTAC_L007_R1_001.fastq.gz"],
    "Lean_356965"  => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-37-i12_CTTGTA_L007_R1_001.fastq.gz"],
    "Lean_500453"  => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-37-i13_AGTCAA_L007_R1_001.fastq.gz"],
    "Lean_319490"  => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-37-i14_AGTTCC_L007_R1_001.fastq.gz"],
    "T2D_55D"      => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-37-i15_ATGTCA_L007_R1_001.fastq.gz"],
    "T2D_33D"      => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-37-i16_CCGTCC_L007_R1_001.fastq.gz"],
    "T2D_21D"      => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-37-i18_GTCCGC_L007_R1_001.fastq.gz"],
    "T2D_75D"      => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-37-i19_GTGAAA_L007_R1_001.fastq.gz"],
    "Lean_309920"  => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-37-i20_GTGGCC_L007_R1_001.fastq.gz"],
    "Obese_87795"  => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-37-i36_CCAACA_L007_R1_001.fastq.gz"],
    "Obese_229220" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-37-i37_CGGAAT_L007_R1_001.fastq.gz"],
    "Obese_149430" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-37-i38_CTAGCT_L007_R1_001.fastq.gz"],
    "Obese_52840"  => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-37-i39_CTATAC_L007_R1_001.fastq.gz"],
    "Obese_505399" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-37-i40_CTCAGA_L007_R1_001.fastq.gz"],
    "T2D_88D"      => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-38-i10_TAGCTT_L008_R1_001.fastq.gz"],
    "T2D_80D"      => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-38-i11_GGCTAC_L008_R1_001.fastq.gz"],
    "T2D_120D"     => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-38-i12_CTTGTA_L008_R1_001.fastq.gz"],
    "T2D_100D"     => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-38-i13_AGTCAA_L008_R1_001.fastq.gz"],
    "T2D_61D"      => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-38-i14_AGTTCC_L008_R1_001.fastq.gz"],
    "Lean_147150"  => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-38-i15_ATGTCA_L008_R1_001.fastq.gz"],
    "Lean_510480"  => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-38-i16_CCGTCC_L008_R1_001.fastq.gz"],
    "Lean_262122"  => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-38-i17_GTAGAG_L008_R1_001.fastq.gz"],
    "Lean_511677"  => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-38-i18_GTCCGC_L008_R1_001.fastq.gz"],
    "Lean_133875"  => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-38-i19_GTGAAA_L008_R1_001.fastq.gz"],
    "Obese_500506" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-38-i31_CACGAT_L008_R1_001.fastq.gz"],
    "Obese_511855" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-38-i32_CACTCA_L008_R1_001.fastq.gz"],
    "Obese_508023" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-38-i33_CAGGCG_L008_R1_001.fastq.gz"],
    "Obese_313320" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-38-i34_CATGGC_L008_R1_001.fastq.gz"],
    "Obese_161490" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-38-i35_CATTTT_L008_R1_001.fastq.gz"],
  },
};

my $config = performSmallRNA_hg19($def, 0);
performTask($config, "bowtie1_miRbase_pm_count");

1;

