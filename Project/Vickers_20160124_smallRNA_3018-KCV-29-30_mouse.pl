#!/usr/bin/perl
use strict;
use warnings;

use CQS::PerformSmallRNA;

my $def = {

  #General options
  task_name  => "KCV-29-30-PaperExample",
  email      => "quanhu.sheng\@vanderbilt.edu",
  target_dir => "/scratch/cqs/shengq1/vickers/20160124_3018-KCV-29-30_PaperExample_mouse",
  max_thread => 8,
  cqstools   => "/home/shengq1/cqstools/CQS.Tools.exe",
    
    #Default software parameter (don't change it except you really know it)
    fastq_remove_N        => 0,
    remove_sequences      => "'CCACGTTCCCGTGG;ACAGTCCGACGATC'",
    search_unmapped_reads => 1,

  #Data
  files => {
    "APOB_WildType1" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-30/3018-KCV-30-i1_ATCACG_L006_R1_001.fastq.gz"],
    "HDL_SRBI_Knockout1" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-30/3018-KCV-30-i10_TAGCTT_L006_R1_001.fastq.gz"],
    "HDL_SRBI_Knockout2" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-30/3018-KCV-30-i11_GGCTAC_L006_R1_001.fastq.gz"],
    "HDL_SRBI_Knockout3" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-30/3018-KCV-30-i12_CTTGTA_L006_R1_001.fastq.gz"],
    "APOB_WildType2" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-30/3018-KCV-30-i2_CGATGT_L006_R1_001.fastq.gz"],
    "APOB_WildType3" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-30/3018-KCV-30-i3_TTAGGC_L006_R1_001.fastq.gz"],
    "APOB_SRBI_Knockout1" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-30/3018-KCV-30-i4_TGACCA_L006_R1_001.fastq.gz"],
    "APOB_SRBI_Knockout2" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-30/3018-KCV-30-i5_ACAGTG_L006_R1_001.fastq.gz"],
    "APOB_SRBI_Knockout3" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-30/3018-KCV-30-i6_GCCAAT_L006_R1_001.fastq.gz"],
    "HDL_WildType1" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-30/3018-KCV-30-i7_CAGATC_L006_R1_001.fastq.gz"],
    "HDL_WildType2" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-30/3018-KCV-30-i8_ACTTGA_L006_R1_001.fastq.gz"],
    "HDL_WildType3" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-30/3018-KCV-30-i9_GATCAG_L006_R1_001.fastq.gz"],
    "LIVER_WildType1" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-29/3018-KCV-29-i13_AGTCAA_L005_R1_001.fastq.gz"],
    "LIVER_WildType2" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-29/3018-KCV-29-i14_AGTTCC_L005_R1_001.fastq.gz"],
    "LIVER_WildType3" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-29/3018-KCV-29-i15_ATGTCA_L005_R1_001.fastq.gz"],
    "LIVER_SRBI_Knockout1" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-29/3018-KCV-29-i16_CCGTCC_L005_R1_001.fastq.gz"],
    "LIVER_SRBI_Knockout2" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-29/3018-KCV-29-i17_GTAGAG_L005_R1_001.fastq.gz"],
    "LIVER_SRBI_Knockout3" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-29/3018-KCV-29-i18_GTCCGC_L005_R1_001.fastq.gz"],
    "BILE_WildType1" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-29/3018-KCV-29-i19_GTGAAA_L005_R1_001.fastq.gz"],
    "BILE_WildType2" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-29/3018-KCV-29-i20_GTGGCC_L005_R1_001.fastq.gz"],
    "BILE_WildType3" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-29/3018-KCV-29-i21_GTTTCG_L005_R1_001.fastq.gz"],
    "BILE_SRBI_Knockout1" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-29/3018-KCV-29-i22_CGTACG_L005_R1_001.fastq.gz"],
    "BILE_SRBI_Knockout2" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-29/3018-KCV-29-i23_GAGTGG_L005_R1_001.fastq.gz"],
    "BILE_SRBI_Knockout3" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-29/3018-KCV-29-i24_GGTAGC_L005_R1_001.fastq.gz"],
    },
    groups => {
        "APOB_WildType"    => [ "APOB_WildType1","APOB_WildType2","APOB_WildType3"],        
        "APOB_SRBI_Knockout"    => [ "APOB_SRBI_Knockout1","APOB_SRBI_Knockout2","APOB_SRBI_Knockout3"],  
        "HDL_WildType"    => [ "HDL_WildType1","HDL_WildType2","HDL_WildType3"],        
        "HDL_SRBI_Knockout"    => [ "HDL_SRBI_Knockout1","HDL_SRBI_Knockout2","HDL_SRBI_Knockout3"],  
        "LIVER_WildType"    => [ "LIVER_WildType1","LIVER_WildType2","LIVER_WildType3"],        
        "LIVER_SRBI_Knockout"    => [ "LIVER_SRBI_Knockout1","LIVER_SRBI_Knockout2","LIVER_SRBI_Knockout3"],  
        "BILE_WildType"    => [ "BILE_WildType1","BILE_WildType2","BILE_WildType3"],        
        "BILE_SRBI_Knockout"    => [ "BILE_SRBI_Knockout1","BILE_SRBI_Knockout2","BILE_SRBI_Knockout3"],        
    },
    pairs => { 
        "APOB_Knockout_VS_WildType" => { groups => [ "APOB_WildType", "APOB_SRBI_Knockout" ], },
        "HDL_Knockout_VS_WildType" => { groups => [ "HDL_WildType", "HDL_SRBI_Knockout" ], },
        "LIVER_Knockout_VS_WildType" => { groups => [ "LIVER_WildType", "LIVER_SRBI_Knockout" ], },
        "BILE_Knockout_VS_WildType" => { groups => [ "BILE_WildType", "BILE_SRBI_Knockout" ], },
         },
};

my $config = performSmallRNA_mm10($def, 0);
performTask($config, "reads_in_tasks");


1;

