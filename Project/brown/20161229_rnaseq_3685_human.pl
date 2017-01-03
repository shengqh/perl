#!/usr/bin/perl
use strict;
use warnings;

use Pipeline::RNASeq;
use CQS::ClassFactory;

my $def = {

  #General options
  task_name                 => "human_3685",
  email                     => "quanhu.sheng\@vanderbilt.edu",
  target_dir                => "/scratch/cqs/shengq1/brown/20161229_rnaseq_3685_human",
  max_thread                => 8,
  cqstools                  => "/home/shengq1/cqstools/cqstools.exe",
  sequencetask_run_time     => 6,
  table_vis_group_text_size => 12,

  use_pearson_in_hca         => 1,
  use_green_red_color_in_hca => 1,
  top25cv_in_hca             => 0,

  transcript_gtf => "/scratch/cqs/shengq1/references/gencode/hg19/gencode.v19.chr_patch_hapl_scaff.annotation.gtf",
  name_map_file  => "/scratch/cqs/shengq1/references/gencode/hg19/gencode.v19.chr_patch_hapl_scaff.annotation.map",
  star_index     => "/scratch/cqs/shengq1/references/gencode/hg19/STAR_index_2.5.2b_gencodeV19_sjdb99",

  perform_rnaseqc      => 0,
  rnaseqc_init_command => "setpkgs -a java",
  rnaseqc_jar          => "/home/shengq1/local/bin/RNA-SeQC_v1.1.8.jar",
  fasta_file           => "/scratch/cqs/shengq1/references/gencode/hg19/GRCh37.p13.genome.fa",
  rrna_fasta           => "/home/shengq1/local/bin/RNA-SeQC-resources/human_all_rRNA.fasta",

  #qc3
  perform_qc3bam => 1,
  qc3_perl       => "/scratch/cqs/shengq1/local/bin/qc3/qc3.pl",

  #Data
  files => {
    "Endothelial_01" => [ "/gpfs21/scratch/vantage_repo/Brown/3685/3685-JDB-1_S1_R1_001.fastq.gz",   "/gpfs21/scratch/vantage_repo/Brown/3685/3685-JDB-1_S1_R2_001.fastq.gz" ],
    "Leukocyte_02"   => [ "/gpfs21/scratch/vantage_repo/Brown/3685/3685-JDB-2_S2_R1_001.fastq.gz",   "/gpfs21/scratch/vantage_repo/Brown/3685/3685-JDB-2_S2_R2_001.fastq.gz" ],
    "Endothelial_03" => [ "/gpfs21/scratch/vantage_repo/Brown/3685/3685-JDB-3_S3_R1_001.fastq.gz",   "/gpfs21/scratch/vantage_repo/Brown/3685/3685-JDB-3_S3_R2_001.fastq.gz" ],
    "Leukocyte_04"   => [ "/gpfs21/scratch/vantage_repo/Brown/3685/3685-JDB-4_S4_R1_001.fastq.gz",   "/gpfs21/scratch/vantage_repo/Brown/3685/3685-JDB-4_S4_R2_001.fastq.gz" ],
    "Endothelial_05" => [ "/gpfs21/scratch/vantage_repo/Brown/3685/3685-JDB-5_S5_R1_001.fastq.gz",   "/gpfs21/scratch/vantage_repo/Brown/3685/3685-JDB-5_S5_R2_001.fastq.gz" ],
    "Leukocyte_06"   => [ "/gpfs21/scratch/vantage_repo/Brown/3685/3685-JDB-6_S6_R1_001.fastq.gz",   "/gpfs21/scratch/vantage_repo/Brown/3685/3685-JDB-6_S6_R2_001.fastq.gz" ],
    "Endothelial_07" => [ "/gpfs21/scratch/vantage_repo/Brown/3685/3685-JDB-7_S7_R1_001.fastq.gz",   "/gpfs21/scratch/vantage_repo/Brown/3685/3685-JDB-7_S7_R2_001.fastq.gz" ],
    "Leukocyte_08"   => [ "/gpfs21/scratch/vantage_repo/Brown/3685/3685-JDB-8_S8_R1_001.fastq.gz",   "/gpfs21/scratch/vantage_repo/Brown/3685/3685-JDB-8_S8_R2_001.fastq.gz" ],
    "Endothelial_09" => [ "/gpfs21/scratch/vantage_repo/Brown/3685/3685-JDB-9_S9_R1_001.fastq.gz",   "/gpfs21/scratch/vantage_repo/Brown/3685/3685-JDB-9_S9_R2_001.fastq.gz" ],
    "HUVEC_11"       => [ "/gpfs21/scratch/vantage_repo/Brown/3685/3685-JDB-11_S11_R1_001.fastq.gz", "/gpfs21/scratch/vantage_repo/Brown/3685/3685-JDB-11_S11_R2_001.fastq.gz" ],
    "HUVEC_12"       => [ "/gpfs21/scratch/vantage_repo/Brown/3685/3685-JDB-12_S12_R1_001.fastq.gz", "/gpfs21/scratch/vantage_repo/Brown/3685/3685-JDB-12_S12_R2_001.fastq.gz" ],
    "HUVEC_13"       => [ "/gpfs21/scratch/vantage_repo/Brown/3685/3685-JDB-13_S13_R1_001.fastq.gz", "/gpfs21/scratch/vantage_repo/Brown/3685/3685-JDB-13_S13_R2_001.fastq.gz" ],
    "HUVEC_14"       => [ "/gpfs21/scratch/vantage_repo/Brown/3685/3685-JDB-14_S14_R1_001.fastq.gz", "/gpfs21/scratch/vantage_repo/Brown/3685/3685-JDB-14_S14_R2_001.fastq.gz" ],
    "HAEC_15"        => [ "/gpfs21/scratch/vantage_repo/Brown/3685/3685-JDB-15_S15_R1_001.fastq.gz", "/gpfs21/scratch/vantage_repo/Brown/3685/3685-JDB-15_S15_R2_001.fastq.gz" ],
    "HAEC_16"        => [ "/gpfs21/scratch/vantage_repo/Brown/3685/3685-JDB-16_S16_R1_001.fastq.gz", "/gpfs21/scratch/vantage_repo/Brown/3685/3685-JDB-16_S16_R2_001.fastq.gz" ],
    "HAEC_17"        => [ "/gpfs21/scratch/vantage_repo/Brown/3685/3685-JDB-17_S17_R1_001.fastq.gz", "/gpfs21/scratch/vantage_repo/Brown/3685/3685-JDB-17_S17_R2_001.fastq.gz" ],
    "HAEC_18"        => [ "/gpfs21/scratch/vantage_repo/Brown/3685/3685-JDB-18_S18_R1_001.fastq.gz", "/gpfs21/scratch/vantage_repo/Brown/3685/3685-JDB-18_S18_R2_001.fastq.gz" ],
  },
  groups => {
    "Endothelial" => [ "Endothelial_01", "Endothelial_03", "Endothelial_05", "Endothelial_07", "Endothelial_09" ],
    "Leukocyte"   => [ "Leukocyte_02",   "Leukocyte_04",   "Leukocyte_06",   "Leukocyte_08" ],
    "HUVEC"       => [ "HUVEC_11",       "HUVEC_12",       "HUVEC_13",       "HUVEC_14" ],
    "HAEC"        => [ "HAEC_15",        "HAEC_16",        "HAEC_17",        "HAEC_18" ]
  },
  pairs => {
    "Leukocyte_vs_Endothelial" => [ "Endothelial", "Leukocyte" ],
    "HUVEC_vs_Endothelial"     => [ "Endothelial", "HUVEC" ],
    "HAEC_vs_Endothelial"      => [ "Endothelial", "HAEC" ],
    "HUVEC_vs_Leukocyte"       => [ "Leukocyte",   "HUVEC" ],
    "HAEC_vs_Leukocyte"        => [ "Leukocyte",   "HAEC" ],
    "HAEC_vs_HUVEC"            => [ "HUVEC",       "HAEC" ],
  },
};

performRNASeq($def);

1;

