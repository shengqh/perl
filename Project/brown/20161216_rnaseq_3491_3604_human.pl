#!/usr/bin/perl
use strict;
use warnings;

use Pipeline::RNASeq;
use CQS::ClassFactory;

my $def = {

  #General options
  task_name                 => "human_3491_3604",
  email                     => "quanhu.sheng\@vanderbilt.edu",
  target_dir                => "/scratch/cqs/shengq1/brown/20161216_rnaseq_3491_3604_human",
  max_thread                => 8,
  cqstools                  => "/home/shengq1/cqstools/cqstools.exe",
  sequencetask_run_time     => 6,
  table_vis_group_text_size => 12,

  use_pearson_in_hca => 1,
  use_green_red_color_in_hca => 1,
  top25cv_in_hca => 1,

  transcript_gtf => "/scratch/cqs/shengq1/references/gencode/hg19/gencode.v25lift37.annotation.gtf",
  name_map_file  => "/scratch/cqs/shengq1/references/gencode/hg19/gencode.v25lift37.annotation.map",
  star_index     => "/scratch/cqs/shengq1/references/gencode/hg19/STAR_index_2.5.2b_gencodeV25_sjdb99",

  #Data
  files => {
    "C1_4"      => [ "/gpfs21/scratch/cqs/shengq1/brown/data/3491/3491-DR-1_1_combined_sequence_trim.txt.gz", "/gpfs21/scratch/cqs/shengq1/brown/data/3491/3491-DR-1_2_combined_sequence_trim.txt.gz" ],
    "C2_1"      => [ "/gpfs21/scratch/cqs/shengq1/brown/data/3491/3491-DR-2_1_combined_sequence_trim.txt.gz", "/gpfs21/scratch/cqs/shengq1/brown/data/3491/3491-DR-2_2_combined_sequence_trim.txt.gz" ],
    "C4_3"      => [ "/gpfs21/scratch/cqs/shengq1/brown/data/3491/3491-DR-3_1_combined_sequence_trim.txt.gz", "/gpfs21/scratch/cqs/shengq1/brown/data/3491/3491-DR-3_2_combined_sequence_trim.txt.gz" ],
    "C6_1"      => [ "/gpfs21/scratch/cqs/shengq1/brown/data/3491/3491-DR-4_1_combined_sequence_trim.txt.gz", "/gpfs21/scratch/cqs/shengq1/brown/data/3491/3491-DR-4_2_combined_sequence_trim.txt.gz" ],
    "BrS026_03" => [ "/gpfs21/scratch/cqs/shengq1/brown/data/3491/3491-DR-5_1_combined_sequence_trim.txt.gz", "/gpfs21/scratch/cqs/shengq1/brown/data/3491/3491-DR-5_2_combined_sequence_trim.txt.gz" ],
    "BrS026_05" => [ "/gpfs21/scratch/cqs/shengq1/brown/data/3491/3491-DR-6_1_combined_sequence_trim.txt.gz", "/gpfs21/scratch/cqs/shengq1/brown/data/3491/3491-DR-6_2_combined_sequence_trim.txt.gz" ],
    "BrS026_10" => [ "/gpfs21/scratch/cqs/shengq1/brown/data/3491/3491-DR-7_1_combined_sequence_trim.txt.gz", "/gpfs21/scratch/cqs/shengq1/brown/data/3491/3491-DR-7_2_combined_sequence_trim.txt.gz" ],
    "BrS026_11" => [ "/gpfs21/scratch/cqs/shengq1/brown/data/3491/3491-DR-8_1_combined_sequence_trim.txt.gz", "/gpfs21/scratch/cqs/shengq1/brown/data/3491/3491-DR-8_2_combined_sequence_trim.txt.gz" ],
    "Isogenic_CR_Y" =>
      [ "/gpfs21/scratch/cqs/shengq1/brown/data/3491/3491-DR-10_1_combined_sequence_trim.txt.gz", "/gpfs21/scratch/cqs/shengq1/brown/data/3491/3491-DR-10_2_combined_sequence_trim.txt.gz" ],
    "Isogenic_CR_37" =>
      [ "/gpfs21/scratch/cqs/shengq1/brown/data/3491/3491-DR-11_1_combined_sequence_trim.txt.gz", "/gpfs21/scratch/cqs/shengq1/brown/data/3491/3491-DR-11_2_combined_sequence_trim.txt.gz" ],
    "Isogenic_S16" => [ "/gpfs21/scratch/cqs/shengq1/brown/data/3491/3604-DR-1_1_sequence.txt.gz", "/gpfs21/scratch/cqs/shengq1/brown/data/3491/3604-DR-1_2_sequence.txt.gz" ],
    "Isogenic_S17" => [ "/gpfs21/scratch/cqs/shengq1/brown/data/3491/3604-DR-2_1_sequence.txt.gz", "/gpfs21/scratch/cqs/shengq1/brown/data/3491/3604-DR-2_2_sequence.txt.gz" ],
    "TBX5_Haplo_Z" =>
      [ "/gpfs21/scratch/cqs/shengq1/brown/data/3491/3491-DR-12_1_combined_sequence_trim.txt.gz", "/gpfs21/scratch/cqs/shengq1/brown/data/3491/3491-DR-12_2_combined_sequence_trim.txt.gz" ],
    "TBX5_Haplo_PP" =>
      [ "/gpfs21/scratch/cqs/shengq1/brown/data/3491/3491-DR-13_1_combined_sequence_trim.txt.gz", "/gpfs21/scratch/cqs/shengq1/brown/data/3491/3491-DR-13_2_combined_sequence_trim.txt.gz" ],
    "TBX5_S18" => [ "/gpfs21/scratch/cqs/shengq1/brown/data/3491/3604-DR-3_1_sequence.txt.gz", "/gpfs21/scratch/cqs/shengq1/brown/data/3491/3604-DR-3_2_sequence.txt.gz" ],
  },
  groups => {
    "Control"  => [ "C1_4",          "C2_1",           "C4_3",         "C6_1" ],
    "Brugada"  => [ "BrS026_03",     "BrS026_05",      "BrS026_10",    "BrS026_11" ],
    "Isogenic" => [ "Isogenic_CR_Y", "Isogenic_CR_37", "Isogenic_S16", "Isogenic_S17" ],
    "TBX5"     => [ "TBX5_Haplo_Z",  "TBX5_Haplo_PP",  "TBX5_S18" ]
  },
  pairs => {
    "Brugada_vs_Control"  => [ "Control",  "Brugada" ],
    "Isogenic_vs_Control" => [ "Control",  "Isogenic" ],
    "TBX5_vs_Control"     => [ "Control",  "TBX5" ],
    "Isogenic_vs_Brugada" => [ "Brugada",  "Isogenic" ],
    "TBX5_vs_Brugada"     => [ "Brugada",  "TBX5" ],
    "TBX5_vs_Isogenic"    => [ "Isogenic", "TBX5" ],
  },
};

performRNASeq($def);

1;

