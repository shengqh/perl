#!/usr/bin/perl
use strict;
use warnings;

use CQS::PerformSmallRNA2;
use CQS::ClassFactory;

my $def = {

  #General options
  task_name                 => "rat_2570",
  email                     => "quanhu.sheng\@vanderbilt.edu",
  target_dir                => "/scratch/cqs/shengq1/vickers/20161005_smallRNA_2570_rat_v2",
  max_thread                => 8,
  cqstools                  => "/home/shengq1/cqstools/cqstools.exe",
  sequencetask_run_time     => 6,
  table_vis_group_text_size => 12,

  #Default software parameter (don't change it except you really know it)
  fastq_remove_N        => 0,
  remove_sequences      => "'CCACGTTCCCGTGG;ACAGTCCGACGATC'",
  search_unmapped_reads => 1,
  blast_unmapped_reads  => 0,
  top_read_number       => 100,
  blast_top_reads       => 0,
  blast_localdb         => "/scratch/cqs/shengq1/references/blastdb",

  #next flex
  fastq_remove_random => 0,

  #Data
  files => {
    "Lean_3182_1" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/2570_Rat_Liver_diabetes/2570-KCV-1_1_sequence.txt.gz"],
    "Lean_3182_2" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/2570_Rat_Liver_diabetes/2570-KCV-2_1_sequence.txt.gz"],
    "Lean_3183_1" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/2570_Rat_Liver_diabetes/2570-KCV-3_1_sequence.txt.gz"],
    "Lean_3184_1" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/2570_Rat_Liver_diabetes/2570-KCV-4_1_sequence.txt.gz"],
    "Lean_3184_2" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/2570_Rat_Liver_diabetes/2570-KCV-5_1_sequence.txt.gz"],
    "Lean_3185_2" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/2570_Rat_Liver_diabetes/2570-KCV-6_1_sequence.txt.gz"],
    "Veh_3271_1"  => ["/gpfs21/scratch/cqs/shengq1/vickers/data/2570_Rat_Liver_diabetes/2570-KCV-7_1_sequence.txt.gz"],
    "Veh_3259_2"  => ["/gpfs21/scratch/cqs/shengq1/vickers/data/2570_Rat_Liver_diabetes/2570-KCV-8_1_sequence.txt.gz"],
    "Veh_3195_1"  => ["/gpfs21/scratch/cqs/shengq1/vickers/data/2570_Rat_Liver_diabetes/2570-KCV-9_1_sequence.txt.gz"],
    "Veh_3209_1"  => ["/gpfs21/scratch/cqs/shengq1/vickers/data/2570_Rat_Liver_diabetes/2570-KCV-10_1_sequence.txt.gz"],
    "Veh_3259_1"  => ["/gpfs21/scratch/cqs/shengq1/vickers/data/2570_Rat_Liver_diabetes/2570-KCV-11_1_sequence.txt.gz"],
    "Veh_3263_1"  => ["/gpfs21/scratch/cqs/shengq1/vickers/data/2570_Rat_Liver_diabetes/2570-KCV-12_1_sequence.txt.gz"],
    "Col_3214_1"  => ["/gpfs21/scratch/cqs/shengq1/vickers/data/2570_Rat_Liver_diabetes/2570-KCV-13_1_sequence.txt.gz"],
    "Col_3251_1"  => ["/gpfs21/scratch/cqs/shengq1/vickers/data/2570_Rat_Liver_diabetes/2570-KCV-14_1_sequence.txt.gz"],
    "Col_3251_2"  => ["/gpfs21/scratch/cqs/shengq1/vickers/data/2570_Rat_Liver_diabetes/2570-KCV-15_1_sequence.txt.gz"],
    "Col_3258_1"  => ["/gpfs21/scratch/cqs/shengq1/vickers/data/2570_Rat_Liver_diabetes/2570-KCV-16_1_sequence.txt.gz"],
    "Col_3258_2"  => ["/gpfs21/scratch/cqs/shengq1/vickers/data/2570_Rat_Liver_diabetes/2570-KCV-17_1_sequence.txt.gz"],
    "Col_3269_1"  => ["/gpfs21/scratch/cqs/shengq1/vickers/data/2570_Rat_Liver_diabetes/2570-KCV-18_1_sequence.txt.gz"],
  },
  groups => {
    "Col"  => [ "Col_3214_1",  "Col_3251_1",  "Col_3251_2",  "Col_3258_1",  "Col_3258_2",  "Col_3269_1" ],
    "Lean" => [ "Lean_3182_1", "Lean_3182_2", "Lean_3183_1", "Lean_3184_1", "Lean_3184_2", "Lean_3185_2" ],
    "Veh"  => [ "Veh_3271_1",  "Veh_3259_2",  "Veh_3195_1",  "Veh_3209_1",  "Veh_3259_1",  "Veh_3263_1" ],
  },
  pairs => {
    "Veh_vs_Lean" => [ "Lean", "Veh" ],
    "Col_vs_Veh"  => [ "Veh",  "Col" ],
    "Col_vs_Lean" => [ "Lean", "Col" ],
    },
};

performSmallRNA_rn5($def);

1;

