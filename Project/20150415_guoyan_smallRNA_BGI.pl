#!/usr/bin/perl
use strict;
use warnings;

use CQS::PerformSmallRNA;

my $def_human = {

  #General options
  task_name       => "BGI_human",
  email           => "quanhu.sheng\@vanderbilt.edu",
  target_dir      => "/scratch/cqs/shengq1/smallRNA/20150415_guoyan_smallRNA_BGI/",
  max_thread      => 8,
  min_read_length => 16,
  cluster         => "slurm",
  fastq_remove_N  => 0,
  run_cutadapt    => 0,

  #Data
  files => {
    "01" => ["/gpfs21/scratch/cqs/guoy1/BGI_miRNA/01clean.fq.gz"],
    "02" => ["/gpfs21/scratch/cqs/guoy1/BGI_miRNA/02clean.fq.gz"],
    "03" => ["/gpfs21/scratch/cqs/guoy1/BGI_miRNA/03clean.fq.gz"],
    "04" => ["/gpfs21/scratch/cqs/guoy1/BGI_miRNA/04clean.fq.gz"],
    "05" => ["/gpfs21/scratch/cqs/guoy1/BGI_miRNA/05clean.fq.gz"],
    "06" => ["/gpfs21/scratch/cqs/guoy1/BGI_miRNA/06clean.fq.gz"],
    "07" => ["/gpfs21/scratch/cqs/guoy1/BGI_miRNA/07clean.fq.gz"],
    "08" => ["/gpfs21/scratch/cqs/guoy1/BGI_miRNA/08clean.fq.gz"],
    "09" => ["/gpfs21/scratch/cqs/guoy1/BGI_miRNA/09clean.fq.gz"],
    "10" => ["/gpfs21/scratch/cqs/guoy1/BGI_miRNA/10clean.fq.gz"],
    "11" => ["/gpfs21/scratch/cqs/guoy1/BGI_miRNA/11clean.fq.gz"],
    "12" => ["/gpfs21/scratch/cqs/guoy1/BGI_miRNA/12clean.fq.gz"],
    "13" => ["/gpfs21/scratch/cqs/guoy1/BGI_miRNA/13clean.fq.gz"],
    "14" => ["/gpfs21/scratch/cqs/guoy1/BGI_miRNA/14clean.fq.gz"],
    "15" => ["/gpfs21/scratch/cqs/guoy1/BGI_miRNA/15clean.fq.gz"],
    "16" => ["/gpfs21/scratch/cqs/guoy1/BGI_miRNA/16clean.fq.gz"],
    "17" => ["/gpfs21/scratch/cqs/guoy1/BGI_miRNA/17clean.fq.gz"],
    "18" => ["/gpfs21/scratch/cqs/guoy1/BGI_miRNA/18clean.fq.gz"],
    "19" => ["/gpfs21/scratch/cqs/guoy1/BGI_miRNA/19clean.fq.gz"],
    "20" => ["/gpfs21/scratch/cqs/guoy1/BGI_miRNA/20clean.fq.gz"],
  }
};

performSmallRNA_hg19($def_human);
1;
