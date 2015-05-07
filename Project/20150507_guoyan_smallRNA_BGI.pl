#!/usr/bin/perl
use strict;
use warnings;

use CQS::PerformSmallRNA;

my $def_human = {

  #General options
  task_name       => "BGI_human",
  email           => "quanhu.sheng\@vanderbilt.edu",
  target_dir      => "/scratch/cqs/shengq1/smallRNA/20150507_guoyan_smallRNA_BGI/",
  max_thread      => 8,
  min_read_length => 16,
  cluster         => "slurm",
  fastq_remove_N  => 0,
  run_cutadapt    => 1,

  #Data
  files => {
    "01" => ["/gpfs21/scratch/cqs/guoy1/BGI_miRNA/raw/01.gz"],
    "02" => ["/gpfs21/scratch/cqs/guoy1/BGI_miRNA/raw/02.gz"],
    "03" => ["/gpfs21/scratch/cqs/guoy1/BGI_miRNA/raw/03.gz"],
    "04" => ["/gpfs21/scratch/cqs/guoy1/BGI_miRNA/raw/04.gz"],
    "05" => ["/gpfs21/scratch/cqs/guoy1/BGI_miRNA/raw/05.gz"],
    "06" => ["/gpfs21/scratch/cqs/guoy1/BGI_miRNA/raw/06.gz"],
    "07" => ["/gpfs21/scratch/cqs/guoy1/BGI_miRNA/raw/07.gz"],
    "08" => ["/gpfs21/scratch/cqs/guoy1/BGI_miRNA/raw/08.gz"],
    "09" => ["/gpfs21/scratch/cqs/guoy1/BGI_miRNA/raw/09.gz"],
    "10" => ["/gpfs21/scratch/cqs/guoy1/BGI_miRNA/raw/10.gz"],
    "11" => ["/gpfs21/scratch/cqs/guoy1/BGI_miRNA/raw/11.gz"],
    "12" => ["/gpfs21/scratch/cqs/guoy1/BGI_miRNA/raw/12.gz"],
    "13" => ["/gpfs21/scratch/cqs/guoy1/BGI_miRNA/raw/13.gz"],
    "14" => ["/gpfs21/scratch/cqs/guoy1/BGI_miRNA/raw/14.gz"],
    "15" => ["/gpfs21/scratch/cqs/guoy1/BGI_miRNA/raw/15.gz"],
    "16" => ["/gpfs21/scratch/cqs/guoy1/BGI_miRNA/raw/16.gz"],
    "17" => ["/gpfs21/scratch/cqs/guoy1/BGI_miRNA/raw/17.gz"],
    "18" => ["/gpfs21/scratch/cqs/guoy1/BGI_miRNA/raw/18.gz"],
    "19" => ["/gpfs21/scratch/cqs/guoy1/BGI_miRNA/raw/19.gz"],
    "20" => ["/gpfs21/scratch/cqs/guoy1/BGI_miRNA/raw/20.gz"],
  }
};

performSmallRNA_hg19($def_human);
1;
