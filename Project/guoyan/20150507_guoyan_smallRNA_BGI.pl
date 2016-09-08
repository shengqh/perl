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
    "01" => ["/gpfs21/scratch/cqs/shengq1/smallRNA/20150507_guoyan_smallRNA_BGI/raw/01.fastq.gz"],
    "02" => ["/gpfs21/scratch/cqs/shengq1/smallRNA/20150507_guoyan_smallRNA_BGI/raw/02.fastq.gz"],
    "03" => ["/gpfs21/scratch/cqs/shengq1/smallRNA/20150507_guoyan_smallRNA_BGI/raw/03.fastq.gz"],
    "04" => ["/gpfs21/scratch/cqs/shengq1/smallRNA/20150507_guoyan_smallRNA_BGI/raw/04.fastq.gz"],
    "05" => ["/gpfs21/scratch/cqs/shengq1/smallRNA/20150507_guoyan_smallRNA_BGI/raw/05.fastq.gz"],
    "06" => ["/gpfs21/scratch/cqs/shengq1/smallRNA/20150507_guoyan_smallRNA_BGI/raw/06.fastq.gz"],
    "07" => ["/gpfs21/scratch/cqs/shengq1/smallRNA/20150507_guoyan_smallRNA_BGI/raw/07.fastq.gz"],
    "08" => ["/gpfs21/scratch/cqs/shengq1/smallRNA/20150507_guoyan_smallRNA_BGI/raw/08.fastq.gz"],
    "09" => ["/gpfs21/scratch/cqs/shengq1/smallRNA/20150507_guoyan_smallRNA_BGI/raw/09.fastq.gz"],
    "10" => ["/gpfs21/scratch/cqs/shengq1/smallRNA/20150507_guoyan_smallRNA_BGI/raw/10.fastq.gz"],
    "11" => ["/gpfs21/scratch/cqs/shengq1/smallRNA/20150507_guoyan_smallRNA_BGI/raw/11.fastq.gz"],
    "12" => ["/gpfs21/scratch/cqs/shengq1/smallRNA/20150507_guoyan_smallRNA_BGI/raw/12.fastq.gz"],
    "13" => ["/gpfs21/scratch/cqs/shengq1/smallRNA/20150507_guoyan_smallRNA_BGI/raw/13.fastq.gz"],
    "14" => ["/gpfs21/scratch/cqs/shengq1/smallRNA/20150507_guoyan_smallRNA_BGI/raw/14.fastq.gz"],
    "15" => ["/gpfs21/scratch/cqs/shengq1/smallRNA/20150507_guoyan_smallRNA_BGI/raw/15.fastq.gz"],
    "16" => ["/gpfs21/scratch/cqs/shengq1/smallRNA/20150507_guoyan_smallRNA_BGI/raw/16.fastq.gz"],
    "17" => ["/gpfs21/scratch/cqs/shengq1/smallRNA/20150507_guoyan_smallRNA_BGI/raw/17.fastq.gz"],
    "18" => ["/gpfs21/scratch/cqs/shengq1/smallRNA/20150507_guoyan_smallRNA_BGI/raw/18.fastq.gz"],
    "19" => ["/gpfs21/scratch/cqs/shengq1/smallRNA/20150507_guoyan_smallRNA_BGI/raw/19.fastq.gz"],
    "20" => ["/gpfs21/scratch/cqs/shengq1/smallRNA/20150507_guoyan_smallRNA_BGI/raw/20.fastq.gz"],
    }
};

performSmallRNA_hg19($def_human);
1;
