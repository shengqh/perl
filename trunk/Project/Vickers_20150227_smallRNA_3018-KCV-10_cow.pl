#!/usr/bin/perl
use strict;
use warnings;

use CQS::PerformSmallRNA;

my $def = {

  #General options
  task_name  => "3018-KCV-10",
  email      => "quanhu.sheng\@vanderbilt.edu",
  target_dir => "/scratch/cqs/shengq1/vickers/20150227_smallRNA_3018-KCV-10_cow",
  max_thread => 8,

  #Default software parameter (don't change it except you really know it)
  bowtie1_option_1mm         => "-a -m 100 --best --strata -v 1 -p 8",
  bowtie1_option_pm          => "-a -m 100 --best --strata -v 0 -p 8",
  mirnacount_option          => "-s",                                         #ignore score
  smallrnacount_option       => "-s --min_overlap 0.5 --length --sequence",
  mirna_overlap_count_option => "-s --min_overlap 0.5 --gtf_key miRNA",

  #Software and miRBase database options
  samtools => "/scratch/cqs/shengq1/local/bin/samtools",
  cqstools => "/home/shengq1/cqstools/CQS.Tools.exe",

  #genome database
  bowtie1_index => "/scratch/cqs/shengq1/references/cow/bowtie_index_1.1.1/Bos_taurus_UMD_3.1",

};

my $cluster = "slurm";

my $config = {
  general => {
    task_name => $def->{task_name},
    cluster   => $cluster,
  },
  cutadapt => {
    "3018-KCV-10-37" => ["/gpfs21/scratch/cqs/shengq1/vickers/20150224_smallRNA_3018-KCV-10_human/cutadapt/result/3018-KCV-10-37_clipped.fastq.gz"],
    "3018-KCV-10-38" => ["/gpfs21/scratch/cqs/shengq1/vickers/20150224_smallRNA_3018-KCV-10_human/cutadapt/result/3018-KCV-10-38_clipped.fastq.gz"],
    "3018-KCV-10-39" => ["/gpfs21/scratch/cqs/shengq1/vickers/20150224_smallRNA_3018-KCV-10_human/cutadapt/result/3018-KCV-10-39_clipped.fastq.gz"],
    "3018-KCV-10-40" => ["/gpfs21/scratch/cqs/shengq1/vickers/20150224_smallRNA_3018-KCV-10_human/cutadapt/result/3018-KCV-10-40_clipped.fastq.gz"],
    "3018-KCV-10-41" => ["/gpfs21/scratch/cqs/shengq1/vickers/20150224_smallRNA_3018-KCV-10_human/cutadapt/result/3018-KCV-10-41_clipped.fastq.gz"],
    "3018-KCV-10-42" => ["/gpfs21/scratch/cqs/shengq1/vickers/20150224_smallRNA_3018-KCV-10_human/cutadapt/result/3018-KCV-10-42_clipped.fastq.gz"],
    "3018-KCV-10-43" => ["/gpfs21/scratch/cqs/shengq1/vickers/20150224_smallRNA_3018-KCV-10_human/cutadapt/result/3018-KCV-10-43_clipped.fastq.gz"],
    "3018-KCV-10-44" => ["/gpfs21/scratch/cqs/shengq1/vickers/20150224_smallRNA_3018-KCV-10_human/cutadapt/result/3018-KCV-10-44_clipped.fastq.gz"],
    "3018-KCV-10-45" => ["/gpfs21/scratch/cqs/shengq1/vickers/20150224_smallRNA_3018-KCV-10_human/cutadapt/result/3018-KCV-10-45_clipped.fastq.gz"],
    "3018-KCV-10-46" => ["/gpfs21/scratch/cqs/shengq1/vickers/20150224_smallRNA_3018-KCV-10_human/cutadapt/result/3018-KCV-10-46_clipped.fastq.gz"],
    "3018-KCV-10-47" => ["/gpfs21/scratch/cqs/shengq1/vickers/20150224_smallRNA_3018-KCV-10_human/cutadapt/result/3018-KCV-10-47_clipped.fastq.gz"],
    "3018-KCV-10-48" => ["/gpfs21/scratch/cqs/shengq1/vickers/20150224_smallRNA_3018-KCV-10_human/cutadapt/result/3018-KCV-10-48_clipped.fastq.gz"],
  },
  identical => {
    "3018-KCV-10-37" => ["/gpfs21/scratch/cqs/shengq1/vickers/20150224_smallRNA_3018-KCV-10_human/identical/result/3018-KCV-10-37_clipped_identical.fastq.gz"],
    "3018-KCV-10-38" => ["/gpfs21/scratch/cqs/shengq1/vickers/20150224_smallRNA_3018-KCV-10_human/identical/result/3018-KCV-10-38_clipped_identical.fastq.gz"],
    "3018-KCV-10-39" => ["/gpfs21/scratch/cqs/shengq1/vickers/20150224_smallRNA_3018-KCV-10_human/identical/result/3018-KCV-10-39_clipped_identical.fastq.gz"],
    "3018-KCV-10-40" => ["/gpfs21/scratch/cqs/shengq1/vickers/20150224_smallRNA_3018-KCV-10_human/identical/result/3018-KCV-10-40_clipped_identical.fastq.gz"],
    "3018-KCV-10-41" => ["/gpfs21/scratch/cqs/shengq1/vickers/20150224_smallRNA_3018-KCV-10_human/identical/result/3018-KCV-10-41_clipped_identical.fastq.gz"],
    "3018-KCV-10-42" => ["/gpfs21/scratch/cqs/shengq1/vickers/20150224_smallRNA_3018-KCV-10_human/identical/result/3018-KCV-10-42_clipped_identical.fastq.gz"],
    "3018-KCV-10-43" => ["/gpfs21/scratch/cqs/shengq1/vickers/20150224_smallRNA_3018-KCV-10_human/identical/result/3018-KCV-10-43_clipped_identical.fastq.gz"],
    "3018-KCV-10-44" => ["/gpfs21/scratch/cqs/shengq1/vickers/20150224_smallRNA_3018-KCV-10_human/identical/result/3018-KCV-10-44_clipped_identical.fastq.gz"],
    "3018-KCV-10-45" => ["/gpfs21/scratch/cqs/shengq1/vickers/20150224_smallRNA_3018-KCV-10_human/identical/result/3018-KCV-10-45_clipped_identical.fastq.gz"],
    "3018-KCV-10-46" => ["/gpfs21/scratch/cqs/shengq1/vickers/20150224_smallRNA_3018-KCV-10_human/identical/result/3018-KCV-10-46_clipped_identical.fastq.gz"],
    "3018-KCV-10-47" => ["/gpfs21/scratch/cqs/shengq1/vickers/20150224_smallRNA_3018-KCV-10_human/identical/result/3018-KCV-10-47_clipped_identical.fastq.gz"],
    "3018-KCV-10-48" => ["/gpfs21/scratch/cqs/shengq1/vickers/20150224_smallRNA_3018-KCV-10_human/identical/result/3018-KCV-10-48_clipped_identical.fastq.gz"],
  },

  #perfect match search to cow only
  bowtie1_genome_1mm_NTA_pmnames => {
    "3018-KCV-10-37" => ["/gpfs21/scratch/cqs/shengq1/vickers/20150224_smallRNA_3018-KCV-10_human/bowtie1_genome_1mm_NTA_pmnames/result/3018-KCV-10-37.pmnames"],
    "3018-KCV-10-38" => ["/gpfs21/scratch/cqs/shengq1/vickers/20150224_smallRNA_3018-KCV-10_human/bowtie1_genome_1mm_NTA_pmnames/result/3018-KCV-10-38.pmnames"],
    "3018-KCV-10-39" => ["/gpfs21/scratch/cqs/shengq1/vickers/20150224_smallRNA_3018-KCV-10_human/bowtie1_genome_1mm_NTA_pmnames/result/3018-KCV-10-39.pmnames"],
    "3018-KCV-10-40" => ["/gpfs21/scratch/cqs/shengq1/vickers/20150224_smallRNA_3018-KCV-10_human/bowtie1_genome_1mm_NTA_pmnames/result/3018-KCV-10-40.pmnames"],
    "3018-KCV-10-41" => ["/gpfs21/scratch/cqs/shengq1/vickers/20150224_smallRNA_3018-KCV-10_human/bowtie1_genome_1mm_NTA_pmnames/result/3018-KCV-10-41.pmnames"],
    "3018-KCV-10-42" => ["/gpfs21/scratch/cqs/shengq1/vickers/20150224_smallRNA_3018-KCV-10_human/bowtie1_genome_1mm_NTA_pmnames/result/3018-KCV-10-42.pmnames"],
    "3018-KCV-10-43" => ["/gpfs21/scratch/cqs/shengq1/vickers/20150224_smallRNA_3018-KCV-10_human/bowtie1_genome_1mm_NTA_pmnames/result/3018-KCV-10-43.pmnames"],
    "3018-KCV-10-44" => ["/gpfs21/scratch/cqs/shengq1/vickers/20150224_smallRNA_3018-KCV-10_human/bowtie1_genome_1mm_NTA_pmnames/result/3018-KCV-10-44.pmnames"],
    "3018-KCV-10-45" => ["/gpfs21/scratch/cqs/shengq1/vickers/20150224_smallRNA_3018-KCV-10_human/bowtie1_genome_1mm_NTA_pmnames/result/3018-KCV-10-45.pmnames"],
    "3018-KCV-10-46" => ["/gpfs21/scratch/cqs/shengq1/vickers/20150224_smallRNA_3018-KCV-10_human/bowtie1_genome_1mm_NTA_pmnames/result/3018-KCV-10-46.pmnames"],
    "3018-KCV-10-47" => ["/gpfs21/scratch/cqs/shengq1/vickers/20150224_smallRNA_3018-KCV-10_human/bowtie1_genome_1mm_NTA_pmnames/result/3018-KCV-10-47.pmnames"],
    "3018-KCV-10-48" => ["/gpfs21/scratch/cqs/shengq1/vickers/20150224_smallRNA_3018-KCV-10_human/bowtie1_genome_1mm_NTA_pmnames/result/3018-KCV-10-48.pmnames"],
  },
  bowtie1_cow_pm => {
    class         => "Alignment::Bowtie1",
    perform       => 1,
    target_dir    => $def->{target_dir} . "/bowtie1_cow_pm",
    option        => $def->{bowtie1_option_pm},
    source_ref    => "cutadapt",
    bowtie1_index => $def->{bowtie1_index},
    samonly       => 0,
    sh_direct     => 0,
    cluster       => $cluster,
    pbs           => {
      "email"    => $def->{email},
      "nodes"    => "1:ppn=" . $def->{max_thread},
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
};

performConfig($config);

1;

