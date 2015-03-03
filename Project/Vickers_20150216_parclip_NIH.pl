#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;
use CQS::PerformSmallRNA;

my $def = {

  #General options
  task_name  => "parclip_NIH",
  email      => "quanhu.sheng\@vanderbilt.edu",
  target_dir => "/scratch/cqs/shengq1/vickers/20150216_parclip_NIH/",
  max_thread => 8,
  cluster    => "slurm",

  #Default software parameter (don't change it except you really know it)
  fastq_remove_N => 0,
  adapter        => "TGGAATTCTCGGGTGCCAAGG",

  binding_db => "/data/cqs/shengq1/reference/targetscan/targetscan_v61_hg19.bed",
  utr3_db    => "/data/cqs/shengq1/reference/utr3/20140612_ucsc_hg19_3UTR.txt",
  samtools   => "/scratch/cqs/shengq1/local/bin/samtools",
  cqstools   => "/home/shengq1/cqstools/CQS.Tools.exe",

  #Data
  files => {
    "Parclip_01" => ["/scratch/cqs/shengq1/vickers/data/201312_parclip_NIH/Vickers_Parclip_1_ATCACG_L002_R1.fastq.gz"],
    "Parclip_02" => ["/scratch/cqs/shengq1/vickers/data/201312_parclip_NIH/Vickers_Parclip_2_CGATGT_L002_R1.fastq.gz"],
    "Parclip_03" => ["/scratch/cqs/shengq1/vickers/data/201312_parclip_NIH/Vickers_Parclip_3_TTAGGC_L002_R1.fastq.gz"],
    "Parclip_04" => ["/scratch/cqs/shengq1/vickers/data/201312_parclip_NIH/Vickers_Parclip_4_TGACCA_L002_R1.fastq.gz"],
    "Parclip_05" => ["/scratch/cqs/shengq1/vickers/data/201312_parclip_NIH/Vickers_Parclip_5_ACAGTG_L002_R1.fastq.gz"],
    "Parclip_06" => ["/scratch/cqs/shengq1/vickers/data/201312_parclip_NIH/Vickers_Parclip_6_GCCAAT_L002_R1.fastq.gz"],
    "Parclip_07" => ["/scratch/cqs/shengq1/vickers/data/201312_parclip_NIH/Vickers_Parclip_7_CAGATC_L002_R1.fastq.gz"],
    "Parclip_08" => ["/scratch/cqs/shengq1/vickers/data/201312_parclip_NIH/Vickers_Parclip_8_ACTTGA_L002_R1.fastq.gz"],
  },
};

#performSmallRNA_hg19($def);

my $config = {
  general => { "task_name" => $def->{task_name}, },
  identificalfiles => {
    "Parclip_01" => ["/gpfs21/scratch/cqs/shengq1/vickers/20150216_parclip_NIH/identical/result/Parclip_01_clipped_identical.fastq.gz"],
    "Parclip_02" => ["/gpfs21/scratch/cqs/shengq1/vickers/20150216_parclip_NIH/identical/result/Parclip_02_clipped_identical.fastq.gz"],
    "Parclip_03" => ["/gpfs21/scratch/cqs/shengq1/vickers/20150216_parclip_NIH/identical/result/Parclip_03_clipped_identical.fastq.gz"],
    "Parclip_04" => ["/gpfs21/scratch/cqs/shengq1/vickers/20150216_parclip_NIH/identical/result/Parclip_04_clipped_identical.fastq.gz"],
    "Parclip_05" => ["/gpfs21/scratch/cqs/shengq1/vickers/20150216_parclip_NIH/identical/result/Parclip_05_clipped_identical.fastq.gz"],
    "Parclip_06" => ["/gpfs21/scratch/cqs/shengq1/vickers/20150216_parclip_NIH/identical/result/Parclip_06_clipped_identical.fastq.gz"],
    "Parclip_07" => ["/gpfs21/scratch/cqs/shengq1/vickers/20150216_parclip_NIH/identical/result/Parclip_07_clipped_identical.fastq.gz"],
    "Parclip_08" => ["/gpfs21/scratch/cqs/shengq1/vickers/20150216_parclip_NIH/identical/result/Parclip_08_clipped_identical.fastq.gz"],
  },
  countfiles => {
    "Parclip_01" => ["/gpfs21/scratch/cqs/shengq1/vickers/20150216_parclip_NIH/identical/result/Parclip_01_clipped_identical.fastq.dupcount"],
    "Parclip_02" => ["/gpfs21/scratch/cqs/shengq1/vickers/20150216_parclip_NIH/identical/result/Parclip_02_clipped_identical.fastq.dupcount"],
    "Parclip_03" => ["/gpfs21/scratch/cqs/shengq1/vickers/20150216_parclip_NIH/identical/result/Parclip_03_clipped_identical.fastq.dupcount"],
    "Parclip_04" => ["/gpfs21/scratch/cqs/shengq1/vickers/20150216_parclip_NIH/identical/result/Parclip_04_clipped_identical.fastq.dupcount"],
    "Parclip_05" => ["/gpfs21/scratch/cqs/shengq1/vickers/20150216_parclip_NIH/identical/result/Parclip_05_clipped_identical.fastq.dupcount"],
    "Parclip_06" => ["/gpfs21/scratch/cqs/shengq1/vickers/20150216_parclip_NIH/identical/result/Parclip_06_clipped_identical.fastq.dupcount"],
    "Parclip_07" => ["/gpfs21/scratch/cqs/shengq1/vickers/20150216_parclip_NIH/identical/result/Parclip_07_clipped_identical.fastq.dupcount"],
    "Parclip_08" => ["/gpfs21/scratch/cqs/shengq1/vickers/20150216_parclip_NIH/identical/result/Parclip_08_clipped_identical.fastq.dupcount"],
  },
  bamfiles => {
    "Parclip_01" => ["/gpfs21/scratch/cqs/shengq1/vickers/20150216_parclip_NIH/bowtie1_genome_1mm_NTA/result/Parclip_01/Parclip_01.bam"],
    "Parclip_02" => ["/gpfs21/scratch/cqs/shengq1/vickers/20150216_parclip_NIH/bowtie1_genome_1mm_NTA/result/Parclip_02/Parclip_02.bam"],
    "Parclip_03" => ["/gpfs21/scratch/cqs/shengq1/vickers/20150216_parclip_NIH/bowtie1_genome_1mm_NTA/result/Parclip_03/Parclip_03.bam"],
    "Parclip_04" => ["/gpfs21/scratch/cqs/shengq1/vickers/20150216_parclip_NIH/bowtie1_genome_1mm_NTA/result/Parclip_04/Parclip_04.bam"],
    "Parclip_05" => ["/gpfs21/scratch/cqs/shengq1/vickers/20150216_parclip_NIH/bowtie1_genome_1mm_NTA/result/Parclip_05/Parclip_05.bam"],
    "Parclip_06" => ["/gpfs21/scratch/cqs/shengq1/vickers/20150216_parclip_NIH/bowtie1_genome_1mm_NTA/result/Parclip_06/Parclip_06.bam"],
    "Parclip_07" => ["/gpfs21/scratch/cqs/shengq1/vickers/20150216_parclip_NIH/bowtie1_genome_1mm_NTA/result/Parclip_07/Parclip_07.bam"],
    "Parclip_08" => ["/gpfs21/scratch/cqs/shengq1/vickers/20150216_parclip_NIH/bowtie1_genome_1mm_NTA/result/Parclip_08/Parclip_08.bam"],
  },
  utr3_count => {
    class           => "CQSMappedCount",
    perform         => 1,
    target_dir      => $def->{target_dir} . "/count_3utr",
    option          => "-m 0",
    source_ref      => "bamfiles",
    fastq_files_ref => "identicalfiles",
    seqcount_ref    => "countfiles",
    cqs_tools       => $def->{cqstools},
    gff_file        => $def->{utr3_db},
    samtools        => $def->{samtools},
    sh_direct       => 1,
    pbs             => {
      "email"    => $def->{email},
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "20gb"
    },
  },
  binding_count => {
    class           => "CQSMappedCount",
    perform         => 1,
    target_dir      => $def->{target_dir} . "/count_binding",
    option          => "-m 0",
    source_ref      => "bowtie1bam",
    fastq_files_ref => "identicalfiles",
    seqcount_ref    => "countfiles",
    cqs_tools       => $def->{cqstools},
    gff_file        => $def->{binding_db},
    samtools        => $def->{samtools},
    sh_direct       => 1,
    pbs             => {
      "email"    => $def->{email},
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "20gb"
    },
  }
};

performConfig($config);

1;
