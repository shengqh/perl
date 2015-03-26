#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;
use CQS::PerformSmallRNA;

my $def = {

  #General options
  task_name  => "parclip_NIH",
  email      => "quanhu.sheng\@vanderbilt.edu",
  target_dir => "/scratch/cqs/shengq1/vickers/20150312_parclip_3018-KCV-15/",
  max_thread => 8,
  cluster    => "slurm",

  #Default software parameter (don't change it except you really know it)
  fastq_remove_N => 0,
  adapter        => "TGGAATTCTCGGGTGCCAAGG",

  binding_db => "/data/cqs/shengq1/reference/targetscan/targetscan_v61_hg19.bed",
  utr3_db    => "/data/cqs/shengq1/reference/utr3/20140612_ucsc_hg19_3UTR.txt",
  samtools   => "/scratch/cqs/shengq1/local/bin/samtools",
  cqstools   => "/home/shengq1/cqstools/CQS.Tools.exe",

  #for parclip target
  fasta_file   => "/gpfs21/scratch/cqs/shengq1/references/hg19_16569_M/hg19_16569_M.fa",
  refgene_file => "/gpfs21/scratch/cqs/shengq1/references/hg19_16569_M/hg19_refgene.tsv",

  #Data
  files => {
    "3018-KCV-15-15" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-15_parclip/3018-KCV-15-15_ATGTCA_L006_R1_001.fastq.gz"],
    "3018-KCV-15-36" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-15_parclip/3018-KCV-15-36_CCAACA_L006_R1_001.fastq.gz"],
    "3018-KCV-15-37" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-15_parclip/3018-KCV-15-37_CGGAAT_L006_R1_001.fastq.gz"],
    "3018-KCV-15-46" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-15_parclip/3018-KCV-15-46_TCCCGA_L006_R1_001.fastq.gz"],
    "3018-KCV-15-47" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-15_parclip/3018-KCV-15-47_TCGAAG_L006_R1_001.fastq.gz"],
  },
};

my $hg19Config = performSmallRNA_hg19( $def, 0 );

my $config = {
  general => { "task_name" => $def->{task_name}, },
  t2c     => {
    class             => "CQS::ParclipT2CFinder",
    perform           => 0,
    target_dir        => $def->{target_dir} . "/t2c_finder",
    option            => "-p 0.05 -e 0.013",
    source_config_ref => [ $hg19Config, "bowtie1_genome_1mm_NTA_smallRNA_count", ".mapped.xml\$" ],
    cqs_tools         => $def->{cqstools},
    sh_direct         => 1,
    pbs               => {
      "email"    => $def->{email},
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "20gb"
    },
  },
  utr3_count => {
    class                  => "CQS::SmallRNACount",
    perform                => 0,
    target_dir             => $def->{target_dir} . "/count_3utr",
    option                 => "-m 0",
    source_config_ref      => [ $hg19Config, "bowtie1_genome_1mm_NTA", ".bam\$" ],
    fastq_files_config_ref => [ $hg19Config, "identical_NTA", "fastq.gz\$" ],
    seqcount_config_ref    => [ $hg19Config, "identical_NTA", ".dupcount\$" ],
    cqs_tools              => $def->{cqstools},
    coordinate_file        => $def->{utr3_db},
    samtools               => $def->{samtools},
    sh_direct              => 1,
    pbs                    => {
      "email"    => $def->{email},
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "20gb"
    },
  },
  mirna_target => {
    class        => "CQS::ParclipMirnaTarget",
    perform      => 1,
    target_dir   => $def->{target_dir} . "/mirna_3utr_target",
    option       => "",
    source_ref   => [ "t2c", ".xml\$" ],
    target_ref   => [ "utr3_count", ".xml\$" ],
    fasta_file   => $def->{fasta_file},
    refgene_file => $def->{refgene_file},
    cqs_tools    => $def->{cqstools},
    sh_direct    => 1,
    pbs          => {
      "email"    => $def->{email},
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "20gb"
    },
  },

  binding_count => {
    class                  => "CQS::SmallRNACount",
    perform                => 0,
    target_dir             => $def->{target_dir} . "/count_binding",
    option                 => "-m 0",
    source_config_ref      => [ $hg19Config, "bowtie1_genome_1mm_NTA", ".bam\$" ],
    fastq_files_config_ref => [ $hg19Config, "identical_NTA", "fastq.gz\$" ],
    seqcount_config_ref    => [ $hg19Config, "identical_NTA", ".dupcount\$" ],
    cqs_tools              => $def->{cqstools},
    coordinate_file        => $def->{binding_db},
    samtools               => $def->{samtools},
    sh_direct              => 1,
    pbs                    => {
      "email"    => $def->{email},
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "20gb"
    },
  },
  sequencetask => {
    class      => "CQS::SequenceTask",
    perform    => 0,
    target_dir => $def->{target_dir} . "/sequencetask2",
    option     => "",
    source     => { count => [ "utr3_count", "binding_count" ], },
    sh_direct  => 1,
    pbs        => {
      "email"    => $def->{email},
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "20gb"
    },
  },
};

performConfig($config);

1;
