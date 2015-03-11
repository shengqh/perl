#!/usr/bin/perl
use strict;
use warnings;

use CQS::PerformSmallRNA;

my $def_human = {

  #General options
  task_name            => "2687_b2",
  email                => "quanhu.sheng\@vanderbilt.edu",
  target_dir           => "/scratch/cqs/shengq1/smallRNA/20150310_guoyan_2687_human_mouse_batch2/human",
  max_thread           => 8,
  min_read_length      => 16,
  cluster              => "slurm",
  smallrnacount_option => "-s",

  #Data
  files => {
    "2687-GG-095" => ["/gpfs21/scratch/cqs/guoy1/2687/2687-GG-95_1_sequence.txt.gz"],
    "2687-GG-096" => ["/gpfs21/scratch/cqs/guoy1/2687/2687-GG-96_1_sequence.txt.gz"],
    "2687-GG-100" => ["/gpfs21/scratch/cqs/guoy1/2687/2687-GG-100_1_sequence.txt.gz"],
    "2687-GG-101" => ["/gpfs21/scratch/cqs/guoy1/2687/2687-GG-101_1_sequence.txt.gz"],
    "2687-GG-102" => ["/gpfs21/scratch/cqs/guoy1/2687/2687-GG-102_1_sequence.txt.gz"],
    "2687-GG-103" => ["/gpfs21/scratch/cqs/guoy1/2687/2687-GG-103_1_sequence.txt.gz"],
    "2687-GG-104" => ["/gpfs21/scratch/cqs/guoy1/2687/2687-GG-104_1_sequence.txt.gz"],
    "2687-GG-105" => ["/gpfs21/scratch/cqs/guoy1/2687/2687-GG-105_1_sequence.txt.gz"],
    "2687-GG-106" => ["/gpfs21/scratch/cqs/guoy1/2687/2687-GG-106_1_sequence.txt.gz"],
    "2687-GG-107" => ["/gpfs21/scratch/cqs/guoy1/2687/2687-GG-107_1_sequence.txt.gz"],
    "2687-GG-113" => ["/gpfs21/scratch/cqs/guoy1/2687/2687-GG-113_1_sequence.txt.gz"],
    "2687-GG-114" => ["/gpfs21/scratch/cqs/guoy1/2687/2687-GG-114_1_sequence.txt.gz"],
    "2687-GG-115" => ["/gpfs21/scratch/cqs/guoy1/2687/2687-GG-115_1_sequence.txt.gz"],
    "2687-GG-116" => ["/gpfs21/scratch/cqs/guoy1/2687/2687-GG-116_1_sequence.txt.gz"],
    "2687-GG-117" => ["/gpfs21/scratch/cqs/guoy1/2687/2687-GG-117_1_sequence.txt.gz"],
    "2687-GG-118" => ["/gpfs21/scratch/cqs/guoy1/2687/2687-GG-118_1_sequence.txt.gz"],
    "2687-GG-119" => ["/gpfs21/scratch/cqs/guoy1/2687/2687-GG-119_1_sequence.txt.gz"],
    "2687-GG-120" => ["/gpfs21/scratch/cqs/guoy1/2687/2687-GG-120_1_sequence.txt.gz"],
  }
};

#my $config_human = performSmallRNA_hg19($def_human);
my $config_human = getDefinition($def_human, hg19_genome());

my $mm10_genome = {

  #genome database
  mirbase_count_option  => "-p mmu",
  coordinate            => "/scratch/cqs/shengq1/references/smallrna/mm10_miRBase21_ucsc-tRNA_ensembl78.bed",
  coordinate_fasta      => "/scratch/cqs/shengq1/references/smallrna/mm10_miRBase21_ucsc-tRNA_ensembl78.bed.fa",
  bowtie1_index         => "/data/cqs/shengq1/reference/mm10/bowtie_index/mm10",
  bowtie1_miRBase_index => "/data/cqs/shengq1/reference/miRBase21/bowtie_index_1.1.1/mature.dna",
};

my $def_mouse = {

  #General options
  task_name            => "2687_b2",
  email                => "quanhu.sheng\@vanderbilt.edu",
  target_dir           => "/scratch/cqs/shengq1/smallRNA/20150310_guoyan_2687_human_mouse_batch2/mouse/",
  max_thread           => 8,
  min_read_length      => 16,
  cluster              => "slurm",
  smallrnacount_option => "-s",

  #Default software parameter (don't change it except you really know it)
  bowtie1_option_1mm => "-a -m 100 --best --strata -v 1 -p 8",
  bowtie1_option_pm  => "-a -m 100 --best --strata -v 0 -p 8",

  #Software and miRBase database options
  samtools => "/scratch/cqs/shengq1/local/bin/samtools",
  cqstools => "/home/shengq1/cqstools/CQS.Tools.exe",

  coordinate       => $mm10_genome->{coordinate},
  coordinate_fasta => $mm10_genome->{coordinate_fasta},
  bowtie1_index    => $mm10_genome->{bowtie1_index},
};

my $config_mouse = {
  general => { task_name => $def_mouse->{task_name}, },

  #not identical, for IGV
  bowtie1_genome_1mm_notidentical => {
    class             => "Bowtie1",
    perform           => 1,
    target_dir        => $def_mouse->{target_dir} . "/bowtie1_genome_1mm_notidentical",
    option            => $def_mouse->{bowtie1_option_1mm},
    source_config_ref => [ $config_human, "cutadapt", ".fastq.gz\$" ],
    bowtie1_index     => $def_mouse->{bowtie1_index},
    samonly           => 0,
    sh_direct         => 0,
    pbs               => {
      "email"    => $def_mouse->{email},
      "nodes"    => "1:ppn=" . $def_mouse->{max_thread},
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },

  #1 mismatch search, NTA
  bowtie1_genome_1mm_NTA => {
    class             => "Bowtie1",
    perform           => 1,
    target_dir        => $def_mouse->{target_dir} . "/bowtie1_genome_1mm_NTA",
    option            => $def_mouse->{bowtie1_option_1mm},
    source_config_ref => [ $config_human, "identical_NTA", ".fastq.gz\$" ],
    bowtie1_index     => $def_mouse->{bowtie1_index},
    samonly           => 0,
    sh_direct         => 1,
    cluster           => $def_mouse->{cluster},
    pbs               => {
      "email"    => $def_mouse->{email},
      "nodes"    => "1:ppn=" . $def_mouse->{max_thread},
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  bowtie1_genome_1mm_NTA_smallRNA_count => {
    class                  => "CQS::SmallRNACount",
    perform                => 1,
    target_dir             => $def_mouse->{target_dir} . "/bowtie1_genome_1mm_NTA_smallRNA_count",
    option                 => $def_mouse->{smallrnacount_option},
    source_ref             => "bowtie1_genome_1mm_NTA",
    fastq_files_config_ref => [ $config_human, "identical_NTA" ],
    seqcount_config_ref    => [ $config_human, "identical_NTA", ".dupcount\$" ],
    cqs_tools              => $def_mouse->{cqstools},
    coordinate_file        => $def_mouse->{coordinate},
    fasta_file             => $def_mouse->{coordinate_fasta},
    samtools               => $def_mouse->{samtools},
    sh_direct              => 1,
    cluster                => $def_mouse->{cluster},
    pbs                    => {
      "email"    => $def_mouse->{email},
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  bowtie1_genome_1mm_NTA_smallRNA_table => {
    class      => "CQS::SmallRNATable",
    perform    => 1,
    target_dir => $def_mouse->{target_dir} . "/bowtie1_genome_1mm_NTA_smallRNA_table",
    option     => "",
    source_ref => [ "bowtie1_genome_1mm_NTA_smallRNA_count", ".mapped.xml" ],
    cqs_tools  => $def_mouse->{cqstools},
    prefix     => "smallRNA_1mm_",
    sh_direct  => 1,
    cluster    => $def_mouse->{cluster},
    pbs        => {
      "email"    => $def_mouse->{email},
      "nodes"    => "1:ppn=1",
      "walltime" => "10",
      "mem"      => "10gb"
    },
  },
  bowtie1_genome_1mm_NTA_smallRNA_category => {
    class      => "CQS::SmallRNACategory",
    perform    => 1,
    target_dir => $def_mouse->{target_dir} . "/bowtie1_genome_1mm_NTA_smallRNA_category",
    option     => "",
    source_ref => [ "bowtie1_genome_1mm_NTA_smallRNA_count", ".info\$" ],
    cqs_tools  => $def_mouse->{cqstools},
    sh_direct  => 1,
    cluster    => $def_mouse->{cluster},
    pbs        => {
      "email"    => $def_mouse->{email},
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  sequencetask => {
    class      => "CQS::SequenceTask",
    perform    => 1,
    target_dir => $def_mouse->{target_dir} . "/sequencetask",
    option     => "",
    source     => {
      one => [ "bowtie1_genome_1mm_NTA", "bowtie1_genome_1mm_NTA_smallRNA_count", "bowtie1_genome_1mm_notidentical" ],
      all => [

        #NTA table
        "bowtie1_genome_1mm_NTA_smallRNA_table",
        "bowtie1_genome_1mm_NTA_smallRNA_category",
      ],
    },
    sh_direct => 0,
    cluster   => $def_mouse->{cluster},
    pbs       => {
      "email"    => $def_mouse->{email},
      "nodes"    => "1:ppn=" . $def_mouse->{max_thread},
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
};

performConfig($config_mouse);

1;

