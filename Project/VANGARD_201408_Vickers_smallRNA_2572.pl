#!/usr/bin/perl
use strict;
use warnings;

use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;

my $target_dir = create_directory_or_die("/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_2572/");
my $cqstools   = "/home/shengq1/cqstools/CQS.Tools.exe";
my $samtools   = "/home/shengq1/local/bin/samtools/samtools";

my $email     = "quanhu.sheng\@vanderbilt.edu";
my $task_name = "2572";

my $raconfig = {
  source => {
    "2572-KCV-1-19" =>
      ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_smallRNA/result/2572-KCV-1-19/2572-KCV-1-19.bam.count.mapped.xml"],
    "2572-KCV-1-20" =>
      ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_smallRNA/result/2572-KCV-1-20/2572-KCV-1-20.bam.count.mapped.xml"],
    "2572-KCV-1-21" =>
      ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_smallRNA/result/2572-KCV-1-21/2572-KCV-1-21.bam.count.mapped.xml"],
    "2572-KCV-1-22" =>
      ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_smallRNA/result/2572-KCV-1-22/2572-KCV-1-22.bam.count.mapped.xml"],
    "2572-KCV-1-23" =>
      ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_smallRNA/result/2572-KCV-1-23/2572-KCV-1-23.bam.count.mapped.xml"],
    "2572-KCV-1-24" =>
      ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_smallRNA/result/2572-KCV-1-24/2572-KCV-1-24.bam.count.mapped.xml"],
  },
  groups => {
    "RA_Control" => [ "2572-KCV-1-19", "2572-KCV-1-20", "2572-KCV-1-21" ],
    "RA_Sample"  => [ "2572-KCV-1-22", "2572-KCV-1-23", "2572-KCV-1-24" ],
  },
  pairs => { "RA" => [ "RA_Control", "RA_Sample" ], },
};

my $sleconfig = {
  source => {
    "2572-KCV-1-25" =>
      ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_smallRNA/result/2572-KCV-1-25/2572-KCV-1-25.bam.count.mapped.xml"],
    "2572-KCV-1-26" =>
      ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_smallRNA/result/2572-KCV-1-26/2572-KCV-1-26.bam.count.mapped.xml"],
    "2572-KCV-1-27" =>
      ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_smallRNA/result/2572-KCV-1-27/2572-KCV-1-27.bam.count.mapped.xml"],
    "2572-KCV-1-28" =>
      ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_smallRNA/result/2572-KCV-1-28/2572-KCV-1-28.bam.count.mapped.xml"],
    "2572-KCV-1-29" =>
      ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_smallRNA/result/2572-KCV-1-29/2572-KCV-1-29.bam.count.mapped.xml"],
    "2572-KCV-1-30" =>
      ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_smallRNA/result/2572-KCV-1-30/2572-KCV-1-30.bam.count.mapped.xml"],
  },
  groups => {
    "SLE_Control" => [ "2572-KCV-1-25", "2572-KCV-1-26", "2572-KCV-1-27" ],
    "SLE_Sample"  => [ "2572-KCV-1-28", "2572-KCV-1-29", "2572-KCV-1-30" ],
  },
  pairs => { "SLE" => [ "SLE_Control", "SLE_Sample" ], },
};

my @defs = ( $raconfig, $sleconfig );

foreach my $def (@defs) {
  my $config = {
    general            => { "task_name" => $task_name },
    smallRNA_1mm_table => {
      class      => "CQSMappedTable",
      perform    => 1,
      target_dir => "${target_dir}/topN_bowtie1_genome_cutadapt_1mm_count_smallRNA_table",
      option     => "",
      source     => $def->{source},
      cqs_tools  => $cqstools,
      sh_direct  => 1,
      pbs        => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "10",
        "mem"      => "10gb"
      },
    },
    deseq2 => {
      class         => "DESeq2",
      perform       => 1,
      target_dir    => "${target_dir}/deseq2",
      option        => "",
      source_ref    => $def->{pairs},
      groups_ref    => $def->{groups},
      countfile_ref => "smallRNA_1mm_table",
      sh_direct     => 1,
      pbs           => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "10",
        "mem"      => "10gb"
      },
    },
  };

  performConfig($config);
}

1;
