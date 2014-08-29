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
  name   => "RA_",
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
  name   => "SLE_",
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
      target_dir => "${target_dir}/" . $def->{name} . "topN_bowtie1_genome_cutadapt_1mm_count_smallRNA_table",
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
      target_dir    => "${target_dir}/" . $def->{name} . "topN_bowtie1_genome_cutadapt_1mm_count_smallRNA_table_deseq2",
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

  #performConfig($config);
}

my $bowtie1_option_pm = "-a -m 100 --best --strata -v 0 -l 12 -p 8";

my $config = {
  general     => { "task_name" => $task_name },
  fastq_files => {
    "2572-KCV-1-19" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/identical/result/2572-KCV-1-19_clipped_identical.fastq.gz"],
    "2572-KCV-1-20" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/identical/result/2572-KCV-1-20_clipped_identical.fastq.gz"],
    "2572-KCV-1-21" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/identical/result/2572-KCV-1-21_clipped_identical.fastq.gz"],
    "2572-KCV-1-22" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/identical/result/2572-KCV-1-22_clipped_identical.fastq.gz"],
    "2572-KCV-1-23" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/identical/result/2572-KCV-1-23_clipped_identical.fastq.gz"],
    "2572-KCV-1-24" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/identical/result/2572-KCV-1-24_clipped_identical.fastq.gz"],
    "2572-KCV-1-25" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/identical/result/2572-KCV-1-25_clipped_identical.fastq.gz"],
    "2572-KCV-1-26" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/identical/result/2572-KCV-1-26_clipped_identical.fastq.gz"],
    "2572-KCV-1-27" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/identical/result/2572-KCV-1-27_clipped_identical.fastq.gz"],
    "2572-KCV-1-28" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/identical/result/2572-KCV-1-28_clipped_identical.fastq.gz"],
    "2572-KCV-1-29" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/identical/result/2572-KCV-1-29_clipped_identical.fastq.gz"],
    "2572-KCV-1-30" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/identical/result/2572-KCV-1-30_clipped_identical.fastq.gz"],
  },
  bam_files => {
    "2572-KCV-1-19" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm/result/2572-KCV-1-19/2572-KCV-1-19.bam"],
    "2572-KCV-1-20" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm/result/2572-KCV-1-20/2572-KCV-1-20.bam"],
    "2572-KCV-1-21" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm/result/2572-KCV-1-21/2572-KCV-1-21.bam"],
    "2572-KCV-1-22" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm/result/2572-KCV-1-22/2572-KCV-1-22.bam"],
    "2572-KCV-1-23" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm/result/2572-KCV-1-23/2572-KCV-1-23.bam"],
    "2572-KCV-1-24" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm/result/2572-KCV-1-24/2572-KCV-1-24.bam"],
    "2572-KCV-1-25" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm/result/2572-KCV-1-25/2572-KCV-1-25.bam"],
    "2572-KCV-1-26" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm/result/2572-KCV-1-26/2572-KCV-1-26.bam"],
    "2572-KCV-1-27" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm/result/2572-KCV-1-27/2572-KCV-1-27.bam"],
    "2572-KCV-1-28" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm/result/2572-KCV-1-28/2572-KCV-1-28.bam"],
    "2572-KCV-1-29" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm/result/2572-KCV-1-29/2572-KCV-1-29.bam"],
    "2572-KCV-1-30" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm/result/2572-KCV-1-30/2572-KCV-1-30.bam"],
  },
  count_files => {
    "2572-KCV-1-19" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/identical/result/2572-KCV-1-19_clipped_identical.fastq.dupcount"],
    "2572-KCV-1-20" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/identical/result/2572-KCV-1-20_clipped_identical.fastq.dupcount"],
    "2572-KCV-1-21" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/identical/result/2572-KCV-1-21_clipped_identical.fastq.dupcount"],
    "2572-KCV-1-22" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/identical/result/2572-KCV-1-22_clipped_identical.fastq.dupcount"],
    "2572-KCV-1-23" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/identical/result/2572-KCV-1-23_clipped_identical.fastq.dupcount"],
    "2572-KCV-1-24" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/identical/result/2572-KCV-1-24_clipped_identical.fastq.dupcount"],
    "2572-KCV-1-25" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/identical/result/2572-KCV-1-25_clipped_identical.fastq.dupcount"],
    "2572-KCV-1-26" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/identical/result/2572-KCV-1-26_clipped_identical.fastq.dupcount"],
    "2572-KCV-1-27" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/identical/result/2572-KCV-1-27_clipped_identical.fastq.dupcount"],
    "2572-KCV-1-28" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/identical/result/2572-KCV-1-28_clipped_identical.fastq.dupcount"],
    "2572-KCV-1-29" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/identical/result/2572-KCV-1-29_clipped_identical.fastq.dupcount"],
    "2572-KCV-1-30" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/identical/result/2572-KCV-1-30_clipped_identical.fastq.dupcount"],
  },

  #2 perfect match search to mirbase only
  bowtie1_genome_cutadapt_topN_genome_pmnames => {
    class      => "Samtools::PerfectMappedReadNames",
    perform    => 0,
    target_dir => "${target_dir}/topN_bowtie1_genome_cutadapt_pmnames",
    option     => "",
    source_ref => "bam_files",
    sh_direct  => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  bowtie1_genome_cutadapt_topN_miRbase_pm => {
    class         => "Alignment::Bowtie1",
    perform       => 0,
    target_dir    => "${target_dir}/topN_bowtie1_genome_cutadapt_miRbase_pm",
    option        => $bowtie1_option_pm,
    source_ref    => "fastq_files",
    bowtie1_index => "/data/cqs/shengq1/reference/miRBase21/bowtie_index_1.0.1/mature.dna",
    samonly       => 0,
    sh_direct     => 1,
    pbs           => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  chromosome_count => {
    class                   => "CQS::CQSChromosomeCount",
    perform                 => 1,
    target_dir              => "${target_dir}/topN_bowtie1_genome_cutadapt_miRbase_pm_count",
    option                  => "-p hsa",
    source_ref              => "bowtie1_genome_cutadapt_topN_miRbase_pm",
    seqcount_ref            => "count_files",
    perfect_mapped_name_ref => "bowtie1_genome_cutadapt_topN_genome_pmnames",
    cqs_tools               => $cqstools,
    samtools                => $samtools,
    sh_direct               => 1,
    pbs                     => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    }
  },
};

performConfig($config);

1;
