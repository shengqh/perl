#!/usr/bin/perl
use strict;
use warnings;

use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;

my $root     = create_directory_or_die("/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2");
my $cqstools = "/home/shengq1/cqstools/CQS.Tools.exe";

my $hg19_mrna_gff     = "/data/cqs/shengq1/reference/miRBase20/hsa.gff3";
my $hg19_trna_bed     = "/data/cqs/guoy1/reference/smallrna/hg19_tRNA_ucsc_ensembl.bed";

my $target_dir = $root;

my $email     = "quanhu.sheng\@vanderbilt.edu";
my $task_name = "VANGARD00055";

my $bowtie1_human_index = "/data/cqs/guoy1/reference/hg19/bowtie_index_hg19_rCRS_1.0.0/hg19_rCRS";
my $human_2bit = "/data/cqs/guoy1/reference/hg19/hg19_rCRS.2bit";

my $parclip_files = {
  "Parclip_01" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/PARCLIP/Vickers_Parclip_1_ATCACG_L002_R1.fastq.gz"],
  "Parclip_02" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/PARCLIP/Vickers_Parclip_2_CGATGT_L002_R1.fastq.gz"],
  "Parclip_03" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/PARCLIP/Vickers_Parclip_3_TTAGGC_L002_R1.fastq.gz"],
  "Parclip_04" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/PARCLIP/Vickers_Parclip_4_TGACCA_L002_R1.fastq.gz"],
  "Parclip_05" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/PARCLIP/Vickers_Parclip_5_ACAGTG_L002_R1.fastq.gz"],
  "Parclip_06" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/PARCLIP/Vickers_Parclip_6_GCCAAT_L002_R1.fastq.gz"],
  "Parclip_07" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/PARCLIP/Vickers_Parclip_7_CAGATC_L002_R1.fastq.gz"],
  "Parclip_08" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/PARCLIP/Vickers_Parclip_8_ACTTGA_L002_R1.fastq.gz"],
};

my $kcv2797 = {
  "2797-KCV-1" => ["/autofs/blue_sequencer/Runs/projects/2797-KCV/2014-02-06/2797-KCV-1_1.fastq.gz"],
};

my $target_parclip_dir = create_directory_or_die( $target_dir . "/parclip_test" );

my $parclip_config = {
  general      => { "task_name" => "parclip", },
  fastqfiles   => $parclip_files,
  cutadapt => {
    class      => "Cutadapt",
    perform    => 1,
    target_dir => "${target_parclip_dir}/cutadapt",
    option     => "-O 10 -e 0.083",
    source_ref => "fastqfiles",
    adaptor    => "TGGAATTCTCGGGTGCCAAGG",
    extension  => "_clipped.fastq",
    sh_direct  => 0,
    gzipped    => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "24",
      "mem"      => "20gb"
    },
  },
  fastqlen => {
    class      => "FastqLen",
    perform    => 1,
    target_dir => "${target_parclip_dir}/fastqlen",
    option     => "",
    source_ref => "cutadapt",
    cqstools   => $cqstools,
    sh_direct  => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "24",
      "mem"      => "20gb"
    },
  },
  bowtie1out => {
    class         => "Bowtie1",
    perform       => 1,
    target_dir    => "${target_parclip_dir}/bowtie1out",
    option        => "-v 2 -m 10 --best --strata",
    source_ref    => [ "cutadapt", "fastq.gz\$" ],
    bowtie1_index => $bowtie1_human_index,
    samformat     => 0,
    samonly       => 0,
    sh_direct     => 0,
    pbs           => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "20gb"
    },
  },
  bowtie1bam => {
    class         => "Bowtie1",
    perform       => 1,
    target_dir    => "${target_parclip_dir}/bowtie1bam",
    option        => "-v 2 -m 10 --best --strata",
    source_ref    => [ "cutadapt", "fastq.gz\$" ],
    bowtie1_index => $bowtie1_human_index,
    samformat     => 1,
    samonly       => 0,
    sh_direct     => 0,
    pbs           => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "20gb"
    },
  },
  PARalyzer => {
    class      => "ParClip::PARalyzer",
    perform    => 1,
    target_dir => "${target_parclip_dir}/paralyzer",
    option     => "",
    source_ref => "bowtie1out",
    genome2bit => $human_2bit,
    mirna_db   => "/data/cqs/shengq1/reference/miRBase20/hsa.mature.dna.db",
    sh_direct  => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "20gb"
    },
  },
  annotation => {
    class            => "CQS::ParalyzerClusterAnnotator",
    perform          => 1,
    target_dir       => "${target_parclip_dir}/paralyzer",
    option           => "-f miRNA",
    source_ref       => [ "PARalyzer", ".cluster.csv" ],
    cqstools         => $cqstools,
    coordinate_files => [ $hg19_mrna_gff, $hg19_trna_bed ],
    sh_direct        => 1,
    pbs              => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "20gb"
    },
  },
  overall => {
    class      => "CQS::SequenceTask",
    perform    => 1,
    target_dir => "${target_parclip_dir}/overall",
    option     => "",
    source     => { 
      individual => [ "cutadapt", "fastqlen", "bowtie1out", "PARalyzer", "bowtie1bam" ], 
      summary => [ "annotation" ], 
    },
    sh_direct  => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  }
};

performConfig($parclip_config);

1;
