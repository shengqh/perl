#!/usr/bin/perl
use strict;
use warnings;

use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;

my $target_dir  = "/scratch/cqs/shengq1/vangard/VANGARD_Vicky/201312_parclip";
my $cqstools    = "/home/shengq1/cqstools/CQS.Tools.exe";
my $samtools    = "/home/shengq1/local/bin/samtools/samtools";
my $mirna_fasta = "/data/cqs/shengq1/reference/miRBase20/mature.dna.fa";
my $email       = "quanhu.sheng\@vanderbilt.edu";

my $mirnacount_option = "-s";    #ignore score

my $dataset = {
  files => {
    "Parclip_01" => ["/scratch/cqs/shengq1/vangard/VANGARD_Vicky/201312_parclip/raw/Vickers_Parclip_1_ATCACG_L002_R1.fastq.gz"],
    "Parclip_02" => ["/scratch/cqs/shengq1/vangard/VANGARD_Vicky/201312_parclip/raw/Vickers_Parclip_2_CGATGT_L002_R1.fastq.gz"],
    "Parclip_03" => ["/scratch/cqs/shengq1/vangard/VANGARD_Vicky/201312_parclip/raw/Vickers_Parclip_3_TTAGGC_L002_R1.fastq.gz"],
    "Parclip_04" => ["/scratch/cqs/shengq1/vangard/VANGARD_Vicky/201312_parclip/raw/Vickers_Parclip_4_TGACCA_L002_R1.fastq.gz"],
    "Parclip_05" => ["/scratch/cqs/shengq1/vangard/VANGARD_Vicky/201312_parclip/raw/Vickers_Parclip_5_ACAGTG_L002_R1.fastq.gz"],
    "Parclip_06" => ["/scratch/cqs/shengq1/vangard/VANGARD_Vicky/201312_parclip/raw/Vickers_Parclip_6_GCCAAT_L002_R1.fastq.gz"],
    "Parclip_07" => ["/scratch/cqs/shengq1/vangard/VANGARD_Vicky/201312_parclip/raw/Vickers_Parclip_7_CAGATC_L002_R1.fastq.gz"],
    "Parclip_08" => ["/scratch/cqs/shengq1/vangard/VANGARD_Vicky/201312_parclip/raw/Vickers_Parclip_8_ACTTGA_L002_R1.fastq.gz"],
  },
  task_name        => "parclip_human",
  mirna_coordinate => "/data/cqs/shengq1/reference/miRBase20/hsa.gff3",
  trna_coordinate  => "/data/cqs/guoy1/reference/smallrna/hg19_tRNA_ucsc_ensembl.bed",
  bowtie1_index    => "/data/cqs/guoy1/reference/hg19/bowtie_index_hg19_rCRS_1.0.0/hg19_rCRS",
  genome_2bit      => "/data/cqs/guoy1/reference/hg19/hg19_rCRS.2bit",
  mirna_db         => "/data/cqs/shengq1/reference/miRBase20/hsa.mature.dna.db",
};

my $parclip_config = {
  general  => { "task_name" => "parclip", },
  cutadapt => {
    class      => "Cutadapt",
    perform    => 1,
    target_dir => "${target_dir}/cutadapt",
    option     => "-m 12 -O 10 -e 0.083",
    source_ref => $dataset->{files},
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
    target_dir => "${target_dir}/fastqlen",
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
    target_dir    => "${target_dir}/bowtie1out",
    option        => "-v 2 -m 10 --best --strata",
    source        => $dataset->{files},
    bowtie1_index => $dataset->{bowtie1_index},
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
  PARalyzer => {
    class      => "ParClip::PARalyzer",
    perform    => 1,
    target_dir => "${target_dir}/paralyzer",
    option     => "",
    source_ref => "bowtie1out",
    genome2bit => $dataset->{genome_2bit},
    mirna_db   => $dataset->{mirna_db},
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
    target_dir       => "${target_dir}/paralyzer",
    option           => "-f miRNA",
    source_ref       => [ "PARalyzer", ".cluster.csv" ],
    cqstools         => $cqstools,
    coordinate_files => [ $dataset->{mirna_coordinate}, $dataset->{trna_coordinate} ],
    sh_direct        => 1,
    pbs              => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "20gb"
    },
  },
  bowtie1bam => {
    class         => "Bowtie1",
    perform       => 1,
    target_dir    => "${target_dir}/bowtie1bam",
    option        => "-v 2 -m 10 --best --strata",
    source        => $dataset->{files},
    bowtie1_index => $dataset->{bowtie1_index},
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
  mirna_count => {
    class      => "MirnaCount",
    perform    => 1,
    target_dir => "${target_dir}/count_miRNA",
    option     => $mirnacount_option,
    source_ref => "bowtie1bam",
    cqs_tools  => $cqstools,
    gff_file   => $dataset->{mirna_coordinate},
    fasta_file => $mirna_fasta,
    samtools   => $samtools,
    sh_direct  => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  sequencetask => {
    class      => "CQS::SequenceTask",
    perform    => 1,
    target_dir => "${target_dir}/sequencetask",
    option     => "",
    source     => {
      T1_individual => [ "cutadapt", "fastqlen", "bowtie1out", "PARalyzer", "bowtie1bam", "mirna_count" ],
      T2_summary    => ["annotation"],
    },
    sh_direct => 0,
    pbs       => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  }
};

performConfig($parclip_config);

1;
