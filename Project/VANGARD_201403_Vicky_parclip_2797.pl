#!/usr/bin/perl
use strict;
use warnings;

use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;

my $root     = create_directory_or_die("/scratch/cqs/shengq1/vangard/VANGARD_201403_Vicky_parclip/");
my $cqstools = "/home/shengq1/cqstools/CQS.Tools.exe";

my $hg19_mrna_gff = "/data/cqs/shengq1/reference/miRBase20/hsa.gff3";
my $hg19_trna_bed = "/data/cqs/guoy1/reference/smallrna/hg19_tRNA_ucsc_ensembl.bed";

my $target_dir = $root;

my $email     = "quanhu.sheng\@vanderbilt.edu";
my $task_name = "VANGARD00055";

my $bowtie1_human_index = "/data/cqs/guoy1/reference/hg19/bowtie_index_hg19_rCRS_1.0.0/hg19_rCRS";
my $human_2bit          = "/data/cqs/guoy1/reference/hg19/hg19_rCRS.2bit";

my $kcv2797 = {
  files => {
    "2797-KCV-1_RPI40_Ago2INS1Huh7"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_201403_Vicky_parclip/2797-KCV-1/raw/2797-KCV-1_RPI40_Ago2INS1Huh7.fastq.gz"],
    "2797-KCV-1_RPI41_Ago3INS1Huh7"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_201403_Vicky_parclip/2797-KCV-1/raw/2797-KCV-1_RPI41_Ago3INS1Huh7.fastq.gz"],
    "2797-KCV-1_RPI42_Ago2INS1HCEAC" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_201403_Vicky_parclip/2797-KCV-1/raw/2797-KCV-1_RPI42_Ago2INS1HCEAC.fastq.gz"],
    "2797-KCV-1_RPI43_Ago3INS1HCEAC" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_201403_Vicky_parclip/2797-KCV-1/raw/2797-KCV-1_RPI43_Ago3INS1HCEAC.fastq.gz"],
    "2797-KCV-1_RPI47_Ago2MIN6Huh7"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_201403_Vicky_parclip/2797-KCV-1/raw/2797-KCV-1_RPI47_Ago2MIN6Huh7.fastq.gz"],
    "2797-KCV-1_RPI48_Ago3MIN6Huh7"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_201403_Vicky_parclip/2797-KCV-1/raw/2797-KCV-1_RPI48_Ago3MIN6Huh7.fastq.gz"],
  },
  task_name => "2797-KCV"
};

my $kcv2795 = {
  files => {
    "2795-KCV-1" => ["/gpfs21/scratch/vantage_repo/Vickers/2795/2795-KCV-1_1_sequence.txt.gz"],
    "2795-KCV-2" => ["/gpfs21/scratch/vantage_repo/Vickers/2795/2795-KCV-2_1_sequence.txt.gz"],
    "2795-KCV-3" => ["/gpfs21/scratch/vantage_repo/Vickers/2795/2795-KCV-3_1_sequence.txt.gz"],
    "2795-KCV-4" => ["/gpfs21/scratch/vantage_repo/Vickers/2795/2795-KCV-4_1_sequence.txt.gz"],
  },
  maps => {
    "2795-KCV-1" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_201403_Vicky_parclip/2795-KCV/2795-KCV-1.map"],
    "2795-KCV-2" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_201403_Vicky_parclip/2795-KCV/2795-KCV-2.map"],
    "2795-KCV-3" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_201403_Vicky_parclip/2795-KCV/2795-KCV-3.map"],
    "2795-KCV-4" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_201403_Vicky_parclip/2795-KCV/2795-KCV-4.map"],
  },
  task_name => "2795-KCV"
};

my @datasets = (

  #$kcv2797,
  $kcv2795
);

foreach my $dataset (@datasets) {
  my $target_parclip_dir = create_directory_or_die( $root . $dataset->{task_name} );
  my $parclip_config     = {
    general        => { "task_name" => "parclip", },
    demultiplexing => {
      class      => "Format::Demultiplexing",
      perform    => 1,
      target_dir => "${target_parclip_dir}/demultiplexing",
      option     => "",
      source     => $dataset->{files},
      maps       => $dataset->{maps},
      sh_direct  => 0,
      pbs        => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "24",
        "mem"      => "20gb"
      },
    },
    cutadapt => {
      class      => "Cutadapt",
      perform    => 1,
      target_dir => "${target_parclip_dir}/cutadapt",
      option     => "-O 10 -e 0.083",
      source_ref => "demultiplexing",
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
        summary    => ["annotation"],
      },
      sh_direct => 0,
      pbs       => {
        "email"    => $email,
        "nodes"    => "1:ppn=8",
        "walltime" => "72",
        "mem"      => "40gb"
      },
    }
  };

  performConfig($parclip_config);
}

1;
