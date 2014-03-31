#!/usr/bin/perl
use strict;
use warnings;

use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;

my $root     = create_directory_or_die("/scratch/cqs/shengq1/vangard/VANGARD_Vicky/201403_parclip_2797/");
my $cqstools = "/home/shengq1/cqstools/CQS.Tools.exe";

my $email = "quanhu.sheng\@vanderbilt.edu";

my $demultiplexing_config = {
  general        => { "task_name" => "parclip", },
  demultiplexing => {
    class      => "Format::Demultiplexing",
    perform    => 1,
    target_dir => "${root}/demultiplexing",
    option     => "",
    source     => { "2797-KCV-1" => ["/autofs/blue_sequencer/Runs/projects/2797-KCV/2014-02-06/2797-KCV-1_1.fastq.gz"], },
    maps       => { "2797-KCV-1" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_201403_Vicky_parclip/2797-KCV/demultiplexing/pbs/2797-KCV-1.map"], },
    cqstools   => $cqstools,
    sh_direct  => 0,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "24",
      "mem"      => "20gb"
    },
  },
};

performConfig($demultiplexing_config);

my $kcv2797human = {
  files            => { "2797-KCV-1" => ["/autofs/blue_sequencer/Runs/projects/2797-KCV/2014-02-06/2797-KCV-1_1.fastq.gz"], },
  maps             => { "2797-KCV-1" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_201403_Vicky_parclip/2797-KCV/demultiplexing/pbs/2797-KCV-1.map"], },
  task_name        => "2797-KCV",
  mirna_coordinate => "/data/cqs/shengq1/reference/miRBase20/hsa.gff3",
  trna_coordinate => "/data/cqs/guoy1/reference/smallrna/hg19_tRNA_ucsc_ensembl.bed",
  bowtie1_index   => "/data/cqs/guoy1/reference/hg19/bowtie_index_hg19_rCRS_1.0.0/hg19_rCRS",
  genome2bit      => "/data/cqs/guoy1/reference/hg19/hg19_rCRS.2bit",
  mirna_db        => "/data/cqs/shengq1/reference/miRBase20/hsa.mature.dna.db",
  target_dir      => $root,
};

my @datasets = ($kcv2797human);

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
      cqstools   => $cqstools,
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
    bowtie1bam => {
      class         => "Bowtie1",
      perform       => 1,
      target_dir    => "${target_parclip_dir}/bowtie1bam",
      option        => "-v 2 -m 10 --best --strata",
      source_ref    => [ "cutadapt", "fastq.gz\$" ],
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
    PARalyzer => {
      class      => "ParClip::PARalyzer",
      perform    => 1,
      target_dir => "${target_parclip_dir}/paralyzer",
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
      target_dir       => "${target_parclip_dir}/paralyzer",
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
    sequencetask => {
      class      => "CQS::SequenceTask",
      perform    => 1,
      target_dir => "${target_parclip_dir}/sequencetask",
      option     => "",
      source     => {
        T0_prepare    => ["demultiplexing"],
        T1_individual => [ "cutadapt", "fastqlen", "bowtie1out", "PARalyzer", "bowtie1bam" ],
        T2_summary    => ["annotation"],
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

  #performConfig($parclip_config);
#  performTask( $parclip_config, "demultiplexing" );
}

1;
