#!/usr/bin/perl
use strict;
use warnings;

use CQS::QC;
use CQS::DNASeq;
use CQS::RNASeq;
use CQS::FileUtils;
use CQS::SystemUtils;

my $target_dir = create_directory_or_die("/scratch/cqs/shengq1/miRNA/VANGARD00055");
my $email      = "quanhu.sheng\@vanderbilt.edu";
my $task_name  = "VANGARD00055";
my $bwa_option = "-q 15 -l 8 -n 2";
my $bwa_option_wholegenome = $bwa_option . " -t 8";

my $option_samse_mirna = "-n 100";

my $novoalign_option = "-l 15 -t 30 -r Random -m";

my $config_rat = {
  general => {
    path_file => "/home/shengq1/local/bin/path.txt",
    task_name => $task_name . "_rat"
  },
  fastqfiles => {
    "2516-01" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-1_1_clipped.fastq"],
    "2516-02" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-2_1_clipped.fastq"],
    "2516-03" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-3_1_clipped.fastq"],
    "2516-04" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-4_1_clipped.fastq"],
    "2516-05" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-5_1_clipped.fastq"],
    "2516-06" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-6_1_clipped.fastq"],
    "2516-07" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-7_1_clipped.fastq"],
    "2516-08" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-8_1_clipped.fastq"],
    "2516-09" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-9_1_clipped.fastq"],
  },
  bwa => {
    target_dir      => "${target_dir}/bwa_genome",
    option          => $bwa_option_wholegenome,
    option_samse    => "",
    source_ref      => "fastqfiles",
    fasta_file      => "/data/cqs/shengq1/reference/rn5/rn5.fa",
    estimate_insert => 1,
    source_ref      => "fastqfiles",
    pbs             => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "24",
      "mem"      => "20gb"
    },
  },
  bwa_mature => {
    target_dir      => "${target_dir}/bwa_miRBase_species",
    option          => $bwa_option,
    option_samse    => $option_samse_mirna,
    source_ref      => "fastqfiles",
    fasta_file      => "/data/cqs/shengq1/reference/miRBase19/rno.mature.dna.fa",
    pbs             => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "24",
      "mem"      => "20gb"
    },
  },
};

my $config_human = {
  general => {
    path_file => "/home/shengq1/local/bin/path.txt",
    task_name => $task_name . "_human"
  },
  fastqfiles => {
    "2516-10" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-10_1_clipped.fastq"],
    "2516-11" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-11_1_clipped.fastq"],
    "2516-12" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-12_1_clipped.fastq"],
    "2516-13" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-13_1_clipped.fastq"],
    "2516-14" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-14_1_clipped.fastq"],
    "2516-15" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-15_1_clipped.fastq"],
    "2516-16" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-16_1_clipped.fastq"],
    "2516-17" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-17_1_clipped.fastq"],
    "2516-18" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-18_1_clipped.fastq"],
    "2516-19" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-19_1_clipped.fastq"],
    "2516-20" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-20_1_clipped.fastq"],
    "2516-21" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-21_1_clipped.fastq"],
    "2516-22" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-22_1_clipped.fastq"],
    "2516-23" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-23_1_clipped.fastq"],
    "2516-24" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-24_1_clipped.fastq"],
    "2516-25" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-25_1_clipped.fastq"],
    "2516-26" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-26_1_clipped.fastq"],
    "2516-27" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-27_1_clipped.fastq"],
    "2516-28" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-28_1_clipped.fastq"],
    "2516-29" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-29_1_clipped.fastq"],
    "2516-30" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-30_1_clipped.fastq"],
    "2516-31" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-31_1_clipped.fastq"],
    "2516-32" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-32_1_clipped.fastq"],
    "2516-33" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-33_1_clipped.fastq"],
    "2516-34" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-34_1_clipped.fastq"],
    "2516-35" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-35_1_clipped.fastq"],
    "2516-36" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-36_1_clipped.fastq"],
    "2516-37" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-37_1_clipped.fastq"],
    "2516-38" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-38_1_clipped.fastq"],
    "2516-39" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-39_1_clipped.fastq"],
    "2516-40" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-40_1_clipped.fastq"],
    "2516-41" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-41_1_clipped.fastq"],
    "2516-42" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-42_1_clipped.fastq"],
    "2516-43" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-43_1_clipped.fastq"],
    "2516-44" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-44_1_clipped.fastq"],
    "2516-45" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-45_1_clipped.fastq"],
    "2516-46" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-46_1_clipped.fastq"],
    "2516-47" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-47_1_clipped.fastq"],
    "2516-48" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-48_1_clipped.fastq"],
  },
  bwa => {
    target_dir      => "${target_dir}/bwa_genome",
    option          => $bwa_option_wholegenome,
    option_samse    => "",
    source_ref      => "fastqfiles",
    fasta_file      => "/data/cqs/shengq1/reference/hg19/hg19_chr.fa",
    estimate_insert => 1,
    source_ref      => "fastqfiles",
    pbs             => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "24",
      "mem"      => "20gb"
    },
  },
  bwa_mature => {
    target_dir      => "${target_dir}/bwa_miRBase_species",
    option          => $bwa_option,
    option_samse    => $option_samse_mirna,
    source_ref      => "fastqfiles",
    fasta_file      => "/data/cqs/shengq1/reference/miRBase19/hsa.mature.dna.fa",
    pbs             => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "24",
      "mem"      => "20gb"
    },
  },
};

my $config_mirna = {
  general => {
    path_file => "/home/shengq1/local/bin/path.txt",
    task_name => $task_name . "_mirna"
  },
  fastqfiles => {
    "2516-01" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-1_1_clipped.fastq"],
    "2516-02" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-2_1_clipped.fastq"],
    "2516-03" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-3_1_clipped.fastq"],
    "2516-04" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-4_1_clipped.fastq"],
    "2516-05" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-5_1_clipped.fastq"],
    "2516-06" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-6_1_clipped.fastq"],
    "2516-07" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-7_1_clipped.fastq"],
    "2516-08" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-8_1_clipped.fastq"],
    "2516-09" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-9_1_clipped.fastq"],
    "2516-10" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-10_1_clipped.fastq"],
    "2516-11" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-11_1_clipped.fastq"],
    "2516-12" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-12_1_clipped.fastq"],
    "2516-13" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-13_1_clipped.fastq"],
    "2516-14" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-14_1_clipped.fastq"],
    "2516-15" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-15_1_clipped.fastq"],
    "2516-16" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-16_1_clipped.fastq"],
    "2516-17" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-17_1_clipped.fastq"],
    "2516-18" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-18_1_clipped.fastq"],
    "2516-19" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-19_1_clipped.fastq"],
    "2516-20" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-20_1_clipped.fastq"],
    "2516-21" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-21_1_clipped.fastq"],
    "2516-22" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-22_1_clipped.fastq"],
    "2516-23" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-23_1_clipped.fastq"],
    "2516-24" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-24_1_clipped.fastq"],
    "2516-25" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-25_1_clipped.fastq"],
    "2516-26" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-26_1_clipped.fastq"],
    "2516-27" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-27_1_clipped.fastq"],
    "2516-28" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-28_1_clipped.fastq"],
    "2516-29" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-29_1_clipped.fastq"],
    "2516-30" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-30_1_clipped.fastq"],
    "2516-31" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-31_1_clipped.fastq"],
    "2516-32" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-32_1_clipped.fastq"],
    "2516-33" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-33_1_clipped.fastq"],
    "2516-34" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-34_1_clipped.fastq"],
    "2516-35" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-35_1_clipped.fastq"],
    "2516-36" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-36_1_clipped.fastq"],
    "2516-37" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-37_1_clipped.fastq"],
    "2516-38" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-38_1_clipped.fastq"],
    "2516-39" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-39_1_clipped.fastq"],
    "2516-40" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-40_1_clipped.fastq"],
    "2516-41" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-41_1_clipped.fastq"],
    "2516-42" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-42_1_clipped.fastq"],
    "2516-43" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-43_1_clipped.fastq"],
    "2516-44" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-44_1_clipped.fastq"],
    "2516-45" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-45_1_clipped.fastq"],
    "2516-46" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-46_1_clipped.fastq"],
    "2516-47" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-47_1_clipped.fastq"],
    "2516-48" => ["/scratch/cqs/shengq1/miRNA/VANGARD00055/clipped/2516-KCV-48_1_clipped.fastq"],
  },
  bwa_mature => {
    target_dir      => "${target_dir}/bwa_miRBase_mature",
    option          => $bwa_option,
    option_samse    => $option_samse_mirna,
    source_ref      => "fastqfiles",
    fasta_file      => "/data/cqs/shengq1/reference/miRBase19/mature.dna.fa",
    pbs             => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "24",
      "mem"      => "20gb"
    },
  },
  bwa_hairpin => {
    target_dir      => "${target_dir}/bwa_miRBase_hairpin",
    option          => $bwa_option,
    option_samse    => $option_samse_mirna,
    source_ref      => "fastqfiles",
    fasta_file      => "/data/cqs/shengq1/reference/miRBase19/hairpin.dna.fa",
    pbs             => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "24",
      "mem"      => "20gb"
    },
  },
  bwa_illumina => {
    target_dir      => "${target_dir}/bwa_illumina_miRNA",
    option          => $bwa_option,
    option_samse    => $option_samse_mirna,
    source_ref      => "fastqfiles",
    fasta_file      => "/data/cqs/shengq1/reference/miRNA_illumina/mir.fa",
    pbs             => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "24",
      "mem"      => "20gb"
    },
  },
  novoalign_mature => {
    target_dir      => "${target_dir}/novoalign_miRBase_mature",
    option          => $novoalign_option,
    option_samse    => "",
    source_ref      => "fastqfiles",
    novoindex      => "/data/cqs/shengq1/reference/miRBase19/mature.dna.nix",
    pbs             => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "24",
      "mem"      => "20gb"
    },
  },
  novoalign_hairpin => {
    target_dir      => "${target_dir}/novoalign_miRBase_hairpin",
    option          => $novoalign_option,
    option_samse    => "",
    source_ref      => "fastqfiles",
    novoindex      => "/data/cqs/shengq1/reference/miRBase19/hairpin.dna.nix",
    pbs             => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "24",
      "mem"      => "20gb"
    },
  },
  novoalign_illumina => {
    target_dir      => "${target_dir}/novoalign_illumina",
    option          => $novoalign_option,
    option_samse    => "",
    source_ref      => "fastqfiles",
    novoindex      => "/data/cqs/shengq1/reference/miRNA_illumina/mir.nix",
    pbs             => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "24",
      "mem"      => "20gb"
    },
  },
  
};

#bwa_by_pbs_single( $config_rat, "bwa" );
#bwa_by_pbs_single( $config_human, "bwa" );
bwa_by_pbs_single( $config_mirna, "bwa_mature" );
bwa_by_pbs_single( $config_mirna, "bwa_hairpin" );
bwa_by_pbs_single( $config_mirna, "bwa_illumina" );
#novoalign( $config_mirna, "novoalign_mature" );
#novoalign( $config_mirna, "novoalign_hairpin" );
#novoalign( $config_mirna, "novoalign_illumina" );

bwa_by_pbs_single( $config_rat, "bwa_mature" );
bwa_by_pbs_single( $config_human, "bwa_mature" );


1;
