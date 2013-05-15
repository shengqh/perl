#!/usr/bin/perl
use strict;
use warnings;

use CQS::QC;
use CQS::DNASeq;
use CQS::RNASeq;
use CQS::FileUtils;
use CQS::SystemUtils;

my $root = "/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna";
my $target_dir             = create_directory_or_die($root);
my $email                  = "quanhu.sheng\@vanderbilt.edu";
my $task_name              = "VANGARD00055";
my $bwa_option             = "-q 15 -l 8 -n 2";
my $bwa_option_wholegenome = $bwa_option . " -t 8";
my $option_samse_mirna     = "-n 100";

my $bowtie2_option             = "-N 0";
my $bowtie2_option_wholegenome = $bowtie2_option . " -t 8";

my $novoalign_option = "-l 15 -t 30 -r Random -m";

my $config_rat = {
  general => {
    path_file => "/home/shengq1/local/bin/path.txt",
    task_name => $task_name . "_rat"
  },
  fastqfiles => {
    "2516-01" => ["${root}/cutadapt/2516-KCV-1_1_clipped.fastq"],
    "2516-02" => ["${root}/cutadapt/2516-KCV-2_1_clipped.fastq"],
    "2516-03" => ["${root}/cutadapt/2516-KCV-3_1_clipped.fastq"],
    "2516-04" => ["${root}/cutadapt/2516-KCV-4_1_clipped.fastq"],
    "2516-05" => ["${root}/cutadapt/2516-KCV-5_1_clipped.fastq"],
    "2516-06" => ["${root}/cutadapt/2516-KCV-6_1_clipped.fastq"],
    "2516-07" => ["${root}/cutadapt/2516-KCV-7_1_clipped.fastq"],
    "2516-08" => ["${root}/cutadapt/2516-KCV-8_1_clipped.fastq"],
    "2516-09" => ["${root}/cutadapt/2516-KCV-9_1_clipped.fastq"],
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
    target_dir   => "${target_dir}/bwa_miRBase_species",
    option       => $bwa_option,
    option_samse => $option_samse_mirna,
    source_ref   => "fastqfiles",
    fasta_file   => "/data/cqs/shengq1/reference/miRBase19/rno.mature.dna.fa",
    pbs          => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "24",
      "mem"      => "20gb"
    },
  },
  bowtie2 => {
    target_dir    => "${target_dir}/bowtie2_genome",
    option        => $bowtie2_option_wholegenome,
    source_ref    => "fastqfiles",
    bowtie2_index => "/data/cqs/shengq1/reference/rn4/rn4",
    fasta_file    => "/data/cqs/shengq1/reference/rn4/rn4.fa",
    pbs           => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
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
    "2516-10" => ["${root}/cutadapt/2516-KCV-10_1_clipped.fastq"],
    "2516-11" => ["${root}/cutadapt/2516-KCV-11_1_clipped.fastq"],
    "2516-12" => ["${root}/cutadapt/2516-KCV-12_1_clipped.fastq"],
    "2516-13" => ["${root}/cutadapt/2516-KCV-13_1_clipped.fastq"],
    "2516-14" => ["${root}/cutadapt/2516-KCV-14_1_clipped.fastq"],
    "2516-15" => ["${root}/cutadapt/2516-KCV-15_1_clipped.fastq"],
    "2516-16" => ["${root}/cutadapt/2516-KCV-16_1_clipped.fastq"],
    "2516-17" => ["${root}/cutadapt/2516-KCV-17_1_clipped.fastq"],
    "2516-18" => ["${root}/cutadapt/2516-KCV-18_1_clipped.fastq"],
    "2516-19" => ["${root}/cutadapt/2516-KCV-19_1_clipped.fastq"],
    "2516-20" => ["${root}/cutadapt/2516-KCV-20_1_clipped.fastq"],
    "2516-21" => ["${root}/cutadapt/2516-KCV-21_1_clipped.fastq"],
    "2516-22" => ["${root}/cutadapt/2516-KCV-22_1_clipped.fastq"],
    "2516-23" => ["${root}/cutadapt/2516-KCV-23_1_clipped.fastq"],
    "2516-24" => ["${root}/cutadapt/2516-KCV-24_1_clipped.fastq"],
    "2516-25" => ["${root}/cutadapt/2516-KCV-25_1_clipped.fastq"],
    "2516-26" => ["${root}/cutadapt/2516-KCV-26_1_clipped.fastq"],
    "2516-27" => ["${root}/cutadapt/2516-KCV-27_1_clipped.fastq"],
    "2516-28" => ["${root}/cutadapt/2516-KCV-28_1_clipped.fastq"],
    "2516-29" => ["${root}/cutadapt/2516-KCV-29_1_clipped.fastq"],
    "2516-30" => ["${root}/cutadapt/2516-KCV-30_1_clipped.fastq"],
    "2516-31" => ["${root}/cutadapt/2516-KCV-31_1_clipped.fastq"],
    "2516-32" => ["${root}/cutadapt/2516-KCV-32_1_clipped.fastq"],
    "2516-33" => ["${root}/cutadapt/2516-KCV-33_1_clipped.fastq"],
    "2516-34" => ["${root}/cutadapt/2516-KCV-34_1_clipped.fastq"],
    "2516-35" => ["${root}/cutadapt/2516-KCV-35_1_clipped.fastq"],
    "2516-36" => ["${root}/cutadapt/2516-KCV-36_1_clipped.fastq"],
    "2516-37" => ["${root}/cutadapt/2516-KCV-37_1_clipped.fastq"],
    "2516-38" => ["${root}/cutadapt/2516-KCV-38_1_clipped.fastq"],
    "2516-39" => ["${root}/cutadapt/2516-KCV-39_1_clipped.fastq"],
    "2516-40" => ["${root}/cutadapt/2516-KCV-40_1_clipped.fastq"],
    "2516-41" => ["${root}/cutadapt/2516-KCV-41_1_clipped.fastq"],
    "2516-42" => ["${root}/cutadapt/2516-KCV-42_1_clipped.fastq"],
    "2516-43" => ["${root}/cutadapt/2516-KCV-43_1_clipped.fastq"],
    "2516-44" => ["${root}/cutadapt/2516-KCV-44_1_clipped.fastq"],
    "2516-45" => ["${root}/cutadapt/2516-KCV-45_1_clipped.fastq"],
    "2516-46" => ["${root}/cutadapt/2516-KCV-46_1_clipped.fastq"],
    "2516-47" => ["${root}/cutadapt/2516-KCV-47_1_clipped.fastq"],
    "2516-48" => ["${root}/cutadapt/2516-KCV-48_1_clipped.fastq"],
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
    target_dir   => "${target_dir}/bwa_miRBase_species",
    option       => $bwa_option,
    option_samse => $option_samse_mirna,
    source_ref   => "fastqfiles",
    fasta_file   => "/data/cqs/shengq1/reference/miRBase19/hsa.mature.dna.fa",
    pbs          => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "24",
      "mem"      => "20gb"
    },
  },
  bowtie2 => {
    target_dir    => "${target_dir}/bowtie2_genome",
    option        => $bowtie2_option_wholegenome,
    source_ref    => "fastqfiles",
    bowtie2_index => "/data/cqs/guoy1/reference/hg19/bowtie2_index/hg19",
    fasta_file    => "/data/cqs/shengq1/reference/hg19/hg19_chr.fa",
    pbs           => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
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
    "2516-01" => ["${root}/cutadapt/2516-KCV-1_1_clipped.fastq"],
    "2516-02" => ["${root}/cutadapt/2516-KCV-2_1_clipped.fastq"],
    "2516-03" => ["${root}/cutadapt/2516-KCV-3_1_clipped.fastq"],
    "2516-04" => ["${root}/cutadapt/2516-KCV-4_1_clipped.fastq"],
    "2516-05" => ["${root}/cutadapt/2516-KCV-5_1_clipped.fastq"],
    "2516-06" => ["${root}/cutadapt/2516-KCV-6_1_clipped.fastq"],
    "2516-07" => ["${root}/cutadapt/2516-KCV-7_1_clipped.fastq"],
    "2516-08" => ["${root}/cutadapt/2516-KCV-8_1_clipped.fastq"],
    "2516-09" => ["${root}/cutadapt/2516-KCV-9_1_clipped.fastq"],
    "2516-10" => ["${root}/cutadapt/2516-KCV-10_1_clipped.fastq"],
    "2516-11" => ["${root}/cutadapt/2516-KCV-11_1_clipped.fastq"],
    "2516-12" => ["${root}/cutadapt/2516-KCV-12_1_clipped.fastq"],
    "2516-13" => ["${root}/cutadapt/2516-KCV-13_1_clipped.fastq"],
    "2516-14" => ["${root}/cutadapt/2516-KCV-14_1_clipped.fastq"],
    "2516-15" => ["${root}/cutadapt/2516-KCV-15_1_clipped.fastq"],
    "2516-16" => ["${root}/cutadapt/2516-KCV-16_1_clipped.fastq"],
    "2516-17" => ["${root}/cutadapt/2516-KCV-17_1_clipped.fastq"],
    "2516-18" => ["${root}/cutadapt/2516-KCV-18_1_clipped.fastq"],
    "2516-19" => ["${root}/cutadapt/2516-KCV-19_1_clipped.fastq"],
    "2516-20" => ["${root}/cutadapt/2516-KCV-20_1_clipped.fastq"],
    "2516-21" => ["${root}/cutadapt/2516-KCV-21_1_clipped.fastq"],
    "2516-22" => ["${root}/cutadapt/2516-KCV-22_1_clipped.fastq"],
    "2516-23" => ["${root}/cutadapt/2516-KCV-23_1_clipped.fastq"],
    "2516-24" => ["${root}/cutadapt/2516-KCV-24_1_clipped.fastq"],
    "2516-25" => ["${root}/cutadapt/2516-KCV-25_1_clipped.fastq"],
    "2516-26" => ["${root}/cutadapt/2516-KCV-26_1_clipped.fastq"],
    "2516-27" => ["${root}/cutadapt/2516-KCV-27_1_clipped.fastq"],
    "2516-28" => ["${root}/cutadapt/2516-KCV-28_1_clipped.fastq"],
    "2516-29" => ["${root}/cutadapt/2516-KCV-29_1_clipped.fastq"],
    "2516-30" => ["${root}/cutadapt/2516-KCV-30_1_clipped.fastq"],
    "2516-31" => ["${root}/cutadapt/2516-KCV-31_1_clipped.fastq"],
    "2516-32" => ["${root}/cutadapt/2516-KCV-32_1_clipped.fastq"],
    "2516-33" => ["${root}/cutadapt/2516-KCV-33_1_clipped.fastq"],
    "2516-34" => ["${root}/cutadapt/2516-KCV-34_1_clipped.fastq"],
    "2516-35" => ["${root}/cutadapt/2516-KCV-35_1_clipped.fastq"],
    "2516-36" => ["${root}/cutadapt/2516-KCV-36_1_clipped.fastq"],
    "2516-37" => ["${root}/cutadapt/2516-KCV-37_1_clipped.fastq"],
    "2516-38" => ["${root}/cutadapt/2516-KCV-38_1_clipped.fastq"],
    "2516-39" => ["${root}/cutadapt/2516-KCV-39_1_clipped.fastq"],
    "2516-40" => ["${root}/cutadapt/2516-KCV-40_1_clipped.fastq"],
    "2516-41" => ["${root}/cutadapt/2516-KCV-41_1_clipped.fastq"],
    "2516-42" => ["${root}/cutadapt/2516-KCV-42_1_clipped.fastq"],
    "2516-43" => ["${root}/cutadapt/2516-KCV-43_1_clipped.fastq"],
    "2516-44" => ["${root}/cutadapt/2516-KCV-44_1_clipped.fastq"],
    "2516-45" => ["${root}/cutadapt/2516-KCV-45_1_clipped.fastq"],
    "2516-46" => ["${root}/cutadapt/2516-KCV-46_1_clipped.fastq"],
    "2516-47" => ["${root}/cutadapt/2516-KCV-47_1_clipped.fastq"],
    "2516-48" => ["${root}/cutadapt/2516-KCV-48_1_clipped.fastq"],
  },
  bwa_mature => {
    target_dir   => "${target_dir}/bwa_miRBase_mature",
    option       => $bwa_option,
    option_samse => $option_samse_mirna,
    source_ref   => "fastqfiles",
    fasta_file   => "/data/cqs/shengq1/reference/miRBase19/mature.dna.fa",
    pbs          => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "24",
      "mem"      => "20gb"
    },
  },
  bwa_hairpin => {
    target_dir   => "${target_dir}/bwa_miRBase_hairpin",
    option       => $bwa_option,
    option_samse => $option_samse_mirna,
    source_ref   => "fastqfiles",
    fasta_file   => "/data/cqs/shengq1/reference/miRBase19/hairpin.dna.fa",
    pbs          => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "24",
      "mem"      => "20gb"
    },
  },
  bwa_illumina => {
    target_dir   => "${target_dir}/bwa_illumina_miRNA",
    option       => $bwa_option,
    option_samse => $option_samse_mirna,
    source_ref   => "fastqfiles",
    fasta_file   => "/data/cqs/shengq1/reference/miRNA_illumina/mir.fa",
    pbs          => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "24",
      "mem"      => "20gb"
    },
  },
  novoalign_mature => {
    target_dir   => "${target_dir}/novoalign_miRBase_mature",
    option       => $novoalign_option,
    option_samse => "",
    source_ref   => "fastqfiles",
    novoindex    => "/data/cqs/shengq1/reference/miRBase19/mature.dna.nix",
    pbs          => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "24",
      "mem"      => "20gb"
    },
  },
  novoalign_hairpin => {
    target_dir   => "${target_dir}/novoalign_miRBase_hairpin",
    option       => $novoalign_option,
    option_samse => "",
    source_ref   => "fastqfiles",
    novoindex    => "/data/cqs/shengq1/reference/miRBase19/hairpin.dna.nix",
    pbs          => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "24",
      "mem"      => "20gb"
    },
  },
  novoalign_illumina => {
    target_dir   => "${target_dir}/novoalign_illumina",
    option       => $novoalign_option,
    option_samse => "",
    source_ref   => "fastqfiles",
    novoindex    => "/data/cqs/shengq1/reference/miRNA_illumina/mir.nix",
    pbs          => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "24",
      "mem"      => "20gb"
    },
  },

};

#bwa_by_pbs_single( $config_rat, "bwa" );
#bwa_by_pbs_single( $config_human, "bwa" );
#bwa_by_pbs_single( $config_mirna, "bwa_mature" );
#bwa_by_pbs_single( $config_mirna, "bwa_hairpin" );
#bwa_by_pbs_single( $config_mirna, "bwa_illumina" );

#novoalign( $config_mirna, "novoalign_mature" );
#novoalign( $config_mirna, "novoalign_hairpin" );
#novoalign( $config_mirna, "novoalign_illumina" );

#bwa_by_pbs_single( $config_rat,   "bwa_mature" );
#bwa_by_pbs_single( $config_human, "bwa_mature" );

bowtie2($config_rat, "bowtie2");
bowtie2($config_human, "bowtie2");

1;
