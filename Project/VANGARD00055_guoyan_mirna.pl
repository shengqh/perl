#!/usr/bin/perl
use strict;
use warnings;

use CQS::QC;
use CQS::DNASeq;
use CQS::RNASeq;
use CQS::CQSTools;
use CQS::FileUtils;
use CQS::SystemUtils;

my $root       = "/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna";
my $target_dir = create_directory_or_die($root);

my $target_rat_dir   = create_directory_or_die( $target_dir . "/rat" );
my $target_human_dir = create_directory_or_die( $target_dir . "/human" );

my $email     = "quanhu.sheng\@vanderbilt.edu";
my $task_name = "VANGARD00055";

my $cqs_tools = "/home/shengq1/cqstools/CQS.Tools.exe";
my $hsa_gffs  = "/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/hsa.gff3";
my $rno_gffs  = "/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/rno.gff3";

my $bwa_option             = "-l 8 -n 1 -o 0";
my $bwa_option_wholegenome = $bwa_option . " -t 8";
my $option_samse_mirna     = "";

my $bowtie2_option             = "-N 0 --phred33";
my $bowtie2_option_wholegenome = $bowtie2_option . " -p 8";

my $bowtie1_option             = "-v 1 --best";
my $bowtie1_option_wholegenome = $bowtie1_option . " -p 8";

my $novoalign_option = "-l 15 -t 30 -r Random -m";

my $shrimp2_option = "-Q -N 8 -n 1 -o 1 --qv-offset 33";

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
  unmappedfiles => {
    "2516-01" => ["${target_rat_dir}/bowtie2_genome/result/2516-01/2516-01.bam.unmapped.fastq"],
    "2516-02" => ["${target_rat_dir}/bowtie2_genome/result/2516-02/2516-02.bam.unmapped.fastq"],
    "2516-03" => ["${target_rat_dir}/bowtie2_genome/result/2516-03/2516-03.bam.unmapped.fastq"],
    "2516-04" => ["${target_rat_dir}/bowtie2_genome/result/2516-04/2516-04.bam.unmapped.fastq"],
    "2516-05" => ["${target_rat_dir}/bowtie2_genome/result/2516-05/2516-05.bam.unmapped.fastq"],
    "2516-06" => ["${target_rat_dir}/bowtie2_genome/result/2516-06/2516-06.bam.unmapped.fastq"],
    "2516-07" => ["${target_rat_dir}/bowtie2_genome/result/2516-07/2516-07.bam.unmapped.fastq"],
    "2516-08" => ["${target_rat_dir}/bowtie2_genome/result/2516-08/2516-08.bam.unmapped.fastq"],
    "2516-09" => ["${target_rat_dir}/bowtie2_genome/result/2516-09/2516-09.bam.unmapped.fastq"],
  },
  bwa_mature => {
    target_dir   => "${target_rat_dir}/bwa_miRBase_species",
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
  bwa => {
    target_dir      => "${target_rat_dir}/bwa_genome",
    option          => $bwa_option_wholegenome,
    option_samse    => "",
    source_ref      => "fastqfiles",
    fasta_file      => "/data/cqs/shengq1/reference/rn4/bwa_index/rn4.fa",
    estimate_insert => 0,
    source_ref      => "fastqfiles",
    pbs             => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "24",
      "mem"      => "20gb"
    },
  },
  bowtie1 => {
    target_dir    => "${target_rat_dir}/bowtie1_genome",
    option        => $bowtie1_option_wholegenome,
    source_ref    => "fastqfiles",
    bowtie1_index => "/data/cqs/shengq1/reference/rn4/bowtie1_index/rn4",
    fasta_file    => "/data/cqs/shengq1/reference/rn4/bowtie1_index/rn4.fa",
    samonly       => 0,
    pbs           => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "24",
      "mem"      => "20gb"
    },
  },
  bowtie2 => {
    target_dir    => "${target_rat_dir}/bowtie2_genome",
    option        => $bowtie2_option_wholegenome,
    source_ref    => "fastqfiles",
    bowtie2_index => "/data/cqs/shengq1/reference/rn4/bowtie2_index/rn4",
    fasta_file    => "/data/cqs/shengq1/reference/rn4/bowtie2_index/rn4.fa",
    samonly       => 0,
    pbs           => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "24",
      "mem"      => "20gb"
    },
  },
  mirna_count_bwa => {
    target_dir   => "${target_rat_dir}/bwa_genome",
    option       => "",
    source_ref   => "bwa",
    cqs_tools    => $cqs_tools,
    gff_file     => $rno_gffs,
    fasta_format => 0,
    pbs          => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "24",
      "mem"      => "40gb"
    },
  },
  mirna_count_bowtie1 => {
    target_dir   => "${target_rat_dir}/bowtie1_genome",
    option       => "",
    source_ref   => "bowtie1",
    cqs_tools    => $cqs_tools,
    gff_file     => $rno_gffs,
    fasta_format => 0,
    pbs          => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "24",
      "mem"      => "40gb"
    },
  },
  mirna_count_bowtie2 => {
    target_dir   => "${target_rat_dir}/bowtie2_genome",
    option       => "",
    source_ref   => "bowtie2",
    cqs_tools    => $cqs_tools,
    gff_file     => $rno_gffs,
    fasta_format => 0,
    pbs          => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "24",
      "mem"      => "40gb"
    },
  },
  shrimp2 => {
    target_dir   => "${target_rat_dir}/shrimp2",
    option       => $shrimp2_option,
    unmapped_ref => "bowtie2",
    genome_index => "/data/cqs/shengq1/reference/rn4/shrimp2_index_ls/rn4-ls",
    is_mirna     => 1,
    output_bam   => 1,
    sh_direct    => 1,
    pbs          => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "24",
      "mem"      => "20gb"
    },
  },
  bowtie2_mature => {
    target_dir    => "${target_rat_dir}/bowtie2_mature",
    option        => $bowtie2_option,
    source_ref    => "unmappedfiles",
    bowtie2_index => "/data/cqs/shengq1/reference/miRBase19/rno.mature.dna",
    fasta_file    => "/data/cqs/shengq1/reference/miRBase19/rno.mature.dna.fa",
    samonly       => 0,
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
    "2516-10"  => ["${root}/cutadapt/2516-KCV-10_1_clipped.fastq"],
    "2516-11"  => ["${root}/cutadapt/2516-KCV-11_1_clipped.fastq"],
    "2516-12"  => ["${root}/cutadapt/2516-KCV-12_1_clipped.fastq"],
    "2516-13"  => ["${root}/cutadapt/2516-KCV-13_1_clipped.fastq"],
    "2516-14"  => ["${root}/cutadapt/2516-KCV-14_1_clipped.fastq"],
    "2516-15"  => ["${root}/cutadapt/2516-KCV-15_1_clipped.fastq"],
    "2516-16"  => ["${root}/cutadapt/2516-KCV-16_1_clipped.fastq"],
    "2516-17"  => ["${root}/cutadapt/2516-KCV-17_1_clipped.fastq"],
    "2516-18"  => ["${root}/cutadapt/2516-KCV-18_1_clipped.fastq"],
    "2516-19"  => ["${root}/cutadapt/2516-KCV-19_1_clipped.fastq"],
    "2516-20"  => ["${root}/cutadapt/2516-KCV-20_1_clipped.fastq"],
    "2516-21"  => ["${root}/cutadapt/2516-KCV-21_1_clipped.fastq"],
    "2516-22"  => ["${root}/cutadapt/2516-KCV-22_1_clipped.fastq"],
    "2516-23"  => ["${root}/cutadapt/2516-KCV-23_1_clipped.fastq"],
    "2516-24"  => ["${root}/cutadapt/2516-KCV-24_1_clipped.fastq"],
    "2516-25"  => ["${root}/cutadapt/2516-KCV-25_1_clipped.fastq"],
    "2516-26"  => ["${root}/cutadapt/2516-KCV-26_1_clipped.fastq"],
    "2516-27"  => ["${root}/cutadapt/2516-KCV-27_1_clipped.fastq"],
    "2516-28"  => ["${root}/cutadapt/2516-KCV-28_1_clipped.fastq"],
    "2516-29"  => ["${root}/cutadapt/2516-KCV-29_1_clipped.fastq"],
    "2516-30"  => ["${root}/cutadapt/2516-KCV-30_1_clipped.fastq"],
    "2516-31"  => ["${root}/cutadapt/2516-KCV-31_1_clipped.fastq"],
    "2516-32"  => ["${root}/cutadapt/2516-KCV-32_1_clipped.fastq"],
    "2516-33"  => ["${root}/cutadapt/2516-KCV-33_1_clipped.fastq"],
    "2516-34"  => ["${root}/cutadapt/2516-KCV-34_1_clipped.fastq"],
    "2516-35"  => ["${root}/cutadapt/2516-KCV-35_1_clipped.fastq"],
    "2516-36"  => ["${root}/cutadapt/2516-KCV-36_1_clipped.fastq"],
    "2516-37"  => ["${root}/cutadapt/2516-KCV-37_1_clipped.fastq"],
    "2516-38"  => ["${root}/cutadapt/2516-KCV-38_1_clipped.fastq"],
    "2516-39"  => ["${root}/cutadapt/2516-KCV-39_1_clipped.fastq"],
    "2516-40"  => ["${root}/cutadapt/2516-KCV-40_1_clipped.fastq"],
    "2516-41"  => ["${root}/cutadapt/2516-KCV-41_1_clipped.fastq"],
    "2516-42"  => ["${root}/cutadapt/2516-KCV-42_1_clipped.fastq"],
    "2516-43"  => ["${root}/cutadapt/2516-KCV-43_1_clipped.fastq"],
    "2516-44"  => ["${root}/cutadapt/2516-KCV-44_1_clipped.fastq"],
    "2516-45"  => ["${root}/cutadapt/2516-KCV-45_1_clipped.fastq"],
    "2516-46"  => ["${root}/cutadapt/2516-KCV-46_1_clipped.fastq"],
    "2516-47"  => ["${root}/cutadapt/2516-KCV-47_1_clipped.fastq"],
    "2516-48"  => ["${root}/cutadapt/2516-KCV-48_1_clipped.fastq"],
    "KCV2_1N2" => ["${root}/cutadapt/Fuchs/KCV2_1N2_GCCAAT_L003_R1_001_clipped.fastq"],
    "KCV2_1N3" => [ "${root}/cutadapt/Fuchs/KCV2_1N3_CAGATC_L003_R1_001_clipped.fastq", "${root}/cutadapt/Fuchs/KCV2_1N3_CAGATC_L003_R1_002_clipped.fastq" ],
    "KCV2_1N4" => [ "${root}/cutadapt/Fuchs/KCV2_1N4_ACTTGA_L003_R1_001_clipped.fastq", "${root}/cutadapt/Fuchs/KCV2_1N4_ACTTGA_L003_R1_002_clipped.fastq" ],
    "KCV2_1N5" => [ "${root}/cutadapt/Fuchs/KCV2_1N5_GATCAG_L003_R1_001_clipped.fastq", "${root}/cutadapt/Fuchs/KCV2_1N5_GATCAG_L003_R1_002_clipped.fastq" ],
    "KCV2_1N6" => [ "${root}/cutadapt/Fuchs/KCV2_1N6_TAGCTT_L003_R1_001_clipped.fastq", "${root}/cutadapt/Fuchs/KCV2_1N6_TAGCTT_L003_R1_002_clipped.fastq" ],
    "KCV2_2N1" => [
      "${root}/cutadapt/Fuchs/KCV2_2N1_GGCTAC_L003_R1_001_clipped.fastq", "${root}/cutadapt/Fuchs/KCV2_2N1_GGCTAC_L003_R1_002_clipped.fastq",
      "${root}/cutadapt/Fuchs/KCV2_2N1_GGCTAC_L003_R1_003_clipped.fastq", "${root}/cutadapt/Fuchs/KCV2_2N1_GGCTAC_L003_R1_004_clipped.fastq",
      "${root}/cutadapt/Fuchs/KCV2_2N1_GGCTAC_L003_R1_005_clipped.fastq"
    ],
    "KCV2_2N2" => ["${root}/cutadapt/Fuchs/KCV2_2N2_GCCGCG_L003_R1_001_clipped.fastq"],
    "KCV2_2N3" => [
      "${root}/cutadapt/Fuchs/KCV2_2N3_CTTGTA_L003_R1_001_clipped.fastq", "${root}/cutadapt/Fuchs/KCV2_2N3_CTTGTA_L003_R1_002_clipped.fastq",
      "${root}/cutadapt/Fuchs/KCV2_2N3_CTTGTA_L003_R1_003_clipped.fastq", "${root}/cutadapt/Fuchs/KCV2_2N3_CTTGTA_L003_R1_004_clipped.fastq",
      "${root}/cutadapt/Fuchs/KCV2_2N3_CTTGTA_L003_R1_005_clipped.fastq", "${root}/cutadapt/Fuchs/KCV2_2N3_CTTGTA_L003_R1_006_clipped.fastq",
      "${root}/cutadapt/Fuchs/KCV2_2N3_CTTGTA_L003_R1_007_clipped.fastq", "${root}/cutadapt/Fuchs/KCV2_2N3_CTTGTA_L003_R1_008_clipped.fastq",
      "${root}/cutadapt/Fuchs/KCV2_2N3_CTTGTA_L003_R1_009_clipped.fastq", "${root}/cutadapt/Fuchs/KCV2_2N3_CTTGTA_L003_R1_010_clipped.fastq"
    ],
    "KCV2_2N4" => ["${root}/cutadapt/Fuchs/KCV2_2N4_GCCTTA_L003_R1_001_clipped.fastq"],
    "KCV2_2N5" => ["${root}/cutadapt/Fuchs/KCV2_2N5_GCTCCA_L003_R1_001_clipped.fastq"],
    "KCV3_1C2" => ["${root}/cutadapt/Fuchs/KCV3_1C2_GGCACA_L004_R1_001_clipped.fastq"],
    "KCV3_1C3" => ["${root}/cutadapt/Fuchs/KCV3_1C3_GGCCTG_L004_R1_001_clipped.fastq"],
    "KCV3_1C4" => ["${root}/cutadapt/Fuchs/KCV3_1C4_TCTACC_L004_R1_001_clipped.fastq"],
    "KCV3_1C5" => ["${root}/cutadapt/Fuchs/KCV3_1C5_TGAAGT_L004_R1_001_clipped.fastq"],
    "KCV3_1C6" => ["${root}/cutadapt/Fuchs/KCV3_1C6_TGCCAT_L004_R1_001_clipped.fastq"],
    "KCV3_2C1" => ["${root}/cutadapt/Fuchs/KCV3_2C1_TGCTGG_L004_R1_001_clipped.fastq"],
    "KCV3_2C2" => ["${root}/cutadapt/Fuchs/KCV3_2C2_TGGCGC_L004_R1_001_clipped.fastq"],
    "KCV3_2C3" => ["${root}/cutadapt/Fuchs/KCV3_2C3_TTCGAA_L004_R1_001_clipped.fastq"],
    "Sample1"  => ["${root}/cutadapt/MT/Sample1_12_clipped.fastq"],
    "Sample2"  => ["${root}/cutadapt/MT/Sample2_12_clipped.fastq"],
    "Sample3"  => ["${root}/cutadapt/MT/Sample3_12_clipped.fastq"],
    "Sample4"  => ["${root}/cutadapt/MT/Sample4_12_clipped.fastq"],
    "Sample5"  => ["${root}/cutadapt/MT/Sample5_12_clipped.fastq"],
  },
  unmappedfiles => {
    "2516-10"  => ["${target_human_dir}/bowtie2_genome/result/2516-10/2516-10.bam.unmapped.fastq"],
    "2516-11"  => ["${target_human_dir}/bowtie2_genome/result/2516-11/2516-11.bam.unmapped.fastq"],
    "2516-12"  => ["${target_human_dir}/bowtie2_genome/result/2516-12/2516-12.bam.unmapped.fastq"],
    "2516-13"  => ["${target_human_dir}/bowtie2_genome/result/2516-13/2516-13.bam.unmapped.fastq"],
    "2516-14"  => ["${target_human_dir}/bowtie2_genome/result/2516-14/2516-14.bam.unmapped.fastq"],
    "2516-15"  => ["${target_human_dir}/bowtie2_genome/result/2516-15/2516-15.bam.unmapped.fastq"],
    "2516-16"  => ["${target_human_dir}/bowtie2_genome/result/2516-16/2516-16.bam.unmapped.fastq"],
    "2516-17"  => ["${target_human_dir}/bowtie2_genome/result/2516-17/2516-17.bam.unmapped.fastq"],
    "2516-18"  => ["${target_human_dir}/bowtie2_genome/result/2516-18/2516-18.bam.unmapped.fastq"],
    "2516-19"  => ["${target_human_dir}/bowtie2_genome/result/2516-19/2516-19.bam.unmapped.fastq"],
    "2516-20"  => ["${target_human_dir}/bowtie2_genome/result/2516-20/2516-20.bam.unmapped.fastq"],
    "2516-21"  => ["${target_human_dir}/bowtie2_genome/result/2516-21/2516-21.bam.unmapped.fastq"],
    "2516-22"  => ["${target_human_dir}/bowtie2_genome/result/2516-22/2516-22.bam.unmapped.fastq"],
    "2516-23"  => ["${target_human_dir}/bowtie2_genome/result/2516-23/2516-23.bam.unmapped.fastq"],
    "2516-24"  => ["${target_human_dir}/bowtie2_genome/result/2516-24/2516-24.bam.unmapped.fastq"],
    "2516-25"  => ["${target_human_dir}/bowtie2_genome/result/2516-25/2516-25.bam.unmapped.fastq"],
    "2516-26"  => ["${target_human_dir}/bowtie2_genome/result/2516-26/2516-26.bam.unmapped.fastq"],
    "2516-27"  => ["${target_human_dir}/bowtie2_genome/result/2516-27/2516-27.bam.unmapped.fastq"],
    "2516-28"  => ["${target_human_dir}/bowtie2_genome/result/2516-28/2516-28.bam.unmapped.fastq"],
    "2516-29"  => ["${target_human_dir}/bowtie2_genome/result/2516-29/2516-29.bam.unmapped.fastq"],
    "2516-30"  => ["${target_human_dir}/bowtie2_genome/result/2516-30/2516-30.bam.unmapped.fastq"],
    "2516-31"  => ["${target_human_dir}/bowtie2_genome/result/2516-31/2516-31.bam.unmapped.fastq"],
    "2516-32"  => ["${target_human_dir}/bowtie2_genome/result/2516-32/2516-32.bam.unmapped.fastq"],
    "2516-33"  => ["${target_human_dir}/bowtie2_genome/result/2516-33/2516-33.bam.unmapped.fastq"],
    "2516-34"  => ["${target_human_dir}/bowtie2_genome/result/2516-34/2516-34.bam.unmapped.fastq"],
    "2516-35"  => ["${target_human_dir}/bowtie2_genome/result/2516-35/2516-35.bam.unmapped.fastq"],
    "2516-36"  => ["${target_human_dir}/bowtie2_genome/result/2516-36/2516-36.bam.unmapped.fastq"],
    "2516-37"  => ["${target_human_dir}/bowtie2_genome/result/2516-37/2516-37.bam.unmapped.fastq"],
    "2516-38"  => ["${target_human_dir}/bowtie2_genome/result/2516-38/2516-38.bam.unmapped.fastq"],
    "2516-39"  => ["${target_human_dir}/bowtie2_genome/result/2516-39/2516-39.bam.unmapped.fastq"],
    "2516-40"  => ["${target_human_dir}/bowtie2_genome/result/2516-40/2516-40.bam.unmapped.fastq"],
    "2516-41"  => ["${target_human_dir}/bowtie2_genome/result/2516-41/2516-41.bam.unmapped.fastq"],
    "2516-42"  => ["${target_human_dir}/bowtie2_genome/result/2516-42/2516-42.bam.unmapped.fastq"],
    "2516-43"  => ["${target_human_dir}/bowtie2_genome/result/2516-43/2516-43.bam.unmapped.fastq"],
    "2516-44"  => ["${target_human_dir}/bowtie2_genome/result/2516-44/2516-44.bam.unmapped.fastq"],
    "2516-45"  => ["${target_human_dir}/bowtie2_genome/result/2516-45/2516-45.bam.unmapped.fastq"],
    "2516-46"  => ["${target_human_dir}/bowtie2_genome/result/2516-46/2516-46.bam.unmapped.fastq"],
    "2516-47"  => ["${target_human_dir}/bowtie2_genome/result/2516-47/2516-47.bam.unmapped.fastq"],
    "2516-48"  => ["${target_human_dir}/bowtie2_genome/result/2516-48/2516-48.bam.unmapped.fastq"],
    "KCV2_1N2" => ["${root}/cutadapt/Fuchs/KCV2_1N2_GCCAAT_L003_R1_001_clipped.fastq"],
    "KCV2_1N3" => [ "${root}/cutadapt/Fuchs/KCV2_1N3_CAGATC_L003_R1_001_clipped.fastq", "${root}/cutadapt/Fuchs/KCV2_1N3_CAGATC_L003_R1_002_clipped.fastq" ],
    "KCV2_1N4" => [ "${root}/cutadapt/Fuchs/KCV2_1N4_ACTTGA_L003_R1_001_clipped.fastq", "${root}/cutadapt/Fuchs/KCV2_1N4_ACTTGA_L003_R1_002_clipped.fastq" ],
    "KCV2_1N5" => [ "${root}/cutadapt/Fuchs/KCV2_1N5_GATCAG_L003_R1_001_clipped.fastq", "${root}/cutadapt/Fuchs/KCV2_1N5_GATCAG_L003_R1_002_clipped.fastq" ],
    "KCV2_1N6" => [ "${root}/cutadapt/Fuchs/KCV2_1N6_TAGCTT_L003_R1_001_clipped.fastq", "${root}/cutadapt/Fuchs/KCV2_1N6_TAGCTT_L003_R1_002_clipped.fastq" ],
    "KCV2_2N1" => [
      "${root}/cutadapt/Fuchs/KCV2_2N1_GGCTAC_L003_R1_001_clipped.fastq", "${root}/cutadapt/Fuchs/KCV2_2N1_GGCTAC_L003_R1_002_clipped.fastq",
      "${root}/cutadapt/Fuchs/KCV2_2N1_GGCTAC_L003_R1_003_clipped.fastq", "${root}/cutadapt/Fuchs/KCV2_2N1_GGCTAC_L003_R1_004_clipped.fastq",
      "${root}/cutadapt/Fuchs/KCV2_2N1_GGCTAC_L003_R1_005_clipped.fastq"
    ],
    "KCV2_2N2" => ["${root}/cutadapt/Fuchs/KCV2_2N2_GCCGCG_L003_R1_001_clipped.fastq"],
    "KCV2_2N3" => [
      "${root}/cutadapt/Fuchs/KCV2_2N3_CTTGTA_L003_R1_001_clipped.fastq", "${root}/cutadapt/Fuchs/KCV2_2N3_CTTGTA_L003_R1_002_clipped.fastq",
      "${root}/cutadapt/Fuchs/KCV2_2N3_CTTGTA_L003_R1_003_clipped.fastq", "${root}/cutadapt/Fuchs/KCV2_2N3_CTTGTA_L003_R1_004_clipped.fastq",
      "${root}/cutadapt/Fuchs/KCV2_2N3_CTTGTA_L003_R1_005_clipped.fastq", "${root}/cutadapt/Fuchs/KCV2_2N3_CTTGTA_L003_R1_006_clipped.fastq",
      "${root}/cutadapt/Fuchs/KCV2_2N3_CTTGTA_L003_R1_007_clipped.fastq", "${root}/cutadapt/Fuchs/KCV2_2N3_CTTGTA_L003_R1_008_clipped.fastq",
      "${root}/cutadapt/Fuchs/KCV2_2N3_CTTGTA_L003_R1_009_clipped.fastq", "${root}/cutadapt/Fuchs/KCV2_2N3_CTTGTA_L003_R1_010_clipped.fastq"
    ],
    "KCV2_2N4" => ["${root}/cutadapt/Fuchs/KCV2_2N4_GCCTTA_L003_R1_001_clipped.fastq"],
    "KCV2_2N5" => ["${root}/cutadapt/Fuchs/KCV2_2N5_GCTCCA_L003_R1_001_clipped.fastq"],
    "KCV3_1C2" => ["${root}/cutadapt/Fuchs/KCV3_1C2_GGCACA_L004_R1_001_clipped.fastq"],
    "KCV3_1C3" => ["${root}/cutadapt/Fuchs/KCV3_1C3_GGCCTG_L004_R1_001_clipped.fastq"],
    "KCV3_1C4" => ["${root}/cutadapt/Fuchs/KCV3_1C4_TCTACC_L004_R1_001_clipped.fastq"],
    "KCV3_1C5" => ["${root}/cutadapt/Fuchs/KCV3_1C5_TGAAGT_L004_R1_001_clipped.fastq"],
    "KCV3_1C6" => ["${root}/cutadapt/Fuchs/KCV3_1C6_TGCCAT_L004_R1_001_clipped.fastq"],
    "KCV3_2C1" => ["${root}/cutadapt/Fuchs/KCV3_2C1_TGCTGG_L004_R1_001_clipped.fastq"],
    "KCV3_2C2" => ["${root}/cutadapt/Fuchs/KCV3_2C2_TGGCGC_L004_R1_001_clipped.fastq"],
    "KCV3_2C3" => ["${root}/cutadapt/Fuchs/KCV3_2C3_TTCGAA_L004_R1_001_clipped.fastq"],
    "Sample1"  => ["${root}/cutadapt/MT/Sample1_12_clipped.fastq"],
    "Sample2"  => ["${root}/cutadapt/MT/Sample2_12_clipped.fastq"],
    "Sample3"  => ["${root}/cutadapt/MT/Sample3_12_clipped.fastq"],
    "Sample4"  => ["${root}/cutadapt/MT/Sample4_12_clipped.fastq"],
    "Sample5"  => ["${root}/cutadapt/MT/Sample5_12_clipped.fastq"],
  },
  bwa_mature => {
    target_dir   => "${target_human_dir}/bwa_miRBase_species",
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
  bwa => {
    target_dir      => "${target_human_dir}/bwa_genome",
    option          => $bwa_option_wholegenome,
    option_samse    => "",
    source_ref      => "fastqfiles",
    fasta_file      => "/data/cqs/shengq1/reference/hg19/hg19_chr.fa",
    estimate_insert => 0,
    source_ref      => "fastqfiles",
    sh_direct       => 1,
    pbs             => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "24",
      "mem"      => "20gb"
    },
  },
  bowtie1 => {
    target_dir    => "${target_human_dir}/bowtie1_genome",
    option        => $bowtie1_option_wholegenome,
    source_ref    => "fastqfiles",
    bowtie1_index => "/data/cqs/guoy1/reference/hg19/bowtie_index/hg19",
    fasta_file    => "/data/cqs/guoy1/reference/hg19/bowtie_index/hg19.fa",
    samonly       => 0,
    sh_direct     => 1,
    pbs           => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "24",
      "mem"      => "20gb"
    },
  },
  bowtie2 => {
    target_dir    => "${target_human_dir}/bowtie2_genome",
    option        => $bowtie2_option_wholegenome,
    source_ref    => "fastqfiles",
    bowtie2_index => "/data/cqs/guoy1/reference/hg19/bowtie2_index/hg19",
    fasta_file    => "/data/cqs/guoy1/reference/hg19/bowtie2_index/hg19.fa",
    samonly       => 0,
    sh_direct     => 1,
    pbs           => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "24",
      "mem"      => "20gb"
    },
  },
  mirna_count_bwa => {
    target_dir   => "${target_human_dir}/bwa_genome",
    option       => "",
    source_ref   => "bwa",
    cqs_tools    => $cqs_tools,
    gff_file     => $hsa_gffs,
    sh_direct    => 1,
    fasta_format => 0,
    pbs          => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "24",
      "mem"      => "20gb"
    },
  },
  mirna_count_bowtie1 => {
    target_dir   => "${target_human_dir}/bowtie1_genome",
    option       => "",
    source_ref   => "bowtie1",
    cqs_tools    => $cqs_tools,
    gff_file     => $hsa_gffs,
    sh_direct    => 1,
    fasta_format => 0,
    pbs          => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "24",
      "mem"      => "20gb"
    },
  },
  mirna_count_bowtie2 => {
    target_dir   => "${target_human_dir}/bowtie2_genome",
    option       => "",
    source_ref   => "bowtie2",
    cqs_tools    => $cqs_tools,
    gff_file     => $hsa_gffs,
    sh_direct    => 1,
    fasta_format => 0,
    pbs          => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "24",
      "mem"      => "20gb"
    },
  },
  shrimp2 => {
    target_dir   => "${target_human_dir}/shrimp2",
    option       => $shrimp2_option,
    unmapped_ref => "bowtie2",
    genome_index => "/data/cqs/guoy1/reference/hg19/shrimp2_index/hg19_chr-ls",
    is_mirna     => 1,
    output_bam   => 1,
    sh_direct    => 1,
    pbs          => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "24",
      "mem"      => "20gb"
    },
  },
  bowtie2_mature => {
    target_dir    => "${target_human_dir}/bowtie2_mature",
    option        => $bowtie2_option,
    source_ref    => "fastqfiles",
    bowtie2_index => "/data/cqs/shengq1/reference/miRBase19/hsa.mature.dna",
    fasta_file    => "/data/cqs/shengq1/reference/miRBase19/hsa.mature.dna.fa",
    samonly       => 0,
    sh_direct     => 1,
    pbs           => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "24",
      "mem"      => "20gb"
    },
  },
  bowtie2_unmapped_mature => {
    target_dir    => "${target_human_dir}/bowtie2_unmapped_mature",
    option        => $bowtie2_option,
    unmapped_ref => "bowtie2",
    bowtie2_index => "/data/cqs/shengq1/reference/miRBase19/hsa.mature.dna",
    fasta_file    => "/data/cqs/shengq1/reference/miRBase19/hsa.mature.dna.fa",
    samonly       => 0,
    sh_direct     => 1,
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

#bwa_by_pbs_single( $config_mirna, "bwa_mature" );
#bwa_by_pbs_single( $config_mirna, "bwa_hairpin" );
#bwa_by_pbs_single( $config_mirna, "bwa_illumina" );

#novoalign( $config_mirna, "novoalign_mature" );
#novoalign( $config_mirna, "novoalign_hairpin" );
#novoalign( $config_mirna, "novoalign_illumina" );

#bwa_by_pbs_single( $config_rat,   "bwa_mature" );
#bwa_by_pbs_single( $config_human, "bwa_mature" );

#bwa_by_pbs_single( $config_rat, "bwa" );
#bwa_by_pbs_single( $config_human, "bwa" );
#
#bowtie2( $config_rat,   "bowtie2" );
#bowtie2( $config_human, "bowtie2" );
#
#bowtie1( $config_rat,   "bowtie1" );
#bowtie1( $config_human, "bowtie1" );
#
#mirna_count( $config_rat,   "mirna_count_bwa" );
#mirna_count( $config_human, "mirna_count_bwa" );

#mirna_count( $config_rat,   "mirna_count_bowtie1" );
#mirna_count( $config_human, "mirna_count_bowtie1" );

#mirna_count( $config_rat,   "mirna_count_bowtie2" );
#mirna_count( $config_human, "mirna_count_bowtie2" );

shrimp2($config_rat, "shrimp2");
#shrimp2($config_human, "shrimp2");

#bowtie2( $config_rat,   "bowtie2_mature" );
#bowtie2( $config_human, "bowtie2_mature" );

#bowtie2( $config_human, "bowtie2_unmapped_mature" );

1;
