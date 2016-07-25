#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;
use CQS::FileUtils;
use Data::Dumper;

my $target_dir = create_directory_or_die("/scratch/cqs/shengq1/brown/20160721_chipseq_3530");
my $task       = "3530";

my $fasta_file   = "/scratch/cqs/shengq1/references/gencode/hg19/bowtie_index_1.1.2/GRCh37.p13.genome.fa";
my $bowtie_index = "/scratch/cqs/shengq1/references/gencode/hg19/bowtie_index_1.1.2/GRCh37.p13.genome";
my $cqstools     = "/home/shengq1/cqstools/cqstools.exe";

my $email = "quanhu.sheng\@vanderbilt.edu";

my $config = {
  general => { task_name => $task },
  files   => {
    "CR_Y_Input"    => ["/gpfs21/scratch/cqs/shengq1/brown/data/3530/3530-KB-4_1_sequence.txt.gz"],
    "CR_S_Input"    => ["/gpfs21/scratch/cqs/shengq1/brown/data/3530/3530-KB-5_1_sequence.txt.gz"],
    "G145R_C_Input" => ["/gpfs21/scratch/cqs/shengq1/brown/data/3530/3530-KB-6_1_sequence.txt.gz"],
    "G145R_P_Input" => ["/gpfs21/scratch/cqs/shengq1/brown/data/3530/3530-KB-7_1_sequence.txt.gz"],
    "C2_1_Input"    => ["/gpfs21/scratch/cqs/shengq1/brown/data/3530/3530-KB-8_1_sequence.txt.gz"],
    "C4_3_Input"    => ["/gpfs21/scratch/cqs/shengq1/brown/data/3530/3530-KB-9_1_sequence.txt.gz"],
    "CR_Y_TBX5"     => ["/gpfs21/scratch/cqs/shengq1/brown/data/3530/3530-KB-12_1_sequence.txt.gz"],
    "CR_S_TBX5"     => ["/gpfs21/scratch/cqs/shengq1/brown/data/3530/3530-KB-13_1_sequence.txt.gz"],
    "G145R_C_TBX5"  => ["/gpfs21/scratch/cqs/shengq1/brown/data/3530/3530-KB-14_1_sequence.txt.gz"],
    "G145R_P_TBX5"  => ["/gpfs21/scratch/cqs/shengq1/brown/data/3530/3530-KB-15_1_sequence.txt.gz"],
    "C2_1_TBX5"     => ["/gpfs21/scratch/cqs/shengq1/brown/data/3530/3530-KB-10_1_sequence.txt.gz"],
    "C4_3_TBX5"     => ["/gpfs21/scratch/cqs/shengq1/brown/data/3530/3530-KB-11_1_sequence.txt.gz"],
  },
  treatments => {
    "CR_Y_TBX5"    => ["CR_Y_TBX5"],
    "CR_S_TBX5"    => ["CR_S_TBX5"],
    "G145R_C_TBX5" => ["G145R_C_TBX5"],
    "G145R_P_TBX5" => ["G145R_P_TBX5"],
    "C2_1_TBX5"    => ["C2_1_TBX5"],
    "C4_3_TBX5"    => ["C4_3_TBX5"],
  },
  controls => {
    "CR_Y_TBX5"    => ["CR_Y_Input"],
    "CR_S_TBX5"    => ["CR_S_Input"],
    "G145R_C_TBX5" => ["G145R_C_Input"],
    "G145R_P_TBX5" => ["G145R_P_Input"],
    "C2_1_TBX5"    => ["C2_1_Input"],
    "C4_3_TBX5"    => ["C4_3_Input"],
  },
  depthgroups => {
    "CR_Y_TBX5"    => [ "CR_Y_TBX5",    "CR_Y_Input" ],
    "CR_S_TBX5"    => [ "CR_S_TBX5",    "CR_S_Input" ],
    "G145R_C_TBX5" => [ "G145R_C_TBX5", "G145R_C_Input" ],
    "G145R_P_TBX5" => [ "G145R_P_TBX5", "G145R_P_Input" ],
    "C2_1_TBX5"    => [ "C2_1_TBX5",    "C2_1_Input" ],
    "C4_3_TBX5"    => [ "C4_3_TBX5",    "C4_3_Input" ],
  },
  fastqc_pre_trim => {
    class      => "QC::FastQC",
    perform    => 1,
    target_dir => "${target_dir}/fastqc_pre_trim",
    option     => "",
    source_ref => ["files"],
    sh_direct  => 0,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "2",
      "mem"      => "10gb"
    },
  },
  fastqc_pre_trim_summary => {
    class      => "QC::FastQCSummary",
    perform    => 1,
    sh_direct  => 1,
    target_dir => "${target_dir}/fastqc_pre_trim",
    cqstools   => $cqstools,
    option     => "",
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "2",
      "mem"      => "10gb"
    },
  },
  cutadapt => {
    class      => "Trimmer::Cutadapt",
    perform    => 1,
    target_dir => "${target_dir}/cutadapt",
    option     => "-m 30 --trim-n",
    source_ref => ["files"],
    adapter    => "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",    #trueseq adapter
    extension  => "_clipped.fastq",
    sh_direct  => 0,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "24",
      "mem"      => "20gb"
    },
  },
  fastqc_post_trim => {
    class      => "QC::FastQC",
    perform    => 1,
    target_dir => "${target_dir}/fastqc_post_trim",
    option     => "",
    sh_direct  => 1,
    source_ref => [ "cutadapt", ".fastq.gz" ],
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "2",
      "mem"      => "10gb"
    },
  },
  fastqc_post_trim_summary => {
    class      => "QC::FastQCSummary",
    perform    => 1,
    sh_direct  => 1,
    target_dir => "${target_dir}/fastqc_post_trim",
    cqstools   => $cqstools,
    option     => "",
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "2",
      "mem"      => "10gb"
    },
  },
  fastq_len => {
    class      => "CQS::FastqLen",
    perform    => 1,
    target_dir => "$target_dir/fastq_len",
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
  bowtie1 => {
    class                   => "Alignment::Bowtie1",
    perform                 => 1,
    target_dir              => "${target_dir}/bowtie1",
    option                  => "-v 1 -m 1 --best --strata",
    fasta_file              => $fasta_file,
    source_ref              => [ "cutadapt", ".fastq.gz\$" ],
    bowtie1_index           => $bowtie_index,
    chromosome_grep_pattern => "\"^chr\"",
    sh_direct               => 0,
    pbs                     => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  macs1callpeak => {
    class        => "Chipseq::MACS",
    perform      => 1,
    target_dir   => "${target_dir}/macs1callpeak",
    option       => "-p 1e-9 -w -S --space=50",
    source_ref   => "bowtie1",
    groups_ref   => "treatments",
    controls_ref => "controls",
    sh_direct    => 0,
    pbs          => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  macs1callpeak_depth => {
    class         => "Visualization::Depth",
    perform       => 1,
    target_dir    => "${target_dir}/macs1callpeak_depth",
    option        => "",
    source_ref    => [ "macs1callpeak", ".name.bed" ],
    groups_ref    => "depthgroups",
    bam_files_ref => "bowtie1",
    sh_direct     => 1,
    pbs           => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "24",
      "mem"      => "10gb"
    },
  },
  macs2callpeak => {
    class        => "Chipseq::MACS2Callpeak",
    perform      => 1,
    target_dir   => "${target_dir}/macs2callpeak",
    option       => "-g hs --broad -B -p 1e-9",
    source_ref   => "bowtie1",
    groups_ref   => "treatments",
    controls_ref => "controls",
    sh_direct    => 0,
    pbs          => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  macs2callpeak_depth => {
    class         => "Visualization::Depth",
    perform       => 1,
    target_dir    => "${target_dir}/macs2callpeak_depth",
    option        => "",
    source_ref    => [ "macs2callpeak", ".bed\$" ],
    groups_ref    => "depthgroups",
    bam_files_ref => "bowtie1",
    sh_direct     => 1,
    pbs           => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "24",
      "mem"      => "10gb"
    },
  },
  sequencetask => {
    class      => "CQS::SequenceTask",
    perform    => 1,
    target_dir => "${target_dir}/sequencetask",
    option     => "",
    source     => {
      step_1 => [ "fastqc_pre_trim",         "cutadapt",                 "fastqc_post_trim", "fastq_len",           "bowtie1" ],
      step_2 => [ "fastqc_pre_trim_summary", "fastqc_post_trim_summary", "macs1callpeak",    "macs1callpeak_depth", "macs2callpeak", "macs2callpeak_depth" ],
    },
    sh_direct => 0,
    pbs       => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
};

performConfig($config);
#performTask( $config, "depth" );

1;
