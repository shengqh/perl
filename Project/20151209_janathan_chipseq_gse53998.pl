#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;
use CQS::FileUtils;

my $target_dir = create_directory_or_die("/scratch/cqs/shengq1/chipseq/20151208_gse53998");

my $fasta_file   = "/scratch/cqs/shengq1/references/hg18_chrM/bowtie_index_1.1.2/hg18_chrM.fa";
my $bowtie_index = "/scratch/cqs/shengq1/references/hg18_chrM/bowtie_index_1.1.2/hg18_chrM";
my $cqstools     = "/home/shengq1/cqstools/CQS.Tools.exe";

my $email = "quanhu.sheng\@vanderbilt.edu";
my $task  = "gse53998";

my $config = {
  general => { task_name => $task },
  files   => {
    "EC_BRD4_CON"         => ["/gpfs21/scratch/cqs/shengq1/chipseq/20151208_gse53998/sra/SRR1106515.sra"],
    "EC_BRD4_TNF"         => ["/gpfs21/scratch/cqs/shengq1/chipseq/20151208_gse53998/sra/SRR1106516.sra"],
    "EC_BRD4_TNF_JQ1"     => ["/gpfs21/scratch/cqs/shengq1/chipseq/20151208_gse53998/sra/SRR1106517.sra"],
    "EC_H3K27AC_CON"      => ["/gpfs21/scratch/cqs/shengq1/chipseq/20151208_gse53998/sra/SRR1106518.sra"],
    "EC_H3K27AC_TNF"      => ["/gpfs21/scratch/cqs/shengq1/chipseq/20151208_gse53998/sra/SRR1106519.sra"],
    "EC_H3K27AC_TNF_JQ1"  => ["/gpfs21/scratch/cqs/shengq1/chipseq/20151208_gse53998/sra/SRR1106520.sra"],
    "EC_H3K4ME3_CON"      => ["/gpfs21/scratch/cqs/shengq1/chipseq/20151208_gse53998/sra/SRR1106521.sra"],
    "EC_H3K4ME3_TNF"      => ["/gpfs21/scratch/cqs/shengq1/chipseq/20151208_gse53998/sra/SRR1106522.sra"],
    "EC_H3K4ME3_TNF_JQ1"  => ["/gpfs21/scratch/cqs/shengq1/chipseq/20151208_gse53998/sra/SRR1106523.sra"],
    "EC_P65_CON"          => ["/gpfs21/scratch/cqs/shengq1/chipseq/20151208_gse53998/sra/SRR1106524.sra"],
    "EC_P65_TNF"          => ["/gpfs21/scratch/cqs/shengq1/chipseq/20151208_gse53998/sra/SRR1106525.sra"],
    "EC_P65_TNF_JQ1"      => ["/gpfs21/scratch/cqs/shengq1/chipseq/20151208_gse53998/sra/SRR1106526.sra"],
    "EC_RNA_POL2_CON"     => ["/gpfs21/scratch/cqs/shengq1/chipseq/20151208_gse53998/sra/SRR1106527.sra"],
    "EC_RNA_POL2_TNF"     => ["/gpfs21/scratch/cqs/shengq1/chipseq/20151208_gse53998/sra/SRR1106528.sra"],
    "EC_RNA_POL2_TNF_JQ1" => ["/gpfs21/scratch/cqs/shengq1/chipseq/20151208_gse53998/sra/SRR1106529.sra"],
    "EC_WCE_CON"          => ["/gpfs21/scratch/cqs/shengq1/chipseq/20151208_gse53998/sra/SRR1106530.sra"],
    "EC_WCE_CON_2"        => ["/gpfs21/scratch/cqs/shengq1/chipseq/20151208_gse53998/sra/SRR1106531.sra"],
    "EC_WCE_TNF"          => ["/gpfs21/scratch/cqs/shengq1/chipseq/20151208_gse53998/sra/SRR1106532.sra"],
    "EC_WCE_TNF_2"        => ["/gpfs21/scratch/cqs/shengq1/chipseq/20151208_gse53998/sra/SRR1106533.sra"],
    "EC_WCE_TNF_JQ1"      => ["/gpfs21/scratch/cqs/shengq1/chipseq/20151208_gse53998/sra/SRR1106534.sra"],
    "EC_WCE_TNF_JQ1_2"    => ["/gpfs21/scratch/cqs/shengq1/chipseq/20151208_gse53998/sra/SRR1106535.sra"],
  },
  sra2fastq => {
    class      => "SRA::FastqDump",
    perform    => 1,
    ispaired   => 0,
    target_dir => "${target_dir}/FastqDump",
    option     => "",
    source_ref => "files",
    sh_direct  => 0,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "10",
      "mem"      => "10gb"
    },
  },
  fastqc_pre_trim => {
    class      => "QC::FastQC",
    perform    => 1,
    target_dir => "${target_dir}/fastqc_pre_trim",
    option     => "",
    source_ref => "sra2fastq",
    sh_direct  => 0,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=2",
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
    class      => "Cutadapt",
    perform    => 1,
    target_dir => "${target_dir}/cutadapt",
    option     => "-O 10 -m 30",
    source_ref => "sra2fastq",
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
    class         => "Alignment::Bowtie1",
    perform       => 1,
    target_dir    => "${target_dir}/bowtie1",
    option        => "-v 1 -m 1 --best --strata",
    fasta_file    => $fasta_file,
    source_ref    => [ "cutadapt", ".fastq.gz\$" ],
    bowtie1_index => $bowtie_index,
    sh_direct     => 0,
    pbs           => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
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
      step_1 => [ "sra2fastq",               "fastqc_pre_trim", "cutadapt", "fastqc_post_trim", "bowtie1" ],
      step_2 => [ "fastqc_pre_trim_summary", "fastqc_post_trim_summary" ],
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

1;
