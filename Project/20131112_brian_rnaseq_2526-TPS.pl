#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;
use CQS::FileUtils;

my $target_dir = create_directory_or_die("/scratch/cqs/shengq1/rnaseq/20131104_brian_rnaseq_SRA");

my $fasta_file           = "/data/cqs/guoy1/reference/hg19/bwa_index_0.7.4/hg19_chr.fa";
my $transcript_gtf       = "/scratch/cqs/shengq1/references/hg19/Homo_sapiens.GRCh37.73.gtf";
my $transcript_gtf_index = "/scratch/cqs/shengq1/gtfindex/hg19_GRCh37_73";

my $hg19_gff = "/scratch/cqs/shengq1/references/hg19/dexseq_gff/Homo_sapiens.GRCh37.73.dexseq.gff";
my $hg19_map = "/scratch/cqs/shengq1/references/hg19/Homo_sapiens.GRCh37.73.map";

my $bowtie2_index = "/data/cqs/guoy1/reference/hg19/bowtie2_index/hg19";

my $annovar_param = "-protocol refGene,snp137,cosmic64,esp6500si_all,1000g2012apr_all -operation g,f,f,f,f --remove --otherinfo";
my $annovar_db    = "/scratch/cqs/shengq1/references/annovar/humandb/";

my $cqstools = "/home/shengq1/cqstools/CQS.Tools.exe";

my $email = "quanhu.sheng\@vanderbilt.edu";
my $task  = "2526-TPS";

my $files = {
  "2526-TPS-01" => [ "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-07/2526-TPS-1_1.fastq.gz",  "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-07/2526-TPS-1_2.fastq.gz" ],
  "2526-TPS-10" => [ "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-07/2526-TPS-10_1.fastq.gz", "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-07/2526-TPS-10_2.fastq.gz" ],
  "2526-TPS-11" => [ "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-07/2526-TPS-11_1.fastq.gz", "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-07/2526-TPS-11_2.fastq.gz" ],
  "2526-TPS-12" => [ "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-07/2526-TPS-12_1.fastq.gz", "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-07/2526-TPS-12_2.fastq.gz" ],
  "2526-TPS-13" => [ "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-07/2526-TPS-13_1.fastq.gz", "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-07/2526-TPS-13_2.fastq.gz" ],
  "2526-TPS-14" => [ "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-07/2526-TPS-14_1.fastq.gz", "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-07/2526-TPS-14_2.fastq.gz" ],
  "2526-TPS-15" => [ "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-07/2526-TPS-15_1.fastq.gz", "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-07/2526-TPS-15_2.fastq.gz" ],
  "2526-TPS-16" => [ "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-07/2526-TPS-16_1.fastq.gz", "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-07/2526-TPS-16_2.fastq.gz" ],
  "2526-TPS-17" => [ "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-07/2526-TPS-17_1.fastq.gz", "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-07/2526-TPS-17_2.fastq.gz" ],
  "2526-TPS-18" => [ "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-07/2526-TPS-18_1.fastq.gz", "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-07/2526-TPS-18_2.fastq.gz" ],
  "2526-TPS-02" => [ "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-07/2526-TPS-2_1.fastq.gz",  "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-07/2526-TPS-2_2.fastq.gz" ],
  "2526-TPS-20" => [ "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-07/2526-TPS-20_1.fastq.gz", "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-07/2526-TPS-20_2.fastq.gz" ],
  "2526-TPS-21" => [ "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-07/2526-TPS-21_1.fastq.gz", "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-07/2526-TPS-21_2.fastq.gz" ],
  "2526-TPS-22" => [ "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-07/2526-TPS-22_1.fastq.gz", "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-07/2526-TPS-22_2.fastq.gz" ],
  "2526-TPS-23" => [ "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-07/2526-TPS-23_1.fastq.gz", "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-07/2526-TPS-23_2.fastq.gz" ],
  "2526-TPS-24" => [ "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-07/2526-TPS-24_1.fastq.gz", "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-07/2526-TPS-24_2.fastq.gz" ],
  "2526-TPS-25" => [ "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-07/2526-TPS-25_1.fastq.gz", "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-07/2526-TPS-25_2.fastq.gz" ],
  "2526-TPS-26" => [ "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-07/2526-TPS-26_1.fastq.gz", "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-07/2526-TPS-26_2.fastq.gz" ],
  "2526-TPS-27" => [ "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-07/2526-TPS-27_1.fastq.gz", "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-07/2526-TPS-27_2.fastq.gz" ],
  "2526-TPS-28" => [ "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-07/2526-TPS-28_1.fastq.gz", "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-07/2526-TPS-28_2.fastq.gz" ],
  "2526-TPS-29" => [ "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-07/2526-TPS-29_1.fastq.gz", "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-07/2526-TPS-29_2.fastq.gz" ],
  "2526-TPS-03" => [ "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-07/2526-TPS-3_1.fastq.gz",  "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-07/2526-TPS-3_2.fastq.gz" ],
  "2526-TPS-30" => [ "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-07/2526-TPS-30_1.fastq.gz", "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-07/2526-TPS-30_2.fastq.gz" ],
  "2526-TPS-31" => [ "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-07/2526-TPS-31_1.fastq.gz", "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-07/2526-TPS-31_2.fastq.gz" ],
  "2526-TPS-32" => [ "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-07/2526-TPS-32_1.fastq.gz", "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-07/2526-TPS-32_2.fastq.gz" ],
  "2526-TPS-33" => [ "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-07/2526-TPS-33_1.fastq.gz", "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-07/2526-TPS-33_2.fastq.gz" ],
  "2526-TPS-34" => [ "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-07/2526-TPS-34_1.fastq.gz", "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-07/2526-TPS-34_2.fastq.gz" ],
  "2526-TPS-35" => [ "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-07/2526-TPS-35_1.fastq.gz", "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-07/2526-TPS-35_2.fastq.gz" ],
  "2526-TPS-36" => [ "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-07/2526-TPS-36_1.fastq.gz", "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-07/2526-TPS-36_2.fastq.gz" ],
  "2526-TPS-37" => [ "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-07/2526-TPS-37_1.fastq.gz", "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-07/2526-TPS-37_2.fastq.gz" ],
  "2526-TPS-38" => [ "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-07/2526-TPS-38_1.fastq.gz", "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-07/2526-TPS-38_2.fastq.gz" ],
  "2526-TPS-39" => [ "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-27/2526-TPS-39_1.fastq.gz", "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-27/2526-TPS-39_2.fastq.gz" ],
  "2526-TPS-40" => [ "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-27/2526-TPS-40_1.fastq.gz", "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-27/2526-TPS-40_2.fastq.gz" ],
  "2526-TPS-41" => [ "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-27/2526-TPS-41_1.fastq.gz", "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-27/2526-TPS-41_2.fastq.gz" ],
  "2526-TPS-42" => [ "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-27/2526-TPS-42_1.fastq.gz", "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-27/2526-TPS-42_2.fastq.gz" ],
  "2526-TPS-43" => [ "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-27/2526-TPS-43_1.fastq.gz", "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-27/2526-TPS-43_2.fastq.gz" ],
  "2526-TPS-44" => [ "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-27/2526-TPS-44_1.fastq.gz", "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-27/2526-TPS-44_2.fastq.gz" ],
  "2526-TPS-45" => [ "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-27/2526-TPS-45_1.fastq.gz", "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-27/2526-TPS-45_2.fastq.gz" ],
  "2526-TPS-46" => [ "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-27/2526-TPS-46_1.fastq.gz", "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-27/2526-TPS-46_2.fastq.gz" ],
  "2526-TPS-47" => [ "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-27/2526-TPS-47_1.fastq.gz", "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-27/2526-TPS-47_2.fastq.gz" ],
  "2526-TPS-48" => [ "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-27/2526-TPS-48_1.fastq.gz", "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-27/2526-TPS-48_2.fastq.gz" ],
  "2526-TPS-49" => [ "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-27/2526-TPS-49_1.fastq.gz", "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-27/2526-TPS-49_2.fastq.gz" ],
  "2526-TPS-05" => [ "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-07/2526-TPS-5_1.fastq.gz",  "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-07/2526-TPS-5_2.fastq.gz" ],
  "2526-TPS-50" => [ "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-27/2526-TPS-50_1.fastq.gz", "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-27/2526-TPS-50_2.fastq.gz" ],
  "2526-TPS-51" => [ "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-27/2526-TPS-51_1.fastq.gz", "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-27/2526-TPS-51_2.fastq.gz" ],
  "2526-TPS-52" => [ "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-27/2526-TPS-52_1.fastq.gz", "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-27/2526-TPS-52_2.fastq.gz" ],
  "2526-TPS-53" => [ "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-27/2526-TPS-53_1.fastq.gz", "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-27/2526-TPS-53_2.fastq.gz" ],
  "2526-TPS-54" => [ "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-27/2526-TPS-54_1.fastq.gz", "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-27/2526-TPS-54_2.fastq.gz" ],
  "2526-TPS-55" => [ "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-27/2526-TPS-55_1.fastq.gz", "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-27/2526-TPS-55_2.fastq.gz" ],
  "2526-TPS-56" => [ "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-27/2526-TPS-56_1.fastq.gz", "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-27/2526-TPS-56_2.fastq.gz" ],
  "2526-TPS-57" => [ "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-27/2526-TPS-57_1.fastq.gz", "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-27/2526-TPS-57_2.fastq.gz" ],
  "2526-TPS-58" => [ "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-27/2526-TPS-58_1.fastq.gz", "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-27/2526-TPS-58_2.fastq.gz" ],
  "2526-TPS-59" => [ "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-27/2526-TPS-59_1.fastq.gz", "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-27/2526-TPS-59_2.fastq.gz" ],
  "2526-TPS-06" => [ "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-07/2526-TPS-6_1.fastq.gz",  "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-07/2526-TPS-6_2.fastq.gz" ],
  "2526-TPS-60" => [ "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-27/2526-TPS-60_1.fastq.gz", "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-27/2526-TPS-60_2.fastq.gz" ],
  "2526-TPS-62" => [ "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-27/2526-TPS-62_1.fastq.gz", "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-27/2526-TPS-62_2.fastq.gz" ],
  "2526-TPS-63" => [ "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-27/2526-TPS-63_1.fastq.gz", "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-27/2526-TPS-63_2.fastq.gz" ],
  "2526-TPS-64" => [ "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-27/2526-TPS-64_1.fastq.gz", "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-27/2526-TPS-64_2.fastq.gz" ],
  "2526-TPS-67" => [ "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-27/2526-TPS-67_1.fastq.gz", "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-27/2526-TPS-67_2.fastq.gz" ],
  "2526-TPS-68" => [ "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-27/2526-TPS-68_1.fastq.gz", "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-27/2526-TPS-68_2.fastq.gz" ],
  "2526-TPS-69" => [ "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-27/2526-TPS-69_1.fastq.gz", "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-27/2526-TPS-69_2.fastq.gz" ],
  "2526-TPS-07" => [ "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-07/2526-TPS-7_1.fastq.gz",  "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-07/2526-TPS-7_2.fastq.gz" ],
  "2526-TPS-70" => [ "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-27/2526-TPS-70_1.fastq.gz", "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-27/2526-TPS-70_2.fastq.gz" ],
  "2526-TPS-71" => [ "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-27/2526-TPS-71_1.fastq.gz", "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-27/2526-TPS-71_2.fastq.gz" ],
  "2526-TPS-08" => [ "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-07/2526-TPS-8_1.fastq.gz",  "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-07/2526-TPS-8_2.fastq.gz" ],
  "2526-TPS-09" => [ "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-07/2526-TPS-9_1.fastq.gz",  "/autofs/blue_sequencer/Runs/projects/2526-TPS/2013-06-07/2526-TPS-9_2.fastq.gz" ],
};

my $config = {
  general => { task_name => $task },
  files   => $files,
  fastqc  => {
    class      => "FastQC",
    perform    => 1,
    target_dir => "${target_dir}/fastqc",
    option     => "",
    source_ref => "files",
    sh_direct  => 0,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=2",
      "walltime" => "2",
      "mem"      => "10gb"
    },
  },
  tophat2 => {
    class                => "Tophat2",
    perform              => 1,
    target_dir           => "${target_dir}/tophat2",
    option               => "--segment-length 25 -r 0 -p 6",
    source_ref           => "files",
    bowtie2_index        => $bowtie2_index,
    transcript_gtf       => $transcript_gtf,
    transcript_gtf_index => $transcript_gtf_index,
    sh_direct            => 0,
    pbs                  => {
      "email"    => $email,
      "nodes"    => "1:ppn=6",
      "walltime" => "72",
      "mem"      => "30gb"
    },
  },
  sortbam => {
    class         => "Sortbam",
    perform       => 1,
    target_dir    => "${target_dir}/sortname",
    option        => "",
    source_ref    => "tophat2",
    sort_by_query => 1,
    sh_direct     => 0,
    pbs           => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "20gb"
    },
  },
  htseqcount => {
    class      => "HTSeqCount",
    perform    => 1,
    target_dir => "${target_dir}/htseqcount",
    option     => "",
    source_ref => "sortbam",
    gff_file   => $transcript_gtf,
    sh_direct  => 0,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  genetable => {
    class         => "CQSDatatable",
    perform       => 1,
    target_dir    => "${target_dir}/genetable",
    option        => "-p ENS --noheader -o ${task}_gene.count",
    source_ref    => "htseqcount",
    name_map_file => $hg19_map,
    cqs_tools     => $cqstools,
    sh_direct     => 1,
    pbs           => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "10",
      "mem"      => "10gb"
    },
  },
  dexseqcount => {
    class        => "DexseqCount",
    perform      => 1,
    target_dir   => "${target_dir}/dexseqcount",
    option       => "",
    source_ref   => "tophat2",
    gff_file     => $hg19_gff,
    dexseq_count => "/home/shengq1/pylibs/bin/dexseq_count.py",
    sh_direct    => 0,
    pbs          => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  exontable => {
    class         => "CQSDatatable",
    perform       => 1,
    target_dir    => "${target_dir}/exontable",
    option        => "-p ENS --noheader -o ${task}_exon.count",
    name_map_file => $hg19_map,
    source_ref    => "dexseqcount",
    cqs_tools     => $cqstools,
    sh_direct     => 1,
    pbs           => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "10",
      "mem"      => "10gb"
    },
  },
  varscan2 => {
    class           => "VarScan2::Mpileup2snp",
    perform         => 1,
    target_dir      => "${target_dir}/varscan2",
    option          => "--min-coverage 10",
    mpileup_options => "-q 20",
    java_option     => "-Xmx40g",
    source_ref      => "tophat2",
    fasta_file      => $fasta_file,
    somatic_p_value => 0.05,
    sh_direct       => 0,
    VarScan2_jar    => "/home/shengq1/local/bin/VarScan.v2.3.5.jar",
    pbs             => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  annovar_varscan2 => {
    class      => "Annovar",
    perform    => 1,
    target_dir => "${target_dir}/varscan2",
    option     => $annovar_param,
    source_ref => [ "varscan2", "\.vcf\$" ],
    annovar_db => $annovar_db,
    buildver   => "hg19",
    sh_direct  => 1,
    isvcf      => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "10gb"
    },
  },
};

performConfig($config);

1;
