#!/usr/bin/perl
use strict;
use warnings;

use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;

my $target_dir = "/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics";
my $cqstools   = "/home/shengq1/cqstools/CQS.Tools.exe";
my $email      = "quanhu.sheng\@vanderbilt.edu";

my $fasta_file           = "/scratch/cqs/shengq1/references/hg20/bowtie2_index_2.2.3/hg20.fa";
my $bowtie2_index        = "/scratch/cqs/shengq1/references/hg20/bowtie2_index_2.2.3/hg20";
my $transcript_gtf       = "/scratch/cqs/shengq1/references/ensembl/Homo_sapiens.GRCh38.78.gtf";
my $transcript_gtf_index = "/scratch/cqs/shengq1/references/ensembl/index/Homo_sapiens.GRCh38.78";
my $name_map_file        = "/scratch/cqs/shengq1/references/ensembl/Homo_sapiens.GRCh38.78.map";

my $config = {
  general => { task_name => "lichen" },

  files => {
    "DNA_273_SRR926701" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/DNA_273_SRR926701.sra"],
    "DNA_273_SRR926707" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/DNA_273_SRR926707.sra"],
    "DNA_273_SRR926711" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/DNA_273_SRR926711.sra"],
    "DNA_273_SRR926713" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/DNA_273_SRR926713.sra"],
    "DNA_273_SRR926719" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/DNA_273_SRR926719.sra"],
    "DNA_273_SRR926721" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/DNA_273_SRR926721.sra"],
    "DNA_273_SRR926727" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/DNA_273_SRR926727.sra"],
    "DNA_273_SRR926729" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/DNA_273_SRR926729.sra"],
    "DNA_283_SRR926702" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/DNA_283_SRR926702.sra"],
    "DNA_283_SRR926706" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/DNA_283_SRR926706.sra"],
    "DNA_283_SRR926710" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/DNA_283_SRR926710.sra"],
    "DNA_283_SRR926714" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/DNA_283_SRR926714.sra"],
    "DNA_283_SRR926718" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/DNA_283_SRR926718.sra"],
    "DNA_283_SRR926722" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/DNA_283_SRR926722.sra"],
    "DNA_283_SRR926726" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/DNA_283_SRR926726.sra"],
    "DNA_283_SRR926730" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/DNA_283_SRR926730.sra"],
    "DNA_311_SRR926703" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/DNA_311_SRR926703.sra"],
    "DNA_311_SRR926705" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/DNA_311_SRR926705.sra"],
    "DNA_311_SRR926709" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/DNA_311_SRR926709.sra"],
    "DNA_311_SRR926715" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/DNA_311_SRR926715.sra"],
    "DNA_311_SRR926717" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/DNA_311_SRR926717.sra"],
    "DNA_311_SRR926723" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/DNA_311_SRR926723.sra"],
    "DNA_311_SRR926725" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/DNA_311_SRR926725.sra"],
    "DNA_311_SRR926731" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/DNA_311_SRR926731.sra"],
    "DNA_321_SRR926704" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/DNA_321_SRR926704.sra"],
    "DNA_321_SRR926708" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/DNA_321_SRR926708.sra"],
    "DNA_321_SRR926712" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/DNA_321_SRR926712.sra"],
    "DNA_321_SRR926716" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/DNA_321_SRR926716.sra"],
    "DNA_321_SRR926720" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/DNA_321_SRR926720.sra"],
    "DNA_321_SRR926724" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/DNA_321_SRR926724.sra"],
    "DNA_321_SRR926728" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/DNA_321_SRR926728.sra"],
    "DNA_321_SRR926732" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/DNA_321_SRR926732.sra"],
    "RNA_273_SRR926693" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/RNA_273_SRR926693.sra"],
    "RNA_273_SRR926694" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/RNA_273_SRR926694.sra"],
    "RNA_273_SRR926695" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/RNA_273_SRR926695.sra"],
    "RNA_273_SRR926696" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/RNA_273_SRR926696.sra"],
    "RNA_283_SRR926685" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/RNA_283_SRR926685.sra"],
    "RNA_283_SRR926686" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/RNA_283_SRR926686.sra"],
    "RNA_283_SRR926687" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/RNA_283_SRR926687.sra"],
    "RNA_283_SRR926688" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/RNA_283_SRR926688.sra"],
    "RNA_311_SRR926697" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/RNA_311_SRR926697.sra"],
    "RNA_311_SRR926698" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/RNA_311_SRR926698.sra"],
    "RNA_311_SRR926699" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/RNA_311_SRR926699.sra"],
    "RNA_311_SRR926700" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/RNA_311_SRR926700.sra"],
    "RNA_321_SRR926689" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/RNA_321_SRR926689.sra"],
    "RNA_321_SRR926690" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/RNA_321_SRR926690.sra"],
    "RNA_321_SRR926691" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/RNA_321_SRR926691.sra"],
    "RNA_321_SRR926692" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/RNA_321_SRR926692.sra"],
  },
  sra2fastq => {
    class      => "SRA::FastqDump",
    perform    => 0,
    ispaired   => 1,
    target_dir => "${target_dir}/FastqDump",
    option     => "",
    source_ref => "files",
    sh_direct  => 0,

    #cluster => "slurm",
    pbs => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "10",
      "mem"      => "10gb"
    },
  }
};

performConfig($config);

my $rna_config = {
  general => { task_name => "lichen" },

  files => {
    "RNA_273" => [
      "/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/FastqDump/result/RNA_273_1.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/FastqDump/result/RNA_273_2.fastq.gz"
    ],
    "RNA_283" => [
      "/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/FastqDump/result/RNA_283_1.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/FastqDump/result/RNA_283_2.fastq.gz"
    ],
    "RNA_311" => [
      "/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/FastqDump/result/RNA_311_1.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/FastqDump/result/RNA_311_2.fastq.gz"
    ],
    "RNA_321" => [
      "/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/FastqDump/result/RNA_321_1.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/FastqDump/result/RNA_321_2.fastq.gz"
    ],
  },
  fastqc => {
    class      => "QC::FastQC",
    perform    => 1,
    target_dir => "${target_dir}/fastqc",
    option     => "",
    source_ref => "files",
    sh_direct  => 0,
    cluster    => "slurm",
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=2",
      "walltime" => "2",
      "mem"      => "40gb"
    },
  },
  tophat2 => {
    class                => "Tophat2",
    perform              => 1,
    target_dir           => "${target_dir}/tophat2",
    option               => "--segment-length 25 -r 0 -p 8",
    source               => "files",
    bowtie2_index        => $bowtie2_index,
    transcript_gtf       => $transcript_gtf,
    transcript_gtf_index => $transcript_gtf_index,
    rename_bam           => 1,
    sh_direct            => 1,
    cluster              => "slurm",
    pbs                  => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },

};

performConfig($rna_config);

1
