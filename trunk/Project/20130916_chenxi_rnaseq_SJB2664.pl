#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;
use CQS::FileUtils;

my $target_dir     = create_directory_or_die("/scratch/cqs/shengq1/rnaseq/20130916_chenxi_rnaseq_SJB2664");
my $transcript_gtf = "/data/cqs/guoy1/reference/annotation2/mm10/Mus_musculus.GRCm38.68.gtf";
my $cqstools       = "/home/shengq1/cqstools/CQS.Tools.exe";
my $fasta_file     = "/data/cqs/guoy1/reference/mm10/bowtie2_index/mm10.fa";

my $email = "quanhu.sheng\@vanderbilt.edu";
my $task  = "smad4";

my $config = {
  general    => { task_name => $task },
  fastqfiles => {
    "2664-SJB-1" => [ "/data/cqs/chenx/project/ELF2/2664-SJB-1_1_sequence.txt.gz", "/data/cqs/chenx/project/ELF2/2664-SJB-1_2_sequence.txt.gz" ],
    "2664-SJB-2" => [ "/data/cqs/chenx/project/ELF2/2664-SJB-2_1_sequence.txt.gz", "/data/cqs/chenx/project/ELF2/2664-SJB-2_2_sequence.txt.gz" ],
  },
  groups => {
    "2664-SJB-1" => ["2664-SJB-1"],
    "2664-SJB-2" => ["2664-SJB-2"],
  },
  pairs  => { "2664-SJB" => [ "2664-SJB-1", "2664-SJB-2" ] },
  fastqc => {
    class      => "FastQC",
    perform    => 0,
    target_dir => "${target_dir}/fastqc",
    option     => "",
    source_ref => "fastqfiles",
    sh_direct  => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=2",
      "walltime" => "2",
      "mem"      => "10gb"
    },
  },
  tophat2 => {
    class         => "Tophat2",
    perform       => 1,
    target_dir    => "${target_dir}/tophat2",
    option        => "--segment-length 25 -r 0 -p 6",
    source_ref    => "fastqfiles",
    bowtie2_index => "/data/cqs/guoy1/reference/mm10/bowtie2_index/mm10",
    sh_direct     => 0,
    pbs           => {
      "email"    => $email,
      "nodes"    => "1:ppn=6",
      "walltime" => "72",
      "mem"      => "30gb"
    },
  },
  cuffdiff => {
    class          => "Cuffdiff",
    perform        => 1,
    target_dir     => "${target_dir}/cuffdiff",
    option         => "-p 8 -u -N",
    transcript_gtf => $transcript_gtf,
    source_ref     => "tophat2",
    groups_ref     => "groups",
    pairs_ref      => "pairs",
    fasta_file     => $fasta_file,
    sh_direct      => 0,
    pbs            => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "720",
      "mem"      => "40gb"
    },
  },
  sortbam => {
    class         => "Sortbam",
    perform       => 1,
    target_dir    => "${target_dir}/sortname",
    option        => "",
    source_ref    => "tophat2",
    sort_by_query => 1,
    sh_direct     => 1,
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
    sh_direct  => 1,
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
    name_map_file => "/data/cqs/shengq1/reference/mm10/mm10.gene.map",
    cqs_tools     => $cqstools,
    sh_direct     => 1,
    pbs           => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "10",
      "mem"      => "10gb"
    },
  },
  deseq2 => {
    class         => "DESeq2",
    perform       => 1,
    target_dir    => "${target_dir}/deseq2",
    option        => "",
    source_ref    => "pairs",
    groups_ref    => "groups",
    countfile_ref => "genetable",
    sh_direct     => 1,
    pbs           => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "10",
      "mem"      => "10gb"
    },
  },
};

performConfig($config);

1;
