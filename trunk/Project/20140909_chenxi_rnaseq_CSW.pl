#!/usr/bin/perl
use strict;
use warnings;

use CQS::QC;
use CQS::RNASeq;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ClassFactory;

my $task_name  = "CSW";
my $target_dir = create_directory_or_die("/scratch/cqs/shengq1/chenxi/20140909_chenxi_rnaseq_CSW");

my $transcript_gtf       = "/data/cqs/shengq1/reference/ensembl_gtf/Homo_sapiens.GRCh37.75.M.gtf";
my $transcript_gtf_map   = "/data/cqs/shengq1/reference/ensembl_gtf/Homo_sapiens.GRCh37.75.map";
my $transcript_gtf_index = "/scratch/cqs/shengq1/gtfindex/hg19_GRCh37_75";
my $bowtie2_index        = "/data/cqs/shengq1/reference/hg19_16569_M/bowtie2_index_2.1.0/hg19_16569_M";
my $fasta_file           = "/data/cqs/shengq1/reference/hg19_16569_M/bowtie2_index_2.1.0/hg19_16569_M.fa";

my $email = "quanhu.sheng\@vanderbilt.edu";

my $config = {
  general    => { task_name => $task_name, },
  fastqfiles => {
    "3009-CSW-1" => ["/gpfs20/data/cqs/chenx/cw2014/3009-CSW-1_1_sequence.txt.gz"],
    "3009-CSW-2" => ["/gpfs20/data/cqs/chenx/cw2014/3009-CSW-2_1_sequence.txt.gz"],
    "3009-CSW-3" => ["/gpfs20/data/cqs/chenx/cw2014/3009-CSW-3_1_sequence.txt.gz"],
    "3009-CSW-4" => ["/gpfs20/data/cqs/chenx/cw2014/3009-CSW-4_1_sequence.txt.gz"],
    "3009-CSW-5" => ["/gpfs20/data/cqs/chenx/cw2014/3009-CSW-5_1_sequence.txt.gz"],
    "3009-CSW-6" => ["/gpfs20/data/cqs/chenx/cw2014/3009-CSW-6_1_sequence.txt.gz"],
    "3009-CSW-7" => ["/gpfs20/data/cqs/chenx/cw2014/3009-CSW-7_1_sequence.txt.gz"],
    "3009-CSW-8" => ["/gpfs20/data/cqs/chenx/cw2014/3009-CSW-8_1_sequence.txt.gz"],
  },
  groups => {
    "Nontarget_shRNA_control" => [ "3009-CSW-1", "3009-CSW-2" ],
    "STK17A_shRNA_construct"  => [ "3009-CSW-3", "3009-CSW-4" ],
    "Nontarget_siRNA_control" => [ "3009-CSW-5", "3009-CSW-6" ],
    "POP3_siRNA"              => [ "3009-CSW-7", "3009-CSW-8" ],
  },
  pairs => {
    "STK17A_shRNA_construct_VS_Nontarget_shRNA_control" => [ "Nontarget_shRNA_control", "STK17A_shRNA_construct" ],
    "POP3_siRNA_VS_Nontarget_siRNA_control"             => [ "Nontarget_siRNA_control", "POP3_siRNA" ]
  },
  fastqc => {
    class      => "FastQC",
    perform    => 1,
    target_dir => "${target_dir}/fastqc",
    option     => "",
    source_ref => "fastqfiles",
    sh_direct  => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "2",
      "mem"      => "10gb"
    },
  },
  tophat2 => {
    class                => "Tophat2",
    perform              => 1,
    target_dir           => "${target_dir}/tophat2",
    option               => "--segment-length 25 -r 0 -p 8",
    source_ref           => "fastqfiles",
    transcript_gtf       => $transcript_gtf,
    transcript_gtf_index => $transcript_gtf_index,
    bowtie2_index        => $bowtie2_index,
    pbs                  => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  rnaseqc => {
    class          => "RNASeQC",
    perform        => 1,
    target_dir     => "${target_dir}/rnaseqc",
    option         => "",
    source_ref     => "tophat2",
    jar            => "/home/shengq1/local/bin/RNA-SeQC_v1.1.7.jar",
    fasta_file     => $fasta_file,
    transcript_gtf => $transcript_gtf,
    pbs            => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "20gb"
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
    option        => "-p ENSG",
    source_ref    => "htseqcount",
    name_map_file => $transcript_gtf_map,
    cqs_tools     => "/home/shengq1/cqstools/CQS.Tools.exe",
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
  sequencetask => {
    class      => "SequenceTask",
    perform    => 1,
    target_dir => "${target_dir}/sequencetask",
    source     => {
      "sample" => [ "fastqc",  "tophat2",   "sortbam", "htseqcount" ],
      "task"   => [ "rnaseqc", "genetable", "deseq2" ],
    },
    sh_direct => 0,
    pbs       => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "30gb"
    }
    }

};

performConfig($config);

1;
