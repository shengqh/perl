#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;
use CQS::FileUtils;

my $target_dir = create_directory_or_die("/scratch/cqs/shengq1/rnaseq/20140204_yuhui_rnaseq_count");

my $fasta_file           = "/scratch/yuh9/software/bowtie2-2.1.0/index/hg19.fa";
my $transcript_gtf       = "/scratch/yuh9/software/bowtie2-2.1.0/index/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf";

my $hg19_map = "/scratch/cqs/shengq1/references/hg19/Homo_sapiens.GRCh37.73.map";

my $cqstools = "/home/shengq1/cqstools/CQS.Tools.exe";

my $email = "quanhu.sheng\@vanderbilt.edu";
my $task  = "20140204_yuhui";

#cqstools file_def -i /data/lehmanbd/insight_seq/TNBC_RNA-seq_phaseI/Run_3 -n \(.+_\)L001_\(.+\)_001
my $files = {
  "mRNA_62812_0_old" => ["/scratch/yuh9/HCC.RNAseq/allBatches.mRNA/oldBAM/mRNA_62812_0.bam"],
  "mRNA_120512_0_old" => ["/scratch/yuh9/HCC.RNAseq/allBatches.mRNA/oldBAM/mRNA_120512_0.bam"],
  "mRNA_070213_0_old" => ["/scratch/yuh9/HCC.RNAseq/allBatches.mRNA/oldBAM/mRNA_070213_0.bam"],
  "mRNA_62812_0_new" => ["/scratch/yuh9/HCC.RNAseq/allBatches.mRNA/BAM/mRNA_62812_0.bam"],
  "30-mRNA_120512_0_new" => ["/scratch/yuh9/HCC.RNAseq/allBatches.mRNA/BAM/mRNA_120512_0.bam"],
  "30-mRNA_070213_0_new" => ["/scratch/yuh9/HCC.RNAseq/allBatches.mRNA/BAM/mRNA_070213_0.bam"],
};

my $config = {
  general => { task_name => $task },
  files   => $files,
  sortbam => {
    class         => "Sortbam",
    perform       => 1,
    target_dir    => "${target_dir}/sortname",
    option        => "",
    source_ref    => "files",
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
};

performConfig($config);

1;
