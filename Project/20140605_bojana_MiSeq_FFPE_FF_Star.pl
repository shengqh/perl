#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;

my $task = "20140605_bojana_FFPE_FF";

my $target_dir =
  "/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star";

my $email = "quanhu.sheng\@vanderbilt.edu";

my $config = {
  general => { task_name => $task },
  files =>{
    "IG-01" => "/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-1/Aligned.out.sorted.bam",
    "IG-02" => "/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-2/Aligned.out.sorted.bam",
    "IG-03" => "/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-3/Aligned.out.sorted.bam",
    "IG-04" => "/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-4/Aligned.out.sorted.bam",
    "IG-05" => "/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-5/Aligned.out.sorted.bam",
    "IG-06" => "/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-6/Aligned.out.sorted.bam",
    "IG-07" => "/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-7/Aligned.out.sorted.bam",
    "IG-08" => "/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-8/Aligned.out.sorted.bam",
    "IG-09" => "/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-9/Aligned.out.sorted.bam",
    "IG-10" => "/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-10/Aligned.out.sorted.bam",
    "IG-11" => "/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-11/Aligned.out.sorted.bam",
    "IG-12" => "/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-12/Aligned.out.sorted.bam",
    "IG-13" => "/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-13/Aligned.out.sorted.bam",
    "IG-14" => "/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-14/Aligned.out.sorted.bam",
    "IG-15" => "/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-15/Aligned.out.sorted.bam",
    "IG-16" => "/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-16/Aligned.out.sorted.bam",
    "IG-17" => "/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-17/Aligned.out.sorted.bam",
    "IG-18" => "/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-18/Aligned.out.sorted.bam",
    "IG-19" => "/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-19/Aligned.out.sorted.bam",
    "IG-20" => "/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-20/Aligned.out.sorted.bam",
    "IG-21" => "/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-21/Aligned.out.sorted.bam",
    "IG-22" => "/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-22/Aligned.out.sorted.bam",
    "IG-23" => "/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-23/Aligned.out.sorted.bam",
    "IG-24" => "/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-24/Aligned.out.sorted.bam",
    "IG-25" => "/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-25/Aligned.out.sorted.bam",
    "IG-26" => "/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-26/Aligned.out.sorted.bam",
    "IG-27" => "/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-27/Aligned.out.sorted.bam",
    "IG-28" => "/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-28/Aligned.out.sorted.bam",
    "IG-29" => "/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-29/Aligned.out.sorted.bam",
    "IG-33" => "/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-33/Aligned.out.sorted.bam",
    "IG-34" => "/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-34/Aligned.out.sorted.bam",
    "IG-39" => "/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-39/Aligned.out.sorted.bam",
    "IG-40" => "/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-40/Aligned.out.sorted.bam",
    "IG-41" => "/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-41/Aligned.out.sorted.bam",
    "IG-42" => "/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-42/Aligned.out.sorted.bam",
    "IG-43" => "/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-43/Aligned.out.sorted.bam",
    "IG-44" => "/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-44/Aligned.out.sorted.bam",
    "IG-45" => "/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-45/Aligned.out.sorted.bam",
    "IG-46" => "/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-46/Aligned.out.sorted.bam",
    "IG-47" => "/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-47/Aligned.out.sorted.bam",
    "IG-49" => "/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-49/Aligned.out.sorted.bam",
    "IG-50" => "/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-50/Aligned.out.sorted.bam",
    "IG-51" => "/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-51/Aligned.out.sorted.bam",
    "IG-52" => "/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-52/Aligned.out.sorted.bam",
    "IG-53" => "/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-53/Aligned.out.sorted.bam",
    "IG-54" => "/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-54/Aligned.out.sorted.bam",
    "IG-55" => "/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-55/Aligned.out.sorted.bam",
    "IG-56" => "/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-56/Aligned.out.sorted.bam",
    "IG-57" => "/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-57/Aligned.out.sorted.bam",
    "IG-58" => "/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-58/Aligned.out.sorted.bam",
    "IG-59" => "/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-59/Aligned.out.sorted.bam",
    "IG-60" => "/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-60/Aligned.out.sorted.bam",
    "IG-61" => "/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-61/Aligned.out.sorted.bam",
  },
};

performConfig($config);

1;
