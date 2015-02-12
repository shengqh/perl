#!/usr/bin/perl
use strict;
use warnings;

use CQS::SmallRNA;

my $def = {

  #General options
  task_name  => "celegans",
  email      => "quanhu.sheng\@vanderbilt.edu",
  target_dir => "/scratch/cqs/shengq1/smallRNA/Zhangp2_20150206_smallRNA_celegans/",
  max_thread => 8,

  #Default software parameter (don't change it except you really know it)
  fastq_remove_N             => 0,
  adapter                    => "CTGTAGGCACCATCAATC",

  #Data
  files => {
    "micro_rep1_daf16"  => ["/gpfs21/scratch/cqs/zhangp2/rawdata/sequencing/microRNA/alldata/micro_rep1_daf16.fq"],
    "micro_rep1_daf2"   => ["/gpfs21/scratch/cqs/zhangp2/rawdata/sequencing/microRNA/alldata/micro_rep1_daf2.fq"],
    "micro_rep1_double" => ["/gpfs21/scratch/cqs/zhangp2/rawdata/sequencing/microRNA/alldata/micro_rep1_double.fq"],
    "micro_rep1_N2"     => ["/gpfs21/scratch/cqs/zhangp2/rawdata/sequencing/microRNA/alldata/micro_rep1_N2.fq"],
    "micro_rep2_daf16"  => ["/gpfs21/scratch/cqs/zhangp2/rawdata/sequencing/microRNA/alldata/micro_rep2_daf16.fq"],
    "micro_rep2_daf2"   => ["/gpfs21/scratch/cqs/zhangp2/rawdata/sequencing/microRNA/alldata/micro_rep2_daf2.fq"],
    "micro_rep2_double" => ["/gpfs21/scratch/cqs/zhangp2/rawdata/sequencing/microRNA/alldata/micro_rep2_double.fq"],
    "micro_rep2_N2"     => ["/gpfs21/scratch/cqs/zhangp2/rawdata/sequencing/microRNA/alldata/micro_rep2_N2.fq"],
    "micro_rep3_daf16"  => ["/gpfs21/scratch/cqs/zhangp2/rawdata/sequencing/microRNA/alldata/micro_rep3_daf16.fq"],
    "micro_rep3_daf2"   => ["/gpfs21/scratch/cqs/zhangp2/rawdata/sequencing/microRNA/alldata/micro_rep3_daf2.fq"],
    "micro_rep3_double" => ["/gpfs21/scratch/cqs/zhangp2/rawdata/sequencing/microRNA/alldata/micro_rep3_double.fq"],
    "micro_rep3_N2"     => ["/gpfs21/scratch/cqs/zhangp2/rawdata/sequencing/microRNA/alldata/micro_rep3_N2.fq"],
  }
};

performSmallRNA_cel235($def);

1;

