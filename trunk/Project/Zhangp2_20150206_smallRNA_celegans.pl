#!/usr/bin/perl
use strict;
use warnings;

use Pipeline::SmallRNA;

my $def = {

  #General options
  task_name  => "celegans",
  email      => "quanhu.sheng\@vanderbilt.edu",
  target_dir => "/scratch/cqs/shengq1/smallRNA/Zhangp2_20150206_smallRNA_celegans/",
  max_thread => 8,

  #Default software parameter (don't change it except you really know it)
  bowtie1_option_1mm         => "-a -m 100 --best --strata -v 1 -p 8",
  bowtie1_option_pm          => "-a -m 100 --best --strata -v 0 -p 8",
  mirnacount_option          => "-s",                                         #ignore score
  smallrnacount_option       => "-s --min_overlap 0.5 --length --sequence",
  mirna_overlap_count_option => "-s --min_overlap 0.5 --gtf_key miRNA",
  min_read_length            => 16,

  #Software and miRBase database options
  samtools              => "/home/shengq1/local/bin/samtools/samtools",
  cqstools              => "/home/shengq1/cqstools/CQS.Tools.exe",
  mirna_fasta           => "/data/cqs/shengq1/reference/miRBase21/mature.dna.fa",
  trna_coordinate       => "/scratch/cqs/shengq1/references/wbcel235/WBcel235_tRNA.bed",
  trna_fasta            => "/scratch/cqs/shengq1/references/wbcel235/WBcel235_tRNA.bed.fa",
  smallrna_coordinate   => "/scratch/cqs/shengq1/references/wbcel235/WBcel235_otherSmallRNA.bed",
  bowtie1_miRBase_index => "/data/cqs/shengq1/reference/miRBase21/bowtie_index_1.0.1/mature.dna",
  mirbase_count_option  => "-p cel",

  #genome database
  mirna_coordinate => "/data/cqs/shengq1/reference/miRBase21/cel.gff3",
  bowtie1_index    => "/scratch/cqs/zhangp2/reference/wormbase/bowtie_index_1.1.0/Caenorhabditis_elegans.WBcel235.dna.toplevel",

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

performSmallRNA($def);

1;

