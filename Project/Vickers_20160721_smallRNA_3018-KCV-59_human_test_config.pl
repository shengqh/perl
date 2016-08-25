#!/usr/bin/perl
use strict;
use warnings;

use CQS::PerformSmallRNA;
use CQS::FileUtils;
use Pipeline::SmallRNAUtils;
use Pipeline::SmallRNA;

my $genome = hg19_genome();

my $root = create_directory_or_die("/scratch/cqs/shengq1/vickers/20160721_smallRNA_3018-KCV-59_human_test_config/");

my $userdef = {

  #General options
  email                     => "quanhu.sheng\@vanderbilt.edu",
  max_thread                => 8,
  cqstools                  => "/home/shengq1/cqstools/cqstools.exe",
  sequencetask_run_time     => 12,
  table_vis_group_text_size => 14,

  #Default software parameter (don't change it except you really know it)
  fastq_remove_N        => 0,
  remove_sequences      => "'CCACGTTCCCGTGG;ACAGTCCGACGATC'",
  search_unmapped_reads => 1,
  blast_unmapped_reads  => 1,

  #General options
  task_name  => "withcontig",
  target_dir => $root . "with_contig",

  #Data
  files => {
    "i1" => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-59/3018-KCV-59-i1_S1_R1_001.fastq.gz'],
    "i2" => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-59/3018-KCV-59-i2_S2_R1_001.fastq.gz'],
    "i3" => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-59/3018-KCV-59-i3_S3_R1_001.fastq.gz'],
    "i4" => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-59/3018-KCV-59-i4_S4_R1_001.fastq.gz'],
  },
  groups => {
    "G1" => [ "i1", "i2" ],
    "G2" => [ "i3", "i4" ]
  },
  pairs => {
    "COMP" => [ "G1", "G2" ]
  }
};

my $defWithContig = getSmallRNADefinition( $userdef, $genome );
performSmallRNA($defWithContig);

#$userdef->{target_dir}   = $root . "without_contig";
#$userdef->{task_name}    = "without_contig";
#$genome->{bowtie1_index} = "/scratch/cqs/shengq1/references/gencode/hg19_slim/bowtie_index_1.1.2/GRCh37.p13.genome.slim";
#
#my $defWithoutContig = getSmallRNADefinition( $userdef, $genome );
#performSmallRNA($defWithoutContig);

1;

