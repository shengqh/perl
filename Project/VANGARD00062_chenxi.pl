#!/usr/bin/perl
use strict;
use warnings;

use CQS::QC;
use CQS::RNASeq;
use CQS::FileUtils;
use CQS::SystemUtils;

my $target_dir = create_directory_or_die("/scratch/cqs/shengq1/rnaseq/VANGARD00062_chenxi");

my $transcript_gtf       = "/data/cqs/guoy1/reference/annotation2/mm10/Mus_musculus.GRCm38.68_chr1-22-X-Y-M.gtf";
my $transcript_gtf_index = "/scratch/cqs/shengq1/gtfindex/mm10_GRCm38.68";

my $email = "quanhu.sheng\@vanderbilt.edu";

my $config = {
  general => {
    bowtie2_index        => "/data/cqs/guoy1/reference/mm10/mm10",
    transcript_gtf       => $transcript_gtf,
    transcript_gtf_index => $transcript_gtf_index,
    path_file            => "/home/shengq1/local/bin/path.txt",
    task_name            => "VANGARD00062"
  },
  fastqfiles => {
    "WT-1" => ["/data/cqs/chenx/project/sepp1/data/2443-CSW-1_1_sequence.txt"],
    "WT-2" => ["/data/cqs/chenx/project/sepp1/data/2443-CSW-2_1_sequence.txt"],
    "WT-3" => ["/data/cqs/chenx/project/sepp1/data/2443-CSW-3_1_sequence.txt"],
    "HET-1" => ["/data/cqs/chenx/project/sepp1/data/2443-CSW-4_1_sequence.txt"],
    "HET-2" => ["/data/cqs/chenx/project/sepp1/data/2443-CSW-5_1_sequence.txt"],
    "HET-3" => ["/data/cqs/chenx/project/sepp1/data/2443-CSW-6_1_sequence.txt"],
    "NULL-1" => ["/data/cqs/chenx/project/sepp1/data/2443-CSW-7_1_sequence.txt"],
    "NULL-2" => ["/data/cqs/chenx/project/sepp1/data/2443-CSW-8_1_sequence.txt"],
    "NULL-3" => ["/data/cqs/chenx/project/sepp1/data/2443-CSW-9_1_sequence.txt"]
  },
  groups => {
    "WT" => [ "WT-1", "WT-2", "WT-3" ],
    "HET" => [ "HET-1", "HET-2", "HET-3" ],
    "NULL" => [ "NULL-1", "NULL-2", "NULL-3" ],
  },
  pairs  => {
     "HET_vs_WT" => [ "HET", "WT" ], 
     "HET_vs_NULL" => [ "HET", "NULL" ], 
     "WT_vs_NULL" => [ "WT", "NULL" ], 
  },
  fastqc => {
    target_dir => "${target_dir}/fastqc",
    option     => "",
    source_ref => "fastqfiles",
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "2",
      "mem"      => "10gb"
    },
  },
  tophat2 => {
    target_dir           => "${target_dir}/tophat2",
    option               => "-p 8",
    batchmode            => 0,
    sortbam              => 1,
    source_ref           => "fastqfiles",
    transcript_gtf       => $transcript_gtf,
    transcript_gtf_index => $transcript_gtf_index,
    pbs                  => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "240",
      "mem"      => "40gb"
    },
  },
  cuffdiff => {
    target_dir     => "${target_dir}/cuffdiff",
    option         => "-p 8 -u -N",
    transcript_gtf => $transcript_gtf,
    source_ref     => "tophat2",
    groups_ref     => "groups",
    pairs_ref      => "pairs",
    pbs            => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "720",
      "mem"      => "40gb"
    },
  },
  rename_diff => {
    target_dir => "${target_dir}/cuffdiff/result/comparison",
    root_dir   => "${target_dir}/cuffdiff/result",
    gene_only  => 1
  },
};

fastqc_by_pbs( $config, "fastqc" );

tophat2_by_pbs( $config, "tophat2" );

cuffdiff_by_pbs( $config, "cuffdiff" );

##copy_and_rename_cuffdiff_file($config, "rename_diff");

1;
