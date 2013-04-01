#!/usr/bin/perl
use strict;
use warnings;

use CQS::QC;
use CQS::RNASeq;
use CQS::FileUtils;
use CQS::SystemUtils;

my $target_dir = create_directory_or_die("/scratch/cqs/shengq1/rnaseq/pierre");

my $transcript_gtf = "/data/cqs/guoy1/reference/annotation2/hg19/Homo_sapiens.GRCh37.68.gtf";

my $email = "quanhu.sheng\@vanderbilt.edu";

my $config = {
  general => {
    bowtie2_index        => "/data/cqs/guoy1/reference/hg19/bowtie2_index/hg19",
    transcript_gtf       => $transcript_gtf,
    transcript_gtf_index => "/scratch/cqs/shengq1/gtfindex/hg19_GRCh37_68",
    path_file            => "/home/shengq1/local/bin/path.txt",
    task_name            => "pierre"
  },
  fastqfiles => {
    "G1S1" => [ "/scratch/cqs/shengq1/rnaseq/pierre/raw/004/s_3_1_sequence.txt.gz", "/scratch/cqs/shengq1/rnaseq/pierre/raw/004/s_3_2_sequence.txt.gz" ],
    "G1S2" => [ "/scratch/cqs/shengq1/rnaseq/pierre/raw/003/s_2_1_sequence.txt.gz", "/scratch/cqs/shengq1/rnaseq/pierre/raw/003/s_2_2_sequence.txt.gz" ],
    "G2S1" => [ "/scratch/cqs/shengq1/rnaseq/pierre/raw/001/s_1_1_sequence.txt.gz", "/scratch/cqs/shengq1/rnaseq/pierre/raw/001/s_1_2_sequence.txt.gz" ],
    "G2S2" => [ "/scratch/cqs/shengq1/rnaseq/pierre/raw/001/s_5_1_sequence.txt.gz", "/scratch/cqs/shengq1/rnaseq/pierre/raw/001/s_5_2_sequence.txt.gz" ],
    "G2S3" => [ "/scratch/cqs/shengq1/rnaseq/pierre/raw/003/s_4_1_sequence.txt.gz", "/scratch/cqs/shengq1/rnaseq/pierre/raw/003/s_4_2_sequence.txt.gz" ],
    "G2S4" => [ "/scratch/cqs/shengq1/rnaseq/pierre/raw/001/s_2_1_sequence.txt.gz", "/scratch/cqs/shengq1/rnaseq/pierre/raw/001/s_2_2_sequence.txt.gz" ],
    "G3S1" => [ "/scratch/cqs/shengq1/rnaseq/pierre/raw/002/s_4_1_sequence.txt.gz", "/scratch/cqs/shengq1/rnaseq/pierre/raw/002/s_4_2_sequence.txt.gz" ],
    "G4S1" => [ "/scratch/cqs/shengq1/rnaseq/pierre/raw/003/s_1_1_sequence.txt.gz", "/scratch/cqs/shengq1/rnaseq/pierre/raw/003/s_1_2_sequence.txt.gz" ],
    "G4S2" => [ "/scratch/cqs/shengq1/rnaseq/pierre/raw/003/s_3_1_sequence.txt.gz", "/scratch/cqs/shengq1/rnaseq/pierre/raw/003/s_3_2_sequence.txt.gz" ],
    "G4S3" => [ "/scratch/cqs/shengq1/rnaseq/pierre/raw/001/s_4_1_sequence.txt.gz", "/scratch/cqs/shengq1/rnaseq/pierre/raw/001/s_4_2_sequence.txt.gz" ],
    "G4S4" => [ "/scratch/cqs/shengq1/rnaseq/pierre/raw/003/s_5_1_sequence.txt.gz", "/scratch/cqs/shengq1/rnaseq/pierre/raw/003/s_5_2_sequence.txt.gz" ],
    "G4S5" => [ "/scratch/cqs/shengq1/rnaseq/pierre/raw/004/s_2_1_sequence.txt.gz", "/scratch/cqs/shengq1/rnaseq/pierre/raw/004/s_2_2_sequence.txt.gz" ],
    "G5S1" => [ "/scratch/cqs/shengq1/rnaseq/pierre/raw/001/s_3_1_sequence.txt.gz", "/scratch/cqs/shengq1/rnaseq/pierre/raw/001/s_3_2_sequence.txt.gz" ],
    "G5S2" => [ "/scratch/cqs/shengq1/rnaseq/pierre/raw/002/s_5_1_sequence.txt.gz", "/scratch/cqs/shengq1/rnaseq/pierre/raw/002/s_5_2_sequence.txt.gz" ],
    "G5S3" => [ "/scratch/cqs/shengq1/rnaseq/pierre/raw/002/s_1_1_sequence.txt.gz", "/scratch/cqs/shengq1/rnaseq/pierre/raw/002/s_1_2_sequence.txt.gz" ],
    "G5S4" => [ "/scratch/cqs/shengq1/rnaseq/pierre/raw/004/s_4_1_sequence.txt.gz", "/scratch/cqs/shengq1/rnaseq/pierre/raw/004/s_4_2_sequence.txt.gz" ],
    "G6S1" => [ "/scratch/cqs/shengq1/rnaseq/pierre/raw/002/s_2_1_sequence.txt.gz", "/scratch/cqs/shengq1/rnaseq/pierre/raw/002/s_2_2_sequence.txt.gz" ],
    "G6S2" => [ "/scratch/cqs/shengq1/rnaseq/pierre/raw/002/s_3_1_sequence.txt.gz", "/scratch/cqs/shengq1/rnaseq/pierre/raw/002/s_3_2_sequence.txt.gz" ],
    "G6S3" => [ "/scratch/cqs/shengq1/rnaseq/pierre/raw/004/s_5_1_sequence.txt.gz", "/scratch/cqs/shengq1/rnaseq/pierre/raw/004/s_5_2_sequence.txt.gz" ],
  },
  groups => {
    "NOCANCER"  => [ "G1S1", "G1S2", "G2S1", "G2S2", "G2S3", "G2S4", "G3S1" ],
    "CANCER"    => [ "G4S1", "G4S2", "G4S3", "G4S4", "G4S5", "G5S1", "G5S2", "G5S3", "G5S4", "G6S1", "G6S2", "G6S3" ],
    "NONSMOKER" => [ "G1S1", "G1S2", "G4S1", "G4S2", "G4S3", "G4S4", "G4S5" ],
    "SMOKER"             => [ "G2S1", "G2S2", "G2S3", "G2S4", "G3S1", "G5S2", "G5S3", "G5S4", "G6S1", "G6S2", "G6S3" ],
    "NONSMOKER_NOCANCER" => [ "G1S1", "G1S2" ],
    "NONSMOKER_CANCER" => [ "G4S1", "G4S2", "G4S3", "G4S4", "G4S5" ],
    "LOWRISK_NOCANCER" => [ "G2S1", "G2S2", "G2S3", "G2S4" ],
    "LOWRISK_CANCER"   => [ "G5S1", "G5S2", "G5S3", "G5S4" ],
    "HIGHRISK_CANCER"  => [ "G6S1", "G6S2", "G6S3" ],
  },
  pairs => {
    "NOCANCER_vs_CANCER"                     => [ "NOCANCER",           "CANCER" ],
    "NONSMOKER_vs_SMOKER"                    => [ "NONSMOKER",          "SMOKER" ],
    "NONSMOKER_NOCANCER_vs_NONSMOKER_CANCER" => [ "NONSMOKER_NOCANCER", "NONSMOKER_CANCER" ],
    "LOWRISK_NOCANCER_vs_LOWRISK_CANCER"     => [ "LOWRISK_NOCANCER",   "LOWRISK_CANCER" ],
    "LOWRISK_CANCER_vs_HIGHRISK_CANCER"      => [ "LOWRISK_CANCER",     "HIGHRISK_CANCER" ]
  },
  fastqc => {
    target_dir => "${target_dir}/fastqc",
    option     => "",
    source_ref => "fastqfiles",
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=2",
      "walltime" => "2",
      "mem"      => "10gb"
    },
  },
  tophat2 => {
    target_dir => "${target_dir}/tophat2",
    option     => "-p 8",
    batchmode  => 0,
    source_ref => "fastqfiles",
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "240",
      "mem"      => "40gb"
    },
  },
  cufflinks => {
    target_dir     => "${target_dir}/cufflinks",
    option         => "-p 8 -u -N",
    transcript_gtf => $transcript_gtf,
    source_ref     => "tophat2",
    pbs            => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "10gb"
    },
  },
  cuffmerge => {
    target_dir => "${target_dir}/cuffmerge",
    option     => "-p 8",
    source_ref => "cufflinks",
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  cufflinks_cuffdiff => {
    target_dir         => "${target_dir}/cufflinks_cuffdiff",
    option             => "-p 8 -u -N",
    transcript_gtf_ref => "cuffmerge",
    source_ref         => "tophat2",
    groups_ref         => "groups",
    pairs_ref          => "pairs",
    pbs                => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "720",
      "mem"      => "40gb"
    },
  },
};

fastqc_by_pbs( $config, "fastqc" );

tophat2_by_pbs( $config, "tophat2" );

cufflinks_by_pbs( $config, "cufflinks" );

cuffmerge_by_pbs( $config, "cuffmerge" );

cuffdiff_by_pbs( $config, "cufflinks_cuffdiff" );

1;
