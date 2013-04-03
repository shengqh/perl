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
    "G1_7071" => [ "/scratch/cqs/shengq1/rnaseq/pierre/raw/004/s_3_1_sequence.txt.gz", "/scratch/cqs/shengq1/rnaseq/pierre/raw/004/s_3_2_sequence.txt.gz" ],
    "G1_7143" => [ "/scratch/cqs/shengq1/rnaseq/pierre/raw/003/s_2_1_sequence.txt.gz", "/scratch/cqs/shengq1/rnaseq/pierre/raw/003/s_2_2_sequence.txt.gz" ],
    "G2_7030" => [ "/scratch/cqs/shengq1/rnaseq/pierre/raw/001/s_1_1_sequence.txt.gz", "/scratch/cqs/shengq1/rnaseq/pierre/raw/001/s_1_2_sequence.txt.gz" ],
    "G2_7178" => [ "/scratch/cqs/shengq1/rnaseq/pierre/raw/001/s_5_1_sequence.txt.gz", "/scratch/cqs/shengq1/rnaseq/pierre/raw/001/s_5_2_sequence.txt.gz" ],
    "G2_7222" => [ "/scratch/cqs/shengq1/rnaseq/pierre/raw/003/s_4_1_sequence.txt.gz", "/scratch/cqs/shengq1/rnaseq/pierre/raw/003/s_4_2_sequence.txt.gz" ],
    "G2_7228" => [ "/scratch/cqs/shengq1/rnaseq/pierre/raw/001/s_2_1_sequence.txt.gz", "/scratch/cqs/shengq1/rnaseq/pierre/raw/001/s_2_2_sequence.txt.gz" ],
    "G3_7089" => [ "/scratch/cqs/shengq1/rnaseq/pierre/raw/002/s_4_1_sequence.txt.gz", "/scratch/cqs/shengq1/rnaseq/pierre/raw/002/s_4_2_sequence.txt.gz" ],
    "G4_7448" => [ "/scratch/cqs/shengq1/rnaseq/pierre/raw/003/s_1_1_sequence.txt.gz", "/scratch/cqs/shengq1/rnaseq/pierre/raw/003/s_1_2_sequence.txt.gz" ],
    "G4_7485" => [ "/scratch/cqs/shengq1/rnaseq/pierre/raw/003/s_3_1_sequence.txt.gz", "/scratch/cqs/shengq1/rnaseq/pierre/raw/003/s_3_2_sequence.txt.gz" ],
    "G4_7522" => [ "/scratch/cqs/shengq1/rnaseq/pierre/raw/001/s_4_1_sequence.txt.gz", "/scratch/cqs/shengq1/rnaseq/pierre/raw/001/s_4_2_sequence.txt.gz" ],
    "G4_7617" => [ "/scratch/cqs/shengq1/rnaseq/pierre/raw/003/s_5_1_sequence.txt.gz", "/scratch/cqs/shengq1/rnaseq/pierre/raw/003/s_5_2_sequence.txt.gz" ],
    "G4_7697" => [ "/scratch/cqs/shengq1/rnaseq/pierre/raw/004/s_2_1_sequence.txt.gz", "/scratch/cqs/shengq1/rnaseq/pierre/raw/004/s_2_2_sequence.txt.gz" ],
    "G5_6820" => [ "/scratch/cqs/shengq1/rnaseq/pierre/raw/001/s_3_1_sequence.txt.gz", "/scratch/cqs/shengq1/rnaseq/pierre/raw/001/s_3_2_sequence.txt.gz" ],
    "G5_7080" => [ "/scratch/cqs/shengq1/rnaseq/pierre/raw/002/s_5_1_sequence.txt.gz", "/scratch/cqs/shengq1/rnaseq/pierre/raw/002/s_5_2_sequence.txt.gz" ],
    "G5_7176" => [ "/scratch/cqs/shengq1/rnaseq/pierre/raw/002/s_1_1_sequence.txt.gz", "/scratch/cqs/shengq1/rnaseq/pierre/raw/002/s_1_2_sequence.txt.gz" ],
    "G5_7182" => [ "/scratch/cqs/shengq1/rnaseq/pierre/raw/004/s_4_1_sequence.txt.gz", "/scratch/cqs/shengq1/rnaseq/pierre/raw/004/s_4_2_sequence.txt.gz" ],
    "G6_7134" => [ "/scratch/cqs/shengq1/rnaseq/pierre/raw/002/s_2_1_sequence.txt.gz", "/scratch/cqs/shengq1/rnaseq/pierre/raw/002/s_2_2_sequence.txt.gz" ],
    "G6_7053" => [ "/scratch/cqs/shengq1/rnaseq/pierre/raw/002/s_3_1_sequence.txt.gz", "/scratch/cqs/shengq1/rnaseq/pierre/raw/002/s_3_2_sequence.txt.gz" ],
    "G6_7116" => [ "/scratch/cqs/shengq1/rnaseq/pierre/raw/004/s_5_1_sequence.txt.gz", "/scratch/cqs/shengq1/rnaseq/pierre/raw/004/s_5_2_sequence.txt.gz" ],
  },
  groups => {
    "NOCANCER"  => [ "G1_7071", "G1_7143", "G2_7030", "G2_7178", "G2_7222", "G2_7228", "G3_7089", ],
    "CANCER"    => [ "G4_7448", "G4_7485", "G4_7522", "G4_7617", "G4_7697", "G5_6820", "G5_7080", "G5_7176", "G5_7182", "G6_7134", "G6_7053", "G6_7116", ],
    "NONSMOKER" => [ "G1_7071", "G1_7143", "G4_7448", "G4_7485", "G4_7522", "G4_7617", "G4_7697", ],
    "SMOKER"    => [ "G2_7030", "G2_7178", "G2_7222", "G2_7228", "G3_7089", "G5_6820", "G5_7080", "G5_7176", "G5_7182", "G6_7134", "G6_7053", "G6_7116", ],
    "NONSMOKER_NOCANCER" => [ "G1_7071", "G1_7143", ],
    "NONSMOKER_CANCER"   => [ "G4_7448", "G4_7485", "G4_7522", "G4_7617", "G4_7697", ],
    "LOWRISK_NOCANCER" => [ "G2_7030", "G2_7178", "G2_7222", "G2_7228", ],
    "LOWRISK_CANCER"   => [ "G5_6820", "G5_7080", "G5_7176", "G5_7182", ],
    "HIGHRISK_CANCER"  => [ "G6_7134", "G6_7053", "G6_7116", ],
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
    target_dir => "${target_dir}/tophat2_RG",
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
  rnaseqc => {
    target_dir     => "${target_dir}/RNASeQC",
    option         => "",
    transcript_gtf => $transcript_gtf,
    genome_fasta   => "/data/cqs/guoy1/reference/hg19/hg19_chr.fa",
    rnaseqc_jar    => "/home/shengq1/local/bin/RNA-SeQC_v1.1.7.jar",
    source_ref     => "tophat2",
    pbs            => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "10gb"
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

#fastqc_by_pbs( $config, "fastqc" );

tophat2_by_pbs( $config, "tophat2" );

#call_RNASeQC($config, "rnaseqc");

#cufflinks_by_pbs( $config, "cufflinks" );

#cuffmerge_by_pbs( $config, "cuffmerge" );

#cuffdiff_by_pbs( $config, "cufflinks_cuffdiff" );

1;
