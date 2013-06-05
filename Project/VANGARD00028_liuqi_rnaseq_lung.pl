#!/usr/bin/perl
use strict;
use warnings;

use CQS::QC;
use CQS::RNASeq;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ClassFactory;

my $target_dir = create_directory_or_die("/scratch/cqs/shengq1/vangard/VANGARD00028_liuqi_rnaseq");

my $bowtie2_index = "/data/cqs/guoy1/reference/hg19/bowtie2_index/hg19";

my $transcript_gtf       = "/data/cqs/guoy1/reference/annotation2/hg19/Homo_sapiens.GRCh37.68.gtf";
my $transcript_gtf_index = "/scratch/cqs/shengq1/gtfindex/hg19_GRCh37_68";

my $email = "quanhu.sheng\@vanderbilt.edu";

my $config = {
  general => {
    bowtie2_index        => "/data/cqs/guoy1/reference/hg19/bowtie2_index/hg19",
    transcript_gtf       => $transcript_gtf,
    transcript_gtf_index => $transcript_gtf_index,
    path_file            => "/home/shengq1/local/bin/path.txt",
    task_name            => "VANGARD00028"
  },
  fastqfiles => {
    "G1_7071" => [ "${target_dir}/raw/004/s_3_1_sequence.txt.gz", "${target_dir}/raw/004/s_3_2_sequence.txt.gz" ],
    "G1_7143" => [ "${target_dir}/raw/003/s_2_1_sequence.txt.gz", "${target_dir}/raw/003/s_2_2_sequence.txt.gz" ],
    "G2_7030" => [ "${target_dir}/raw/001/s_1_1_sequence.txt.gz", "${target_dir}/raw/001/s_1_2_sequence.txt.gz" ],
    "G2_7178" => [ "${target_dir}/raw/001/s_5_1_sequence.txt.gz", "${target_dir}/raw/001/s_5_2_sequence.txt.gz" ],
    "G2_7222" => [ "${target_dir}/raw/003/s_4_1_sequence.txt.gz", "${target_dir}/raw/003/s_4_2_sequence.txt.gz" ],
    "G2_7228" => [ "${target_dir}/raw/001/s_2_1_sequence.txt.gz", "${target_dir}/raw/001/s_2_2_sequence.txt.gz" ],
    "G3_7089" => [ "${target_dir}/raw/002/s_4_1_sequence.txt.gz", "${target_dir}/raw/002/s_4_2_sequence.txt.gz" ],
    "G4_7448" => [ "${target_dir}/raw/003/s_1_1_sequence.txt.gz", "${target_dir}/raw/003/s_1_2_sequence.txt.gz" ],
    "G4_7485" => [ "${target_dir}/raw/003/s_3_1_sequence.txt.gz", "${target_dir}/raw/003/s_3_2_sequence.txt.gz" ],
    "G4_7522" => [ "${target_dir}/raw/001/s_4_1_sequence.txt.gz", "${target_dir}/raw/001/s_4_2_sequence.txt.gz" ],
    "G4_7617" => [ "${target_dir}/raw/003/s_5_1_sequence.txt.gz", "${target_dir}/raw/003/s_5_2_sequence.txt.gz" ],
    "G4_7697" => [ "${target_dir}/raw/004/s_2_1_sequence.txt.gz", "${target_dir}/raw/004/s_2_2_sequence.txt.gz" ],
    "G5_6820" => [ "${target_dir}/raw/001/s_3_1_sequence.txt.gz", "${target_dir}/raw/001/s_3_2_sequence.txt.gz" ],
    "G5_7080" => [ "${target_dir}/raw/002/s_5_1_sequence.txt.gz", "${target_dir}/raw/002/s_5_2_sequence.txt.gz" ],
    "G5_7176" => [ "${target_dir}/raw/002/s_1_1_sequence.txt.gz", "${target_dir}/raw/002/s_1_2_sequence.txt.gz" ],
    "G5_7182" => [ "${target_dir}/raw/004/s_4_1_sequence.txt.gz", "${target_dir}/raw/004/s_4_2_sequence.txt.gz" ],
    "G6_7134" => [ "${target_dir}/raw/002/s_2_1_sequence.txt.gz", "${target_dir}/raw/002/s_2_2_sequence.txt.gz" ],
    "G6_7053" => [ "${target_dir}/raw/002/s_3_1_sequence.txt.gz", "${target_dir}/raw/002/s_3_2_sequence.txt.gz" ],
    "G6_7116" => [ "${target_dir}/raw/004/s_5_1_sequence.txt.gz", "${target_dir}/raw/004/s_5_2_sequence.txt.gz" ],
  },
  groups => {
    "NOCANCER"  => [ "G1_7071", "G1_7143", "G2_7030", "G2_7178", "G2_7222", "G2_7228", "G3_7089", ],
    "CANCER"    => [ "G4_7448", "G4_7485", "G4_7522", "G4_7617", "G4_7697", "G5_6820", "G5_7080", "G5_7176", "G5_7182", "G6_7134", "G6_7053", "G6_7116", ],
    "NONSMOKER" => [ "G1_7071", "G1_7143", "G4_7448", "G4_7485", "G4_7522", "G4_7617", "G4_7697", ],
    "SMOKER"    => [ "G2_7030", "G2_7178", "G2_7222", "G2_7228", "G3_7089", "G5_6820", "G5_7080", "G5_7176", "G5_7182", "G6_7134", "G6_7053", "G6_7116", ],

    "G1" => [ "G1_7071", "G1_7143", ],                          #NONSMOKER_NOCANCER
    "G2" => [ "G2_7030", "G2_7178", "G2_7222", "G2_7228", ],    #LOWRISK_NOCANCER
    "G3" => [ "G3_7089", ],                                                #HIGHRISK_NOCANCER
    "G4" => [ "G4_7448", "G4_7485", "G4_7522", "G4_7617", "G4_7697", ],    #NONSMOKER_CANCER
    "G5" => [ "G5_6820", "G5_7080", "G5_7176", "G5_7182", ],               #LOWRISK_CANCER
    "G6" => [ "G6_7134", "G6_7053", "G6_7116", ],                          #HIGHRISK_CANCER
  },
  pairs => {

    #    "NOCANCER_vs_CANCER"                     => [ "NOCANCER",  "CANCER" ],
    #    "NONSMOKER_vs_SMOKER"                    => [ "NONSMOKER", "SMOKER" ],
    #    "NONSMOKER_NOCANCER_vs_NONSMOKER_CANCER" => [ "G1",        "G4" ],
    #    "LOWRISK_NOCANCER_vs_LOWRISK_CANCER"     => [ "G2",        "G5" ],
    #    "LOWRISK_CANCER_vs_HIGHRISK_CANCER"      => [ "G5",        "G6" ],
    #    "G1236"                                  => [ "G1",        "G2", "G3", "G6" ],
    "G13" => [ "G1", "G3", ]
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
    target_dir           => "${target_dir}/tophat2",
    option               => "-p 8",
    batchmode            => 0,
    sortbam              => 0,
    indexbam             => 1,
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
  tophat2_class => {
    class                => "Tophat2",
    bowtie2_index        => $bowtie2_index,
    target_dir           => "${target_dir}/tophat2_class",
    option               => "-p 8",
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
  rename_diff => {
    target_dir => "${target_dir}/cufflinks_cuffdiff/result/comparison",
    root_dir   => "${target_dir}/cufflinks_cuffdiff/result",
    gene_only  => 1
  },
};

#fastqc_by_pbs( $config, "fastqc" );

performTask( $config, "tophat2_class" );

#call_RNASeQC($config, "rnaseqc");

#cufflinks_by_pbs( $config, "cufflinks" );

#cuffmerge_by_pbs( $config, "cuffmerge" );

#cuffdiff_by_pbs( $config, "cufflinks_cuffdiff" );

#copy_and_rename_cuffdiff_file( $config, "rename_diff" );

1;
