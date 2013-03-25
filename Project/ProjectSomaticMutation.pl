#!/usr/bin/perl
use strict;
use warnings;

use CQS::QC;
use CQS::RNASeq;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::SomaticMutation;

my $target_dir = "/scratch/cqs/shengq1/somaticmutation";

my $email = "quanhu.sheng\@vanderbilt.edu";

my $transcript_gtf = "/data/cqs/guoy1/reference/annotation2/hg19/Homo_sapiens.GRCh37.68.gtf";

my $config = {
  general => {
    bowtie2_index        => "/data/cqs/guoy1/reference/hg19/bowtie2_index/hg19",
    transcript_gtf       => $transcript_gtf,
    transcript_gtf_index => "/scratch/cqs/shengq1/gtfindex/hg19_GRCh37_68",
    path_file            => "/home/shengq1/bin/path.txt",
    task_name            => "somaticmutation"
  },
  fastqfiles => {
    "TCGA-A7-A0D9-NT" => ["/scratch/cqs/shengq1/somaticmutation/raw/TCGA-A7-A0D9-RNA_NT_sorted.fastq"],
    "TCGA-A7-A0D9-TP" => ["/scratch/cqs/shengq1/somaticmutation/raw/TCGA-A7-A0D9-RNA_TP_sorted.fastq"],
  },
  bam_single => {
    "TCGA-A7-A0D9-TOPHAT2-SINGLE" =>
      [ "/scratch/cqs/shengq1/somaticmutation/tophat2/result/TCGA-A7-A0D9-NT/accepted_hits.bam", "/scratch/cqs/shengq1/somaticmutation/tophat2/result/TCGA-A7-A0D9-TP/accepted_hits.bam" ],
    "TCGA-A7-A0D9-SINGLE" => [ "/scratch/cqs/shengq1/somaticmutation/raw/TCGA-A7-A0D9-RNA_NT_sorted.bam", "/scratch/cqs/shengq1/somaticmutation/raw/TCGA-A7-A0D9-RNA_TP_sorted.bam" ]
  },
  bamlocal => {
    "TCGA-A7-A0D9-TOPHAT2-THREAD" =>
      [ "/scratch/cqs/shengq1/somaticmutation/tophat2/result/TCGA-A7-A0D9-NT/accepted_hits.bam", "/scratch/cqs/shengq1/somaticmutation/tophat2/result/TCGA-A7-A0D9-TP/accepted_hits.bam", ],
    "TCGA-A7-A0D9-THREAD" => [ "/scratch/cqs/shengq1/somaticmutation/raw/TCGA-A7-A0D9-RNA_NT_sorted.bam", "/scratch/cqs/shengq1/somaticmutation/raw/TCGA-A7-A0D9-RNA_TP_sorted.bam", ],

    #    "TCGA-BH-A0B3" => [
    #      "/scratch/cqs/shengq1/somaticmutation/raw/TCGA-BH-A0B3-RNA_NT_sorted.bam"
    #      "/scratch/cqs/shengq1/somaticmutation/raw/TCGA-BH-A0B3-RNA_TP_sorted.bam",
    #    ],
    #    "TCGA-BH-A0B8" => [
    #      "/scratch/cqs/shengq1/somaticmutation/raw/TCGA-BH-A0B8-RNA_NT_sorted.bam"
    #      "/scratch/cqs/shengq1/somaticmutation/raw/TCGA-BH-A0B8-RNA_TP_sorted.bam",
    #    ],
    #    "TCGA-BH-A0BJ" => [
    #      "/scratch/cqs/shengq1/somaticmutation/raw/TCGA-BH-A0BJ-RNA_NT_sorted.bam"
    #      "/scratch/cqs/shengq1/somaticmutation/raw/TCGA-BH-A0BJ-RNA_TP_sorted.bam",
    #    ],
    #    "TCGA-BH-A0BM" => [
    #      "/scratch/cqs/shengq1/somaticmutation/raw/TCGA-BH-A0BM-RNA_NT_sorted.bam"
    #      "/scratch/cqs/shengq1/somaticmutation/raw/TCGA-BH-A0BM-RNA_TP_sorted.bam",
    #    ],
    #    "TCGA-BH-A0C0" => [
    #      "/scratch/cqs/shengq1/somaticmutation/raw/TCGA-BH-A0C0-RNA_NT_sorted.bam"
    #      "/scratch/cqs/shengq1/somaticmutation/raw/TCGA-BH-A0C0-RNA_TP_sorted.bam",
    #    ],
    #    "TCGA-BH-A0DK" => [
    #      "/scratch/cqs/shengq1/somaticmutation/raw/TCGA-BH-A0DK-RNA_NT_sorted.bam"
    #      "/scratch/cqs/shengq1/somaticmutation/raw/TCGA-BH-A0DK-RNA_TP_sorted.bam",
    #    ],
    #    "TCGA-BH-A0DP" => [
    #      "/scratch/cqs/shengq1/somaticmutation/raw/TCGA-BH-A0DP-RNA_NT_sorted.bam"
    #      "/scratch/cqs/shengq1/somaticmutation/raw/TCGA-BH-A0DP-RNA_TP_sorted.bam",
    #    ],
    #    "TCGA-BH-A0E0" => [
    #      "/scratch/cqs/shengq1/somaticmutation/raw/TCGA-BH-A0E0-RNA_NT_sorted.bam"
    #      "/scratch/cqs/shengq1/somaticmutation/raw/TCGA-BH-A0E0-RNA_TP_sorted.bam",
    #    ],
    #    "TCGA-BH-A0H7" => [
    #      "/scratch/cqs/shengq1/somaticmutation/raw/TCGA-BH-A0H7-RNA_NT_sorted.bam"
    #      "/scratch/cqs/shengq1/somaticmutation/raw/TCGA-BH-A0H7-RNA_TP_sorted.bam",
    #    ],
  },
  bamfiles => {
    "TCGA-A7-A0D9" => [
      "/workspace/guoy1/GeneTorrent/lij17/processed/BRCA/TCGA-A7-A0D9/TCGA-A7-A0D9-RNA_NT_sorted.bam",
      "/workspace/guoy1/GeneTorrent/lij17/processed/BRCA/TCGA-A7-A0D9/TCGA-A7-A0D9-RNA_TP_sorted.bam",
    ],
    "TCGA-BH-A0B3" => [
      "/workspace/guoy1/GeneTorrent/lij17/processed/BRCA/TCGA-BH-A0B3/TCGA-BH-A0B3-RNA_NT_sorted.bam",
      "/workspace/guoy1/GeneTorrent/lij17/processed/BRCA/TCGA-BH-A0B3/TCGA-BH-A0B3-RNA_TP_sorted.bam",
    ],
    "TCGA-BH-A0B8" => [
      "/workspace/guoy1/GeneTorrent/lij17/processed/BRCA/TCGA-BH-A0B8/TCGA-BH-A0B8-RNA_NT_sorted.bam",
      "/workspace/guoy1/GeneTorrent/lij17/processed/BRCA/TCGA-BH-A0B8/TCGA-BH-A0B8-RNA_TP_sorted.bam",
    ],
    "TCGA-BH-A0BJ" => [
      "/workspace/guoy1/GeneTorrent/lij17/processed/BRCA/TCGA-BH-A0BJ/TCGA-BH-A0BJ-RNA_NT_sorted.bam",
      "/workspace/guoy1/GeneTorrent/lij17/processed/BRCA/TCGA-BH-A0BJ/TCGA-BH-A0BJ-RNA_TP_sorted.bam",
    ],
    "TCGA-BH-A0BM" => [
      "/workspace/guoy1/GeneTorrent/lij17/processed/BRCA/TCGA-BH-A0BM/TCGA-BH-A0BM-RNA_NT_sorted.bam",
      "/workspace/guoy1/GeneTorrent/lij17/processed/BRCA/TCGA-BH-A0BM/TCGA-BH-A0BM-RNA_TP_sorted.bam",
    ],
    "TCGA-BH-A0C0" => [
      "/workspace/guoy1/GeneTorrent/lij17/processed/BRCA/TCGA-BH-A0C0/TCGA-BH-A0C0-RNA_NT_sorted.bam",
      "/workspace/guoy1/GeneTorrent/lij17/processed/BRCA/TCGA-BH-A0C0/TCGA-BH-A0C0-RNA_TP_sorted.bam",
    ],
    "TCGA-BH-A0DK" => [
      "/workspace/guoy1/GeneTorrent/lij17/processed/BRCA/TCGA-BH-A0DK/TCGA-BH-A0DK-RNA_NT_sorted.bam",
      "/workspace/guoy1/GeneTorrent/lij17/processed/BRCA/TCGA-BH-A0DK/TCGA-BH-A0DK-RNA_TP_sorted.bam",
    ],
    "TCGA-BH-A0DP" => [
      "/workspace/guoy1/GeneTorrent/lij17/processed/BRCA/TCGA-BH-A0DP/TCGA-BH-A0DP-RNA_NT_sorted.bam",
      "/workspace/guoy1/GeneTorrent/lij17/processed/BRCA/TCGA-BH-A0DP/TCGA-BH-A0DP-RNA_TP_sorted.bam",
    ],
    "TCGA-BH-A0E0" => [
      "/workspace/guoy1/GeneTorrent/lij17/processed/BRCA/TCGA-BH-A0E0/TCGA-BH-A0E0-RNA_NT_sorted.bam",
      "/workspace/guoy1/GeneTorrent/lij17/processed/BRCA/TCGA-BH-A0E0/TCGA-BH-A0E0-RNA_TP_sorted.bam",
    ],
    "TCGA-BH-A0H7" => [
      "/workspace/guoy1/GeneTorrent/lij17/processed/BRCA/TCGA-BH-A0H7/TCGA-BH-A0H7-RNA_NT_sorted.bam",
      "/workspace/guoy1/GeneTorrent/lij17/processed/BRCA/TCGA-BH-A0H7/TCGA-BH-A0H7-RNA_TP_sorted.bam",
    ],
  },
  tophat2 => {
    target_dir => "${target_dir}/tophat2",
    option     => "--segment-length 25 -r 0 -p 8",
    batchmode  => 0,
    source_ref => "fastqfiles",
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "240",
      "mem"      => "40gb"
    },
  },
  wsmdetector => {
    target_dir                => "${target_dir}/wsmdetector",
    option                    => "-c 8",                                                   #thread mode
    source_ref                => "bamlocal",
    source_type               => "bam",                                                    #source_type can be bam/mpileup
    mpileup_sequence          => "/data/cqs/guoy1/reference/hg19/bowtie2_index/hg19.fa",
    annovar_buildver          => "hg19",
    annovar_database_location => "/data/cqs/guoy1/reference/annovar",

    #mpileup_option   => "-q 20",
    execute_file => "/home/shengq1/wsmdetector/wsmdetector.exe",
    pbs          => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "10gb"
    },
  },
  wsmdetector_single => {
    target_dir                => "${target_dir}/wsmdetector",
    option                    => "",                                                       #thread mode
    source_ref                => "bam_single",
    source_type               => "bam",                                                    #source_type can be bam/mpileup
    mpileup_sequence          => "/data/cqs/guoy1/reference/hg19/bowtie2_index/hg19.fa",
    annovar_buildver          => "hg19",
    annovar_database_location => "/data/cqs/guoy1/reference/annovar",
    #mpileup_option   => "-q 20",
    execute_file => "/home/shengq1/wsmdetector/wsmdetector.exe",
    pbs          => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "10gb"
    },
  },
};

#call_tophat2($config, "tophat2");

call_wsmdetector( $config, "wsmdetector" );

1;
