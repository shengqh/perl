#!/usr/bin/perl

use strict;
use warnings;

my $files = {

  #cqstools file_def -i /data/lehmanbd/insight_seq/TNBC_RNA-seq_phaseI/Run_3 -n \(.+_\)L001_\(.+\)_001
  "run3" => {
    "30-PA_S1_R1"  => ["/data/lehmanbd/insight_seq/TNBC_RNA-seq_phaseI/Run_3/30-PA_S1_L001_R1_001.fastq.gz"],
    "30-PA_S1_R2"  => ["/data/lehmanbd/insight_seq/TNBC_RNA-seq_phaseI/Run_3/30-PA_S1_L001_R2_001.fastq.gz"],
    "30-SB_S2_R1"  => ["/data/lehmanbd/insight_seq/TNBC_RNA-seq_phaseI/Run_3/30-SB_S2_L001_R1_001.fastq.gz"],
    "30-SB_S2_R2"  => ["/data/lehmanbd/insight_seq/TNBC_RNA-seq_phaseI/Run_3/30-SB_S2_L001_R2_001.fastq.gz"],
    "30-TA_S3_R1"  => ["/data/lehmanbd/insight_seq/TNBC_RNA-seq_phaseI/Run_3/30-TA_S3_L001_R1_001.fastq.gz"],
    "30-TA_S3_R2"  => ["/data/lehmanbd/insight_seq/TNBC_RNA-seq_phaseI/Run_3/30-TA_S3_L001_R2_001.fastq.gz"],
    "40-PA_S4_R1"  => ["/data/lehmanbd/insight_seq/TNBC_RNA-seq_phaseI/Run_3/40-PA_S4_L001_R1_001.fastq.gz"],
    "40-PA_S4_R2"  => ["/data/lehmanbd/insight_seq/TNBC_RNA-seq_phaseI/Run_3/40-PA_S4_L001_R2_001.fastq.gz"],
    "40-TXA_S5_R1" => ["/data/lehmanbd/insight_seq/TNBC_RNA-seq_phaseI/Run_3/40-TXA_S5_L001_R1_001.fastq.gz"],
    "40-TXA_S5_R2" => ["/data/lehmanbd/insight_seq/TNBC_RNA-seq_phaseI/Run_3/40-TXA_S5_L001_R2_001.fastq.gz"],
    "42-PC_S6_R1"  => ["/data/lehmanbd/insight_seq/TNBC_RNA-seq_phaseI/Run_3/42-PC_S6_L001_R1_001.fastq.gz"],
    "42-PC_S6_R2"  => ["/data/lehmanbd/insight_seq/TNBC_RNA-seq_phaseI/Run_3/42-PC_S6_L001_R2_001.fastq.gz"],
  },

  # cqstools file_def -i /data/lehmanbd/insight_seq/TNBC_RNA-seq_phaseI/Run_4/ -r -n \(.+_\)L001_\(.+\)_001
  "run4" => {
    "30-SE_S1_R1" => ["/data/lehmanbd/insight_seq/TNBC_RNA-seq_phaseI/Run_4/30-SE/30-SE_S1_L001_R1_001.fastq.gz"],
    "30-SE_S1_R2" => ["/data/lehmanbd/insight_seq/TNBC_RNA-seq_phaseI/Run_4/30-SE/30-SE_S1_L001_R2_001.fastq.gz"],
    "32-TE_S2_R1" => ["/data/lehmanbd/insight_seq/TNBC_RNA-seq_phaseI/Run_4/32-TE/32-TE_S2_L001_R1_001.fastq.gz"],
    "32-TE_S2_R2" => ["/data/lehmanbd/insight_seq/TNBC_RNA-seq_phaseI/Run_4/32-TE/32-TE_S2_L001_R2_001.fastq.gz"],
    "40-P_S3_R1"  => ["/data/lehmanbd/insight_seq/TNBC_RNA-seq_phaseI/Run_4/40-P/40-P_S3_L001_R1_001.fastq.gz"],
    "40-P_S3_R2"  => ["/data/lehmanbd/insight_seq/TNBC_RNA-seq_phaseI/Run_4/40-P/40-P_S3_L001_R2_001.fastq.gz"],
    "40-TX_S4_R1" => ["/data/lehmanbd/insight_seq/TNBC_RNA-seq_phaseI/Run_4/40-TX/40-TX_S4_L001_R1_001.fastq.gz"],
    "40-TX_S4_R2" => ["/data/lehmanbd/insight_seq/TNBC_RNA-seq_phaseI/Run_4/40-TX/40-TX_S4_L001_R2_001.fastq.gz"],
  },
};

my @runs = ( #"run3", 
"run4");

foreach my $run (@runs) {
  print ($files->{$run});
}