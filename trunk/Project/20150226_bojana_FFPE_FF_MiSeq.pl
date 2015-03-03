#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;

my $task = "20150226_bojana_FFPE_FF";

my $target_dir = "/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq";

my $transcript_gtf       = "/scratch/cqs/shengq1/references/ensembl_gtf/v75/Homo_sapiens.GRCh37.75.M.gtf";
my $name_map_file        = "/scratch/cqs/shengq1/references/ensembl_gtf/v75/Homo_sapiens.GRCh37.75.M.map";
my $transcript_gtf_index = "/scratch/cqs/shengq1/references/ensembl_gtf/v75/gtfindex/Homo_sapiens.GRCh37.75.M";
my $fasta_file_16569_M   = "/data/cqs/guoy1/reference/hg19/bowtie2_index/hg19.fa";
my $bowtie2_index        = "/data/cqs/guoy1/reference/hg19/bowtie2_index/hg19";
my $cqstools             = "/home/shengq1/cqstools/CQS.Tools.exe";

my $email = "quanhu.sheng\@vanderbilt.edu";

my $config = {
  general => { task_name => $task },
  files   => {

#Unknown => ["/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/30-TA_S3_L001_R1_001.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/30-TA_S3_L001_R2_001.fastq.gz"],
#"IG-12_S3" => ["/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/IG-12_S3_L001_R1_001.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/IG-12_S3_L001_R2_001.fastq.gz"],
#"IG-29B_S2" => ["/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/IG-29B_S2_L001_R1_001.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/IG-29B_S2_L001_R2_001.fastq.gz"],
#repeat "IG-39-2nd_S3" => ["/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/IG-39-2nd_S3_L001_R1_001.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/IG-39-2nd_S3_L001_R2_001.fastq.gz"],
#repeat "IG-049" => ["/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/IG-49_S1_L001_R1_001.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/IG-49_S1_L001_R2_001.fastq.gz"],
    "IG-001" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/30-PA_S1_L001_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/30-PA_S1_L001_R2_001.fastq.gz"
    ],
    "IG-003" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/30-TC_S1_L001_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/30-TC_S1_L001_R2_001.fastq.gz"
    ],
    "IG-004" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/30-TXE_S1_L001_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/30-TXE_S1_L001_R2_001.fastq.gz"
    ],
    "IG-005" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/30-SB_S2_L001_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/30-SB_S2_L001_R2_001.fastq.gz"
    ],
    "IG-006" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/30-SE_S1_L001_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/30-SE_S1_L001_R2_001.fastq.gz"
    ],
    "IG-007" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/32-PD_S1_L001_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/32-PD_S1_L001_R2_001.fastq.gz"
    ],
    "IG-008" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/32-PE_S2_L001_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/32-PE_S2_L001_R2_001.fastq.gz"
    ],
    "IG-009" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/32-TA_S3_L001_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/32-TA_S3_L001_R2_001.fastq.gz"
    ],
    "IG-010" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/32-TE_S2_L001_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/32-TE_S2_L001_R2_001.fastq.gz"
    ],
    "IG-011" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/33-PA_S2_L001_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/33-PA_S2_L001_R2_001.fastq.gz"
    ],
    "IG-012" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/33-PI_S3_L001_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/33-PI_S3_L001_R2_001.fastq.gz"
    ],
    "IG-013" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/33-TA_S2_L001_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/33-TA_S2_L001_R2_001.fastq.gz"
    ],
    "IG-014" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/33-TG_S1_L001_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/33-TG_S1_L001_R2_001.fastq.gz"
    ],
    "IG-015" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/33-SN_S2_L001_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/33-SN_S2_L001_R2_001.fastq.gz"
    ],
    "IG-016" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/40-PA_S4_L001_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/40-PA_S4_L001_R2_001.fastq.gz"
    ],
    "IG-017" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/40-P_S3_L001_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/40-P_S3_L001_R2_001.fastq.gz"
    ],
    "IG-018" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/40-T_S4_L001_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/40-T_S4_L001_R2_001.fastq.gz"
    ],
    "IG-019" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/40-TXA_S5_L001_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/40-TXA_S5_L001_R2_001.fastq.gz"
    ],
    "IG-020" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/40-TX_S4_L001_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/40-TX_S4_L001_R2_001.fastq.gz"
    ],
    "IG-021" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/42-PC_S6_L001_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/42-PC_S6_L001_R2_001.fastq.gz"
    ],
    "IG-022" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/42-P_S2_L001_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/42-P_S2_L001_R2_001.fastq.gz"
    ],
    "IG-023" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/B42-TD_S3_L001_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/B42-TD_S3_L001_R2_001.fastq.gz"
    ],
    "IG-024" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/42-T_S3_L001_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/42-T_S3_L001_R2_001.fastq.gz"
    ],
    "IG-025" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/402750_S4_L001_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/402750_S4_L001_R2_001.fastq.gz"
    ],
    "IG-026" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/402781_S1_L001_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/402781_S1_L001_R2_001.fastq.gz"
    ],
    "IG-027" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/402830_S2_L001_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/402830_S2_L001_R2_001.fastq.gz"
    ],
    "IG-028" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/405375_S3_L001_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/405375_S3_L001_R2_001.fastq.gz"
    ],
    "IG-029" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/407730_S4_L001_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/407730_S4_L001_R2_001.fastq.gz"
    ],
    "IG-033" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/408637_S4_L001_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/408637_S4_L001_R2_001.fastq.gz"
    ],
    "IG-034" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/IG-34_S3_L001_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/IG-34_S3_L001_R2_001.fastq.gz"
    ],
    "IG-039" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/408648_S5_L001_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/408648_S5_L001_R2_001.fastq.gz"
    ],
    "IG-040" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/IG-40_S1_L001_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/IG-40_S1_L001_R2_001.fastq.gz"
    ],
    "IG-041" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/IG-41_S1_L001_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/IG-41_S1_L001_R2_001.fastq.gz"
    ],
    "IG-042" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/IG-42_S1_L001_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/IG-42_S1_L001_R2_001.fastq.gz"
    ],
    "IG-043" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/IG-43_S2_L001_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/IG-43_S2_L001_R2_001.fastq.gz"
    ],
    "IG-044" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/IG-44_S2_L001_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/IG-44_S2_L001_R2_001.fastq.gz"
    ],
    "IG-045" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/IG-45_S1_L001_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/IG-45_S1_L001_R2_001.fastq.gz"
    ],
    "IG-046" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/IG-46_S2_L001_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/IG-46_S2_L001_R2_001.fastq.gz"
    ],
    "IG-047" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/IG-47_S3_L001_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/IG-47_S3_L001_R2_001.fastq.gz"
    ],
    "IG-048" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/IG-48_S2_L001_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/IG-48_S2_L001_R2_001.fastq.gz"
    ],
    "IG-049" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/IG-49_S3_L001_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/IG-49_S3_L001_R2_001.fastq.gz"
    ],
    "IG-050" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/IG-50_S2_L001_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/IG-50_S2_L001_R2_001.fastq.gz"
    ],
    "IG-051" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/IG-51_S1_L001_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/IG-51_S1_L001_R2_001.fastq.gz"
    ],
    "IG-052" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/IG-52_S2_L001_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/IG-52_S2_L001_R2_001.fastq.gz"
    ],
    "IG-053" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/IG-53_S3_L001_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/IG-53_S3_L001_R2_001.fastq.gz"
    ],
    "IG-054" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/IG-54_S2_L001_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/IG-54_S2_L001_R2_001.fastq.gz"
    ],
    "IG-055" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/IG-55_S3_L001_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/IG-55_S3_L001_R2_001.fastq.gz"
    ],
    "IG-056" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/IG-56_S4_L001_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/IG-56_S4_L001_R2_001.fastq.gz"
    ],
    "IG-057" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/IG-57_S1_L001_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/IG-57_S1_L001_R2_001.fastq.gz"
    ],
    "IG-058" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/IG-58_S3_L001_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/IG-58_S3_L001_R2_001.fastq.gz"
    ],
    "IG-059" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/IG-59_S4_L001_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/IG-59_S4_L001_R2_001.fastq.gz"
    ],
    "IG-060" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/IG-60_S3_L001_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/IG-60_S3_L001_R2_001.fastq.gz"
    ],
    "IG-061" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/IG-61_S4_L001_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/IG-61_S4_L001_R2_001.fastq.gz"
    ],
  },
  groups => {

    "MiSeq_FF"       => [ "IG-001", "IG-007", "IG-011", "IG-016", "IG-021", "IG-042", "IG-043", "IG-060", "IG-061" ],
    "MiSeq_FFPE"     => [ "IG-002", "IG-008", "IG-012", "IG-017", "IG-022", "IG-051", "IG-052", "IG-059", "IG-058" ],
    "MiSeq_FF_OLD"   => [ "IG-001", "IG-007", "IG-011", "IG-016", "IG-021" ],
    "MiSeq_FFPE_OLD" => [ "IG-002", "IG-008", "IG-012", "IG-017", "IG-022" ],
    "MiSeq_FF_NEW"   => [ "IG-042", "IG-043", "IG-060", "IG-061" ],
    "MiSeq_FFPE_NEW" => [ "IG-051", "IG-052", "IG-059", "IG-058" ],
  },
  pairs => {

    "MiSeq_FFPE_VS_FF" => {
      groups => [ "MiSeq_FFPE", "MiSeq_FF" ],
      paired => [ "B30A",       "B32A", "B33A", "B40A", "B42", "P06", "P07", "P14", "P13" ]
    },
    "MiSeq_FFPE_VS_FF_OLD" => {
      groups => [ "MiSeq_FFPE_OLD", "MiSeq_FF_OLD" ],
      paired => [ "B30A",           "B32A", "B33A", "B40A", "B42" ]
    },
    "MiSeq_FFPE_VS_FF_NEW" => {
      groups => [ "MiSeq_FFPE_NEW", "MiSeq_FF_NEW" ],
      paired => [ "P06",            "P07", "P14", "P13" ]
    },
  },
  trimmer => {
    class      => "CQS::FastqTrimmer",
    perform    => 1,
    target_dir => "${target_dir}/FastqTrimmer",
    option     => "-n -z",
    extension  => "_trim.fastq.gz",
    source_ref => "files",
    cqstools   => $cqstools,
    sh_direct  => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "2",
      "mem"      => "10gb"
    },
  },
  fastqc => {
    class      => "QC::FastQC",
    perform    => 1,
    target_dir => "${target_dir}/fastqc",
    option     => "",
    source_ref => "files",
    sh_direct  => 0,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=2",
      "walltime" => "2",
      "mem"      => "10gb"
    },
  },
  tophat2 => {
    class                => "Alignment::Tophat2",
    perform              => 1,
    target_dir           => "${target_dir}/tophat2",
    option               => "--segment-length 25 -r 0 -p 8",
    source_ref           => "trimmer",
    bowtie2_index        => $bowtie2_index,
    transcript_gtf       => $transcript_gtf,
    transcript_gtf_index => $transcript_gtf_index,
    rename_bam           => 1,
    sh_direct            => 0,
    pbs                  => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "30gb"
    },
  },
  sortbam => {
    class         => "Samtools::Sort",
    perform       => 1,
    target_dir    => "${target_dir}/sortname",
    option        => "",
    source_ref    => "tophat2",
    sort_by_query => 1,
    sh_direct     => 0,
    pbs           => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "20gb"
    },
  },
  htseqcount => {
    class      => "Count::HTSeqCount",
    perform    => 1,
    target_dir => "${target_dir}/htseqcount",
    option     => "",
    source_ref => "sortbam",
    gff_file   => $transcript_gtf,
    sh_direct  => 0,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  genetable => {
    class         => "CQS::CQSDatatable",
    perform       => 1,
    target_dir    => "${target_dir}/genetable",
    option        => "-p ENS --noheader -o ${task}_gene.count",
    source_ref    => "htseqcount",
    name_map_file => $name_map_file,
    cqs_tools     => $cqstools,
    sh_direct     => 1,
    pbs           => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "10",
      "mem"      => "10gb"
    },
  },
  deseq2 => {
    class         => "Comparison::DESeq2",
    perform       => 1,
    target_dir    => "${target_dir}/deseq2",
    option        => "",
    source_ref    => "pairs",
    groups_ref    => "groups",
    countfile_ref => "genetable",
    sh_direct     => 1,
    pbs           => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "10",
      "mem"      => "10gb"
    },
  },
  sequencetask => {
    class      => "CQS::SequenceTask",
    perform    => 1,
    target_dir => "${target_dir}/sequencetask",
    option     => "",
    source     => {
      individual => [ "trimmer",   "fastqc", "tophat2", "sortbam", "htseqcount" ],
      summary    => [ "genetable", "deseq2" ],
    },
    sh_direct => 0,
    pbs       => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
};

performConfig($config);

1;
