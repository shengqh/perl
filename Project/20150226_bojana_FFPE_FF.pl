#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;

my $task = "20150226_bojana_FFPE_FF";

my $target_dir = "/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF";

my $transcript_gtf       = "/scratch/cqs/shengq1/references/ensembl_gtf/v75/Homo_sapiens.GRCh37.75.M.gtf";
my $name_map_file        = "/scratch/cqs/shengq1/references/ensembl_gtf/v75/Homo_sapiens.GRCh37.75.M.map";
my $transcript_gtf_index = "/scratch/cqs/shengq1/references/ensembl_gtf/v75/Homo_sapiens.GRCh37.75.M.index/Homo_sapiens.GRCh37.75.M";
my $fasta_file_16569_M   = "/scratch/cqs/shengq1/references/hg19_16569_M/hg19_16569_M.fa";
my $bowtie2_index        = "/scratch/cqs/shengq1/references/hg19_16569_M/bowtie2_index_2.2.4/hg19_16569_M";
my $cqstools             = "/home/shengq1/cqstools/CQS.Tools.exe";

my $email = "quanhu.sheng\@vanderbilt.edu";

my $config = {
  general    => { task_name => $task },
  fastqfiles => {
#    "IG-001" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-1_1.fastq.gz",  "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-1_2.fastq.gz" ],
#    "IG-002" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-2_1.fastq.gz",  "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-2_2.fastq.gz" ],
#    "IG-003" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-3_1.fastq.gz",  "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-3_2.fastq.gz" ],
#    "IG-004" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-4_1.fastq.gz",  "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-4_2.fastq.gz" ],
#    "IG-005" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-5_1.fastq.gz",  "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-5_2.fastq.gz" ],
#    "IG-006" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-6_1.fastq.gz",  "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-6_2.fastq.gz" ],
#    "IG-007" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-7_1.fastq.gz",  "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-7_2.fastq.gz" ],
#    "IG-008" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-8_1.fastq.gz",  "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-8_2.fastq.gz" ],
#    "IG-009" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-9_1.fastq.gz",  "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-9_2.fastq.gz" ],
#    "IG-010" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-10_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-10_2.fastq.gz" ],
#    "IG-011" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-11_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-11_2.fastq.gz" ],
#    "IG-012" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-12_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-12_2.fastq.gz" ],
#    "IG-013" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-13_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-13_2.fastq.gz" ],
#    "IG-014" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-14_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-14_2.fastq.gz" ],
#    "IG-016" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-16_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-16_2.fastq.gz" ],
#    "IG-017" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-17_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-17_2.fastq.gz" ],
#    "IG-019" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-19_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-19_2.fastq.gz" ],
#    "IG-020" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-20_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-20_2.fastq.gz" ],
#    "IG-021" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-21_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-21_2.fastq.gz" ],
#    "IG-022" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-22_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-22_2.fastq.gz" ],
#    "IG-023" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-23_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-23_2.fastq.gz" ],
#    "IG-024" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-24_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-24_2.fastq.gz" ],
#    #three unpaired data begin
#    "IG-033" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-33_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-33_2.fastq.gz" ],
#    "IG-034" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-34_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-34_2.fastq.gz" ],
#    "IG-039" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-39_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-39_2.fastq.gz" ],
#    #three unpaired data end
#    "IG-040" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-40_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-40_2.fastq.gz" ],
#    "IG-041" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-41_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-41_2.fastq.gz" ],
#    "IG-042" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-42_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-42_2.fastq.gz" ],
#    "IG-043" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-43_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-43_2.fastq.gz" ],
#    "IG-044" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-44_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-44_2.fastq.gz" ],
#    "IG-045" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-45_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-45_2.fastq.gz" ],
#    "IG-046" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-46_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-46_2.fastq.gz" ],
#    "IG-047" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-47_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-47_2.fastq.gz" ],
#    "IG-049" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-49_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-49_2.fastq.gz" ],
#    "IG-050" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-50_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-50_2.fastq.gz" ],
#    "IG-051" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-51_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-51_2.fastq.gz" ],
#    "IG-052" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-52_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-52_2.fastq.gz" ],
#    "IG-053" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-53_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-53_2.fastq.gz" ],
#    "IG-054" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-54_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-54_2.fastq.gz" ],
#    "IG-055" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-55_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-55_2.fastq.gz" ],
#    "IG-056" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-56_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-56_2.fastq.gz" ],
#    "IG-057" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-57_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-57_2.fastq.gz" ],
#    "IG-058" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-58_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-58_2.fastq.gz" ],
#    "IG-059" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-59_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-59_2.fastq.gz" ],
#    "IG-060" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-60_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-60_2.fastq.gz" ],
#    "IG-061" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-61_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/IG-61_2.fastq.gz" ],
    "IG-062" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/2059-JP-10-1_AGTCAA_L001_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/2059-JP-10-1_AGTCAA_L001_R2_001.fastq.gz"
    ],
    "IG-063" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/2059-JP-12-1_CAGGCG_L003_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/2059-JP-12-1_CAGGCG_L003_R2_001.fastq.gz"
    ],
    "IG-064" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/2059-JP-10-2_ATGTCA_L001_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/2059-JP-10-2_ATGTCA_L001_R2_001.fastq.gz"
    ],
    "IG-065" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/2059-JP-12-2_CATTTT_L003_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/2059-JP-12-2_CATTTT_L003_R2_001.fastq.gz"
    ],
    "IG-066" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/2059-JP-16-1_CGTACG_L007_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/2059-JP-16-1_CGTACG_L007_R2_001.fastq.gz"
    ],
    "IG-067" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/2059-JP-12-3_CGGAAT_L003_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/2059-JP-12-3_CGGAAT_L003_R2_001.fastq.gz"
    ],
    "IG-068" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/2059-JP-16-2_CACGAT_L007_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/2059-JP-16-2_CACGAT_L007_R2_001.fastq.gz"
    ],
    "IG-069" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/2059-JP-12-4_CTATAC_L003_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/2059-JP-12-4_CTATAC_L003_R2_001.fastq.gz"
    ],
    "IG-070" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/2059-JP-17-1_GCCAAT_L008_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/2059-JP-17-1_GCCAAT_L008_R2_001.fastq.gz"
    ],
    "IG-071" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/2059-JP-12-5_TCATTC_L003_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/2059-JP-12-5_TCATTC_L003_R2_001.fastq.gz"
    ],
    "IG-072" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/2059-JP-17-2_CTTGTA_L008_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/2059-JP-17-2_CTTGTA_L008_R2_001.fastq.gz"
    ],
    "IG-073" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/2059-JP-13-1_GTTTCG_L004_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/2059-JP-13-1_GTTTCG_L004_R2_001.fastq.gz"
    ],
    "IG-074" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/2059-JP-10-3_GTGGCC_L001_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/2059-JP-10-3_GTGGCC_L001_R2_001.fastq.gz"
    ],
    "IG-075" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/2059-JP-13-2_CAAAAG_L004_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/2059-JP-13-2_CAAAAG_L004_R2_001.fastq.gz"
    ],
    "IG-076" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/2059-JP-10-4_GAGTGG_L001_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/2059-JP-10-4_GAGTGG_L001_R2_001.fastq.gz"
    ],
    "IG-077" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/2059-JP-13-3_CACTCA_L004_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/2059-JP-13-3_CACTCA_L004_R2_001.fastq.gz"
    ],
    "IG-078" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/2059-JP-10-5_CACCGG_L001_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/2059-JP-10-5_CACCGG_L001_R2_001.fastq.gz"
    ],
    "IG-079" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/2059-JP-13-4_CATGGC_L004_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/2059-JP-13-4_CATGGC_L004_R2_001.fastq.gz"
    ],
    "IG-080" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/2059-JP-17-3_GAGTGG_L008_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/2059-JP-17-3_GAGTGG_L008_R2_001.fastq.gz"
    ],
    "IG-081" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/2059-JP-13-5_CTCAGA_L004_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/2059-JP-13-5_CTCAGA_L004_R2_001.fastq.gz"
    ],
    "IG-082" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/2059-JP-17-4_TCCCGA_L008_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/2059-JP-17-4_TCCCGA_L008_R2_001.fastq.gz"
    ],
    "IG-083" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/2059-JP-14-1_GGTAGC_L005_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/2059-JP-14-1_GGTAGC_L005_R2_001.fastq.gz"
    ],
    "IG-085" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/2059-JP-14-2_CCAACA_L005_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/2059-JP-14-2_CCAACA_L005_R2_001.fastq.gz"
    ],
    "IG-086" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/2059-JP-11-1_CCGTCC_L002_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/2059-JP-11-1_CCGTCC_L002_R2_001.fastq.gz"
    ],
    "IG-087" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/2059-JP-14-3_CTAGCT_L005_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/2059-JP-14-3_CTAGCT_L005_R2_001.fastq.gz"
    ],
    "IG-088" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/2059-JP-11-2_GTAGAG_L002_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/2059-JP-11-2_GTAGAG_L002_R2_001.fastq.gz"
    ],
    "IG-089" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/2059-JP-14-4_TATAAT_L005_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/2059-JP-14-4_TATAAT_L005_R2_001.fastq.gz"
    ],
    "IG-090" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/2059-JP-11-3_GTGAAA_L002_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/2059-JP-11-3_GTGAAA_L002_R2_001.fastq.gz"
    ],
    "IG-091" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/2059-JP-14-5_TCGGCA_L005_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/2059-JP-14-5_TCGGCA_L005_R2_001.fastq.gz"
    ],
    "IG-092" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/2059-JP-16-3_ATTCCT_L007_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/2059-JP-16-3_ATTCCT_L007_R2_001.fastq.gz"
    ],
    "IG-093" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/2059-JP-15-1_AGTTCC_L006_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/2059-JP-15-1_AGTTCC_L006_R2_001.fastq.gz"
    ],
    "IG-094" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/2059-JP-16-4_CAACTA_L007_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/2059-JP-16-4_CAACTA_L007_R2_001.fastq.gz"
    ],
    "IG-095" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/2059-JP-15-2_GTCCGC_L006_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/2059-JP-15-2_GTCCGC_L006_R2_001.fastq.gz"
    ],
    "IG-096" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/2059-JP-16-5_GACGAC_L007_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/2059-JP-16-5_GACGAC_L007_R2_001.fastq.gz"
    ],
    "IG-097" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/2059-JP-15-3_ACTGAT_L006_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/2059-JP-15-3_ACTGAT_L006_R2_001.fastq.gz"
    ],
    "IG-098" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/2059-JP-15-4_TACAGC_L006_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/2059-JP-15-4_TACAGC_L006_R2_001.fastq.gz"
    ],
    "IG-099" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/2059-JP-11-4_ATGAGC_L002_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/2059-JP-11-4_ATGAGC_L002_R2_001.fastq.gz"
    ],
    "IG-100" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/2059-JP-15-5_TCGAAG_L006_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/2059-JP-15-5_TCGAAG_L006_R2_001.fastq.gz"
    ],
    "IG-101" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/2059-JP-11-5_TAATCG_L002_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/2059-JP-11-5_TAATCG_L002_R2_001.fastq.gz"
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

  fastqc => {
    class      => "QC::FastQC",
    perform    => 1,
    target_dir => "${target_dir}/fastqc",
    option     => "",
    source_ref => "fastqfiles",
    sh_direct  => 1,
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
    source_ref           => "fastqfiles",
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
  htseqcount => {
    class      => "Count::HTSeqCount",
    perform    => 1,
    target_dir => "${target_dir}/htseqcount",
    option     => "-r pos",
    source_ref => "tophat2",
    gff_file   => $transcript_gtf,
    sh_direct  => 1,
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
    source     => { individual => [ "fastqc", "tophat2", "htseqcount" ] },
    sh_direct  => 0,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
};

performConfig($config);

1;
