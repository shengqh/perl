#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;

my $task = "FFPE_FF_MiSeq";

my $target_dir = "/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq";

#my $target_dir = "e:/temp";

my $transcript_gtf       = "/scratch/cqs/shengq1/references/ensembl_gtf/v75/Homo_sapiens.GRCh37.75.M.gtf";
my $name_map_file        = "/scratch/cqs/shengq1/references/ensembl_gtf/v75/Homo_sapiens.GRCh37.75.M.map";
my $transcript_gtf_index = "/scratch/cqs/shengq1/references/ensembl_gtf/v75/gtfindex/Homo_sapiens.GRCh37.75.M";
my $fasta_file_16569_M   = "/scratch/cqs/shengq1/references/hg19_16569_M/hg19_16569_M.fa";
my $bowtie2_index        = "/scratch/cqs/shengq1/references/hg19_16569_M/bowtie2_index_2.2.4/hg19_16569_M";
my $cqstools             = "/home/shengq1/cqstools/CQS.Tools.exe";
my $dbsnp                = "/data/cqs/shengq1/reference/dbsnp/human_GRCh37_v141_16569_M.vcf";
my $gatk_jar             = "/home/shengq1/local/bin/GATK/GenomeAnalysisTK.jar";
my $picard_jar           = "/scratch/cqs/shengq1/local/bin/picard/picard.jar";
my $star_index           = "/scratch/cqs/shengq1/references/hg19_16569_M/STAR_index_v37.75_2.4.0j_sjdb75";
my $annovar_param        = "-protocol refGene,snp138,cosmic70 -operation g,f,f --remove";
my $annovar_db           = "/scratch/cqs/shengq1/references/annovar/humandb/";
my $email                = "quanhu.sheng\@vanderbilt.edu";

my $config = {
  general => { task_name => $task },
  files   => {
    "IG-001" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/30-PA_S1_L001_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/30-PA_S1_L001_R2_001.fastq.gz"
    ],
    "IG-002" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/30-TA_S3_L001_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/30-TA_S3_L001_R2_001.fastq.gz"
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
    "IG-012-2nd" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/IG-12_S3_L001_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/IG-12_S3_L001_R2_001.fastq.gz"
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
    "IG-029-2nd" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/IG-29B_S2_L001_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/IG-29B_S2_L001_R2_001.fastq.gz"
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
    "IG-039-2nd" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/IG-39-2nd_S3_L001_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/IG-39-2nd_S3_L001_R2_001.fastq.gz"
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
    "IG-049-2nd" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/IG-49_S1_L001_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/IG-49_S1_L001_R2_001.fastq.gz"
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
    "MiSeq_FF_All" => [
      "IG-001",     "IG-003", "IG-005", "IG-007", "IG-009", "IG-011", "IG-013", "IG-016", "IG-019", "IG-021", "IG-023", "IG-033", "IG-034", "IG-039",
      "IG-039-2nd", "IG-040", "IG-041", "IG-042", "IG-043", "IG-044", "IG-045", "IG-046", "IG-047", "IG-048", "IG-060", "IG-061"
    ],
    "MiSeq_FFPE_All" => [
      "IG-002", "IG-004", "IG-006",     "IG-008", "IG-010",     "IG-012", "IG-012-2nd", "IG-014", "IG-015", "IG-017", "IG-018", "IG-020", "IG-022", "IG-024", "IG-025", "IG-026", "IG-027",
      "IG-028", "IG-029", "IG-029-2nd", "IG-049", "IG-049-2nd", "IG-050", "IG-051",     "IG-052", "IG-053", "IG-054", "IG-055", "IG-056", "IG-057", "IG-058", "IG-059"
    ],

    "MiSeq_FF_Pilot"       => [ "IG-001", "IG-007", "IG-011", "IG-016", "IG-021", "IG-042", "IG-043", "IG-060", "IG-061" ],
    "MiSeq_FFPE_Pilot"     => [ "IG-002", "IG-008", "IG-012", "IG-017", "IG-022", "IG-051", "IG-052", "IG-059", "IG-058" ],
    "MiSeq_FF_Pilot_OLD"   => [ "IG-001", "IG-007", "IG-011", "IG-016", "IG-021" ],
    "MiSeq_FFPE_Pilot_OLD" => [ "IG-002", "IG-008", "IG-012", "IG-017", "IG-022" ],
    "MiSeq_FF_Pilot_NEW"   => [ "IG-042", "IG-043", "IG-060", "IG-061" ],
    "MiSeq_FFPE_Pilot_NEW" => [ "IG-051", "IG-052", "IG-059", "IG-058" ],
    "MiSeq_FF_ALL_NEW"   => [ "IG-033", "IG-034",     "IG-039", "IG-039-2nd", "IG-040", "IG-041", "IG-042", "IG-043", "IG-044", "IG-045", "IG-046", "IG-047", "IG-048", "IG-060", "IG-061" ],
    "MiSeq_FFPE_ALL_NEW" => [ "IG-049", "IG-049-2nd", "IG-050", "IG-051",     "IG-052", "IG-053", "IG-054", "IG-055", "IG-056", "IG-057", "IG-058", "IG-059" ],
  },
  pairs => {
    "MiSeq_FFPE_VS_FF_All"     => { groups => [ "MiSeq_FFPE_All",     "MiSeq_FF_All" ], },
    "MiSeq_FFPE_VS_FF_ALL_NEW" => { groups => [ "MiSeq_FFPE_ALL_NEW", "MiSeq_FF_ALL_NEW" ], },
    "MiSeq_FFPE_VS_FF_Pilot"   => {
      groups => [ "MiSeq_FFPE_Pilot", "MiSeq_FF_Pilot" ],
      paired => [ "B30A",             "B32A", "B33A", "B40A", "B42", "P06", "P07", "P14", "P13" ]
    },
    "MiSeq_FFPE_VS_FF_Pilot_OLD" => {
      groups => [ "MiSeq_FFPE_Pilot_OLD", "MiSeq_FF_Pilot_OLD" ],
      paired => [ "B30A",                 "B32A", "B33A", "B40A", "B42" ]
    },
    "MiSeq_FFPE_VS_FF_Pilot_NEW" => {
      groups => [ "MiSeq_FFPE_Pilot_NEW", "MiSeq_FF_Pilot_NEW" ],
      paired => [ "P06",                  "P07", "P14", "P13" ]
    },
  },
  trimmer => {
    class      => "CQS::FastqTrimmer",
    perform    => 0,
    target_dir => "${target_dir}/trim_terminalN",
    option     => "-n -z -m 30",
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
  fastqlen => {
    class      => "FastqLen",
    perform    => 0,
    target_dir => "${target_dir}/fastqlen",
    option     => "",
    source_ref => "trimmer",
    cqstools   => $cqstools,
    sh_direct  => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "24",
      "mem"      => "20gb"
    },
  },
  fastqc => {
    class      => "QC::FastQC",
    perform    => 0,
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
  star => {
    class      => "Alignment::STAR",
    perform    => 0,
    target_dir => "${target_dir}/star",
    option     => "",
    source_ref => "trimmer",
    genome_dir => $star_index,
    sh_direct  => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "30gb"
    },
  },
  star_index => {
    class          => "Alignment::STARIndex",
    perform        => 0,
    target_dir     => "${target_dir}/star_index",
    option         => "--sjdbOverhang 75",
    source_ref     => [ "star", "tab\$" ],
    fasta_file     => $fasta_file_16569_M,
    transcript_gtf => $transcript_gtf,
    sh_direct      => 1,
    pbs            => {
      "email"    => $email,
      "nodes"    => "1:ppn=24",
      "walltime" => "72",
      "mem"      => "30gb"
    },
  },
  star_2nd_pass => {
    class           => "Alignment::STAR",
    perform         => 0,
    target_dir      => "${target_dir}/star_2nd_pass",
    option          => "",
    source_ref      => "trimmer",
    genome_dir_ref  => "star_index",
    output_unsorted => 1,
    sh_direct       => 0,
    pbs             => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "30gb"
    },
  },
  star_htseqcount => {
    class      => "Count::HTSeqCount",
    perform    => 0,
    target_dir => "${target_dir}/star_htseqcount",
    option     => "",
    source_ref => [ "star_2nd_pass", "_Aligned.out.bam" ],
    gff_file   => $transcript_gtf,
    ispairend  => 1,
    sh_direct  => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  star_genetable => {
    class         => "CQS::CQSDatatable",
    perform       => 1,
    target_dir    => "${target_dir}/star_genetable",
    option        => "-p ENS --noheader -o ${task}_gene.count",
    source_ref    => "star_htseqcount",
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
  star_deseq2 => {
    class         => "Comparison::DESeq2",
    perform       => 1,
    target_dir    => "${target_dir}/star_deseq2",
    option        => "",
    source_ref    => "pairs",
    groups_ref    => "groups",
    countfile_ref => "star_genetable",
    sh_direct     => 1,
    pbs           => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "10",
      "mem"      => "10gb"
    },
  },
  star_2nd_pass_refine => {
    class      => "GATK::RNASeqRefine",
    perform    => 0,
    target_dir => "${target_dir}/star_2nd_pass_refine",
    option     => "-Xmx40g",
    fasta_file => $fasta_file_16569_M,
    source_ref => "star_2nd_pass",
    vcf_files  => [$dbsnp],
    gatk_jar   => $gatk_jar,
    picard_jar => $picard_jar,
    sorted     => 0,
    sh_direct  => 0,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  star_2nd_pass_refine_SNPindel => {
    class       => "GATK::SNPIndel",
    perform     => 0,
    target_dir  => "${target_dir}/star_2nd_pass_refine_SNPindel",
    option      => "-l INFO -G Standard -stand_call_conf 50.0 -stand_emit_conf 10.0 -dcov 200",
    source_ref  => "star_2nd_pass_refine",
    java_option => "",
    fasta_file  => $fasta_file_16569_M,
    vcf_files   => [$dbsnp],
    gatk_jar    => $gatk_jar,
    pbs         => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  star_2nd_pass_refine_SNPindel_annovar => {
    class      => "Annovar",
    perform    => 0,
    target_dir => "${target_dir}/star_2nd_pass_refine_SNPindel_annovar",
    source_ref => "star_2nd_pass_refine_SNPindel",
    option     => $annovar_param,
    annovar_db => $annovar_db,
    buildver   => "hg19",
    sh_direct  => 1,
    isvcf      => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "10gb"
    },
  },

  tophat2 => {
    class                => "Alignment::Tophat2",
    perform              => 0,
    target_dir           => "${target_dir}/tophat2",
    option               => "-p 8",
    source_ref           => "trimmer",
    bowtie2_index        => $bowtie2_index,
    transcript_gtf       => $transcript_gtf,
    transcript_gtf_index => $transcript_gtf_index,
    rename_bam           => 1,
    sh_direct            => 1,
    pbs                  => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "30gb"
    },
  },
  tophat2_sortbam => {
    class         => "Samtools::Sort",
    perform       => 0,
    target_dir    => "${target_dir}/tophat2_sortname",
    option        => "",
    source_ref    => "tophat2",
    sort_by_query => 1,
    sh_direct     => 1,
    pbs           => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "20gb"
    },
  },
  tophat2_htseqcount => {
    class      => "Count::HTSeqCount",
    perform    => 0,
    target_dir => "${target_dir}/tophat2_htseqcount",
    option     => "",
    source_ref => "tophat2_sortbam",
    gff_file   => $transcript_gtf,
    ispairend  => 1,
    sh_direct  => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  tophat2_genetable => {
    class         => "CQS::CQSDatatable",
    perform       => 1,
    target_dir    => "${target_dir}/tophat2_genetable",
    option        => "-p ENS --noheader -o ${task}_gene.count",
    source_ref    => "tophat2_htseqcount",
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
  tophat2_deseq2 => {
    class         => "Comparison::DESeq2",
    perform       => 1,
    target_dir    => "${target_dir}/tophat2_deseq2",
    option        => "",
    source_ref    => "pairs",
    groups_ref    => "groups",
    countfile_ref => "tophat2_genetable",
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
    perform    => 0,
    target_dir => "${target_dir}/sequencetask",
    option     => "",
    source     => {
      step_1 => [ "trimmer", "fastqlen", "fastqc", "star" ],
      step_2 => ["star_index"],
      step_3 => [ "star_2nd_pass",                 "star_htseqcount", "star_2nd_pass_refine" ],
      step_4 => [ "star_genetable",                "star_deseq2", ],
      step_5 => [ "star_2nd_pass_refine_SNPindel", "star_2nd_pass_refine_SNPindel_annovar" ]
    },
    sh_direct => 1,
    pbs       => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
};

performConfig($config);

#performTask( $config, "star_genetable" );
#performTask( $config, "tophat2_genetable" );

1;

