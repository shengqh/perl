#!/usr/bin/perl
use strict;
use warnings;

use CQS::FileUtils;
use CQS::ClassFactory;

my $target_dir = create_directory_or_die("/scratch/cqs/shengq1/somaticmutation_comparison");

my $email         = "quanhu.sheng\@vanderbilt.edu";
my $cqstools      = "/home/shengq1/cqstools/CQS.Tools.exe";
my $samtools      = "/home/shengq1/local/bin/samtools/samtools";
my $rnaediting_db = "/data/cqs/shengq1/reference/rnaediting/hg19.txt";
my $rsmc          = "/home/shengq1/rsmc/rsmc.exe";

##hg19.16569.MT###
my $fasta_file_16569_MT    = "/data/cqs/shengq1/reference/hg19_16569_MT/bwa_index_0.7.8/hg19_16569_MT.fa";
my $bowtie2_index_16569_MT = "/data/cqs/shengq1/reference/hg19_16569_MT/bowtie2_index_2.1.0/hg19_16569_MT";
my $cosmic_file_16569_MT   = "/data/cqs/shengq1/reference/cosmic/cosmic_v69_hg19_16569_MT.vcf";
my $snp_file_16569_MT      = "/data/cqs/shengq1/reference/dbsnp/human_GRCh37_v141_16569_MT.vcf";

##hg19.16569.M###
my $fasta_file_16569_M    = "/data/cqs/shengq1/reference/hg19_16569_M/hg19_16569_M.fa";
my $bowtie2_index_16569_M = "/data/cqs/shengq1/reference/hg19_16569_M/bowtie2_index_2.1.0/hg19_16569_M";
my $cosmic_file_16569_M   = "/data/cqs/shengq1/reference/cosmic/cosmic_v69_hg19_16569_M.vcf";
my $snp_file_16569_M      = "/data/cqs/shengq1/reference/dbsnp/human_GRCh37_v141_16569_M.vcf";

my $transcript_gtf       = "/data/cqs/shengq1/reference/ensembl_gtf/Homo_sapiens.GRCh37.75.MT.gtf";
my $transcript_gtf_index = "/scratch/cqs/shengq1/gtfindex/Homo_sapiens.GRCh37.75.MT";
my $hg19_map             = "/data/cqs/shengq1/reference/ensembl_gtf/Homo_sapiens.GRCh37.75.map";

my $annovar_param = "-protocol refGene,snp138,cosmic70 -operation g,f,f --remove";
my $annovar_db    = "/scratch/cqs/shengq1/references/annovar/humandb/";

my $mutect        = "/home/shengq1/local/bin/mutect-1.1.7.jar";
my $varscan2      = "/home/shengq1/local/bin/VarScan.v2.3.7.jar";
my $affy_file     = "/data/cqs/shengq1/reference/affy/HG-U133_Plus_2.na33.annot.csv";
my $picard        = "/home/shengq1/local/bin/picard/";
my $rnaediting_db = "/data/cqs/shengq1/reference/rnaediting/hg19.txt";
my $rsmc          = "/home/shengq1/rsmc/rsmc.exe";

my $tcga = {
  dna => {
    "TCGA-A7-A0D9-DNA-NB" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/DNA_NB/TCGA-A7-A0D9-10A-01W-A071-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0B3-DNA-NB" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/DNA_NB/TCGA-BH-A0B3-10A-01W-A071-09_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0B8-DNA-NB" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/DNA_NB/TCGA-BH-A0B8-10A-01W-A071-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0BJ-DNA-NB" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/DNA_NB/TCGA-BH-A0BJ-10A-01W-A071-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0BM-DNA-NB" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/DNA_NB/TCGA-BH-A0BM-10A-01W-A071-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0C0-DNA-NB" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/DNA_NB/TCGA-BH-A0C0-10A-01W-A071-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0DK-DNA-NB" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/DNA_NB/TCGA-BH-A0DK-10A-01W-A071-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0DP-DNA-NB" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/DNA_NB/TCGA-BH-A0DP-10A-01W-A071-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0E0-DNA-NB" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/DNA_NB/TCGA-BH-A0E0-10A-01W-A071-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0H7-DNA-NB" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/DNA_NB/TCGA-BH-A0H7-10A-01W-A071-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-A7-A0D9-DNA-NT" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/DNA_NT/TCGA-A7-A0D9-11A-53W-A100-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0B3-DNA-NT" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/DNA_NT/TCGA-BH-A0B3-11B-21W-A100-09_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0B8-DNA-NT" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/DNA_NT/TCGA-BH-A0B8-11A-41W-A10F-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0BJ-DNA-NT" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/DNA_NT/TCGA-BH-A0BJ-11A-23W-A100-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0BM-DNA-NT" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/DNA_NT/TCGA-BH-A0BM-11A-12W-A10F-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0C0-DNA-NT" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/DNA_NT/TCGA-BH-A0C0-11A-21W-A100-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0DK-DNA-NT" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/DNA_NT/TCGA-BH-A0DK-11A-13W-A100-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0DP-DNA-NT" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/DNA_NT/TCGA-BH-A0DP-11A-12W-A100-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0E0-DNA-NT" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/DNA_NT/TCGA-BH-A0E0-11A-13W-A10F-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0H7-DNA-NT" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/DNA_NT/TCGA-BH-A0H7-11A-13W-A100-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-A7-A0D9-DNA-TP" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/DNA_TP/TCGA-A7-A0D9-01A-31W-A071-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0B3-DNA-TP" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/DNA_TP/TCGA-BH-A0B3-01A-11W-A071-09_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0B8-DNA-TP" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/DNA_TP/TCGA-BH-A0B8-01A-21W-A071-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0BJ-DNA-TP" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/DNA_TP/TCGA-BH-A0BJ-01A-11W-A071-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0BM-DNA-TP" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/DNA_TP/TCGA-BH-A0BM-01A-11W-A071-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0C0-DNA-TP" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/DNA_TP/TCGA-BH-A0C0-01A-21W-A071-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0DK-DNA-TP" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/DNA_TP/TCGA-BH-A0DK-01A-21W-A071-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0DP-DNA-TP" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/DNA_TP/TCGA-BH-A0DP-01A-21W-A071-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0E0-DNA-TP" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/DNA_TP/TCGA-BH-A0E0-01A-11W-A071-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0H7-DNA-TP" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/DNA_TP/TCGA-BH-A0H7-01A-13W-A071-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
  },
  rna => {
    "TCGA-A7-A0D9-RNA-NT" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/RNA_NT/TCGA-A7-A0D9-RNA_NT_sorted.bam"],
    "TCGA-BH-A0B3-RNA-NT" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/RNA_NT/TCGA-BH-A0B3-RNA_NT_sorted.bam"],
    "TCGA-BH-A0B8-RNA-NT" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/RNA_NT/TCGA-BH-A0B8-RNA_NT_sorted.bam"],
    "TCGA-BH-A0BJ-RNA-NT" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/RNA_NT/TCGA-BH-A0BJ-RNA_NT_sorted.bam"],
    "TCGA-BH-A0BM-RNA-NT" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/RNA_NT/TCGA-BH-A0BM-RNA_NT_sorted.bam"],
    "TCGA-BH-A0C0-RNA-NT" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/RNA_NT/TCGA-BH-A0C0-RNA_NT_sorted.bam"],
    "TCGA-BH-A0DK-RNA-NT" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/RNA_NT/TCGA-BH-A0DK-RNA_NT_sorted.bam"],
    "TCGA-BH-A0DP-RNA-NT" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/RNA_NT/TCGA-BH-A0DP-RNA_NT_sorted.bam"],
    "TCGA-BH-A0E0-RNA-NT" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/RNA_NT/TCGA-BH-A0E0-RNA_NT_sorted.bam"],
    "TCGA-BH-A0H7-RNA-NT" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/RNA_NT/TCGA-BH-A0H7-RNA_NT_sorted.bam"],
    "TCGA-A7-A0D9-RNA-TP" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/RNA_TP/TCGA-A7-A0D9-RNA_TP_sorted.bam"],
    "TCGA-BH-A0B3-RNA-TP" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/RNA_TP/TCGA-BH-A0B3-RNA_TP_sorted.bam"],
    "TCGA-BH-A0B8-RNA-TP" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/RNA_TP/TCGA-BH-A0B8-RNA_TP_sorted.bam"],
    "TCGA-BH-A0BJ-RNA-TP" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/RNA_TP/TCGA-BH-A0BJ-RNA_TP_sorted.bam"],
    "TCGA-BH-A0BM-RNA-TP" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/RNA_TP/TCGA-BH-A0BM-RNA_TP_sorted.bam"],
    "TCGA-BH-A0C0-RNA-TP" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/RNA_TP/TCGA-BH-A0C0-RNA_TP_sorted.bam"],
    "TCGA-BH-A0DK-RNA-TP" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/RNA_TP/TCGA-BH-A0DK-RNA_TP_sorted.bam"],
    "TCGA-BH-A0DP-RNA-TP" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/RNA_TP/TCGA-BH-A0DP-RNA_TP_sorted.bam"],
    "TCGA-BH-A0E0-RNA-TP" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/RNA_TP/TCGA-BH-A0E0-RNA_TP_sorted.bam"],
    "TCGA-BH-A0H7-RNA-TP" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/RNA_TP/TCGA-BH-A0H7-RNA_TP_sorted.bam"],
  },
  rna_dna_groups => {
    "TCGA-A7-A0D9-RNA-TP-DNA-NB" => [ "TCGA-A7-A0D9-DNA-NB", "TCGA-A7-A0D9-RNA-TP" ],
    "TCGA-BH-A0B3-RNA-TP-DNA-NB" => [ "TCGA-BH-A0B3-DNA-NB", "TCGA-BH-A0B3-RNA-TP" ],
    "TCGA-BH-A0B8-RNA-TP-DNA-NB" => [ "TCGA-BH-A0B8-DNA-NB", "TCGA-BH-A0B8-RNA-TP" ],
    "TCGA-BH-A0BJ-RNA-TP-DNA-NB" => [ "TCGA-BH-A0BJ-DNA-NB", "TCGA-BH-A0BJ-RNA-TP" ],
    "TCGA-BH-A0BM-RNA-TP-DNA-NB" => [ "TCGA-BH-A0BM-DNA-NB", "TCGA-BH-A0BM-RNA-TP" ],
    "TCGA-BH-A0C0-RNA-TP-DNA-NB" => [ "TCGA-BH-A0C0-DNA-NB", "TCGA-BH-A0C0-RNA-TP" ],
    "TCGA-BH-A0DK-RNA-TP-DNA-NB" => [ "TCGA-BH-A0DK-DNA-NB", "TCGA-BH-A0DK-RNA-TP" ],
    "TCGA-BH-A0DP-RNA-TP-DNA-NB" => [ "TCGA-BH-A0DP-DNA-NB", "TCGA-BH-A0DP-RNA-TP" ],
    "TCGA-BH-A0E0-RNA-TP-DNA-NB" => [ "TCGA-BH-A0E0-DNA-NB", "TCGA-BH-A0E0-RNA-TP" ],
    "TCGA-BH-A0H7-RNA-TP-DNA-NB" => [ "TCGA-BH-A0H7-DNA-NB", "TCGA-BH-A0H7-RNA-TP" ],

    "TCGA-A7-A0D9-RNA-TP-DNA-NT" => [ "TCGA-A7-A0D9-DNA-NT", "TCGA-A7-A0D9-RNA-TP" ],
    "TCGA-BH-A0B3-RNA-TP-DNA-NT" => [ "TCGA-BH-A0B3-DNA-NT", "TCGA-BH-A0B3-RNA-TP" ],
    "TCGA-BH-A0B8-RNA-TP-DNA-NT" => [ "TCGA-BH-A0B8-DNA-NT", "TCGA-BH-A0B8-RNA-TP" ],
    "TCGA-BH-A0BJ-RNA-TP-DNA-NT" => [ "TCGA-BH-A0BJ-DNA-NT", "TCGA-BH-A0BJ-RNA-TP" ],
    "TCGA-BH-A0BM-RNA-TP-DNA-NT" => [ "TCGA-BH-A0BM-DNA-NT", "TCGA-BH-A0BM-RNA-TP" ],
    "TCGA-BH-A0C0-RNA-TP-DNA-NT" => [ "TCGA-BH-A0C0-DNA-NT", "TCGA-BH-A0C0-RNA-TP" ],
    "TCGA-BH-A0DK-RNA-TP-DNA-NT" => [ "TCGA-BH-A0DK-DNA-NT", "TCGA-BH-A0DK-RNA-TP" ],
    "TCGA-BH-A0DP-RNA-TP-DNA-NT" => [ "TCGA-BH-A0DP-DNA-NT", "TCGA-BH-A0DP-RNA-TP" ],
    "TCGA-BH-A0E0-RNA-TP-DNA-NT" => [ "TCGA-BH-A0E0-DNA-NT", "TCGA-BH-A0E0-RNA-TP" ],
    "TCGA-BH-A0H7-RNA-TP-DNA-NT" => [ "TCGA-BH-A0H7-DNA-NT", "TCGA-BH-A0H7-RNA-TP" ],

    "TCGA-A7-A0D9-RNA-TP-DNA-TP" => [ "TCGA-A7-A0D9-DNA-TP", "TCGA-A7-A0D9-RNA-TP" ],
    "TCGA-BH-A0B3-RNA-TP-DNA-TP" => [ "TCGA-BH-A0B3-DNA-TP", "TCGA-BH-A0B3-RNA-TP" ],
    "TCGA-BH-A0B8-RNA-TP-DNA-TP" => [ "TCGA-BH-A0B8-DNA-TP", "TCGA-BH-A0B8-RNA-TP" ],
    "TCGA-BH-A0BJ-RNA-TP-DNA-TP" => [ "TCGA-BH-A0BJ-DNA-TP", "TCGA-BH-A0BJ-RNA-TP" ],
    "TCGA-BH-A0BM-RNA-TP-DNA-TP" => [ "TCGA-BH-A0BM-DNA-TP", "TCGA-BH-A0BM-RNA-TP" ],
    "TCGA-BH-A0C0-RNA-TP-DNA-TP" => [ "TCGA-BH-A0C0-DNA-TP", "TCGA-BH-A0C0-RNA-TP" ],
    "TCGA-BH-A0DK-RNA-TP-DNA-TP" => [ "TCGA-BH-A0DK-DNA-TP", "TCGA-BH-A0DK-RNA-TP" ],
    "TCGA-BH-A0DP-RNA-TP-DNA-TP" => [ "TCGA-BH-A0DP-DNA-TP", "TCGA-BH-A0DP-RNA-TP" ],
    "TCGA-BH-A0E0-RNA-TP-DNA-TP" => [ "TCGA-BH-A0E0-DNA-TP", "TCGA-BH-A0E0-RNA-TP" ],
    "TCGA-BH-A0H7-RNA-TP-DNA-TP" => [ "TCGA-BH-A0H7-DNA-TP", "TCGA-BH-A0H7-RNA-TP" ],

    "TCGA-A7-A0D9-RNA-NT-DNA-NB" => [ "TCGA-A7-A0D9-DNA-NB", "TCGA-A7-A0D9-RNA-NT" ],
    "TCGA-BH-A0B3-RNA-NT-DNA-NB" => [ "TCGA-BH-A0B3-DNA-NB", "TCGA-BH-A0B3-RNA-NT" ],
    "TCGA-BH-A0B8-RNA-NT-DNA-NB" => [ "TCGA-BH-A0B8-DNA-NB", "TCGA-BH-A0B8-RNA-NT" ],
    "TCGA-BH-A0BJ-RNA-NT-DNA-NB" => [ "TCGA-BH-A0BJ-DNA-NB", "TCGA-BH-A0BJ-RNA-NT" ],
    "TCGA-BH-A0BM-RNA-NT-DNA-NB" => [ "TCGA-BH-A0BM-DNA-NB", "TCGA-BH-A0BM-RNA-NT" ],
    "TCGA-BH-A0C0-RNA-NT-DNA-NB" => [ "TCGA-BH-A0C0-DNA-NB", "TCGA-BH-A0C0-RNA-NT" ],
    "TCGA-BH-A0DK-RNA-NT-DNA-NB" => [ "TCGA-BH-A0DK-DNA-NB", "TCGA-BH-A0DK-RNA-NT" ],
    "TCGA-BH-A0DP-RNA-NT-DNA-NB" => [ "TCGA-BH-A0DP-DNA-NB", "TCGA-BH-A0DP-RNA-NT" ],
    "TCGA-BH-A0E0-RNA-NT-DNA-NB" => [ "TCGA-BH-A0E0-DNA-NB", "TCGA-BH-A0E0-RNA-NT" ],
    "TCGA-BH-A0H7-RNA-NT-DNA-NB" => [ "TCGA-BH-A0H7-DNA-NB", "TCGA-BH-A0H7-RNA-NT" ],

    "TCGA-A7-A0D9-RNA-NT-DNA-NT" => [ "TCGA-A7-A0D9-DNA-NT", "TCGA-A7-A0D9-RNA-NT" ],
    "TCGA-BH-A0B3-RNA-NT-DNA-NT" => [ "TCGA-BH-A0B3-DNA-NT", "TCGA-BH-A0B3-RNA-NT" ],
    "TCGA-BH-A0B8-RNA-NT-DNA-NT" => [ "TCGA-BH-A0B8-DNA-NT", "TCGA-BH-A0B8-RNA-NT" ],
    "TCGA-BH-A0BJ-RNA-NT-DNA-NT" => [ "TCGA-BH-A0BJ-DNA-NT", "TCGA-BH-A0BJ-RNA-NT" ],
    "TCGA-BH-A0BM-RNA-NT-DNA-NT" => [ "TCGA-BH-A0BM-DNA-NT", "TCGA-BH-A0BM-RNA-NT" ],
    "TCGA-BH-A0C0-RNA-NT-DNA-NT" => [ "TCGA-BH-A0C0-DNA-NT", "TCGA-BH-A0C0-RNA-NT" ],
    "TCGA-BH-A0DK-RNA-NT-DNA-NT" => [ "TCGA-BH-A0DK-DNA-NT", "TCGA-BH-A0DK-RNA-NT" ],
    "TCGA-BH-A0DP-RNA-NT-DNA-NT" => [ "TCGA-BH-A0DP-DNA-NT", "TCGA-BH-A0DP-RNA-NT" ],
    "TCGA-BH-A0E0-RNA-NT-DNA-NT" => [ "TCGA-BH-A0E0-DNA-NT", "TCGA-BH-A0E0-RNA-NT" ],
    "TCGA-BH-A0H7-RNA-NT-DNA-NT" => [ "TCGA-BH-A0H7-DNA-NT", "TCGA-BH-A0H7-RNA-NT" ],

    "TCGA-A7-A0D9-DNA-TP-RNA-NT" => [ "TCGA-A7-A0D9-RNA-NT", "TCGA-A7-A0D9-DNA-TP" ],
    "TCGA-BH-A0B3-DNA-TP-RNA-NT" => [ "TCGA-BH-A0B3-RNA-NT", "TCGA-BH-A0B3-DNA-TP" ],
    "TCGA-BH-A0B8-DNA-TP-RNA-NT" => [ "TCGA-BH-A0B8-RNA-NT", "TCGA-BH-A0B8-DNA-TP" ],
    "TCGA-BH-A0BJ-DNA-TP-RNA-NT" => [ "TCGA-BH-A0BJ-RNA-NT", "TCGA-BH-A0BJ-DNA-TP" ],
    "TCGA-BH-A0BM-DNA-TP-RNA-NT" => [ "TCGA-BH-A0BM-RNA-NT", "TCGA-BH-A0BM-DNA-TP" ],
    "TCGA-BH-A0C0-DNA-TP-RNA-NT" => [ "TCGA-BH-A0C0-RNA-NT", "TCGA-BH-A0C0-DNA-TP" ],
    "TCGA-BH-A0DK-DNA-TP-RNA-NT" => [ "TCGA-BH-A0DK-RNA-NT", "TCGA-BH-A0DK-DNA-TP" ],
    "TCGA-BH-A0DP-DNA-TP-RNA-NT" => [ "TCGA-BH-A0DP-RNA-NT", "TCGA-BH-A0DP-DNA-TP" ],
    "TCGA-BH-A0E0-DNA-TP-RNA-NT" => [ "TCGA-BH-A0E0-RNA-NT", "TCGA-BH-A0E0-DNA-TP" ],
    "TCGA-BH-A0H7-DNA-TP-RNA-NT" => [ "TCGA-BH-A0H7-RNA-NT", "TCGA-BH-A0H7-DNA-TP" ],
  },
  dna_groups => {
    "TCGA-A7-A0D9-DNA-TP-NB" => [ "TCGA-A7-A0D9-DNA-NB", "TCGA-A7-A0D9-DNA-TP" ],
    "TCGA-BH-A0B3-DNA-TP-NB" => [ "TCGA-BH-A0B3-DNA-NB", "TCGA-BH-A0B3-DNA-TP" ],
    "TCGA-BH-A0B8-DNA-TP-NB" => [ "TCGA-BH-A0B8-DNA-NB", "TCGA-BH-A0B8-DNA-TP" ],
    "TCGA-BH-A0BJ-DNA-TP-NB" => [ "TCGA-BH-A0BJ-DNA-NB", "TCGA-BH-A0BJ-DNA-TP" ],
    "TCGA-BH-A0BM-DNA-TP-NB" => [ "TCGA-BH-A0BM-DNA-NB", "TCGA-BH-A0BM-DNA-TP" ],
    "TCGA-BH-A0C0-DNA-TP-NB" => [ "TCGA-BH-A0C0-DNA-NB", "TCGA-BH-A0C0-DNA-TP" ],
    "TCGA-BH-A0DK-DNA-TP-NB" => [ "TCGA-BH-A0DK-DNA-NB", "TCGA-BH-A0DK-DNA-TP" ],
    "TCGA-BH-A0DP-DNA-TP-NB" => [ "TCGA-BH-A0DP-DNA-NB", "TCGA-BH-A0DP-DNA-TP" ],
    "TCGA-BH-A0E0-DNA-TP-NB" => [ "TCGA-BH-A0E0-DNA-NB", "TCGA-BH-A0E0-DNA-TP" ],
    "TCGA-BH-A0H7-DNA-TP-NB" => [ "TCGA-BH-A0H7-DNA-NB", "TCGA-BH-A0H7-DNA-TP" ],

    #
    #    "TCGA-A7-A0D9-DNA-TP-NT" => [ "TCGA-A7-A0D9-DNA-NT", "TCGA-A7-A0D9-DNA-TP" ],
    #    "TCGA-BH-A0B3-DNA-TP-NT" => [ "TCGA-BH-A0B3-DNA-NT", "TCGA-BH-A0B3-DNA-TP" ],
    #    "TCGA-BH-A0B8-DNA-TP-NT" => [ "TCGA-BH-A0B8-DNA-NT", "TCGA-BH-A0B8-DNA-TP" ],
    #    "TCGA-BH-A0BJ-DNA-TP-NT" => [ "TCGA-BH-A0BJ-DNA-NT", "TCGA-BH-A0BJ-DNA-TP" ],
    #    "TCGA-BH-A0BM-DNA-TP-NT" => [ "TCGA-BH-A0BM-DNA-NT", "TCGA-BH-A0BM-DNA-TP" ],
    #    "TCGA-BH-A0C0-DNA-TP-NT" => [ "TCGA-BH-A0C0-DNA-NT", "TCGA-BH-A0C0-DNA-TP" ],
    #    "TCGA-BH-A0DK-DNA-TP-NT" => [ "TCGA-BH-A0DK-DNA-NT", "TCGA-BH-A0DK-DNA-TP" ],
    #    "TCGA-BH-A0DP-DNA-TP-NT" => [ "TCGA-BH-A0DP-DNA-NT", "TCGA-BH-A0DP-DNA-TP" ],
    #    "TCGA-BH-A0E0-DNA-TP-NT" => [ "TCGA-BH-A0E0-DNA-NT", "TCGA-BH-A0E0-DNA-TP" ],
    #    "TCGA-BH-A0H7-DNA-TP-NT" => [ "TCGA-BH-A0H7-DNA-NT", "TCGA-BH-A0H7-DNA-TP" ],
    #
    #    "TCGA-A7-A0D9-DNA-NT-NB" => [ "TCGA-A7-A0D9-DNA-NB", "TCGA-A7-A0D9-DNA-NT" ],
    #    "TCGA-BH-A0B3-DNA-NT-NB" => [ "TCGA-BH-A0B3-DNA-NB", "TCGA-BH-A0B3-DNA-NT" ],
    #    "TCGA-BH-A0B8-DNA-NT-NB" => [ "TCGA-BH-A0B8-DNA-NB", "TCGA-BH-A0B8-DNA-NT" ],
    #    "TCGA-BH-A0BJ-DNA-NT-NB" => [ "TCGA-BH-A0BJ-DNA-NB", "TCGA-BH-A0BJ-DNA-NT" ],
    #    "TCGA-BH-A0BM-DNA-NT-NB" => [ "TCGA-BH-A0BM-DNA-NB", "TCGA-BH-A0BM-DNA-NT" ],
    #    "TCGA-BH-A0C0-DNA-NT-NB" => [ "TCGA-BH-A0C0-DNA-NB", "TCGA-BH-A0C0-DNA-NT" ],
    #    "TCGA-BH-A0DK-DNA-NT-NB" => [ "TCGA-BH-A0DK-DNA-NB", "TCGA-BH-A0DK-DNA-NT" ],
    #    "TCGA-BH-A0DP-DNA-NT-NB" => [ "TCGA-BH-A0DP-DNA-NB", "TCGA-BH-A0DP-DNA-NT" ],
    #    "TCGA-BH-A0E0-DNA-NT-NB" => [ "TCGA-BH-A0E0-DNA-NB", "TCGA-BH-A0E0-DNA-NT" ],
    #    "TCGA-BH-A0H7-DNA-NT-NB" => [ "TCGA-BH-A0H7-DNA-NB", "TCGA-BH-A0H7-DNA-NT" ],
  },
  rna_groups => {
    "TCGA-A7-A0D9-RNA-TP-NT" => [ "TCGA-A7-A0D9-RNA-NT", "TCGA-A7-A0D9-RNA-TP" ],
    "TCGA-BH-A0B3-RNA-TP-NT" => [ "TCGA-BH-A0B3-RNA-NT", "TCGA-BH-A0B3-RNA-TP" ],
    "TCGA-BH-A0B8-RNA-TP-NT" => [ "TCGA-BH-A0B8-RNA-NT", "TCGA-BH-A0B8-RNA-TP" ],
    "TCGA-BH-A0BJ-RNA-TP-NT" => [ "TCGA-BH-A0BJ-RNA-NT", "TCGA-BH-A0BJ-RNA-TP" ],
    "TCGA-BH-A0BM-RNA-TP-NT" => [ "TCGA-BH-A0BM-RNA-NT", "TCGA-BH-A0BM-RNA-TP" ],
    "TCGA-BH-A0C0-RNA-TP-NT" => [ "TCGA-BH-A0C0-RNA-NT", "TCGA-BH-A0C0-RNA-TP" ],
    "TCGA-BH-A0DK-RNA-TP-NT" => [ "TCGA-BH-A0DK-RNA-NT", "TCGA-BH-A0DK-RNA-TP" ],
    "TCGA-BH-A0DP-RNA-TP-NT" => [ "TCGA-BH-A0DP-RNA-NT", "TCGA-BH-A0DP-RNA-TP" ],
    "TCGA-BH-A0E0-RNA-TP-NT" => [ "TCGA-BH-A0E0-RNA-NT", "TCGA-BH-A0E0-RNA-TP" ],
    "TCGA-BH-A0H7-RNA-TP-NT" => [ "TCGA-BH-A0H7-RNA-NT", "TCGA-BH-A0H7-RNA-TP" ],
  },
  sample_groups => {
    "TCGA-A7-A0D9" => [ "TCGA-A7-A0D9-DNA-NT", "TCGA-A7-A0D9-DNA-TP", "TCGA-A7-A0D9-RNA-NT", "TCGA-A7-A0D9-RNA-TP" ],
    "TCGA-BH-A0B3" => [ "TCGA-BH-A0B3-DNA-NT", "TCGA-BH-A0B3-DNA-TP", "TCGA-BH-A0B3-RNA-NT", "TCGA-BH-A0B3-RNA-TP" ],
    "TCGA-BH-A0B8" => [ "TCGA-BH-A0B8-DNA-NB", "TCGA-BH-A0B8-DNA-TP", "TCGA-BH-A0B8-RNA-NT", "TCGA-BH-A0B8-RNA-TP" ],
    "TCGA-BH-A0BJ" => [ "TCGA-BH-A0BJ-DNA-NB", "TCGA-BH-A0BJ-DNA-TP", "TCGA-BH-A0BJ-RNA-NT", "TCGA-BH-A0BJ-RNA-TP" ],
    "TCGA-BH-A0BM" => [ "TCGA-BH-A0BM-DNA-NB", "TCGA-BH-A0BM-DNA-TP", "TCGA-BH-A0BM-RNA-NT", "TCGA-BH-A0BM-RNA-TP" ],
    "TCGA-BH-A0C0" => [ "TCGA-BH-A0C0-DNA-NB", "TCGA-BH-A0C0-DNA-TP", "TCGA-BH-A0C0-RNA-NT", "TCGA-BH-A0C0-RNA-TP" ],
    "TCGA-BH-A0DK" => [ "TCGA-BH-A0DK-DNA-NB", "TCGA-BH-A0DK-DNA-TP", "TCGA-BH-A0DK-RNA-NT", "TCGA-BH-A0DK-RNA-TP" ],
    "TCGA-BH-A0DP" => [ "TCGA-BH-A0DP-DNA-NB", "TCGA-BH-A0DP-DNA-TP", "TCGA-BH-A0DP-RNA-NT", "TCGA-BH-A0DP-RNA-TP" ],
    "TCGA-BH-A0E0" => [ "TCGA-BH-A0E0-DNA-NB", "TCGA-BH-A0E0-DNA-TP", "TCGA-BH-A0E0-RNA-NT", "TCGA-BH-A0E0-RNA-TP" ],
    "TCGA-BH-A0H7" => [ "TCGA-BH-A0H7-DNA-NB", "TCGA-BH-A0H7-DNA-TP", "TCGA-BH-A0H7-RNA-NT", "TCGA-BH-A0H7-RNA-TP" ],
  }
};

my $preparation = {
  general     => { task_name => "preparation" },
  dna         => $tcga->{dna},
  rna         => $tcga->{rna},
  sortbam_dna => {
    class         => "Sortbam",
    perform       => 1,
    target_dir    => "${target_dir}/preparation_dna_sortname",
    option        => "",
    source_ref    => "dna",
    sort_by_query => 1,
    sh_direct     => 0,
    pbs           => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "20gb"
    },
  },
  bam2fastq_dna => {
    class               => "Bam2Fastq",
    perform             => 1,
    target_dir          => "${target_dir}/preparation_dna_bam2fastq",
    option              => "",
    source_ref          => "sortbam_dna",
    cqstools            => $cqstools,
    ispaired            => 1,
    sort_before_convert => 0,
    sort_thread         => 12,
    sh_direct           => 1,
    pbs                 => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "20gb"
    },
  },
  dna_bwa => {
    class                      => "Alignment::BWA",
    perform                    => 1,
    target_dir                 => "${target_dir}/preparation_dna_bwa",
    option                     => "-T 15 -t 8",
    fasta_file                 => $fasta_file_16569_MT,
    source_ref                 => "bam2fastq_dna",
    addOrReplaceReadGroups_jar => "/home/shengq1/local/bin/picard/AddOrReplaceReadGroups.jar",
    sh_direct                  => 0,
    pbs                        => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  dna_bwa_refine => {
    class              => "GATK::Refine",
    perform            => 1,
    target_dir         => "${target_dir}/preparation_dna_bwa_refine",
    option             => "-Xmx40g",
    fasta_file         => $fasta_file_16569_MT,
    source_ref         => "dna_bwa",
    thread_count       => 8,
    vcf_files          => [$snp_file_16569_MT],
    gatk_jar           => "/home/shengq1/local/bin/GATK/GenomeAnalysisTK.jar",
    markDuplicates_jar => "/home/shengq1/local/bin/picard/MarkDuplicates.jar",
    sh_direct          => 0,
    pbs                => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  bam2fastq_rna => {
    class               => "Bam2Fastq",
    perform             => 1,
    target_dir          => "${target_dir}/preparation_rna_bam2fastq",
    option              => "",
    source_ref          => "rna",
    cqstools            => $cqstools,
    ispaired            => 1,
    sort_before_convert => 0,
    sort_thread         => 12,
    sh_direct           => 1,
    pbs                 => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  tophat2_rna => {
    class                => "Alignment::Tophat2",
    perform              => 1,
    target_dir           => "${target_dir}/preparation_rna_tophat2",
    option               => "--segment-length 25 -r 0 -p 8",
    source_ref           => "bam2fastq_rna",
    bowtie2_index        => $bowtie2_index_16569_MT,
    transcript_gtf       => $transcript_gtf,
    transcript_gtf_index => $transcript_gtf_index,
    rename_bam           => 1,
    sh_direct            => 0,
    pbs                  => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  tophat2_rna_removeduplicates => {
    class              => "Picard::MarkDuplicates",
    perform            => 1,
    target_dir         => "${target_dir}/preparation_rna_tophat2_redup",
    option             => "-Xmx20g",
    source_ref         => "tophat2_rna",
    markDuplicates_jar => $picard . "MarkDuplicates.jar",
    sh_direct          => 0,
    pbs                => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "20gb"
    },
  },
  sortbam_tophat2_rna => {
    class         => "Samtools::Sort",
    perform       => 1,
    target_dir    => "${target_dir}/preparation_rna_tophat2_sortname",
    option        => "",
    source_ref    => "tophat2_rna",
    sort_by_query => 1,
    sh_direct     => 0,
    pbs           => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "20gb"
    },
  },
  htseqcount_rna => {
    class      => "Count::HTSeqCount",
    perform    => 1,
    target_dir => "${target_dir}/preparation_rna_htseqcount",
    option     => "",
    source_ref => "sortbam_tophat2_rna",
    gff_file   => $transcript_gtf,
    sh_direct  => 0,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  genetable_rna => {
    class         => "CQS::CQSDatatable",
    perform       => 1,
    target_dir    => "${target_dir}/preparation_rna_genetable",
    option        => "-p ENS --noheader -o TCGA_rna_gene.count",
    source_ref    => "htseqcount_rna",
    name_map_file => $hg19_map,
    cqs_tools     => $cqstools,
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
    target_dir => "${target_dir}/preparation_sequencetask",
    option     => "",
    source     => {
      preparation => [

        #data preparation
        # "sortbam_dna", "bam2fastq_dna", "bam2fastq_rna",
        "dna_bwa", "dna_bwa_refine", "tophat2_rna", "tophat2_rna_removeduplicates", "sortbam_tophat2_rna", "htseqcount_rna",
      ],
      table => [ "genetable_rna", ],
    },
    sh_direct => 0,
    pbs       => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  }
};

performConfig($preparation);

my $tcga_dna = {
  general     => { task_name => "tcga_dna" },
  files       => $tcga->{dna},
  groups      => $tcga->{dna_groups},
  fasta_file  => $fasta_file_16569_MT,
  cosmic_file => $cosmic_file_16569_MT,
  dbsnp_file  => $snp_file_16569_MT,
};

my $tcga_rna = {
  general     => { task_name => "tcga_rna" },
  files       => $tcga->{rna},
  groups      => $tcga->{rna_groups},
  fasta_file  => $fasta_file_16569_M,
  cosmic_file => $cosmic_file_16569_M,
  dbsnp_file  => $snp_file_16569_M,
};

my $realign_dna = {
  general     => { task_name => "realign_dna" },
  files       => $tcga->{dna},
  fasta_file  => $fasta_file_16569_MT,
  cosmic_file => $cosmic_file_16569_MT,
  dbsnp_file  => $snp_file_16569_MT,
  groups      => $tcga->{dna_groups},
};

my $realign_rna = {
  general     => { task_name => "realign_rna" },
  files       => $tcga->{rna},
  fasta_file  => $fasta_file_16569_MT,
  cosmic_file => $cosmic_file_16569_MT,
  dbsnp_file  => $snp_file_16569_MT,
  groups      => $tcga->{rna_groups},
};

my @cfgs = ( $tcga_dna, $tcga_rna, $realign_dna, $realign_rna );

for my $cfg (@cfgs) {
  my $task_name = $cfg->{general}{task_name};
  my $def       = {
    general => { task_name => $task_name },
    muTect  => {
      class        => "GATK::MuTect",
      perform      => 1,
      target_dir   => "${target_dir}/${task_name}_muTect",
      option       => "--min_qscore 20 --filter_reads_with_N_cigar",
      java_option  => "-Xmx40g",
      source_ref   => $cfg->{files},
      groups_ref   => $cfg->{groups},
      fasta_file   => $cfg->{fasta_file},
      cosmic_file  => $cfg->{cosmic_file},
      dbsnp_file   => $cfg->{dbsnp_file},
      bychromosome => 0,
      sh_direct    => 0,
      muTect_jar   => $mutect,
      pbs          => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "240",
        "mem"      => "40gb"
      },
    },
    annovar_muTect => {
      class      => "Annotation::Annovar",
      perform    => 0,
      target_dir => "${target_dir}/${task_name}_muTect",
      option     => $annovar_param,
      source_ref => [ "muTect", ".pass.vcf\$" ],
      annovar_db => $annovar_db,
      buildver   => "hg19",
      cqstools   => $cqstools,
      affy_file  => $affy_file,
      sh_direct  => 1,
      isvcf      => 1,
      pbs        => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "72",
        "mem"      => "10gb"
      },
    },
    varscan2 => {
      class           => "VarScan2::Somatic",
      perform         => 0,
      target_dir      => "${target_dir}/${task_name}_varscan2",
      option          => "--min-coverage 10",
      mpileup_options => "-A -q 20 -Q 20",
      java_option     => "-Xmx40g",
      source_ref      => $cfg->{files},
      groups_ref      => $cfg->{groups},
      fasta_file      => $cfg->{fasta_file},
      somatic_p_value => 0.05,
      sh_direct       => 0,
      VarScan2_jar    => $varscan2,
      pbs             => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "72",
        "mem"      => "40gb"
      },
    },
    annovar_varscan2 => {
      class      => "Annotation::Annovar",
      perform    => 1,
      target_dir => "${target_dir}/${task_name}_varscan2",
      option     => $annovar_param,
      source_ref => [ "varscan2", "snp.Somatic.hc.vcf\$" ],
      annovar_db => $annovar_db,
      buildver   => "hg19",
      cqstools   => $cqstools,
      affy_file  => $affy_file,
      sh_direct  => 1,
      isvcf      => 1,
      pbs        => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "72",
        "mem"      => "10gb"
      },
    },
    rsmc => {
      class            => "RSMC",
      perform          => 0,
      target_dir       => "${target_dir}/${task_name}_rsmc",
      option           => "-c 8",                              #thread mode
      source_type      => "BAM",                               #source_type can be BAM/Mpileup
      source_ref       => $cfg->{files},
      groups_ref       => $cfg->{groups},
      fasta_file       => $cfg->{fasta_file},
      annovar_buildver => "hg19",
      rnaediting_db    => $rnaediting_db,
      sh_direct        => 0,
      execute_file     => $rsmc,
      pbs              => {
        "email"    => $email,
        "nodes"    => "1:ppn=8",
        "walltime" => "72",
        "mem"      => "20gb"
      },
    },
    sequencetask => {
      class      => "CQS::SequenceTask",
      perform    => 1,
      target_dir => "${target_dir}/${task_name}_sequencetask",
      option     => "",
      source     => { individual => [ "muTect", "annovar_muTect", "varscan2", "annovar_varscan2", "rsmc" ] },
      sh_direct  => 0,
      pbs        => {
        "email"    => $email,
        "nodes"    => "1:ppn=8",
        "walltime" => "72",
        "mem"      => "40gb"
      },
    },
  };

#  performConfig($def);
}

1;

