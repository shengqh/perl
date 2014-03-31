#!/usr/bin/perl
use strict;
use warnings;

use CQS::FileUtils;
use CQS::ClassFactory;

my $target_dir = create_directory_or_die("/scratch/cqs/shengq1/somaticmutation_comparison");

my $email    = "quanhu.sheng\@vanderbilt.edu";
my $cqstools = "/home/shengq1/cqstools/CQS.Tools.exe";
my $samtools = "/home/shengq1/local/bin/samtools/samtools";

##hg19.16569###
my $fasta_file_16569    = "/data/cqs/shengq1/reference/hg19.16569/bowtie2_index_2.1.0/hg19_rCRS.fa";
my $cosmic_file_16569   = "/data/cqs/shengq1/reference/cosmic/cosmic_v66_20130725.hg19.16569.vcf";
my $snp_file_16569      = "/data/cqs/shengq1/reference/snp137/hg19.16569/dbsnp_137.b37.vcf";
my $bowtie2_index_16569 = "/data/cqs/shengq1/reference/hg19.16569/bowtie2_index_2.1.0/hg19_rCRS";

##hg19.16571###
my $fasta_file    = "/data/cqs/guoy1/reference/hg19/bwa_index_0.7.4/hg19_chr.fa";
my $cosmic_file   = "/data/cqs/shengq1/reference/cosmic/cosmic_v65_28052013.hg19.16571.vcf";
my $snp_file      = "/data/cqs/shengq1/reference/snp137/hg19.16571/00-All.vcf";
my $bowtie2_index = "/data/cqs/guoy1/reference/hg19/bowtie2_index/hg19";

my $transcript_gtf       = "/scratch/cqs/shengq1/references/hg19/Homo_sapiens.GRCh37.73.gtf";
my $transcript_gtf_index = "/scratch/cqs/shengq1/gtfindex/hg19_GRCh37_73";
my $hg19_map             = "/scratch/cqs/shengq1/references/hg19/Homo_sapiens.GRCh37.73.map";

my $annovar_param = "-protocol refGene,snp137,cosmic64 -operation g,f,f --remove";
my $annovar_db    = "/scratch/cqs/shengq1/references/annovar/humandb/";

my $config = {
  general => { task_name => "somaticmutation" },
  dna     => {
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
  sortbam_dna => {
    class         => "Sortbam",
    perform       => 0,
    target_dir    => "${target_dir}/dna_sortname",
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
    perform             => 0,
    target_dir          => "${target_dir}/dna_bam2fastq",
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
    class      => "BWA",
    perform    => 0,
    target_dir => "${target_dir}/dna_bwa",
    option     => "-T 15 -t 8",
    fasta_file => $fasta_file,
    source_ref => "bam2fastq_dna",
    sh_direct  => 0,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  dna_bwa_refine => {
    class              => "GATKRefine",
    perform            => 0,
    target_dir         => "${target_dir}/dna_bwa_refine",
    option             => "-Xmx40g",
    fasta_file         => $fasta_file,
    source_ref         => "dna_bwa",
    thread_count       => 8,
    vcf_files          => [$snp_file],
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
  dna_bowtie2 => {
    class         => "Bowtie2",
    perform       => 0,
    target_dir    => "${target_dir}/dna_bowtie2",
    option        => "-p 8",
    source_ref    => "bam2fastq_dna",
    bowtie2_index => $bowtie2_index,
    samonly       => 0,
    sh_direct     => 0,
    pbs           => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "24",
      "mem"      => "20gb"
    },
  },
  dna_bowtie2_refine => {
    class              => "GATKRefine",
    perform            => 0,
    target_dir         => "${target_dir}/dna_bowtie2_refine",
    option             => "-Xmx40g",
    fasta_file         => $fasta_file,
    source_ref         => "dna_bowtie2",
    thread_count       => 8,
    vcf_files          => [$snp_file],
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
  sortbam_rna => {
    class         => "Sortbam",
    perform       => 0,
    target_dir    => "${target_dir}/rna_sortname",
    option        => "",
    source_ref    => "rna",
    sort_by_query => 1,
    sh_direct     => 0,
    pbs           => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "20gb"
    },
  },
  bam2fastq_rna => {
    class               => "Bam2Fastq",
    perform             => 0,
    target_dir          => "${target_dir}/rna_bam2fastq",
    option              => "",
    source_ref          => "sortbam_rna",
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
    class                => "Tophat2",
    perform              => 0,
    target_dir           => "${target_dir}/rna_tophat2",
    option               => "--segment-length 25 -r 0 -p 8",
    source_ref           => "bam2fastq_rna",
    bowtie2_index        => $bowtie2_index,
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
    perform            => 0,
    target_dir         => "${target_dir}/rna_tophat2_redup",
    option             => "-Xmx20g",
    source_ref         => "tophat2_rna",
    markDuplicates_jar => "/home/shengq1/local/bin/picard/MarkDuplicates.jar",
    sh_direct          => 0,
    pbs                => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "20gb"
    },
  },
  sortbam_tophat2_rna => {
    class         => "Sortbam",
    perform       => 0,
    target_dir    => "${target_dir}/rna_tophat2_sortname",
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
    class      => "HTSeqCount",
    perform    => 0,
    target_dir => "${target_dir}/rna_htseqcount",
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
    class         => "CQSDatatable",
    perform       => 0,
    target_dir    => "${target_dir}/rna_genetable",
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
  muTect => {
    class        => "GATK::MuTect",
    perform      => 0,
    target_dir   => "${target_dir}/all_muTect",
    option       => "--min_qscore 20",
    java_option  => "-Xmx40g",
    source_ref   => [ "dna_bwa_refine", "tophat2_rna_removeduplicates" ],
    groups_ref   => [ "dna_groups", "rna_groups" ],
    fasta_file   => $fasta_file,
    cosmic_file  => $cosmic_file,
    dbsnp_file   => $snp_file,
    bychromosome => 0,
    sh_direct    => 0,
    muTect_jar   => "/home/shengq1/local/bin/muTect-1.1.4.jar",
    pbs          => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "240",
      "mem"      => "40gb"
    },
  },
  varscan2 => {
    class           => "VarScan2::Somatic",
    perform         => 1,
    target_dir      => "${target_dir}/all_varscan2",
    option          => "--min-coverage 10",
    mpileup_options => "-A -q 20 -Q 20",
    java_option     => "-Xmx40g",
    source_ref      => [ "dna_bwa_refine", "tophat2_rna_removeduplicates" ],
    groups_ref      => [ "dna_groups", "rna_groups" ],
    fasta_file      => $fasta_file,
    somatic_p_value => 0.05,
    sh_direct       => 0,
    VarScan2_jar    => "/home/shengq1/local/bin/VarScan.v2.3.5.jar",
    pbs             => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  rsmc => {
    class            => "RSMC",
    perform          => 1,
    target_dir       => "${target_dir}/all_rsmc",
    option           => "",                                                     #thread mode
    source_ref       => [ "dna_bwa_refine", "tophat2_rna_removeduplicates" ],
    groups_ref       => [ "dna_groups", "rna_groups" ],
    source_type      => "BAM",                                                  #source_type can be BAM/Mpileup
    fasta_file       => $fasta_file,
    annovar_buildver => "hg19",
    rnaediting_db    => "/data/cqs/shengq1/reference/rnaediting/hg19.txt",
    sh_direct        => 0,
    execute_file     => "/home/shengq1/rsmc/rsmc.exe",
    pbs              => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "20gb"
    },
  },
  muTect_16569 => {
    class        => "GATK::MuTect",
    perform      => 0,
    target_dir   => "${target_dir}/16569_muTect",
    option       => "--min_qscore 20",
    java_option  => "-Xmx40g",
    source_ref   => [ "dna", "rna" ],
    groups_ref   => [ "dna_groups", "rna_groups" ],
    fasta_file   => $fasta_file_16569,
    cosmic_file  => $cosmic_file_16569,
    dbsnp_file   => $snp_file_16569,
    bychromosome => 0,
    sh_direct    => 0,
    muTect_jar   => "/home/shengq1/local/bin/muTect-1.1.4.jar",
    pbs          => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "240",
      "mem"      => "40gb"
    },
  },
  varscan2_16569 => {
    class           => "VarScan2::Somatic",
    perform         => 1,
    target_dir      => "${target_dir}/16569_varscan2",
    option          => "--min-coverage 10",
    mpileup_options => "-A -q 20 -Q 20",
    java_option     => "-Xmx40g",
    source_ref      => [ "dna", "rna" ],
    groups_ref      => [ "dna_groups", "rna_groups" ],
    fasta_file      => $fasta_file_16569,
    somatic_p_value => 0.05,
    sh_direct       => 0,
    VarScan2_jar    => "/home/shengq1/local/bin/VarScan.v2.3.5.jar",
    pbs             => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  rsmc_16569 => {
    class            => "RSMC",
    perform          => 1,
    target_dir       => "${target_dir}/16569_rsmc",
    option           => "",
    source_ref       => [ "dna", "rna" ],
    groups_ref       => ["rna_groups"],
    source_type      => "BAM",                                               #source_type can be BAM/Mpileup
    fasta_file       => $fasta_file_16569,
    annovar_buildver => "hg19",
    rnaediting_db    => "/data/cqs/shengq1/reference/rnaediting/hg19.txt",
    sh_direct        => 0,
    execute_file     => "/home/shengq1/rsmc/rsmc.exe",
    pbs              => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "20gb"
    },
  },
};

performConfig($config);
#performTask($config, "muTect2");


1;
