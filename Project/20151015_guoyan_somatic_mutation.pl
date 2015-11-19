#!/usr/bin/perl
use strict;
use warnings;

use CQS::FileUtils;
use CQS::ClassFactory;

my $target_dir = create_directory_or_die("/scratch/cqs/shengq1/variants/20151015_guoyan_somatic_mutation");

my $email             = "quanhu.sheng\@vanderbilt.edu";
my $cqstools          = "/home/shengq1/cqstools/CQS.Tools.exe";
my $samtools          = "/home/shengq1/local/bin/samtools/samtools";
my $glmvc             = "/home/shengq1/glmvc/glmvc.exe";
my $mutect            = "/home/shengq1/local/bin/mutect-1.1.7.jar";
my $varscan2          = "/home/shengq1/local/bin/VarScan.v2.3.7.jar";
my $gatk_jar          = "/home/shengq1/local/bin/GATK/GenomeAnalysisTK.jar";
my $picard_jar        = "/scratch/cqs/shengq1/local/bin/picard/picard.jar";
my $hg19_map          = "/data/cqs/shengq1/reference/ensembl_gtf/Homo_sapiens.GRCh37.75.map";
my $affy_file         = "/data/cqs/shengq1/reference/affy/HG-U133_Plus_2.na33.annot.csv";
my $rnaediting_db     = "/data/cqs/shengq1/reference/rnaediting/hg19.txt";
my $annovar_protocol  = "refGene,snp138,cosmic70";
my $annovar_operation = "g,f,f";
my $annovar_param     = "-protocol ${annovar_protocol} -operation ${annovar_operation} --remove";
my $annovar_db        = "/scratch/cqs/shengq1/references/annovar/humandb/";

my $qc3_perl = "/scratch/cqs/shengq1/local/bin/qc3/qc3.pl";

##hg19.16569.MT###
my $fasta_file_16569_MT  = "/scratch/cqs/shengq1/references/hg19_16569_MT/hg19_16569_MT.fa";
my $cosmic_file_16569_MT = "/scratch/cqs/shengq1/references/cosmic/cosmic_v71_hg19_16569_MT.vcf";
my $snp_file_16569_MT    = "/scratch/cqs/shengq1/references/dbsnp/human_GRCh37_v142_16569_MT.vcf";
my $gtf_file_16569_MT    = "/scratch/cqs/shengq1/references/ensembl_gtf/v75/Homo_sapiens.GRCh37.75.MT.gtf";

##hg19.16569.M###
my $fasta_file_16569_M  = "/scratch/cqs/shengq1/references/hg19_16569_M/hg19_16569_M.fa";
my $bwa_index_16569_M   = "/scratch/cqs/shengq1/references/hg19_16569_M/bwa_index_0.7.12/hg19_16569_M.fa";
my $cosmic_file_16569_M = "/scratch/cqs/shengq1/references/cosmic/cosmic_v71_hg19_16569_M.vcf";
my $snp_file_16569_M    = "/scratch/cqs/shengq1/references/dbsnp/human_GRCh37_v142_16569_M.vcf";
my $gtf_file_16569_M    = "/scratch/cqs/shengq1/references/ensembl_gtf/v75/Homo_sapiens.GRCh37.75.M.gtf";
my $star_index_16569_M  = "/scratch/cqs/shengq1/references/hg19_16569_M/STAR_index_v37.75_2.4.0j_sjdb50";

##hg19.tcga.dna###
my $fasta_file_tcga_dna = "/scratch/cqs/shengq1/references/tcga/dna/GRCh37-lite.fa";

my $tcga = {
  dna => {
    "TCGA-A7-A0D9-DNA-NB" => ["/scratch/cqs/shengq1/variants/tcga/bam/DNA_NB/TCGA-A7-A0D9-10A-01W-A071-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0B3-DNA-NB" => ["/scratch/cqs/shengq1/variants/tcga/bam/DNA_NB/TCGA-BH-A0B3-10A-01W-A071-09_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0B8-DNA-NB" => ["/scratch/cqs/shengq1/variants/tcga/bam/DNA_NB/TCGA-BH-A0B8-10A-01W-A071-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0BJ-DNA-NB" => ["/scratch/cqs/shengq1/variants/tcga/bam/DNA_NB/TCGA-BH-A0BJ-10A-01W-A071-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0BM-DNA-NB" => ["/scratch/cqs/shengq1/variants/tcga/bam/DNA_NB/TCGA-BH-A0BM-10A-01W-A071-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0C0-DNA-NB" => ["/scratch/cqs/shengq1/variants/tcga/bam/DNA_NB/TCGA-BH-A0C0-10A-01W-A071-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0DK-DNA-NB" => ["/scratch/cqs/shengq1/variants/tcga/bam/DNA_NB/TCGA-BH-A0DK-10A-01W-A071-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0DP-DNA-NB" => ["/scratch/cqs/shengq1/variants/tcga/bam/DNA_NB/TCGA-BH-A0DP-10A-01W-A071-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0E0-DNA-NB" => ["/scratch/cqs/shengq1/variants/tcga/bam/DNA_NB/TCGA-BH-A0E0-10A-01W-A071-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0H7-DNA-NB" => ["/scratch/cqs/shengq1/variants/tcga/bam/DNA_NB/TCGA-BH-A0H7-10A-01W-A071-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-A7-A0D9-DNA-NT" => ["/scratch/cqs/shengq1/variants/tcga/bam/DNA_NT/TCGA-A7-A0D9-11A-53W-A100-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0B3-DNA-NT" => ["/scratch/cqs/shengq1/variants/tcga/bam/DNA_NT/TCGA-BH-A0B3-11B-21W-A100-09_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0B8-DNA-NT" => ["/scratch/cqs/shengq1/variants/tcga/bam/DNA_NT/TCGA-BH-A0B8-11A-41W-A10F-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0BJ-DNA-NT" => ["/scratch/cqs/shengq1/variants/tcga/bam/DNA_NT/TCGA-BH-A0BJ-11A-23W-A100-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0BM-DNA-NT" => ["/scratch/cqs/shengq1/variants/tcga/bam/DNA_NT/TCGA-BH-A0BM-11A-12W-A10F-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0C0-DNA-NT" => ["/scratch/cqs/shengq1/variants/tcga/bam/DNA_NT/TCGA-BH-A0C0-11A-21W-A100-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0DK-DNA-NT" => ["/scratch/cqs/shengq1/variants/tcga/bam/DNA_NT/TCGA-BH-A0DK-11A-13W-A100-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0DP-DNA-NT" => ["/scratch/cqs/shengq1/variants/tcga/bam/DNA_NT/TCGA-BH-A0DP-11A-12W-A100-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0E0-DNA-NT" => ["/scratch/cqs/shengq1/variants/tcga/bam/DNA_NT/TCGA-BH-A0E0-11A-13W-A10F-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0H7-DNA-NT" => ["/scratch/cqs/shengq1/variants/tcga/bam/DNA_NT/TCGA-BH-A0H7-11A-13W-A100-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-A7-A0D9-DNA-TP" => ["/scratch/cqs/shengq1/variants/tcga/bam/DNA_TP/TCGA-A7-A0D9-01A-31W-A071-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0B3-DNA-TP" => ["/scratch/cqs/shengq1/variants/tcga/bam/DNA_TP/TCGA-BH-A0B3-01A-11W-A071-09_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0B8-DNA-TP" => ["/scratch/cqs/shengq1/variants/tcga/bam/DNA_TP/TCGA-BH-A0B8-01A-21W-A071-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0BJ-DNA-TP" => ["/scratch/cqs/shengq1/variants/tcga/bam/DNA_TP/TCGA-BH-A0BJ-01A-11W-A071-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0BM-DNA-TP" => ["/scratch/cqs/shengq1/variants/tcga/bam/DNA_TP/TCGA-BH-A0BM-01A-11W-A071-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0C0-DNA-TP" => ["/scratch/cqs/shengq1/variants/tcga/bam/DNA_TP/TCGA-BH-A0C0-01A-21W-A071-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0DK-DNA-TP" => ["/scratch/cqs/shengq1/variants/tcga/bam/DNA_TP/TCGA-BH-A0DK-01A-21W-A071-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0DP-DNA-TP" => ["/scratch/cqs/shengq1/variants/tcga/bam/DNA_TP/TCGA-BH-A0DP-01A-21W-A071-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0E0-DNA-TP" => ["/scratch/cqs/shengq1/variants/tcga/bam/DNA_TP/TCGA-BH-A0E0-01A-11W-A071-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0H7-DNA-TP" => ["/scratch/cqs/shengq1/variants/tcga/bam/DNA_TP/TCGA-BH-A0H7-01A-13W-A071-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
  },
  rna => {
    "TCGA-A7-A0D9-RNA-NT" => ["/scratch/cqs/shengq1/variants/tcga/bam/RNA_NT/TCGA-A7-A0D9-RNA_NT_sorted.bam"],
    "TCGA-BH-A0B3-RNA-NT" => ["/scratch/cqs/shengq1/variants/tcga/bam/RNA_NT/TCGA-BH-A0B3-RNA_NT_sorted.bam"],
    "TCGA-BH-A0B8-RNA-NT" => ["/scratch/cqs/shengq1/variants/tcga/bam/RNA_NT/TCGA-BH-A0B8-RNA_NT_sorted.bam"],
    "TCGA-BH-A0BJ-RNA-NT" => ["/scratch/cqs/shengq1/variants/tcga/bam/RNA_NT/TCGA-BH-A0BJ-RNA_NT_sorted.bam"],
    "TCGA-BH-A0BM-RNA-NT" => ["/scratch/cqs/shengq1/variants/tcga/bam/RNA_NT/TCGA-BH-A0BM-RNA_NT_sorted.bam"],
    "TCGA-BH-A0C0-RNA-NT" => ["/scratch/cqs/shengq1/variants/tcga/bam/RNA_NT/TCGA-BH-A0C0-RNA_NT_sorted.bam"],
    "TCGA-BH-A0DK-RNA-NT" => ["/scratch/cqs/shengq1/variants/tcga/bam/RNA_NT/TCGA-BH-A0DK-RNA_NT_sorted.bam"],
    "TCGA-BH-A0DP-RNA-NT" => ["/scratch/cqs/shengq1/variants/tcga/bam/RNA_NT/TCGA-BH-A0DP-RNA_NT_sorted.bam"],
    "TCGA-BH-A0E0-RNA-NT" => ["/scratch/cqs/shengq1/variants/tcga/bam/RNA_NT/TCGA-BH-A0E0-RNA_NT_sorted.bam"],
    "TCGA-BH-A0H7-RNA-NT" => ["/scratch/cqs/shengq1/variants/tcga/bam/RNA_NT/TCGA-BH-A0H7-RNA_NT_sorted.bam"],
    "TCGA-A7-A0D9-RNA-TP" => ["/scratch/cqs/shengq1/variants/tcga/bam/RNA_TP/TCGA-A7-A0D9-RNA_TP_sorted.bam"],
    "TCGA-BH-A0B3-RNA-TP" => ["/scratch/cqs/shengq1/variants/tcga/bam/RNA_TP/TCGA-BH-A0B3-RNA_TP_sorted.bam"],
    "TCGA-BH-A0B8-RNA-TP" => ["/scratch/cqs/shengq1/variants/tcga/bam/RNA_TP/TCGA-BH-A0B8-RNA_TP_sorted.bam"],
    "TCGA-BH-A0BJ-RNA-TP" => ["/scratch/cqs/shengq1/variants/tcga/bam/RNA_TP/TCGA-BH-A0BJ-RNA_TP_sorted.bam"],
    "TCGA-BH-A0BM-RNA-TP" => ["/scratch/cqs/shengq1/variants/tcga/bam/RNA_TP/TCGA-BH-A0BM-RNA_TP_sorted.bam"],
    "TCGA-BH-A0C0-RNA-TP" => ["/scratch/cqs/shengq1/variants/tcga/bam/RNA_TP/TCGA-BH-A0C0-RNA_TP_sorted.bam"],
    "TCGA-BH-A0DK-RNA-TP" => ["/scratch/cqs/shengq1/variants/tcga/bam/RNA_TP/TCGA-BH-A0DK-RNA_TP_sorted.bam"],
    "TCGA-BH-A0DP-RNA-TP" => ["/scratch/cqs/shengq1/variants/tcga/bam/RNA_TP/TCGA-BH-A0DP-RNA_TP_sorted.bam"],
    "TCGA-BH-A0E0-RNA-TP" => ["/scratch/cqs/shengq1/variants/tcga/bam/RNA_TP/TCGA-BH-A0E0-RNA_TP_sorted.bam"],
    "TCGA-BH-A0H7-RNA-TP" => ["/scratch/cqs/shengq1/variants/tcga/bam/RNA_TP/TCGA-BH-A0H7-RNA_TP_sorted.bam"],
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
  dna_nb_groups => {

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
  },
  dna_nt_groups => {
    "TCGA-A7-A0D9-DNA-TP-NT" => [ "TCGA-A7-A0D9-DNA-NT", "TCGA-A7-A0D9-DNA-TP" ],
    "TCGA-BH-A0B3-DNA-TP-NT" => [ "TCGA-BH-A0B3-DNA-NT", "TCGA-BH-A0B3-DNA-TP" ],
    "TCGA-BH-A0B8-DNA-TP-NT" => [ "TCGA-BH-A0B8-DNA-NT", "TCGA-BH-A0B8-DNA-TP" ],
    "TCGA-BH-A0BJ-DNA-TP-NT" => [ "TCGA-BH-A0BJ-DNA-NT", "TCGA-BH-A0BJ-DNA-TP" ],
    "TCGA-BH-A0BM-DNA-TP-NT" => [ "TCGA-BH-A0BM-DNA-NT", "TCGA-BH-A0BM-DNA-TP" ],
    "TCGA-BH-A0C0-DNA-TP-NT" => [ "TCGA-BH-A0C0-DNA-NT", "TCGA-BH-A0C0-DNA-TP" ],
    "TCGA-BH-A0DK-DNA-TP-NT" => [ "TCGA-BH-A0DK-DNA-NT", "TCGA-BH-A0DK-DNA-TP" ],
    "TCGA-BH-A0DP-DNA-TP-NT" => [ "TCGA-BH-A0DP-DNA-NT", "TCGA-BH-A0DP-DNA-TP" ],
    "TCGA-BH-A0E0-DNA-TP-NT" => [ "TCGA-BH-A0E0-DNA-NT", "TCGA-BH-A0E0-DNA-TP" ],
    "TCGA-BH-A0H7-DNA-TP-NT" => [ "TCGA-BH-A0H7-DNA-NT", "TCGA-BH-A0H7-DNA-TP" ],

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

  tcga_rna_files => {
    "TCGA-A7-A0D9-RNA-TP-NT" => ["/gpfs21/scratch/cqs/shengq1/variants/tcga/tcga_validation/TCGA-A7-A0D9.tsv"],
    "TCGA-BH-A0B3-RNA-TP-NT" => ["/gpfs21/scratch/cqs/shengq1/variants/tcga/tcga_validation/TCGA-BH-A0B3.tsv"],
    "TCGA-BH-A0B8-RNA-TP-NT" => ["/gpfs21/scratch/cqs/shengq1/variants/tcga/tcga_validation/TCGA-BH-A0B8.tsv"],
    "TCGA-BH-A0BJ-RNA-TP-NT" => ["/gpfs21/scratch/cqs/shengq1/variants/tcga/tcga_validation/TCGA-BH-A0BJ.tsv"],
    "TCGA-BH-A0BM-RNA-TP-NT" => ["/gpfs21/scratch/cqs/shengq1/variants/tcga/tcga_validation/TCGA-BH-A0BM.tsv"],
    "TCGA-BH-A0C0-RNA-TP-NT" => ["/gpfs21/scratch/cqs/shengq1/variants/tcga/tcga_validation/TCGA-BH-A0C0.tsv"],
    "TCGA-BH-A0DK-RNA-TP-NT" => ["/gpfs21/scratch/cqs/shengq1/variants/tcga/tcga_validation/TCGA-BH-A0DK.tsv"],
    "TCGA-BH-A0DP-RNA-TP-NT" => ["/gpfs21/scratch/cqs/shengq1/variants/tcga/tcga_validation/TCGA-BH-A0DP.tsv"],
    "TCGA-BH-A0E0-RNA-TP-NT" => ["/gpfs21/scratch/cqs/shengq1/variants/tcga/tcga_validation/TCGA-BH-A0E0.tsv"],
    "TCGA-BH-A0H7-RNA-TP-NT" => ["/gpfs21/scratch/cqs/shengq1/variants/tcga/tcga_validation/TCGA-BH-A0H7.tsv"]
  },

  tcga_dna_files => {
    "TCGA-A7-A0D9-DNA-TP-NB" => ["/gpfs21/scratch/cqs/shengq1/variants/tcga/tcga_validation/TCGA-A7-A0D9.tsv"],
    "TCGA-BH-A0B3-DNA-TP-NB" => ["/gpfs21/scratch/cqs/shengq1/variants/tcga/tcga_validation/TCGA-BH-A0B3.tsv"],
    "TCGA-BH-A0B8-DNA-TP-NB" => ["/gpfs21/scratch/cqs/shengq1/variants/tcga/tcga_validation/TCGA-BH-A0B8.tsv"],
    "TCGA-BH-A0BJ-DNA-TP-NB" => ["/gpfs21/scratch/cqs/shengq1/variants/tcga/tcga_validation/TCGA-BH-A0BJ.tsv"],
    "TCGA-BH-A0BM-DNA-TP-NB" => ["/gpfs21/scratch/cqs/shengq1/variants/tcga/tcga_validation/TCGA-BH-A0BM.tsv"],
    "TCGA-BH-A0C0-DNA-TP-NB" => ["/gpfs21/scratch/cqs/shengq1/variants/tcga/tcga_validation/TCGA-BH-A0C0.tsv"],
    "TCGA-BH-A0DK-DNA-TP-NB" => ["/gpfs21/scratch/cqs/shengq1/variants/tcga/tcga_validation/TCGA-BH-A0DK.tsv"],
    "TCGA-BH-A0DP-DNA-TP-NB" => ["/gpfs21/scratch/cqs/shengq1/variants/tcga/tcga_validation/TCGA-BH-A0DP.tsv"],
    "TCGA-BH-A0E0-DNA-TP-NB" => ["/gpfs21/scratch/cqs/shengq1/variants/tcga/tcga_validation/TCGA-BH-A0E0.tsv"],
    "TCGA-BH-A0H7-DNA-TP-NB" => ["/gpfs21/scratch/cqs/shengq1/variants/tcga/tcga_validation/TCGA-BH-A0H7.tsv"]
  },

  tcga_dna_nt_files => {
    "TCGA-A7-A0D9-DNA-TP-NT" => ["/gpfs21/scratch/cqs/shengq1/variants/tcga/tcga_validation/TCGA-A7-A0D9.tsv"],
    "TCGA-BH-A0B3-DNA-TP-NT" => ["/gpfs21/scratch/cqs/shengq1/variants/tcga/tcga_validation/TCGA-BH-A0B3.tsv"],
    "TCGA-BH-A0B8-DNA-TP-NT" => ["/gpfs21/scratch/cqs/shengq1/variants/tcga/tcga_validation/TCGA-BH-A0B8.tsv"],
    "TCGA-BH-A0BJ-DNA-TP-NT" => ["/gpfs21/scratch/cqs/shengq1/variants/tcga/tcga_validation/TCGA-BH-A0BJ.tsv"],
    "TCGA-BH-A0BM-DNA-TP-NT" => ["/gpfs21/scratch/cqs/shengq1/variants/tcga/tcga_validation/TCGA-BH-A0BM.tsv"],
    "TCGA-BH-A0C0-DNA-TP-NT" => ["/gpfs21/scratch/cqs/shengq1/variants/tcga/tcga_validation/TCGA-BH-A0C0.tsv"],
    "TCGA-BH-A0DK-DNA-TP-NT" => ["/gpfs21/scratch/cqs/shengq1/variants/tcga/tcga_validation/TCGA-BH-A0DK.tsv"],
    "TCGA-BH-A0DP-DNA-TP-NT" => ["/gpfs21/scratch/cqs/shengq1/variants/tcga/tcga_validation/TCGA-BH-A0DP.tsv"],
    "TCGA-BH-A0E0-DNA-TP-NT" => ["/gpfs21/scratch/cqs/shengq1/variants/tcga/tcga_validation/TCGA-BH-A0E0.tsv"],
    "TCGA-BH-A0H7-DNA-TP-NT" => ["/gpfs21/scratch/cqs/shengq1/variants/tcga/tcga_validation/TCGA-BH-A0H7.tsv"]
  },

  tcga_files => {
    "TCGA-A7-A0D9" => ["/gpfs21/scratch/cqs/shengq1/variants/tcga/tcga_validation/TCGA-A7-A0D9.tsv"],
    "TCGA-BH-A0B3" => ["/gpfs21/scratch/cqs/shengq1/variants/tcga/tcga_validation/TCGA-BH-A0B3.tsv"],
    "TCGA-BH-A0B8" => ["/gpfs21/scratch/cqs/shengq1/variants/tcga/tcga_validation/TCGA-BH-A0B8.tsv"],
    "TCGA-BH-A0BJ" => ["/gpfs21/scratch/cqs/shengq1/variants/tcga/tcga_validation/TCGA-BH-A0BJ.tsv"],
    "TCGA-BH-A0BM" => ["/gpfs21/scratch/cqs/shengq1/variants/tcga/tcga_validation/TCGA-BH-A0BM.tsv"],
    "TCGA-BH-A0C0" => ["/gpfs21/scratch/cqs/shengq1/variants/tcga/tcga_validation/TCGA-BH-A0C0.tsv"],
    "TCGA-BH-A0DK" => ["/gpfs21/scratch/cqs/shengq1/variants/tcga/tcga_validation/TCGA-BH-A0DK.tsv"],
    "TCGA-BH-A0DP" => ["/gpfs21/scratch/cqs/shengq1/variants/tcga/tcga_validation/TCGA-BH-A0DP.tsv"],
    "TCGA-BH-A0E0" => ["/gpfs21/scratch/cqs/shengq1/variants/tcga/tcga_validation/TCGA-BH-A0E0.tsv"],
    "TCGA-BH-A0H7" => ["/gpfs21/scratch/cqs/shengq1/variants/tcga/tcga_validation/TCGA-BH-A0H7.tsv"],
  },

  sample_groups => {
    "TCGA-A7-A0D9" => [ "TCGA-A7-A0D9-DNA-NB", "TCGA-A7-A0D9-DNA-TP", "TCGA-A7-A0D9-RNA-NT", "TCGA-A7-A0D9-RNA-TP" ],
    "TCGA-BH-A0B3" => [ "TCGA-BH-A0B3-DNA-NB", "TCGA-BH-A0B3-DNA-TP", "TCGA-BH-A0B3-RNA-NT", "TCGA-BH-A0B3-RNA-TP" ],
    "TCGA-BH-A0B8" => [ "TCGA-BH-A0B8-DNA-NB", "TCGA-BH-A0B8-DNA-TP", "TCGA-BH-A0B8-RNA-NT", "TCGA-BH-A0B8-RNA-TP" ],
    "TCGA-BH-A0BJ" => [ "TCGA-BH-A0BJ-DNA-NB", "TCGA-BH-A0BJ-DNA-TP", "TCGA-BH-A0BJ-RNA-NT", "TCGA-BH-A0BJ-RNA-TP" ],
    "TCGA-BH-A0BM" => [ "TCGA-BH-A0BM-DNA-NB", "TCGA-BH-A0BM-DNA-TP", "TCGA-BH-A0BM-RNA-NT", "TCGA-BH-A0BM-RNA-TP" ],
    "TCGA-BH-A0C0" => [ "TCGA-BH-A0C0-DNA-NB", "TCGA-BH-A0C0-DNA-TP", "TCGA-BH-A0C0-RNA-NT", "TCGA-BH-A0C0-RNA-TP" ],
    "TCGA-BH-A0DK" => [ "TCGA-BH-A0DK-DNA-NB", "TCGA-BH-A0DK-DNA-TP", "TCGA-BH-A0DK-RNA-NT", "TCGA-BH-A0DK-RNA-TP" ],
    "TCGA-BH-A0DP" => [ "TCGA-BH-A0DP-DNA-NB", "TCGA-BH-A0DP-DNA-TP", "TCGA-BH-A0DP-RNA-NT", "TCGA-BH-A0DP-RNA-TP" ],
    "TCGA-BH-A0E0" => [ "TCGA-BH-A0E0-DNA-NB", "TCGA-BH-A0E0-DNA-TP", "TCGA-BH-A0E0-RNA-NT", "TCGA-BH-A0E0-RNA-TP" ],
    "TCGA-BH-A0H7" => [ "TCGA-BH-A0H7-DNA-NB", "TCGA-BH-A0H7-DNA-TP", "TCGA-BH-A0H7-RNA-NT", "TCGA-BH-A0H7-RNA-TP" ],
  },

  all_sample_groups => {
    "TCGA-A7-A0D9" => [ "TCGA-A7-A0D9-DNA-NB", "TCGA-A7-A0D9-DNA-NT", "TCGA-A7-A0D9-DNA-TP", "TCGA-A7-A0D9-RNA-NT", "TCGA-A7-A0D9-RNA-TP" ],
    "TCGA-BH-A0B3" => [ "TCGA-BH-A0B3-DNA-NB", "TCGA-BH-A0B3-DNA-NT", "TCGA-BH-A0B3-DNA-TP", "TCGA-BH-A0B3-RNA-NT", "TCGA-BH-A0B3-RNA-TP" ],
    "TCGA-BH-A0B8" => [ "TCGA-BH-A0B8-DNA-NB", "TCGA-BH-A0B8-DNA-NT", "TCGA-BH-A0B8-DNA-TP", "TCGA-BH-A0B8-RNA-NT", "TCGA-BH-A0B8-RNA-TP" ],
    "TCGA-BH-A0BJ" => [ "TCGA-BH-A0BJ-DNA-NB", "TCGA-BH-A0BJ-DNA-NT", "TCGA-BH-A0BJ-DNA-TP", "TCGA-BH-A0BJ-RNA-NT", "TCGA-BH-A0BJ-RNA-TP" ],
    "TCGA-BH-A0BM" => [ "TCGA-BH-A0BM-DNA-NB", "TCGA-BH-A0BM-DNA-NT", "TCGA-BH-A0BM-DNA-TP", "TCGA-BH-A0BM-RNA-NT", "TCGA-BH-A0BM-RNA-TP" ],
    "TCGA-BH-A0C0" => [ "TCGA-BH-A0C0-DNA-NB", "TCGA-BH-A0C0-DNA-NT", "TCGA-BH-A0C0-DNA-TP", "TCGA-BH-A0C0-RNA-NT", "TCGA-BH-A0C0-RNA-TP" ],
    "TCGA-BH-A0DK" => [ "TCGA-BH-A0DK-DNA-NB", "TCGA-BH-A0DK-DNA-NT", "TCGA-BH-A0DK-DNA-TP", "TCGA-BH-A0DK-RNA-NT", "TCGA-BH-A0DK-RNA-TP" ],
    "TCGA-BH-A0DP" => [ "TCGA-BH-A0DP-DNA-NB", "TCGA-BH-A0DP-DNA-NT", "TCGA-BH-A0DP-DNA-TP", "TCGA-BH-A0DP-RNA-NT", "TCGA-BH-A0DP-RNA-TP" ],
    "TCGA-BH-A0E0" => [ "TCGA-BH-A0E0-DNA-NB", "TCGA-BH-A0E0-DNA-NT", "TCGA-BH-A0E0-DNA-TP", "TCGA-BH-A0E0-RNA-NT", "TCGA-BH-A0E0-RNA-TP" ],
    "TCGA-BH-A0H7" => [ "TCGA-BH-A0H7-DNA-NB", "TCGA-BH-A0H7-DNA-NT", "TCGA-BH-A0H7-DNA-TP", "TCGA-BH-A0H7-RNA-NT", "TCGA-BH-A0H7-RNA-TP" ],
  }
};

my $preparation = {
  general    => { task_name => "preparation" },
  dna        => $tcga->{dna},
  rna        => $tcga->{rna},
  dna_refine => {
    class        => "GATK::Refine",
    perform      => 0,
    target_dir   => "${target_dir}/preparation_tcga_dna_refine",
    option       => "-Xmx40g",
    fasta_file   => $fasta_file_tcga_dna,
    source_ref   => "dna",
    thread_count => 8,
    vcf_files    => [$snp_file_16569_M],
    gatk_jar     => $gatk_jar,
    picard_jar   => $picard_jar,
    sorted       => 1,
    sh_direct    => 0,
    pbs          => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  rna_refine => {
    class              => "GATK::RNASeqRefine",
    perform            => 0,
    target_dir         => "${target_dir}/preparation_tcga_rna_refine",
    option             => "-Xmx40g",
    fasta_file         => $fasta_file_16569_M,
    source_ref         => "rna",
    vcf_files          => [$snp_file_16569_M],
    gatk_jar           => $gatk_jar,
    picard_jar         => $picard_jar,
    sorted             => 1,
    sh_direct          => 0,
    fixMisencodedQuals => 0,
    pbs                => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  depth => {
    class             => "Samtools::Depth",
    perform           => 0,
    target_dir        => "${target_dir}/depth",
    option            => "",
    minimum_depth     => 10,
    source_ref        => [ "dna_refine", "rna_refine" ],
    groups_config_ref => [ $tcga, "dna_groups", $tcga, "rna_groups", $tcga, "sample_groups" ],
    cqstools          => $cqstools,
    sh_direct         => 0,
    pbs               => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  GlmvcExtract => {
    class             => "Variants::GlmvcExtract",
    perform           => 0,
    target_dir        => "${target_dir}/tcga_glmvc_extract",
    option            => "",
    source_config_ref => [ $tcga, "tcga_files" ],
    bam_files_ref     => [ "dna_refine", "rna_refine" ],
    groups_ref        => $tcga->{all_sample_groups},
    fasta_file        => $fasta_file_tcga_dna,
    sh_direct         => 0,
    execute_file      => $glmvc,
    pbs               => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
};

#performConfig($preparation);
#performTask( $preparation, "depth" );

my $tcga_dna = {
  general => { task_name => "tcga_dna" },
  files_config_ref => [ $preparation, "dna_refine" ],
  groups           => $tcga->{dna_nb_groups},
  fasta_file       => $fasta_file_tcga_dna,
  cosmic_file      => $cosmic_file_16569_MT,
  dbsnp_file       => $snp_file_16569_MT,
  gtf_file         => $gtf_file_16569_MT,
  tcga_file        => $tcga->{tcga_dna_files},
  thread           => 8,
  glm_pvalue       => "0.1"
};

my $tcga_dna_nt = {
  general => { task_name => "tcga_dna_nt" },
  files_config_ref => [ $preparation, "dna_refine" ],
  groups           => $tcga->{dna_nt_groups},
  fasta_file       => $fasta_file_tcga_dna,
  cosmic_file      => $cosmic_file_16569_MT,
  dbsnp_file       => $snp_file_16569_MT,
  gtf_file         => $gtf_file_16569_MT,
  tcga_file        => $tcga->{tcga_dna_nt_files},
  thread           => 8,
  glm_pvalue       => "0.1"
};

my $tcga_rna = {
  general => { task_name => "tcga_rna" },
  files_config_ref => [ $preparation, "rna_refine" ],
  groups           => $tcga->{rna_groups},
  fasta_file       => $fasta_file_16569_M,
  cosmic_file      => $cosmic_file_16569_M,
  dbsnp_file       => $snp_file_16569_M,
  gtf_file         => $gtf_file_16569_M,
  tcga_file        => $tcga->{tcga_rna_files},
  thread           => 8,
  glm_pvalue       => "0.1"
};

#my @cfgs = ( $tcga_dna, $tcga_rna );
#my @nps = ( 0.01, 0.02 );
#my @gps = ( 0.01, 0.05, 0.1 );

my @cfgs = ($tcga_dna_nt);
my @nps  = (0.02);
my @gps  = (0.1);

for my $cfg (@cfgs) {
  my $task_name = $cfg->{general}{task_name};
  my $def       = {
    general => { task_name => $task_name },

    muTect => {
      class             => "GATK::MuTect",
      perform           => 0,
      target_dir        => "${target_dir}/${task_name}_muTect",
      option            => "--min_qscore 20 --filter_reads_with_N_cigar",
      java_option       => "-Xmx40g",
      source_config_ref => $cfg->{files_config_ref},
      groups_ref        => $cfg->{groups},
      fasta_file        => $cfg->{fasta_file},
      cosmic_file       => $cfg->{cosmic_file},
      dbsnp_file        => $cfg->{dbsnp_file},
      bychromosome      => 0,
      sh_direct         => 0,
      muTect_jar        => $mutect,
      pbs               => {
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
      class             => "VarScan2::Somatic",
      perform           => 0,
      target_dir        => "${target_dir}/${task_name}_varscan2",
      option            => "--min-coverage 10",
      mpileup_options   => "-A -q 20 -Q 20",
      java_option       => "-Xmx40g",
      source_config_ref => $cfg->{files_config_ref},
      groups_ref        => $cfg->{groups},
      fasta_file        => $cfg->{fasta_file},
      somatic_p_value   => 0.05,
      sh_direct         => 0,
      VarScan2_jar      => $varscan2,
      pbs               => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "72",
        "mem"      => "40gb"
      },
    },
    annovar_varscan2 => {
      class      => "Annotation::Annovar",
      perform    => 0,
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
    GlmvcValidation => {
      class             => "Variants::GlmvcValidate",
      perform           => 0,
      target_dir        => "${target_dir}/${task_name}_glmvc_validation",
      option            => "--glm_pvalue " . $cfg->{glm_pvalue},
      source_type       => "BAM",
      source_config_ref => $cfg->{files_config_ref},
      groups_ref        => $cfg->{groups},
      validation_files  => $cfg->{tcga_file},
      fasta_file        => $cfg->{fasta_file},
      annovar_buildver  => "hg19",
      rnaediting_db     => $rnaediting_db,
      distance_exon_gtf => $cfg->{gtf_file},
      sh_direct         => 1,
      execute_file      => $glmvc,
      pbs               => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "72",
        "mem"      => "40gb"
      },
    },
  };

  #my @individual = ( "muTect", "annovar_muTect", "varscan2", "annovar_varscan2", "GlmvcValidation" );
  my @individual = ( "GlmvcValidation", "GlmvcValidation_varscan2" );

  for my $np (@nps) {
    for my $gp (@gps) {
      $def->{"${task_name}_glmvc_np${np}_g${gp}"} = {
        class             => "Variants::GlmvcCall",
        perform           => 1,
        target_dir        => "${target_dir}/${task_name}_glmvc_np${np}_g${gp}_thread8",
        option            => "--max_normal_percentage ${np} --glm_pvalue ${gp}",
        source_type       => "BAM",
        source_config_ref => $cfg->{files_config_ref},
        groups_ref        => $cfg->{groups},
        fasta_file        => $cfg->{fasta_file},
        annovar_buildver  => "hg19",
        annovar_protocol  => $annovar_protocol,
        annovar_operation => $annovar_operation,
        annovar_db        => $annovar_db,
        rnaediting_db     => $rnaediting_db,
        distance_exon_gtf => $cfg->{gtf_file},
        sh_direct         => 0,
        execute_file      => $glmvc,
        pbs               => {
          "email"    => $email,
          "nodes"    => "1:ppn=" . $cfg->{thread},
          "walltime" => "72",
          "mem"      => "40gb"
        },
      };

      $def->{"${task_name}_glmvc_np${np}_g${gp}_annotation"} = {
        class             => "Variants::GlmvcAnnotation",
        perform           => 0,
        target_dir        => "${target_dir}/${task_name}_glmvc_np${np}_g${gp}",
        option            => "",
        source_ref        => "${task_name}_glmvc_np${np}_g${gp}",
        annovar_buildver  => "hg19",
        annovar_protocol  => $annovar_protocol,
        annovar_operation => $annovar_operation,
        annovar_db        => $annovar_db,
        rnaediting_db     => $rnaediting_db,
        distance_exon_gtf => $cfg->{gtf_file},
        sh_direct         => 1,
        execute_file      => $glmvc,
        pbs               => {
          "email"    => $email,
          "nodes"    => "1:ppn=8",
          "walltime" => "72",
          "mem"      => "40gb"
        },
      };

      push( @individual, "${task_name}_glmvc_np${np}_g${gp}" );
    }
  }

  $def->{sequencetask} = {
    class      => "CQS::SequenceTask",
    perform    => 0,
    target_dir => "${target_dir}/${task_name}_sequencetask",
    option     => "",
    source     => { one => \@individual },
    sh_direct  => 0,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  };

  performConfig($def);

  #performTask( $def, "annovar_muTect" );
  #performTask( $def, "annovar_varscan2" );
}

my $annotation = {
  general => { task_name => "ann" },
  annovar => {
    class      => "Annotation::Annovar",
    perform    => 0,
    target_dir => "${target_dir}/detected_annotation",
    option     => $annovar_param,
    source     => { "detected" => ["/gpfs21/scratch/cqs/shengq1/variants/20151015_guoyan_somatic_mutation/detected_annotation/data/detected_sites.tsv"], },
    annovar_db => $annovar_db,
    buildver   => "hg19",
    sh_direct  => 1,
    isvcf      => 0,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "10gb"
    },
  },
  GlmvcExtract => {
    class      => "Variants::GlmvcExtract",
    perform    => 0,
    target_dir => "${target_dir}/detected_extract",
    option     => "",
    source     => {
      "TCGA-A7-A0D9" => ["/gpfs21/scratch/cqs/shengq1/variants/20151015_guoyan_somatic_mutation/detected_extract/data/TCGA-A7-A0D9.tsv"],
      "TCGA-BH-A0B3" => ["/gpfs21/scratch/cqs/shengq1/variants/20151015_guoyan_somatic_mutation/detected_extract/data/TCGA-BH-A0B3.tsv"],
      "TCGA-BH-A0B8" => ["/gpfs21/scratch/cqs/shengq1/variants/20151015_guoyan_somatic_mutation/detected_extract/data/TCGA-BH-A0B8.tsv"],
      "TCGA-BH-A0BJ" => ["/gpfs21/scratch/cqs/shengq1/variants/20151015_guoyan_somatic_mutation/detected_extract/data/TCGA-BH-A0BJ.tsv"],
      "TCGA-BH-A0BM" => ["/gpfs21/scratch/cqs/shengq1/variants/20151015_guoyan_somatic_mutation/detected_extract/data/TCGA-BH-A0BM.tsv"],
      "TCGA-BH-A0C0" => ["/gpfs21/scratch/cqs/shengq1/variants/20151015_guoyan_somatic_mutation/detected_extract/data/TCGA-BH-A0C0.tsv"],
      "TCGA-BH-A0DK" => ["/gpfs21/scratch/cqs/shengq1/variants/20151015_guoyan_somatic_mutation/detected_extract/data/TCGA-BH-A0DK.tsv"],
      "TCGA-BH-A0DP" => ["/gpfs21/scratch/cqs/shengq1/variants/20151015_guoyan_somatic_mutation/detected_extract/data/TCGA-BH-A0DP.tsv"],
      "TCGA-BH-A0E0" => ["/gpfs21/scratch/cqs/shengq1/variants/20151015_guoyan_somatic_mutation/detected_extract/data/TCGA-BH-A0E0.tsv"],
      "TCGA-BH-A0H7" => ["/gpfs21/scratch/cqs/shengq1/variants/20151015_guoyan_somatic_mutation/detected_extract/data/TCGA-BH-A0H7.tsv"],
    },
    bam_files_config_ref => [ $preparation, "dna_refine", $preparation, "rna_refine" ],
    groups               => $tcga->{all_sample_groups},
    fasta_file           => $fasta_file_16569_MT,
    sh_direct            => 0,
    execute_file         => $glmvc,
    pbs                  => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
};

#performConfig($annotation);
1;

