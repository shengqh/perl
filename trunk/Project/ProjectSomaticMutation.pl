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
my $fasta_file    = "/data/cqs/shengq1/reference/hg19.16569/bowtie2_index_2.1.0/hg19_rCRS.fa";
my $cosmic_file   = "/data/cqs/shengq1/reference/cosmic/cosmic_v66_20130725.hg19.16569.vcf";
my $snp_file      = "/data/cqs/shengq1/reference/snp137/hg19.16569/dbsnp_137.b37.vcf";
my $bowtie2_index = "/data/cqs/shengq1/reference/hg19.16569/bowtie2_index_2.1.0/hg19_rCRS";
my $annovar_param = "-protocol refGene,snp137,cosmic64 -operation g,f,f --remove";
my $annovar_db    = "/scratch/cqs/shengq1/references/annovar/humandb/";

my $config = {
  general => { task_name => "somaticmutation" },
  files   => {
    "TCGA-A7-A0D9-DNA-NB" => ["/workspace/guoy1/GeneTorrent/lij17/download/BRCA/DNA_NB/b5c34def-778d-4bec-919e-e8917b9430bb/TCGA-A7-A0D9-10A-01W-A071-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0B3-DNA-NB" => ["/workspace/guoy1/GeneTorrent/lij17/download/BRCA/DNA_NB/745a64f1-5a91-405f-848e-4883462fc865/TCGA-BH-A0B3-10A-01W-A071-09_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0B8-DNA-NB" => ["/workspace/guoy1/GeneTorrent/lij17/download/BRCA/DNA_NB/6cf80125-08a7-43bb-b951-731b9f579932/TCGA-BH-A0B8-10A-01W-A071-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0BJ-DNA-NB" => ["/workspace/guoy1/GeneTorrent/lij17/download/BRCA/DNA_NB/ee0101c9-5d60-42de-a5ae-7ffa470c3097/TCGA-BH-A0BJ-10A-01W-A071-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0BM-DNA-NB" => ["/workspace/guoy1/GeneTorrent/lij17/download/BRCA/DNA_NB/c8ca19f5-746d-4108-b0be-16a41343eaa8/TCGA-BH-A0BM-10A-01W-A071-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0C0-DNA-NB" => ["/workspace/guoy1/GeneTorrent/lij17/download/BRCA/DNA_NB/424eb414-d7a3-4a34-99e2-29fb74d9dcb4/TCGA-BH-A0C0-10A-01W-A071-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0DK-DNA-NB" => ["/workspace/guoy1/GeneTorrent/lij17/download/BRCA/DNA_NB/ceb42cb6-59e1-4f5b-8f5c-c25ebf3cfd03/TCGA-BH-A0DK-10A-01W-A071-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0DP-DNA-NB" => ["/workspace/guoy1/GeneTorrent/lij17/download/BRCA/DNA_NB/2136bd75-82ef-4418-b905-3fbb7e5f77ef/TCGA-BH-A0DP-10A-01W-A071-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0E0-DNA-NB" => ["/workspace/guoy1/GeneTorrent/lij17/download/BRCA/DNA_NB/dd2c7052-bb6c-466b-b6a9-554f4958efa2/TCGA-BH-A0E0-10A-01W-A071-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0H7-DNA-NB" => ["/workspace/guoy1/GeneTorrent/lij17/download/BRCA/DNA_NB/eb05c499-7b27-4b8e-a7ec-692bc431eb76/TCGA-BH-A0H7-10A-01W-A071-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-A7-A0D9-DNA-NT" => ["/workspace/guoy1/GeneTorrent/lij17/download/BRCA/DNA_NT/fb943802-658b-4139-8539-c4a58f22129b/TCGA-A7-A0D9-11A-53W-A100-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0B3-DNA-NT" => ["/workspace/guoy1/GeneTorrent/lij17/download/BRCA/DNA_NT/e395bbf3-b0d3-450d-b6c8-e604a260d249/TCGA-BH-A0B3-11B-21W-A100-09_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0B8-DNA-NT" => ["/workspace/guoy1/GeneTorrent/lij17/download/BRCA/DNA_NT/2f7e3df4-0dc7-442f-a6a1-381c58f92629/TCGA-BH-A0B8-11A-41W-A10F-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0BJ-DNA-NT" => ["/workspace/guoy1/GeneTorrent/lij17/download/BRCA/DNA_NT/e4f46e6d-50e7-45e0-aae1-0c968762254e/TCGA-BH-A0BJ-11A-23W-A100-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0BM-DNA-NT" => ["/workspace/guoy1/GeneTorrent/lij17/download/BRCA/DNA_NT/babafabd-e085-40eb-9799-86a77a716317/TCGA-BH-A0BM-11A-12W-A10F-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0C0-DNA-NT" => ["/workspace/guoy1/GeneTorrent/lij17/download/BRCA/DNA_NT/b67574ff-0b97-488e-8e25-be9e5bfdd309/TCGA-BH-A0C0-11A-21W-A100-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0DK-DNA-NT" => ["/workspace/guoy1/GeneTorrent/lij17/download/BRCA/DNA_NT/773de6cd-45be-476f-a030-40d9db08d739/TCGA-BH-A0DK-11A-13W-A100-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0DP-DNA-NT" => ["/workspace/guoy1/GeneTorrent/lij17/download/BRCA/DNA_NT/bc943237-8872-4f9f-ba3b-9d0a560d6f10/TCGA-BH-A0DP-11A-12W-A100-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0E0-DNA-NT" => ["/workspace/guoy1/GeneTorrent/lij17/download/BRCA/DNA_NT/ff8d3854-beb2-4952-83f9-3fb805776e37/TCGA-BH-A0E0-11A-13W-A10F-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0H7-DNA-NT" => ["/workspace/guoy1/GeneTorrent/lij17/download/BRCA/DNA_NT/75d410cd-67b2-4cd1-b069-ccd9a23b8f45/TCGA-BH-A0H7-11A-13W-A100-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-A7-A0D9-DNA-TP" => ["/workspace/guoy1/GeneTorrent/lij17/download/BRCA/DNA_TP/4ef90332-fd3e-482d-9a03-a29e48ed3ded/TCGA-A7-A0D9-01A-31W-A071-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0B3-DNA-TP" => ["/workspace/guoy1/GeneTorrent/lij17/download/BRCA/DNA_TP/f43160ae-c1ea-4151-8e1c-c7a888250234/TCGA-BH-A0B3-01A-11W-A071-09_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0B8-DNA-TP" => ["/workspace/guoy1/GeneTorrent/lij17/download/BRCA/DNA_TP/06437b2b-97a7-47ce-bbed-8344099d59e2/TCGA-BH-A0B8-01A-21W-A071-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0BJ-DNA-TP" => ["/workspace/guoy1/GeneTorrent/lij17/download/BRCA/DNA_TP/e08cb2e2-e46b-48c3-990e-c40d54ea5928/TCGA-BH-A0BJ-01A-11W-A071-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0BM-DNA-TP" => ["/workspace/guoy1/GeneTorrent/lij17/download/BRCA/DNA_TP/1dc540e9-212e-493b-9448-31374323942e/TCGA-BH-A0BM-01A-11W-A071-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0C0-DNA-TP" => ["/workspace/guoy1/GeneTorrent/lij17/download/BRCA/DNA_TP/d2db8412-8cac-4d81-bc4f-b4f9abd5bba7/TCGA-BH-A0C0-01A-21W-A071-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0DK-DNA-TP" => ["/workspace/guoy1/GeneTorrent/lij17/download/BRCA/DNA_TP/8c749fc4-a755-41d6-ac87-010f30eab784/TCGA-BH-A0DK-01A-21W-A071-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0DP-DNA-TP" => ["/workspace/guoy1/GeneTorrent/lij17/download/BRCA/DNA_TP/b7c7d82d-2635-4ae4-a5fd-6e8cd4fb9ec8/TCGA-BH-A0DP-01A-21W-A071-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0E0-DNA-TP" => ["/workspace/guoy1/GeneTorrent/lij17/download/BRCA/DNA_TP/6051b0aa-5543-424e-91b1-a30209a2a267/TCGA-BH-A0E0-01A-11W-A071-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0H7-DNA-TP" => ["/workspace/guoy1/GeneTorrent/lij17/download/BRCA/DNA_TP/95cf0f38-99ff-4529-8830-39c0801de38d/TCGA-BH-A0H7-01A-13W-A071-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-A7-A0D9-RNA-NT" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/TCGA-A7-A0D9-RNA_NT_sorted.bam"],
    "TCGA-BH-A0B3-RNA-NT" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/TCGA-BH-A0B3-RNA_NT_sorted.bam"],
    "TCGA-BH-A0B8-RNA-NT" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/TCGA-BH-A0B8-RNA_NT_sorted.bam"],
    "TCGA-BH-A0BJ-RNA-NT" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/TCGA-BH-A0BJ-RNA_NT_sorted.bam"],
    "TCGA-BH-A0BM-RNA-NT" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/TCGA-BH-A0BM-RNA_NT_sorted.bam"],
    "TCGA-BH-A0C0-RNA-NT" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/TCGA-BH-A0C0-RNA_NT_sorted.bam"],
    "TCGA-BH-A0DK-RNA-NT" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/TCGA-BH-A0DK-RNA_NT_sorted.bam"],
    "TCGA-BH-A0DP-RNA-NT" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/TCGA-BH-A0DP-RNA_NT_sorted.bam"],
    "TCGA-BH-A0E0-RNA-NT" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/TCGA-BH-A0E0-RNA_NT_sorted.bam"],
    "TCGA-BH-A0H7-RNA-NT" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/TCGA-BH-A0H7-RNA_NT_sorted.bam"],
    "TCGA-A7-A0D9-RNA-TP" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/TCGA-A7-A0D9-RNA_TP_sorted.bam"],
    "TCGA-BH-A0B3-RNA-TP" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/TCGA-BH-A0B3-RNA_TP_sorted.bam"],
    "TCGA-BH-A0B8-RNA-TP" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/TCGA-BH-A0B8-RNA_TP_sorted.bam"],
    "TCGA-BH-A0BJ-RNA-TP" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/TCGA-BH-A0BJ-RNA_TP_sorted.bam"],
    "TCGA-BH-A0BM-RNA-TP" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/TCGA-BH-A0BM-RNA_TP_sorted.bam"],
    "TCGA-BH-A0C0-RNA-TP" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/TCGA-BH-A0C0-RNA_TP_sorted.bam"],
    "TCGA-BH-A0DK-RNA-TP" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/TCGA-BH-A0DK-RNA_TP_sorted.bam"],
    "TCGA-BH-A0DP-RNA-TP" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/TCGA-BH-A0DP-RNA_TP_sorted.bam"],
    "TCGA-BH-A0E0-RNA-TP" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/TCGA-BH-A0E0-RNA_TP_sorted.bam"],
    "TCGA-BH-A0H7-RNA-TP" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/TCGA-BH-A0H7-RNA_TP_sorted.bam"],
  },
  groups => {
    "TCGA-A7-A0D9-DNA-TP-NB"     => [ "TCGA-A7-A0D9-DNA-NB", "TCGA-A7-A0D9-DNA-TP" ],
    "TCGA-BH-A0B3-DNA-TP-NB"     => [ "TCGA-BH-A0B3-DNA-NB", "TCGA-BH-A0B3-DNA-TP" ],
    "TCGA-BH-A0B8-DNA-TP-NB"     => [ "TCGA-BH-A0B8-DNA-NB", "TCGA-BH-A0B8-DNA-TP" ],
    "TCGA-BH-A0BJ-DNA-TP-NB"     => [ "TCGA-BH-A0BJ-DNA-NB", "TCGA-BH-A0BJ-DNA-TP" ],
    "TCGA-BH-A0BM-DNA-TP-NB"     => [ "TCGA-BH-A0BM-DNA-NB", "TCGA-BH-A0BM-DNA-TP" ],
    "TCGA-BH-A0C0-DNA-TP-NB"     => [ "TCGA-BH-A0C0-DNA-NB", "TCGA-BH-A0C0-DNA-TP" ],
    "TCGA-BH-A0DK-DNA-TP-NB"     => [ "TCGA-BH-A0DK-DNA-NB", "TCGA-BH-A0DK-DNA-TP" ],
    "TCGA-BH-A0DP-DNA-TP-NB"     => [ "TCGA-BH-A0DP-DNA-NB", "TCGA-BH-A0DP-DNA-TP" ],
    "TCGA-BH-A0E0-DNA-TP-NB"     => [ "TCGA-BH-A0E0-DNA-NB", "TCGA-BH-A0E0-DNA-TP" ],
    "TCGA-BH-A0H7-DNA-TP-NB"     => [ "TCGA-BH-A0H7-DNA-NB", "TCGA-BH-A0H7-DNA-TP" ],
    "TCGA-A7-A0D9-DNA-TP-NT"     => [ "TCGA-A7-A0D9-DNA-NT", "TCGA-A7-A0D9-DNA-TP" ],
    "TCGA-BH-A0B3-DNA-TP-NT"     => [ "TCGA-BH-A0B3-DNA-NT", "TCGA-BH-A0B3-DNA-TP" ],
    "TCGA-BH-A0B8-DNA-TP-NT"     => [ "TCGA-BH-A0B8-DNA-NT", "TCGA-BH-A0B8-DNA-TP" ],
    "TCGA-BH-A0BJ-DNA-TP-NT"     => [ "TCGA-BH-A0BJ-DNA-NT", "TCGA-BH-A0BJ-DNA-TP" ],
    "TCGA-BH-A0BM-DNA-TP-NT"     => [ "TCGA-BH-A0BM-DNA-NT", "TCGA-BH-A0BM-DNA-TP" ],
    "TCGA-BH-A0C0-DNA-TP-NT"     => [ "TCGA-BH-A0C0-DNA-NT", "TCGA-BH-A0C0-DNA-TP" ],
    "TCGA-BH-A0DK-DNA-TP-NT"     => [ "TCGA-BH-A0DK-DNA-NT", "TCGA-BH-A0DK-DNA-TP" ],
    "TCGA-BH-A0DP-DNA-TP-NT"     => [ "TCGA-BH-A0DP-DNA-NT", "TCGA-BH-A0DP-DNA-TP" ],
    "TCGA-BH-A0E0-DNA-TP-NT"     => [ "TCGA-BH-A0E0-DNA-NT", "TCGA-BH-A0E0-DNA-TP" ],
    "TCGA-BH-A0H7-DNA-TP-NT"     => [ "TCGA-BH-A0H7-DNA-NT", "TCGA-BH-A0H7-DNA-TP" ],
    "TCGA-A7-A0D9-RNA-TP-NT"     => [ "TCGA-A7-A0D9-RNA-NT", "TCGA-A7-A0D9-RNA-TP" ],
    "TCGA-BH-A0B3-RNA-TP-NT"     => [ "TCGA-BH-A0B3-RNA-NT", "TCGA-BH-A0B3-RNA-TP" ],
    "TCGA-BH-A0B8-RNA-TP-NT"     => [ "TCGA-BH-A0B8-RNA-NT", "TCGA-BH-A0B8-RNA-TP" ],
    "TCGA-BH-A0BJ-RNA-TP-NT"     => [ "TCGA-BH-A0BJ-RNA-NT", "TCGA-BH-A0BJ-RNA-TP" ],
    "TCGA-BH-A0BM-RNA-TP-NT"     => [ "TCGA-BH-A0BM-RNA-NT", "TCGA-BH-A0BM-RNA-TP" ],
    "TCGA-BH-A0C0-RNA-TP-NT"     => [ "TCGA-BH-A0C0-RNA-NT", "TCGA-BH-A0C0-RNA-TP" ],
    "TCGA-BH-A0DK-RNA-TP-NT"     => [ "TCGA-BH-A0DK-RNA-NT", "TCGA-BH-A0DK-RNA-TP" ],
    "TCGA-BH-A0DP-RNA-TP-NT"     => [ "TCGA-BH-A0DP-RNA-NT", "TCGA-BH-A0DP-RNA-TP" ],
    "TCGA-BH-A0E0-RNA-TP-NT"     => [ "TCGA-BH-A0E0-RNA-NT", "TCGA-BH-A0E0-RNA-TP" ],
    "TCGA-BH-A0H7-RNA-TP-NT"     => [ "TCGA-BH-A0H7-RNA-NT", "TCGA-BH-A0H7-RNA-TP" ],
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
  },
  rsmc_nofilter => {
    class            => "RSMC",
    perform          => 1,
    target_dir       => "${target_dir}/rsmc_nofilter",
    option           => "-c 12",                                             #thread mode
    source_ref       => "files",
    groups_ref       => "groups",
    source_type      => "BAM",                                               #source_type can be BAM/Mpileup
    fasta_file       => $fasta_file,
    annovar_buildver => "hg19",
    rnaediting_db    => "/data/cqs/shengq1/reference/rnaediting/hg19.txt",
    sh_direct        => 1,
    execute_file     => "/home/shengq1/rsmc/rsmc.exe",
    pbs              => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "20gb"
    },
  },
  muTect => {
    class        => "MuTect",
    perform      => 1,
    target_dir   => "${target_dir}/muTect",
    option       => "",
    java_option  => "-Xmx40g",
    source_ref   => "files",
    groups_ref   => "groups",
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
  annovar_mutect => {
    class      => "Annovar",
    perform    => 0,
    target_dir => "${target_dir}/muTect",
    option     => $annovar_param,
    source_ref => [ "muTect", ".pass.vcf\$" ],
    annovar_db => $annovar_db,
    buildver   => "hg19",
    cqstools   => $cqstools,
    affy_file  => "/data/cqs/shengq1/reference/affy/HG-U133_Plus_2.na33.annot.csv",
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
    class           => "VarScan2",
    perform         => 1,
    target_dir      => "${target_dir}/varscan2",
    option          => "--min-coverage 10",
    mpileup_options => "-q 20",
    java_option     => "-Xmx40g",
    source_ref      => "files",
    groups_ref      => "groups",
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
  annovar_varscan2 => {
    class      => "Annovar",
    perform    => 0,
    target_dir => "${target_dir}/varscan2",
    option     => $annovar_param,
    source_ref => [ "varscan2", ".Somatic.hc\$" ],
    annovar_db => $annovar_db,
    buildver   => "hg19",
    cqstools   => $cqstools,
    affy_file  => "/data/cqs/shengq1/reference/affy/HG-U133_Plus_2.na33.annot.csv",
    sh_direct  => 1,
    isvcf      => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "10gb"
    },
  },
};

performConfig($config);

1;
