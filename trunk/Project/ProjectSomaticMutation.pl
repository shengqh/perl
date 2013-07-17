#!/usr/bin/perl
use strict;
use warnings;

use CQS::QC;
use CQS::RNASeq;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::SomaticMutation;
use CQS::ClassFactory;

my $target_dir = create_directory_or_die("/scratch/cqs/shengq1/somaticmutation_2");

my $email = "quanhu.sheng\@vanderbilt.edu";

my $transcript_gtf = "/data/cqs/guoy1/reference/annotation2/hg19/Homo_sapiens.GRCh37.68.gtf";
my $fasta_file = "/data/cqs/shengq1/reference/hg19_rCRS/hg19_rCRS.fa";

my $config = {
  general => {
    path_file => "/home/shengq1/local/bin/path.txt",
    task_name => "somaticmutation"
  },
  fastqfiles => {
    "TCGA-A7-A0D9-NT" => ["/scratch/cqs/shengq1/somaticmutation/raw/TCGA-A7-A0D9-RNA_NT_sorted.fastq"],
    "TCGA-A7-A0D9-TP" => ["/scratch/cqs/shengq1/somaticmutation/raw/TCGA-A7-A0D9-RNA_TP_sorted.fastq"],
  },
  rnafiles => {
    "TCGA-A7-A0D9-NT" => ["/workspace/guoy1/GeneTorrent/lij17/processed/BRCA/TCGA-A7-A0D9/TCGA-A7-A0D9-RNA_NT_sorted.bam"],
    "TCGA-A7-A0D9-TP" => ["/workspace/guoy1/GeneTorrent/lij17/processed/BRCA/TCGA-A7-A0D9/TCGA-A7-A0D9-RNA_TP_sorted.bam"],
    "TCGA-BH-A0B3-NT" => ["/workspace/guoy1/GeneTorrent/lij17/processed/BRCA/TCGA-BH-A0B3/TCGA-BH-A0B3-RNA_NT_sorted.bam"],
    "TCGA-BH-A0B3-TP" => ["/workspace/guoy1/GeneTorrent/lij17/processed/BRCA/TCGA-BH-A0B3/TCGA-BH-A0B3-RNA_TP_sorted.bam"],
    "TCGA-BH-A0B8-NT" => ["/workspace/guoy1/GeneTorrent/lij17/processed/BRCA/TCGA-BH-A0B8/TCGA-BH-A0B8-RNA_NT_sorted.bam"],
    "TCGA-BH-A0B8-TP" => ["/workspace/guoy1/GeneTorrent/lij17/processed/BRCA/TCGA-BH-A0B8/TCGA-BH-A0B8-RNA_TP_sorted.bam"],
    "TCGA-BH-A0BJ-NT" => ["/workspace/guoy1/GeneTorrent/lij17/processed/BRCA/TCGA-BH-A0BJ/TCGA-BH-A0BJ-RNA_NT_sorted.bam"],
    "TCGA-BH-A0BJ-TP" => ["/workspace/guoy1/GeneTorrent/lij17/processed/BRCA/TCGA-BH-A0BJ/TCGA-BH-A0BJ-RNA_TP_sorted.bam"],
    "TCGA-BH-A0BM-NT" => ["/workspace/guoy1/GeneTorrent/lij17/processed/BRCA/TCGA-BH-A0BM/TCGA-BH-A0BM-RNA_NT_sorted.bam"],
    "TCGA-BH-A0BM-TP" => ["/workspace/guoy1/GeneTorrent/lij17/processed/BRCA/TCGA-BH-A0BM/TCGA-BH-A0BM-RNA_TP_sorted.bam"],
    "TCGA-BH-A0C0-NT" => ["/workspace/guoy1/GeneTorrent/lij17/processed/BRCA/TCGA-BH-A0C0/TCGA-BH-A0C0-RNA_NT_sorted.bam"],
    "TCGA-BH-A0C0-TP" => ["/workspace/guoy1/GeneTorrent/lij17/processed/BRCA/TCGA-BH-A0C0/TCGA-BH-A0C0-RNA_TP_sorted.bam"],
    "TCGA-BH-A0DK-NT" => ["/workspace/guoy1/GeneTorrent/lij17/processed/BRCA/TCGA-BH-A0DK/TCGA-BH-A0DK-RNA_NT_sorted.bam"],
    "TCGA-BH-A0DK-TP" => ["/workspace/guoy1/GeneTorrent/lij17/processed/BRCA/TCGA-BH-A0DK/TCGA-BH-A0DK-RNA_TP_sorted.bam"],
    "TCGA-BH-A0DP-NT" => ["/workspace/guoy1/GeneTorrent/lij17/processed/BRCA/TCGA-BH-A0DP/TCGA-BH-A0DP-RNA_NT_sorted.bam"],
    "TCGA-BH-A0DP-TP" => ["/workspace/guoy1/GeneTorrent/lij17/processed/BRCA/TCGA-BH-A0DP/TCGA-BH-A0DP-RNA_TP_sorted.bam"],
    "TCGA-BH-A0E0-NT" => ["/workspace/guoy1/GeneTorrent/lij17/processed/BRCA/TCGA-BH-A0E0/TCGA-BH-A0E0-RNA_NT_sorted.bam"],
    "TCGA-BH-A0E0-TP" => ["/workspace/guoy1/GeneTorrent/lij17/processed/BRCA/TCGA-BH-A0E0/TCGA-BH-A0E0-RNA_TP_sorted.bam"],
    "TCGA-BH-A0H7-NT" => ["/workspace/guoy1/GeneTorrent/lij17/processed/BRCA/TCGA-BH-A0H7/TCGA-BH-A0H7-RNA_NT_sorted.bam"],
    "TCGA-BH-A0H7-TP" => ["/workspace/guoy1/GeneTorrent/lij17/processed/BRCA/TCGA-BH-A0H7/TCGA-BH-A0H7-RNA_TP_sorted.bam"],
  },
  rnagroups => {
    "TCGA-A7-A0D9" => [ "TCGA-A7-A0D9-NT", "TCGA-A7-A0D9-TP" ],
    "TCGA-BH-A0B3" => [ "TCGA-BH-A0B3-NT", "TCGA-BH-A0B3-TP" ],
    "TCGA-BH-A0B8" => [ "TCGA-BH-A0B8-NT", "TCGA-BH-A0B8-TP" ],
    "TCGA-BH-A0BJ" => [ "TCGA-BH-A0BJ-NT", "TCGA-BH-A0BJ-TP" ],
    "TCGA-BH-A0BM" => [ "TCGA-BH-A0BM-NT", "TCGA-BH-A0BM-TP" ],
    "TCGA-BH-A0C0" => [ "TCGA-BH-A0C0-NT", "TCGA-BH-A0C0-TP" ],
    "TCGA-BH-A0DK" => [ "TCGA-BH-A0DK-NT", "TCGA-BH-A0DK-TP" ],
    "TCGA-BH-A0DP" => [ "TCGA-BH-A0DP-NT", "TCGA-BH-A0DP-TP" ],
    "TCGA-BH-A0E0" => [ "TCGA-BH-A0E0-NT", "TCGA-BH-A0E0-TP" ],
    "TCGA-BH-A0H7" => [ "TCGA-BH-A0H7-NT", "TCGA-BH-A0H7-TP" ],
  },
  tophat2 => {
    class      => "Tophat2",
    perform    => 0,
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
  rsmc => {
    target_dir       => "${target_dir}/rsmc",
    option           => "-c 8",                                                 #thread mode
    source_ref       => "bamfiles",
    source_type      => "bam",                                                  #source_type can be bam/mpileup
    fasta_file       => "/data/cqs/shengq1/reference/hg19_rCRS/hg19_rCRS.fa",
    annovar_buildver => "hg19",
    rnaediting_db    => "/data/cqs/shengq1/reference/rnaediting/hg19.txt",
    sh_direct        => 1,
    execute_file     => "/home/shengq1/rsmc/rsmc.exe",
    pbs              => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "10gb"
    },
  },
  muTect => {
    class       => "MuTect",
    perform     => 1,
    target_dir  => "${target_dir}/muTect",
    option      => "-nt 8",
    source_ref  => "rnafiles",
    groups_ref  => "rnagroups",
    java_option => "-Xmx40g",
    fasta_file  => $fasta_file,
    cosmic_file => "/data/cqs/shengq1/reference/cosmic/cosmic_v65_28052013.hg19.16569.vcf",
    dbsnp_file  => "/data/cqs/shengq1/reference/snp137/human_b37/dbsnp_137.b37.vcf",
    sh_direct   => 1,
    muTect_jar  => "/home/shengq1/local/bin/muTect-1.1.4.jar",
    pbs         => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  varscan2 => {
    class           => "VarScan2",
    perform         => 1,
    target_dir      => "${target_dir}/varscan2",
    option          => "",
    source_ref      => "rnafiles",
    groups_ref      => "rnagroups",
    fasta_file      => $fasta_file,
    min_coverage    => 10,
    somatic_p_value => 0.01,
    sh_direct       => 1,
    VarScan2_jar    => "/home/shengq1/local/bin/VarScan.v2.3.5.jar",
    pbs             => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
};

performConfig($config);
#rsmc( $config, "rsmc" );

1;
