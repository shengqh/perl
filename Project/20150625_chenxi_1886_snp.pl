#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;

my $task = "P1886";

my $target_dir = "/scratch/cqs/shengq1/rnaseq/20150625_chenxi_1886_snp";

#my $target_dir = "e:/temp";

my $transcript_gtf = "/scratch/cqs/shengq1/references/ensembl_gtf/v67/Mus_musculus.NCBIM37.67.gtf";
my $name_map_file  = "/scratch/cqs/shengq1/references/ensembl_gtf/v67/Mus_musculus.NCBIM37.67.map";
my $fasta_file     = "/data/cqs/guoy1/reference/mm9/mm9.fa";
my $cqstools       = "/home/shengq1/cqstools/CQS.Tools.exe";
my $dbsnp          = "/scratch/cqs/shengq1/references/dbsnp/mouse/mm9_dbsnp128_sorted.bed";
my $gatk_jar       = "/home/shengq1/local/bin/GATK/GenomeAnalysisTK.jar";
my $picard_jar     = "/scratch/cqs/shengq1/local/bin/picard/picard.jar";
my $email          = "quanhu.sheng\@vanderbilt.edu";
my $qc3_perl       = "/scratch/cqs/shengq1/local/bin/qc3/qc3.pl";

my $config = {
  general => { task_name => $task },
  files   => {
    "CSW-1" => ["/gpfs21/scratch/cqs/shengq1/rnaseq/20150625_chenxi_1886_snp/tophat/tophat_1886-CSW-1_out.bam"],
    "CSW-2" => ["/gpfs21/scratch/cqs/shengq1/rnaseq/20150625_chenxi_1886_snp/tophat/tophat_1886-CSW-2_out.bam"],
    "CSW-3" => ["/gpfs21/scratch/cqs/shengq1/rnaseq/20150625_chenxi_1886_snp/tophat/tophat_1886-CSW-3_out.bam"],
    "CSW-4" => ["/gpfs21/scratch/cqs/shengq1/rnaseq/20150625_chenxi_1886_snp/tophat/tophat_1886-CSW-4_out.bam"],
    "CSW-5" => ["/gpfs21/scratch/cqs/shengq1/rnaseq/20150625_chenxi_1886_snp/tophat/tophat_1886-CSW-5_out.bam"],
    "CSW-6" => ["/gpfs21/scratch/cqs/shengq1/rnaseq/20150625_chenxi_1886_snp/tophat/tophat_1886-CSW-6_out.bam"],
  },
  groups => {
    "Mtg16" => [ "CSW-1", "CSW-2", "CSW-3" ],
    "WT"    => [ "CSW-4", "CSW-5", "CSW-6" ],
  },
  qc3bam => {
    class          => "QC::QC3bam",
    perform        => 1,
    target_dir     => "${target_dir}/QC3bam",
    option         => "",
    transcript_gtf => $transcript_gtf,
    qc3_perl       => $qc3_perl,
    source_ref     => "files",
    pbs            => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  tophat_refine => {
    class              => "GATK::RNASeqRefine",
    perform            => 1,
    target_dir         => "${target_dir}/tophat_refine",
    option             => "-Xmx40g",
    fasta_file         => $fasta_file,
    source_ref         => "files",
    vcf_files          => [$dbsnp],
    gatk_jar           => $gatk_jar,
    picard_jar         => $picard_jar,
    sorted             => 0,
    replace_read_group => 1,
    reorder_chromosome => 1,
    fixMisencodedQuals => 0,
    sh_direct          => 0,
    pbs                => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  tophat_refine_SNPindel => {
    class       => "GATK::SNPIndel",
    perform     => 1,
    target_dir  => "${target_dir}/tophat_refine_SNPindel",
    option      => "",
    source_ref  => "tophat_refine",
    java_option => "",
    fasta_file  => $fasta_file,
    dbsnp_vcf   => $dbsnp,
    gatk_jar    => $gatk_jar,
    is_rna      => 1,
    pbs         => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
};

performConfig($config);

1;
