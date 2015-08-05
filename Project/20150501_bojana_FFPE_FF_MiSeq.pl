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
my $dexseq_gff           = "/scratch/cqs/shengq1/references/ensembl_gtf/v75/Homo_sapiens.GRCh37.75.M.dexseq.gff";
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
my $rnaseqqc_jar         = "/scratch/cqs/shengq1/local/bin/RNA-SeQC_v1.1.8.jar";
my $rnaseqqc_gtf         = "/scratch/cqs/shengq1/references/ensembl_gtf/v75/Homo_sapiens.GRCh37.75.M.RNASeQC.gtf";
my $qc3_perl             = "/scratch/cqs/shengq1/local/bin/qc3/qc3.pl";

my $config = {
  general => { task_name => $task },
  files   => {
    "IG-007" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/32-PD_S1_L001_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/32-PD_S1_L001_R2_001.fastq.gz"
    ],
    "IG-008" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/32-PE_S2_L001_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/32-PE_S2_L001_R2_001.fastq.gz"
    ],
    "IG-011" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/33-PA_S2_L001_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/33-PA_S2_L001_R2_001.fastq.gz"
    ],
    "IG-012" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/33-PI_S3_L001_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/33-PI_S3_L001_R2_001.fastq.gz"
    ],
    "IG-016" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/40-PA_S4_L001_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/40-PA_S4_L001_R2_001.fastq.gz"
    ],
    "IG-017" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/40-P_S3_L001_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/40-P_S3_L001_R2_001.fastq.gz"
    ],
    "IG-021" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/42-PC_S6_L001_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/42-PC_S6_L001_R2_001.fastq.gz"
    ],
    "IG-022" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/42-P_S2_L001_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/42-P_S2_L001_R2_001.fastq.gz"
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
    "IG-048" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/IG-48_S2_L001_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/miseq/rawdata/IG-48_S2_L001_R2_001.fastq.gz"
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
    "MiSeq_FF"       => [ "IG-007", "IG-011", "IG-016", "IG-021", "IG-042", "IG-043", "IG-044", "IG-048", "IG-061", "IG-060" ],
    "MiSeq_FFPE"     => [ "IG-008", "IG-012", "IG-017", "IG-022", "IG-051", "IG-052", "IG-053", "IG-057", "IG-058", "IG-059" ],
    "MiSeq_FF_OLD"   => [ "IG-007", "IG-011", "IG-016", "IG-021" ],
    "MiSeq_FFPE_OLD" => [ "IG-008", "IG-012", "IG-017", "IG-022" ],
    "MiSeq_FF_NEW"   => [ "IG-042", "IG-043", "IG-044", "IG-048", "IG-061", "IG-060" ],
    "MiSeq_FFPE_NEW" => [ "IG-051", "IG-052", "IG-053", "IG-057", "IG-058", "IG-059" ],
  },
  pairs => {
    "MiSeq_FFPE_VS_FF" => {
      groups => [ "MiSeq_FFPE", "MiSeq_FF" ],
      paired => [ "B32A",       "B33A", "B40A", "B42", "P06", "P07", "P08", "P12", "P13", "P14" ]
    },
    "MiSeq_FFPE_VS_FF_OLD" => {
      groups => [ "MiSeq_FFPE_OLD", "MiSeq_FF_OLD" ],
      paired => [ "B32A",           "B33A", "B40A", "B42" ]
    },
    "MiSeq_FFPE_VS_FF_NEW" => {
      groups => [ "MiSeq_FFPE_NEW", "MiSeq_FF_NEW" ],
      paired => [ "P06",            "P07", "P08", "P12", "P13", "P14" ]
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
  qc3fastq => {
    class      => "QC::QC3fastq",
    perform    => 1,
    target_dir => "${target_dir}/QC3fastq",
    option     => "",
    qc3_perl   => $qc3_perl,
    source_ref => "trimmer",
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
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
    perform    => 1,
    target_dir => "${target_dir}/fastqc",
    option     => "",
    source_ref => "trimmer",
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
  star_2nd_pass_sort => {
    class         => "Samtools::Sort",
    perform       => 1,
    target_dir    => "${target_dir}/star_2nd_pass_sort",
    option        => "",
    source_ref    => [ "star_2nd_pass", "_Aligned.out.bam" ],
    sort_by_query => 0,
    sh_direct     => 0,
    pbs           => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "20gb"
    },
  },
  rnaseqc => {
    class          => "QC::RNASeQC",
    perform        => 1,
    target_dir     => "${target_dir}/RNASeQC",
    option         => "",
    transcript_gtf => $rnaseqqc_gtf,
    fasta_file     => $fasta_file_16569_M,
    jar            => $rnaseqqc_jar,
    source_ref     => "star_2nd_pass_sort",
    pbs            => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  qc3bam => {
    class          => "QC::QC3bam",
    perform        => 1,
    target_dir     => "${target_dir}/QC3bam",
    option         => "",
    transcript_gtf => $rnaseqqc_gtf,
    qc3_perl       => $qc3_perl,
    source_ref     => "star_2nd_pass_sort",
    pbs            => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
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
    perform       => 0,
    target_dir    => "${target_dir}/star_genetable",
    option        => "-p ENS --noheader -e -o ${task}_gene.count",
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
  star_dexseqcount => {
    class        => "Count::DexseqCount",
    perform      => 1,
    target_dir   => "${target_dir}/star_dexseqcount",
    option       => "",
    source_ref   => [ "star_2nd_pass", "_Aligned.out.bam" ],
    gff_file     => $dexseq_gff,
    dexseq_count => "/home/shengq1/pylibs/bin/dexseq_count.py",
    sh_direct    => 0,
    pbs          => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  star_exontable => {
    class         => "CQS::CQSDatatable",
    perform       => 1,
    target_dir    => "${target_dir}/star_exontable",
    option        => "-p ENS --noheader -o ${task}_exon.count",
    source_ref    => "star_dexseqcount",
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
    class                => "Comparison::DESeq2",
    perform              => 1,
    target_dir           => "${target_dir}/star_deseq2",
    option               => "",
    source_ref           => "pairs",
    groups_ref           => "groups",
    countfile_ref        => "star_genetable",
    sh_direct            => 1,
    show_DE_gene_cluster => 1,
    pvalue               => 0.05,
    fold_change          => 2.0,
    min_median_read      => 5,
    pbs                  => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "10",
      "mem"      => "10gb"
    },
  },
  star_deseq2_strict_criteria => {
    class                => "Comparison::DESeq2",
    perform              => 0,
    target_dir           => "${target_dir}/star_deseq2_strict_criteria",
    option               => "",
    source_ref           => "pairs",
    groups_ref           => "groups",
    countfile_ref        => "star_genetable",
    sh_direct            => 1,
    show_DE_gene_cluster => 1,
    pvalue               => 0.01,
    fold_change          => 2.0,
    pbs                  => {
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
    option      => "",
    source_ref  => "star_2nd_pass_refine",
    java_option => "",
    fasta_file  => $fasta_file_16569_M,
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
  sequencetask => {
    class      => "CQS::SequenceTask",
    perform    => 0,
    target_dir => "${target_dir}/sequencetask",
    option     => "",
    source     => {
      step_1 => [ "trimmer", "fastqlen", "fastqc", "star" ],
      step_2 => ["star_index"],
      step_3 => [ "star_2nd_pass", "star_htseqcount", "star_2nd_pass_refine" ],
      step_4 => [ "star_genetable", ],
      step_5 => [ "star_2nd_pass_refine_SNPindel", "star_2nd_pass_refine_SNPindel_annovar" ],
      step_6 => [ "star_deseq2",                   "star_deseq2_strict_criteria" ]
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

#performConfig($config);
#performTask( $config, "star_dexseqcount" );
performTask( $config, "star_exontable" );

1;

