#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;

my $task = "JP3533";

my $target_dir = create_directory_or_die("/scratch/cqs/shengq1/rnaseq/20160816_johanna_resistant_cells");

#my $target_dir = "e:/temp";

my $transcript_gtf = "/scratch/cqs/shengq1/references/gencode/hg19/gencode.v23lift37.annotation.gtf";
my $name_map_file  = "/scratch/cqs/shengq1/references/gencode/hg19/gencode.v23lift37.annotation.map";
my $star_index     = "/scratch/cqs/shengq1/references/gencode/hg19/STAR_index_2.5.2a";
my $fasta_file     = "/scratch/cqs/shengq1/references/gencode/hg19/STAR_index_2.5.2a/GRCh37.p13.genome.fa";
my $cqstools       = "/home/shengq1/cqstools/cqstools.exe";
my $email          = "quanhu.sheng\@vanderbilt.edu";
my $qc3_perl       = "/scratch/cqs/shengq1/local/bin/qc3/qc3.pl";
my $gatk_jar       = "/home/shengq1/local/bin/GATK/GenomeAnalysisTK.jar";
my $picard_jar     = "/scratch/cqs/shengq1/local/bin/picard/picard.jar";
my $annovar_param  = "-protocol refGene,snp138,cosmic70 -operation g,f,f --remove";
my $annovar_db     = "/scratch/cqs/shengq1/references/annovar/humandb/";
my $dbsnp          = "/data/cqs/shengq1/reference/dbsnp/human_9606_b147_GRCh37p13.vcf";
my $config         = {
  general => { task_name => $task },
  files   => {
    "P_CAL148" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20160816_johanna_resistant_cells/fastq/3533-JP-9_1.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20160816_johanna_resistant_cells/fastq/3533-JP-9_2.fastq.gz"
    ],
    "R_CAL148" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20160816_johanna_resistant_cells/fastq/3533-JP-10_1.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20160816_johanna_resistant_cells/fastq/3533-JP-10_2.fastq.gz"
    ],
    "P_MDA" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20160816_johanna_resistant_cells/fastq/3533-JP-11_1.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20160816_johanna_resistant_cells/fastq/3533-JP-11_2.fastq.gz"
    ],
    "R_MDA" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20160816_johanna_resistant_cells/fastq/3533-JP-12_1.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20160816_johanna_resistant_cells/fastq/3533-JP-12_2.fastq.gz"
    ],
    "P_CAL51_1" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20160816_johanna_resistant_cells/fastq/3533-JP-13_1.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20160816_johanna_resistant_cells/fastq/3533-JP-13_2.fastq.gz"
    ],
    "P_CAL51_2" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20160816_johanna_resistant_cells/fastq/3533-JP-14_1.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20160816_johanna_resistant_cells/fastq/3533-JP-14_2.fastq.gz"
    ],
    "R_CAL51_1" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20160816_johanna_resistant_cells/fastq/3533-JP-15_1.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20160816_johanna_resistant_cells/fastq/3533-JP-15_2.fastq.gz"
    ],
    "R_CAL51_2" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20160816_johanna_resistant_cells/fastq/3533-JP-16_1.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20160816_johanna_resistant_cells/fastq/3533-JP-16_2.fastq.gz"
    ],

  },
  groups => {
    "Parental"  => [ "P_CAL148", "P_MDA", "P_CAL51_1", "P_CAL51_2" ],
    "Resistant" => [ "R_CAL148", "R_MDA", "R_CAL51_1", "R_CAL51_2" ],
  },
  pairs => {
    "Parental_vs_Resistant" => {
      groups    => [ "Resistant", "Parental" ],
      celllines => [ "CAL148",    "MDA", "CAL51", "CAL51", "CAL148", "MDA", "CAL51", "CAL51" ],
    },
  },
  fastqc => {
    class      => "QC::FastQC",
    perform    => 1,
    target_dir => "${target_dir}/fastqc",
    option     => "",
    source_ref => "files",
    sh_direct  => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=2",
      "walltime" => "2",
      "mem"      => "10gb"
    },
  },
  fastqc_summary => {
    class      => "QC::FastQCSummary",
    perform    => 1,
    target_dir => "${target_dir}/fastqc",
    option     => "",
    source_ref => "fastqc",
    cqstools   => $cqstools,
    sh_direct  => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "2",
      "mem"      => "10gb"
    },
  },
  star => {
    class                     => "Alignment::STAR",
    perform                   => 1,
    target_dir                => "${target_dir}/star",
    option                    => "--twopassMode Basic",
    source_ref                => "files",
    genome_dir                => $star_index,
    output_sort_by_coordinate => 1,
    sh_direct                 => 0,
    pbs                       => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "30gb"
    },
  },
  star_qc3 => {
    class          => "QC::QC3bam",
    perform        => 1,
    target_dir     => "${target_dir}/star_qc3",
    option         => "",
    transcript_gtf => $transcript_gtf,
    qc3_perl       => $qc3_perl,
    source_ref     => [ "star", "_Aligned.sortedByCoord.out.bam" ],
    pbs            => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  star_featurecount => {
    class      => "Count::FeatureCounts",
    perform    => 1,
    target_dir => "${target_dir}/star_featurecount",
    option     => "-g gene_id -t exon",
    source_ref => [ "star", "_Aligned.sortedByCoord.out.bam" ],
    gff_file   => $transcript_gtf,
    ispairend  => 1,
    sh_direct  => 0,
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
    option        => "-k 0 -v 6 -e -o ${task}_gene.count",
    source_ref    => "star_featurecount",
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
  star_genetable_deseq2 => {
    class                => "Comparison::DESeq2",
    perform              => 1,
    target_dir           => "${target_dir}/star_genetable_deseq2",
    option               => "",
    source_ref           => "pairs",
    groups_ref           => "groups",
    countfile_ref        => "star_genetable",
    sh_direct            => 1,
    show_DE_gene_cluster => 0,
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
  star_genetable_correlation => {
    class                    => "CQS::UniqueR",
    perform                  => 1,
    target_dir               => "${target_dir}/star_genetable_correlation",
    rtemplate                => "countTableVisFunctions.R,countTableCorrelation.R",
    output_file              => "parameterSampleFile1",
    output_file_ext          => ".Correlation.png",
    parameterSampleFile1_ref => [ "star_genetable", ".count\$" ],
    parameterSampleFile2_ref => "groups",
    sh_direct                => 1,
    pbs                      => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "1",
      "mem"      => "10gb"
    },
  },
  star_refine => {
    class                    => "GATK::RNASeqRefine",
    perform                  => 1,
    target_dir               => "${target_dir}/star_refine",
    option                   => "-Xmx40g",
    fasta_file               => $fasta_file,
    source_ref               => "star",
    vcf_files                => [$dbsnp],
    gatk_jar                 => $gatk_jar,
    picard_jar               => $picard_jar,
    replace_read_group       => 0,
    reorder_chromosome       => 0,
    fixMisencodedQuals       => 0,
    samtools_baq_calibration => 0,
    sorted                   => 1,
    sh_direct                => 0,
    pbs                      => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  star_refine_SNPindel => {
    class       => "GATK::SNPIndel",
    perform     => 1,
    target_dir  => "${target_dir}/star_refine_SNPindel",
    option      => "",
    source_ref  => "star_refine",
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
  star_refine_SNPindel_annovar => {
    class      => "Annotation::Annovar",
    perform    => 1,
    target_dir => "${target_dir}/star_refine_SNPindel_annovar",
    source_ref => "star_refine_SNPindel",
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
    perform    => 1,
    target_dir => "${target_dir}/sequencetask",
    option     => "",
    source     => {
      step1 => [ "fastqc",         "star",           "star_featurecount",          "star_refine", "star_refine_SNPindel", "star_refine_SNPindel_annovar" ],
      step2 => [ "fastqc_summary", "star_genetable", "star_genetable_correlation", "star_genetable_deseq2" ],
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

#performTask( $config, "star_qc3" );

1;
