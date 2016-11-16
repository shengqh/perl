#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;

my $task = "P3582";

my $target_dir = create_directory_or_die("/scratch/cqs/shengq1/haake/20161114_human_celllines_rnaseq");

#my $target_dir = "e:/temp";

my $transcript_gtf = "/scratch/cqs/shengq1/references/gencode/hg19/gencode.v23lift37.annotation.gtf";
my $name_map_file  = "/scratch/cqs/shengq1/references/gencode/hg19/gencode.v23lift37.annotation.map";
my $star_index     = "/scratch/cqs/shengq1/references/gencode/hg19/STAR_index_2.5.2a";
my $fasta_file     = "/scratch/cqs/shengq1/references/gencode/hg19/STAR_index_2.5.2a/GRCh37.p13.genome.fa";

my $gatk_jar   = "/home/shengq1/local/bin/GATK/GenomeAnalysisTK.jar";
my $picard_jar = "/scratch/cqs/shengq1/local/bin/picard/picard.jar";

my $dbsnp             = "/scratch/cqs/shengq1/references/dbsnp/human_GRCh37_v142_16569_M_chr.vcf";
my $annovar_protocol  = "refGene,avsnp147,cosmic70";
my $annovar_operation = "g,f,f";
my $annovar_param     = "-protocol ${annovar_protocol} -operation ${annovar_operation} --remove";
my $annovar_db        = "/scratch/cqs/shengq1/references/annovar/humandb/";
my $rnaediting_db     = "/data/cqs/shengq1/reference/rnaediting/hg19.txt";

my $cqstools = "/home/shengq1/cqstools/cqstools.exe";
my $glmvc    = "/home/shengq1/glmvc/glmvc.exe";
my $email    = "quanhu.sheng\@vanderbilt.edu";

my $config = {
  general => { task_name => $task },
  files   => {
    "3577-WKR-1" => [
      "/gpfs21/scratch/cqs/shengq1/proteomics/20161019_proteogenomics/data/rnaseq/3577-WKR-1_S1_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20161019_proteogenomics/data/rnaseq/3577-WKR-1_S1_R2_001.fastq.gz"
    ],
    "3577-WKR-2" => [
      "/gpfs21/scratch/cqs/shengq1/proteomics/20161019_proteogenomics/data/rnaseq/3577-WKR-2_S2_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20161019_proteogenomics/data/rnaseq/3577-WKR-2_S2_R2_001.fastq.gz"
    ],
    "3577-WKR-3" => [
      "/gpfs21/scratch/cqs/shengq1/proteomics/20161019_proteogenomics/data/rnaseq/3577-WKR-3_S3_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20161019_proteogenomics/data/rnaseq/3577-WKR-3_S3_R2_001.fastq.gz"
    ],
    "3577-WKR-4" => [
      "/gpfs21/scratch/cqs/shengq1/proteomics/20161019_proteogenomics/data/rnaseq/3577-WKR-4_S4_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20161019_proteogenomics/data/rnaseq/3577-WKR-4_S4_R2_001.fastq.gz"
    ],
    "3577-WKR-5" => [
      "/gpfs21/scratch/cqs/shengq1/proteomics/20161019_proteogenomics/data/rnaseq/3577-WKR-5_S5_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20161019_proteogenomics/data/rnaseq/3577-WKR-5_S5_R2_001.fastq.gz"
    ],
    "3577-WKR-6" => [
      "/gpfs21/scratch/cqs/shengq1/proteomics/20161019_proteogenomics/data/rnaseq/3577-WKR-6_S6_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20161019_proteogenomics/data/rnaseq/3577-WKR-6_S6_R2_001.fastq.gz"
    ],
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
  star_featurecount => {
    class      => "Count::FeatureCounts",
    perform    => 1,
    target_dir => "${target_dir}/star_featurecount",
    option     => "-g gene_id -t exon -M",
    source_ref => [ "star", "_Aligned.sortedByCoord.out.bam" ],
    gff_file   => $transcript_gtf,
    ispairend  => 1,
    sh_direct  => 0,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
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
  star_refine => {
    class            => "GATK::RNASeqRefine",
    perform          => 1,
    target_dir       => "${target_dir}/star_refine",
    option           => "-Xmx40g",
    fasta_file       => $fasta_file,
    source_ref       => "star",
    vcf_files        => [$dbsnp],
    gatk_jar         => $gatk_jar,
    picard_jar       => $picard_jar,
    slim_print_reads => 1,
    sorted           => 1,
    sh_direct        => 1,
    pbs              => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  star_refine_hc => {
    class         => "GATK::HaplotypeCaller",
    perform       => 1,
    target_dir    => "${target_dir}/star_refine_hc",
    option        => "",
    source_ref    => "star_refine",
    java_option   => "",
    fasta_file    => $fasta_file,
    gatk_jar      => $gatk_jar,
    extension     => ".vcf",
    by_chromosome => 0,                               #since we have the bed file, we cannot use by_chromosome.
    gvcf          => 0,                               #http://gatkforums.broadinstitute.org/gatk/discussion/3891/calling-variants-in-rnaseq
    sh_direct     => 0,
    pbs           => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  star_refine_hc_filter => {
    class       => "GATK::VariantFilter",
    perform     => 1,
    target_dir  => "${target_dir}/star_refine_hc_filter",
    option      => "",
    vqsr_mode   => 0,
    source_ref  => "star_refine_hc",
    java_option => "",
    fasta_file  => $fasta_file,
    dbsnp_vcf   => $dbsnp,
    gatk_jar    => $gatk_jar,
    cqstools    => $cqstools,
    is_rna      => 1,
    sh_direct   => 1,
    pbs         => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  star_refine_hc_filter_annovar => {
    class      => "Annotation::Annovar",
    perform    => 1,
    target_dir => "${target_dir}/star_refine_hc_filter_annovar",
    source_ref => "star_refine_hc_filter",
    option     => $annovar_param,
    annovar_db => $annovar_db,
    buildver   => "hg19",
    sh_direct  => 1,
    isvcf      => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "2",
      "mem"      => "10gb"
    },
  },
  sequencetask => {
    class      => "CQS::SequenceTask",
    perform    => 1,
    target_dir => "${target_dir}/sequencetask",
    option     => "",
    source     => {
      step1 => [ "fastqc",         "star",           "star_featurecount",    "star_refine", "star_refine_hc" ],
      step2 => [ "fastqc_summary", "star_genetable", "star_refine_hc_filter", "star_refine_hc_filter_annovar" ],
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
1;

