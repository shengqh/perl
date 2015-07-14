#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;

my $transcript_gtf       = "/scratch/cqs/shengq1/references/ensembl_gtf/v75/Homo_sapiens.GRCh37.75.MT.gtf";
my $name_map_file        = "/scratch/cqs/shengq1/references/ensembl_gtf/v75/Homo_sapiens.GRCh37.75.MT.map";
my $transcript_gtf_index = "/scratch/cqs/shengq1/references/ensembl_gtf/v75/gtfindex/Homo_sapiens.GRCh37.75.MT";
my $fasta_file_16569_MT  = "/scratch/cqs/shengq1/references/hg19_16569_MT/hg19_16569_MT.fa";
my $cqstools             = "/home/shengq1/cqstools/CQS.Tools.exe";
my $gatk_jar             = "/home/shengq1/local/bin/GATK/GenomeAnalysisTK.jar";
my $picard_jar           = "/scratch/cqs/shengq1/local/bin/picard/picard.jar";
my $star_index           = "/scratch/cqs/shengq1/references/hg19_16569_MT/STAR_index_v37.75_2.4.0j_sjdb100";
my $annovar_param        = "-protocol refGene,snp138,cosmic70 -operation g,f,f --remove";
my $annovar_db           = "/scratch/cqs/shengq1/references/annovar/humandb/";
my $qc3_perl             = "/scratch/cqs/shengq1/local/bin/qc3/qc3.pl";

my $dbsnp  = "/scratch/cqs/shengq1/references/dbsnp/human_GRCh37_v142_16569_MT.vcf";
my $hapmap = "/scratch/cqs/shengq1/references/gatk/b37/hapmap_3.3.b37.vcf";
my $omni   = "/scratch/cqs/shengq1/references/gatk/b37/1000G_omni2.5.b37.vcf";
my $g1000  = "/scratch/cqs/shengq1/references/gatk/b37/1000G_phase1.snps.high_confidence.b37.vcf";
my $mills  = "/scratch/cqs/shengq1/references/gatk/b37/Mills_and_1000G_gold_standard.indels.b37.vcf";

my $task       = "20150630_bojana_tnbc";
my $target_dir = "/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc";
my $email      = "quanhu.sheng\@vanderbilt.edu";

my $download = {
  general => { task_name => $task },
  files   => { "samples" => ["/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/fastq.list"], },
  wget    => {
    class      => "Data::Wget",
    perform    => 1,
    target_dir => "${target_dir}/raw",
    option     => "",
    source_ref => "files",
    sh_direct  => 0,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "10gb"
    },
  },
};

#performConfig($download);

my $config = {
  general => { task_name => $task },
  files   => {
    "3193-BJ-0001" =>
      [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL112956_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL112956_2.fastq.gz" ],
    "3193-BJ-0002" =>
      [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL112957_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL112957_2.fastq.gz" ],
    "3193-BJ-0003" =>
      [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL112958_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL112958_2.fastq.gz" ],
    "3193-BJ-0004" =>
      [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL112959_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL112959_2.fastq.gz" ],
    "3193-BJ-0005" =>
      [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL112960_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL112960_2.fastq.gz" ],
    "3193-BJ-0006" =>
      [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL112961_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL112961_2.fastq.gz" ],
    "3193-BJ-0007" =>
      [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL112962_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL112962_2.fastq.gz" ],
    "3193-BJ-0008" =>
      [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL112963_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL112963_2.fastq.gz" ],
    "3193-BJ-0009" =>
      [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL112964_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL112964_2.fastq.gz" ],
    "3193-BJ-0010" =>
      [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL112965_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL112965_2.fastq.gz" ],
    "3193-BJ-0011" =>
      [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL112966_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL112966_2.fastq.gz" ],
    "3193-BJ-0012" =>
      [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL112967_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL112967_2.fastq.gz" ],
    "3193-BJ-0013" =>
      [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL112968_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL112968_2.fastq.gz" ],
    "3193-BJ-0014" =>
      [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL112969_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL112969_2.fastq.gz" ],
    "3193-BJ-0015" =>
      [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL112970_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL112970_2.fastq.gz" ],
    "3193-BJ-0016" =>
      [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL112971_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL112971_2.fastq.gz" ],
    "3193-BJ-0017" =>
      [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL112972_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL112972_2.fastq.gz" ],
    "3193-BJ-0018" =>
      [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL112973_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL112973_2.fastq.gz" ],
    "3193-BJ-0019" =>
      [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL112974_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL112974_2.fastq.gz" ],
    "3193-BJ-0020" =>
      [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL112975_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL112975_2.fastq.gz" ],
    "3193-BJ-0021" =>
      [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL112976_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL112976_2.fastq.gz" ],
    "3193-BJ-0022" =>
      [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL112977_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL112977_2.fastq.gz" ],
    "3193-BJ-0023" =>
      [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL112978_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL112978_2.fastq.gz" ],
    "3193-BJ-0024" =>
      [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL112979_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL112979_2.fastq.gz" ],
    "3193-BJ-0025" =>
      [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL112980_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL112980_2.fastq.gz" ],
    "3193-BJ-0026" =>
      [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL112981_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL112981_2.fastq.gz" ],
    "3193-BJ-0027" =>
      [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL112982_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL112982_2.fastq.gz" ],
    "3193-BJ-0028" =>
      [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL112983_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL112983_2.fastq.gz" ],
    "3193-BJ-0029" =>
      [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL112984_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL112984_2.fastq.gz" ],
    "3193-BJ-0030" =>
      [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL112985_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL112985_2.fastq.gz" ],
    "3193-BJ-0031" =>
      [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL112986_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL112986_2.fastq.gz" ],
    "3193-BJ-0032" =>
      [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL112987_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL112987_2.fastq.gz" ],
    "3193-BJ-0033" =>
      [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL112988_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL112988_2.fastq.gz" ],
    "3193-BJ-0034" =>
      [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL112989_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL112989_2.fastq.gz" ],
    "3193-BJ-0035" =>
      [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL112990_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL112990_2.fastq.gz" ],
    "3193-BJ-0036" =>
      [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL112991_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL112991_2.fastq.gz" ],
    "3193-BJ-0037" =>
      [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL112992_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL112992_2.fastq.gz" ],
    "3193-BJ-0038" =>
      [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL112993_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL112993_2.fastq.gz" ],
    "3193-BJ-0039" =>
      [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL112994_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL112994_2.fastq.gz" ],
    "3193-BJ-0040" =>
      [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL112995_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL112995_2.fastq.gz" ],
    "3193-BJ-0041" =>
      [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL112996_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL112996_2.fastq.gz" ],
    "3193-BJ-0042" =>
      [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL112997_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL112997_2.fastq.gz" ],
    "3193-BJ-0043" =>
      [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL112998_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL112998_2.fastq.gz" ],
    "3193-BJ-0044" =>
      [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL112999_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL112999_2.fastq.gz" ],
    "3193-BJ-0045" =>
      [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL113000_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/SL113000_2.fastq.gz" ],
  },

  groups => {
    "Group1" => [
      "3193-BJ-0008", "3193-BJ-0009", "3193-BJ-0011", "3193-BJ-0012", "3193-BJ-0016", "3193-BJ-0017", "3193-BJ-0020", "3193-BJ-0021", "3193-BJ-0023", "3193-BJ-0026",
      "3193-BJ-0027", "3193-BJ-0028", "3193-BJ-0035", "3193-BJ-0037", "3193-BJ-0038", "3193-BJ-0041", "3193-BJ-0042", "3193-BJ-0043", "3193-BJ-0045"
    ],
    "Group2" => [
      "3193-BJ-0003", "3193-BJ-0006", "3193-BJ-0013", "3193-BJ-0019", "3193-BJ-0024", "3193-BJ-0029", "3193-BJ-0031", "3193-BJ-0036", "3193-BJ-0040"
    ],
    "Group12" => [
      "3193-BJ-0008", "3193-BJ-0009", "3193-BJ-0011", "3193-BJ-0012", "3193-BJ-0016", "3193-BJ-0017", "3193-BJ-0020", "3193-BJ-0021", "3193-BJ-0023", "3193-BJ-0026",
      "3193-BJ-0027", "3193-BJ-0028", "3193-BJ-0035", "3193-BJ-0037", "3193-BJ-0038", "3193-BJ-0041", "3193-BJ-0042", "3193-BJ-0043", "3193-BJ-0045",
      "3193-BJ-0003", "3193-BJ-0006", "3193-BJ-0013", "3193-BJ-0019", "3193-BJ-0024", "3193-BJ-0029", "3193-BJ-0031", "3193-BJ-0036", "3193-BJ-0040"
    ],
     "Group3" => [
      "3193-BJ-0001", "3193-BJ-0002", "3193-BJ-0004", "3193-BJ-0005", "3193-BJ-0007", "3193-BJ-0010", "3193-BJ-0014", "3193-BJ-0015", "3193-BJ-0018", "3193-BJ-0022",
      "3193-BJ-0025", "3193-BJ-0030", "3193-BJ-0032", "3193-BJ-0033", "3193-BJ-0034", "3193-BJ-0039", "3193-BJ-0044"
    ],
  },
  pairs => {
    "ClinicalResponse"     => { groups => [ "Group12",     "Group3" ], },
  },
  fastqc => {
    class      => "QC::FastQC",
    perform    => 1,
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
  fastqc_summary => {
    class      => "QC::FastQCSummary",
    perform    => 1,
    target_dir => "${target_dir}/fastqc",
    option     => "",
    source_ref => "fastqc",
    cqstools   => $cqstools,
    sh_direct  => 0,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "2",
      "mem"      => "10gb"
    },
  },
  star => {
    class      => "Alignment::STAR",
    perform    => 1,
    target_dir => "${target_dir}/star",
    option     => "",
    source_ref => "files",
    genome_dir => $star_index,
    sh_direct  => 0,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "30gb"
    },
  },
  star_index => {
    class          => "Alignment::STARIndex",
    perform        => 1,
    target_dir     => "${target_dir}/star_index",
    option         => "--sjdbOverhang 100",
    source_ref     => [ "star", "tab\$" ],
    fasta_file     => $fasta_file_16569_MT,
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
    perform         => 1,
    target_dir      => "${target_dir}/star_2nd_pass",
    option          => "",
    source_ref      => "files",
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
  qc3 => {
    class          => "QC::QC3bam",
    perform        => 1,
    target_dir     => "${target_dir}/QC3bam",
    option         => "",
    transcript_gtf => $transcript_gtf,
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
    perform    => 1,
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
    pbs                  => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "10",
      "mem"      => "10gb"
    },
  },
  star_2nd_pass_refine => {
    class              => "GATK::RNASeqRefine",
    perform            => 1,
    target_dir         => "${target_dir}/star_2nd_pass_refine",
    option             => "-Xmx40g",
    fasta_file         => $fasta_file_16569_MT,
    source_ref         => "star_2nd_pass",
    vcf_files          => [$dbsnp],
    gatk_jar           => $gatk_jar,
    picard_jar         => $picard_jar,
    fixMisencodedQuals => 0,
    replace_read_group => 0,
    reorderChromosome  => 0,
    sorted             => 0,
    sh_direct          => 0,
    pbs                => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  hc_gvcf => {
    class       => "GATK::HaplotypeCallerGVCF",
    perform     => 1,
    target_dir  => "${target_dir}/hc_gvcf",
    option      => "",
    source_ref  => "star_2nd_pass_refine",
    java_option => "",
    fasta_file  => $fasta_file_16569_MT,
    dbsnp_vcf   => $dbsnp,
    gatk_jar    => $gatk_jar,
    sh_direct   => 0,
    pbs         => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  hc_gvcf_vqsr => {
    class       => "GATK::VariantFilterVQSR",
    perform     => 1,
    target_dir  => "${target_dir}/hc_gvcf_vqsr",
    option      => "",
    source_ref  => "hc_gvcf",
    java_option => "",
    fasta_file  => $fasta_file_16569_MT,
    dbsnp_vcf   => $dbsnp,
    hapmap_vcf  => $hapmap,
    omni_vcf    => $omni,
    g1000_vcf   => $g1000,
    mills_vcf   => $mills,
    gatk_jar    => $gatk_jar,
    sh_direct   => 1,
    pbs         => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "100gb"
    },
  },
  hc_gvcf_vqsr_annovar => {
    class      => "Annotation::Annovar",
    perform    => 1,
    target_dir => "${target_dir}/hc_gvcf_vqsr_annovar",
    source_ref => "hc_gvcf_vqsr",
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
      step1 => [ "fastqc",         "star" ],
      step2 => [ "fastqc_summary", "star_index" ],
      step3 => [ "star_2nd_pass",  "star_htseqcount", "star_2nd_pass_refine", "hc_gvcf", "star_2nd_pass_sort" ],
      step4 => [ "qc3", "star_genetable", "hc_gvcf_vqsr", "hc_gvcf_vqsr_annovar" ],
    },
    sh_direct => 0,
    pbs       => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },

  #  star_2nd_pass_refine_SNPindel_annovar => {
  #    class      => "Annovar",
  #    perform    => 0,
  #    target_dir => "${target_dir}/star_2nd_pass_refine_SNPindel_annovar",
  #    source_ref => "star_2nd_pass_refine_SNPindel",
  #    option     => $annovar_param,
  #    annovar_db => $annovar_db,
  #    buildver   => "hg19",
  #    sh_direct  => 1,
  #    isvcf      => 1,
  #    pbs        => {
  #      "email"    => $email,
  #      "nodes"    => "1:ppn=1",
  #      "walltime" => "72",
  #      "mem"      => "10gb"
  #    },
  #  },
};

performConfig($config);
#performTask( $config, "hc_gvcf_vqsr_annovar" );

1;
