#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;

my $task = "FFPE_FF_HiSeq";

my $target_dir = "/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq";

my $transcript_gtf       = "/scratch/cqs/shengq1/references/ensembl_gtf/v75/Homo_sapiens.GRCh37.75.M.gtf";
my $name_map_file        = "/scratch/cqs/shengq1/references/ensembl_gtf/v75/Homo_sapiens.GRCh37.75.M.map";
my $transcript_gtf_index = "/scratch/cqs/shengq1/references/ensembl_gtf/v75/gtfindex/Homo_sapiens.GRCh37.75.M";
my $fasta_file_16569_M   = "/data/cqs/guoy1/reference/hg19/bowtie2_index/hg19.fa";
my $bowtie2_index        = "/data/cqs/guoy1/reference/hg19/bowtie2_index/hg19";
my $cqstools             = "/home/shengq1/cqstools/CQS.Tools.exe";

my $email = "quanhu.sheng\@vanderbilt.edu";

my $config = {
  general    => { task_name => $task },
  files => {
    "IG-062" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/rawdata/2059-JP-10-1_AGTCAA_L001_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/rawdata/2059-JP-10-1_AGTCAA_L001_R2_001.fastq.gz"
    ],
    "IG-063" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/rawdata/2059-JP-12-1_CAGGCG_L003_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/rawdata/2059-JP-12-1_CAGGCG_L003_R2_001.fastq.gz"
    ],
    "IG-064" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/rawdata/2059-JP-10-2_ATGTCA_L001_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/rawdata/2059-JP-10-2_ATGTCA_L001_R2_001.fastq.gz"
    ],
    "IG-065" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/rawdata/2059-JP-12-2_CATTTT_L003_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/rawdata/2059-JP-12-2_CATTTT_L003_R2_001.fastq.gz"
    ],
    "IG-066" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/rawdata/2059-JP-16-1_CGTACG_L007_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/rawdata/2059-JP-16-1_CGTACG_L007_R2_001.fastq.gz"
    ],
    "IG-067" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/rawdata/2059-JP-12-3_CGGAAT_L003_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/rawdata/2059-JP-12-3_CGGAAT_L003_R2_001.fastq.gz"
    ],
    "IG-068" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/rawdata/2059-JP-16-2_CACGAT_L007_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/rawdata/2059-JP-16-2_CACGAT_L007_R2_001.fastq.gz"
    ],
    "IG-069" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/rawdata/2059-JP-12-4_CTATAC_L003_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/rawdata/2059-JP-12-4_CTATAC_L003_R2_001.fastq.gz"
    ],
    "IG-070" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/rawdata/2059-JP-17-1_GCCAAT_L008_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/rawdata/2059-JP-17-1_GCCAAT_L008_R2_001.fastq.gz"
    ],
    "IG-071" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/rawdata/2059-JP-12-5_TCATTC_L003_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/rawdata/2059-JP-12-5_TCATTC_L003_R2_001.fastq.gz"
    ],
    "IG-072" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/rawdata/2059-JP-17-2_CTTGTA_L008_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/rawdata/2059-JP-17-2_CTTGTA_L008_R2_001.fastq.gz"
    ],
    "IG-073" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/rawdata/2059-JP-13-1_GTTTCG_L004_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/rawdata/2059-JP-13-1_GTTTCG_L004_R2_001.fastq.gz"
    ],
    "IG-074" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/rawdata/2059-JP-10-3_GTGGCC_L001_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/rawdata/2059-JP-10-3_GTGGCC_L001_R2_001.fastq.gz"
    ],
    "IG-075" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/rawdata/2059-JP-13-2_CAAAAG_L004_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/rawdata/2059-JP-13-2_CAAAAG_L004_R2_001.fastq.gz"
    ],
    "IG-076" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/rawdata/2059-JP-10-4_GAGTGG_L001_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/rawdata/2059-JP-10-4_GAGTGG_L001_R2_001.fastq.gz"
    ],
    "IG-077" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/rawdata/2059-JP-13-3_CACTCA_L004_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/rawdata/2059-JP-13-3_CACTCA_L004_R2_001.fastq.gz"
    ],
    "IG-078" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/rawdata/2059-JP-10-5_CACCGG_L001_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/rawdata/2059-JP-10-5_CACCGG_L001_R2_001.fastq.gz"
    ],
    "IG-079" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/rawdata/2059-JP-13-4_CATGGC_L004_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/rawdata/2059-JP-13-4_CATGGC_L004_R2_001.fastq.gz"
    ],
    "IG-080" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/rawdata/2059-JP-17-3_GAGTGG_L008_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/rawdata/2059-JP-17-3_GAGTGG_L008_R2_001.fastq.gz"
    ],
    "IG-081" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/rawdata/2059-JP-13-5_CTCAGA_L004_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/rawdata/2059-JP-13-5_CTCAGA_L004_R2_001.fastq.gz"
    ],
    "IG-082" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/rawdata/2059-JP-17-4_TCCCGA_L008_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/rawdata/2059-JP-17-4_TCCCGA_L008_R2_001.fastq.gz"
    ],
    "IG-083" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/rawdata/2059-JP-14-1_GGTAGC_L005_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/rawdata/2059-JP-14-1_GGTAGC_L005_R2_001.fastq.gz"
    ],
    "IG-085" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/rawdata/2059-JP-14-2_CCAACA_L005_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/rawdata/2059-JP-14-2_CCAACA_L005_R2_001.fastq.gz"
    ],
    "IG-086" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/rawdata/2059-JP-11-1_CCGTCC_L002_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/rawdata/2059-JP-11-1_CCGTCC_L002_R2_001.fastq.gz"
    ],
    "IG-087" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/rawdata/2059-JP-14-3_CTAGCT_L005_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/rawdata/2059-JP-14-3_CTAGCT_L005_R2_001.fastq.gz"
    ],
    "IG-088" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/rawdata/2059-JP-11-2_GTAGAG_L002_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/rawdata/2059-JP-11-2_GTAGAG_L002_R2_001.fastq.gz"
    ],
    "IG-089" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/rawdata/2059-JP-14-4_TATAAT_L005_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/rawdata/2059-JP-14-4_TATAAT_L005_R2_001.fastq.gz"
    ],
    "IG-090" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/rawdata/2059-JP-11-3_GTGAAA_L002_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/rawdata/2059-JP-11-3_GTGAAA_L002_R2_001.fastq.gz"
    ],
    "IG-091" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/rawdata/2059-JP-14-5_TCGGCA_L005_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/rawdata/2059-JP-14-5_TCGGCA_L005_R2_001.fastq.gz"
    ],
    "IG-092" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/rawdata/2059-JP-16-3_ATTCCT_L007_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/rawdata/2059-JP-16-3_ATTCCT_L007_R2_001.fastq.gz"
    ],
    "IG-093" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/rawdata/2059-JP-15-1_AGTTCC_L006_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/rawdata/2059-JP-15-1_AGTTCC_L006_R2_001.fastq.gz"
    ],
    "IG-094" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/rawdata/2059-JP-16-4_CAACTA_L007_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/rawdata/2059-JP-16-4_CAACTA_L007_R2_001.fastq.gz"
    ],
    "IG-095" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/rawdata/2059-JP-15-2_GTCCGC_L006_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/rawdata/2059-JP-15-2_GTCCGC_L006_R2_001.fastq.gz"
    ],
    "IG-096" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/rawdata/2059-JP-16-5_GACGAC_L007_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/rawdata/2059-JP-16-5_GACGAC_L007_R2_001.fastq.gz"
    ],
    "IG-097" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/rawdata/2059-JP-15-3_ACTGAT_L006_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/rawdata/2059-JP-15-3_ACTGAT_L006_R2_001.fastq.gz"
    ],
    "IG-098" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/rawdata/2059-JP-15-4_TACAGC_L006_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/rawdata/2059-JP-15-4_TACAGC_L006_R2_001.fastq.gz"
    ],
    "IG-099" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/rawdata/2059-JP-11-4_ATGAGC_L002_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/rawdata/2059-JP-11-4_ATGAGC_L002_R2_001.fastq.gz"
    ],
    "IG-100" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/rawdata/2059-JP-15-5_TCGAAG_L006_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/rawdata/2059-JP-15-5_TCGAAG_L006_R2_001.fastq.gz"
    ],
    "IG-101" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/rawdata/2059-JP-11-5_TAATCG_L002_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/rawdata/2059-JP-11-5_TAATCG_L002_R2_001.fastq.gz"
    ],
  },
  groups => {
    "HiSeq_FF" =>
      [ "IG-062", "IG-064", "IG-066", "IG-068", "IG-070", "IG-072", "IG-074", "IG-076", "IG-078", "IG-080", "IG-082", "IG-086", "IG-088", "IG-090", "IG-092", "IG-094", "IG-096", "IG-098", "IG-100", ],
    "HiSeq_FFPE" =>
      [ "IG-063", "IG-065", "IG-067", "IG-069", "IG-071", "IG-073", "IG-075", "IG-077", "IG-079", "IG-081", "IG-083", "IG-087", "IG-089", "IG-091", "IG-093", "IG-095", "IG-097", "IG-099", "IG-101", ],
    "HiSeq_FF2"   => [ "IG-062", "IG-064", "IG-066", "IG-068", "IG-070", "IG-072", "IG-074", "IG-076", "IG-078", "IG-080", "IG-082", "IG-086", "IG-088", "IG-090" ],
    "HiSeq_FFPE2" => [ "IG-063", "IG-065", "IG-067", "IG-069", "IG-071", "IG-073", "IG-075", "IG-077", "IG-079", "IG-081", "IG-083", "IG-087", "IG-089", "IG-091" ],
    "HiSeq_FF3"   => [ "IG-092", "IG-094", "IG-096", "IG-098", "IG-100", ],
    "HiSeq_FFPE3" => [ "IG-093", "IG-095", "IG-097", "IG-099", "IG-101", ],
  },
  pairs => {
    "HiSeq_FFPE_VS_FF" => {
      groups => [ "HiSeq_FF", "HiSeq_FFPE" ],
      paired => [ "P01",      "P02", "P03", "P04", "P05", "P06", "P07", "P08", "P09", "P10", "P11", "P13", "P14", "P15", "P16", "P17", "P18", "P19", "P20" ]
    },
    "HiSeq_FFPE2_VS_FF2" => {
      groups => [ "HiSeq_FF2", "HiSeq_FFPE2" ],
      paired => [ "P01",       "P02", "P03", "P04", "P05", "P06", "P07", "P08", "P09", "P10", "P11", "P13", "P14", "P15" ]
    },
    "HiSeq_FFPE3_VS_FF3" => {
      groups => [ "HiSeq_FF3", "HiSeq_FFPE3" ],
      paired => [ "P16",       "P17", "P18", "P19", "P20" ]
    },
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
  cutadapt => {
    class      => "Cutadapt",
    perform    => 1,
    target_dir => "${target_dir}/cutadapt",
    option     => "-O 12 -m 30",
    source_ref => "files",
    adapter    => "AGATCGGAAGAG",
    extension  => "_clipped.fastq",
    sh_direct  => 0,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "24",
      "mem"      => "20gb"
    },
  },
  tophat2 => {
    class                => "Alignment::Tophat2",
    perform              => 1,
    target_dir           => "${target_dir}/tophat2",
    option               => "--segment-length 25 -r 0 -p 8",
    source_ref           => "cutadapt",
    bowtie2_index        => $bowtie2_index,
    transcript_gtf       => $transcript_gtf,
    transcript_gtf_index => $transcript_gtf_index,
    rename_bam           => 1,
    sh_direct            => 0,
    pbs                  => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "30gb"
    },
  },
  sortbam => {
    class         => "Samtools::Sort",
    perform       => 1,
    target_dir    => "${target_dir}/sortname",
    option        => "",
    source_ref    => "tophat2",
    sort_by_query => 1,
    sh_direct     => 0,
    pbs           => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "20gb"
    },
  },
  htseqcount => {
    class      => "Count::HTSeqCount",
    perform    => 1,
    target_dir => "${target_dir}/htseqcount",
    option     => "",
    source_ref => "sortbam",
    gff_file   => $transcript_gtf,
    sh_direct  => 0,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  genetable => {
    class         => "CQS::CQSDatatable",
    perform       => 1,
    target_dir    => "${target_dir}/genetable",
    option        => "-p ENS --noheader -o ${task}_gene.count",
    source_ref    => "htseqcount",
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
  deseq2 => {
    class         => "Comparison::DESeq2",
    perform       => 1,
    target_dir    => "${target_dir}/deseq2",
    option        => "",
    source_ref    => "pairs",
    groups_ref    => "groups",
    countfile_ref => "genetable",
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
    target_dir => "${target_dir}/sequencetask",
    option     => "",
    source     => {
      individual => [ "cutadapt", "tophat2", "sortbam", "htseqcount", "fastqc" ],
      summary    => [ "genetable", "deseq2" ],
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
