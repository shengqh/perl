#!/usr/bin/perl
use strict;
use warnings;

use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ClassFactory;
use CQS::CNV;

my $target_dir = create_directory_or_die("/scratch/cqs/shengq1/chipseq/20140224_jennifer_chipseq_clayton");

my $fasta_file = "/data/cqs/guoy1/reference/hg19/bwa_index_0.7.4/hg19_chr.fa";
my $dbsnp      = "/data/cqs/shengq1/reference/snp137/human/00-All.vcf";
my $gatk       = "/home/shengq1/local/bin/GATK/GenomeAnalysisTK.jar";
my $picard_dir = "/home/shengq1/local/bin/picard/";

my $email = "quanhu.sheng\@vanderbilt.edu";

my $config = {
  general    => { task_name => "chipseq_clayton" },
  fastqfiles => {
    "2653-JP-34_S1" => [
      "/gpfs21/scratch/cqs/shengq1/chipseq/20140224_jennifer_chipseq_clayton/data/analysis_5049046/2653-JP-34_S1_L001_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/chipseq/20140224_jennifer_chipseq_clayton/data/analysis_5049046/2653-JP-34_S1_L001_R2_001.fastq.gz"
    ],
    "2653-JP-35_S2" => [
      "/gpfs21/scratch/cqs/shengq1/chipseq/20140224_jennifer_chipseq_clayton/data/analysis_5049046/2653-JP-35_S2_L001_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/chipseq/20140224_jennifer_chipseq_clayton/data/analysis_5049046/2653-JP-35_S2_L001_R2_001.fastq.gz"
    ],
    "1806-p63IP_S2" => [
      "/gpfs21/scratch/cqs/shengq1/chipseq/20140224_jennifer_chipseq_clayton/data/analysis_5011011/1806-p63IP_S2_L001_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/chipseq/20140224_jennifer_chipseq_clayton/data/analysis_5011011/1806-p63IP_S2_L001_R2_001.fastq.gz"
    ],
    "1806-p73IP_S1" => [
      "/gpfs21/scratch/cqs/shengq1/chipseq/20140224_jennifer_chipseq_clayton/data/analysis_5011011/1806-p73IP_S1_L001_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/chipseq/20140224_jennifer_chipseq_clayton/data/analysis_5011011/1806-p73IP_S1_L001_R2_001.fastq.gz"
    ],
  },
  groups => {
    "Control" => ["2653-JP-34_S1"],
    "POL2IP"  => ["2653-JP-35_S2"],
    "p63IP"   => ["1806-p63IP_S2"],
    "p73IP"   => ["1806-p73IP_S1"],
  },
  pairs => {
    "POL2IP_vs_Control" => { "Control", "POL2IP" },
    "p73IP_vs_Control"  => { "Control", "p63IP" },
    "p63IP_vs_Control"  => { "Control", "p73IP" },
  },
  pretrim_fastqc => {
    class      => "FastQC",
    perform    => 0,
    target_dir => "${target_dir}/pretrim_fastqc",
    option     => "",
    source_ref => "fastqfiles",
    sh_direct  => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "2",
      "mem"      => "10gb"
    },
  },
  pretrim_bwa => {
    class      => "BWA",
    perform    => 0,
    target_dir => "${target_dir}/pretrim_bwa",
    option     => "-t 8",
    fasta_file => $fasta_file,
    source_ref => "fastqfiles",
    sh_direct  => 0,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  pretrim_markdup => {
    class              => "Picard::MarkDuplicates",
    perform            => 0,
    target_dir         => "${target_dir}/pretrim_markdup",
    option             => "-Xmx40g",
    source_ref         => "pretrim_bwa",
    thread_count       => 8,
    markDuplicates_jar => "${picard_dir}/MarkDuplicates.jar",
    sh_direct          => 0,
    pbs                => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  pretrim_refine => {
    class              => "GATKRefine",
    perform            => 0,
    target_dir         => "${target_dir}/pretrim_refine",
    option             => "-Xmx40g",
    fasta_file         => $fasta_file,
    source_ref         => "pretrim_bwa",
    thread_count       => 8,
    vcf_files          => [$dbsnp],
    gatk_jar           => $gatk,
    markDuplicates_jar => "${picard_dir}/MarkDuplicates.jar",
    sh_direct          => 0,
    pbs                => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  trim_scythe => {
    class        => "Trimmer::Scythe",
    perform      => 1,
    target_dir   => "${target_dir}/trim_scythe",
    option       => "-p 0.1",
    source_ref   => "fastqfiles",
    adapter_file => "/scratch/cqs/shengq1/local/bin/illumina_adapters.fa",
    sh_direct    => 1,
    pbs          => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "24",
      "mem"      => "20gb"
    },
  },
  trim_sickle => {
    class      => "Trimmer::Sickle",
    perform    => 1,
    target_dir => "${target_dir}/trim_sickle",
    option     => "",
    qual_type  => "sanger",                      #Type of quality values (solexa (CASAVA < 1.3), illumina (CASAVA 1.3 to 1.7), sanger (which is CASAVA >= 1.8))
    source_ref => "trim_scythe",
    sh_direct  => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "24",
      "mem"      => "20gb"
    },
  },
  posttrim_fastqc => {
    class      => "FastQC",
    perform    => 1,
    target_dir => "${target_dir}/posttrim_fastqc",
    option     => "",
    source_ref => "trim_sickle",
    sh_direct  => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "2",
      "mem"      => "10gb"
    },
  },
  posttrim_bwa => {
    class      => "BWA",
    perform    => 1,
    target_dir => "${target_dir}/posttrim_bwa",
    option     => "-t 8",
    fasta_file => $fasta_file,
    source_ref => "trim_sickle",
    sh_direct  => 0,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  posttrim_markdup => {
    class              => "Picard::MarkDuplicates",
    perform            => 1,
    target_dir         => "${target_dir}/posttrim_markdup",
    option             => "-Xmx40g",
    source_ref         => "posttrim_bwa",
    thread_count       => 8,
    markDuplicates_jar => "${picard_dir}/MarkDuplicates.jar",
    sh_direct          => 0,
    pbs                => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  posttrim_refine => {
    class              => "GATKRefine",
    perform            => 1,
    target_dir         => "${target_dir}/posttrim_refine",
    option             => "-Xmx40g",
    fasta_file         => $fasta_file,
    source_ref         => "posttrim_bwa",
    thread_count       => 8,
    vcf_files          => [$dbsnp],
    gatk_jar           => $gatk,
    markDuplicates_jar => "${picard_dir}/MarkDuplicates.jar",
    sh_direct          => 0,
    pbs                => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  homerMakeTagDirectory => {
    class      => "Homer::MakeTagDirectory",
    perform    => 1,
    target_dir => "${target_dir}/homerMakeTagDirectory",
    option     => "",
    source_ref => "pretrim_markdup",
    groups_ref => "groups",
    sh_direct  => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  homerFindPeaks => {
    class      => "Homer::FindPeaks",
    perform    => 1,
    target_dir => "${target_dir}/homerFindPeaks",
    option     => "",
    source_ref => "homerMakeTagDirectory",
    pairs_ref => "pairs",
    sh_direct  => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  overall => {
    class      => "CQS::SequenceTask",
    perform    => 1,
    target_dir => "${target_dir}/overall",
    option     => "",
    source     => {
      individual => [ "pretrim_fastqc", "pretrim_bwa", "pretrim_markdup", "pretrim_refine", "trim_scythe", "trim_sickle", "posttrim_fastqc", "posttrim_bwa", "posttrim_markdup", "posttrim_refine" ]
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

#performConfig($config);
#performTask( $config, "homerMakeTagDirectory" );
performTask( $config, "homerFindPeaks" );

1;
