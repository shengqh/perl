#!/usr/bin/perl
use strict;
use warnings;

use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;
use Hash::Merge qw( merge );

my $target_dir = create_directory_or_die("/scratch/cqs/shengq1/brown/20161003_atacseq_3585_678");
my $email      = "quanhu.sheng\@vanderbilt.edu";

my $cqstools   = "/home/shengq1/cqstools/CQS.Tools.exe";
my $picard_jar = "/scratch/cqs/shengq1/local/bin/picard/picard.jar";
my $gatk_jar   = "/home/shengq1/local/bin/GATK/GenomeAnalysisTK.jar";
my $dbsnp      = "/scratch/cqs/shengq1/references/dbsnp/human_GRCh37_v142_16569_M.vcf";

my $macs2call_option_qvalue = "-f BEDPE --broad -g hs -B -q 0.01 --broad-cutoff 0.01 --nomodel --slocal 20000 --llocal 20000 --keep-dup all";

my $bwa_fasta = "/scratch/cqs/shengq1/references/gencode/hg19/bwa_index_0.7.12/GRCh37.p13.genome.fa";

my $config = {
  general    => { task_name => "3585_678" },
  files      => {
    "3585-JDB-6" => [ "/gpfs21/scratch/cqs/shengq1/brown/data/3585/3585-JDB-6_1_sequence.txt.gz", "/gpfs21/scratch/cqs/shengq1/brown/data/3585/3585-JDB-6_2_sequence.txt.gz" ],
    "3585-JDB-7" => [ "/gpfs21/scratch/cqs/shengq1/brown/data/3585/3585-JDB-7_1_sequence.txt.gz", "/gpfs21/scratch/cqs/shengq1/brown/data/3585/3585-JDB-7_2_sequence.txt.gz" ],
    "3585-JDB-8" => [ "/gpfs21/scratch/cqs/shengq1/brown/data/3585/3585-JDB-8_1_sequence.txt.gz", "/gpfs21/scratch/cqs/shengq1/brown/data/3585/3585-JDB-8_2_sequence.txt.gz" ],
  },
  treatments => {
    "3585-JDB-6" => "3585-JDB-6",
    "3585-JDB-7" => "3585-JDB-7",
    "3585-JDB-8" => "3585-JDB-8",
  },
  fastqc_raw => {
    class      => "QC::FastQC",
    perform    => 1,
    target_dir => "${target_dir}/fastqc_raw",
    option     => "",
    source_ref => "files",
    sh_direct  => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=2",
      "walltime" => "2",
      "mem"      => "40gb"
    },
  },
  fastqc_raw_summary => {
    class      => "QC::FastQCSummary",
    perform    => 1,
    target_dir => "${target_dir}/fastqc_raw",
    option     => "",
    cqstools   => $cqstools,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "2",
      "mem"      => "10gb"
    },
  },
  cutadapt => {
    class      => "Trimmer::Cutadapt",
    perform    => 1,
    target_dir => "${target_dir}/cutadapt",
    option     => "-m 30 --trim-n",
    source_ref => "files",
    adapter    => "CTGTCTCTTATA",
    extension  => "_clipped.fastq.gz",
    pairend    => 1,
    sh_direct  => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "24",
      "mem"      => "20gb"
    },
  },
  fastqc_trimmed => {
    class      => "QC::FastQC",
    perform    => 1,
    target_dir => "${target_dir}/fastqc_trimmed",
    option     => "",
    source_ref => [ "cutadapt", ".fastq.gz" ],
    sh_direct  => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=2",
      "walltime" => "2",
      "mem"      => "40gb"
    },
  },
  fastqc_trimmed_summary => {
    class      => "QC::FastQCSummary",
    perform    => 1,
    target_dir => "${target_dir}/fastqc_trimmed",
    option     => "",
    cqstools   => $cqstools,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "2",
      "mem"      => "10gb"
    },
  },
  fastqlen => {
    class      => "CQS::FastqLen",
    perform    => 1,
    target_dir => "${target_dir}/fastqlen",
    option     => "",
    source_ref => "cutadapt",
    cqstools   => $cqstools,
    sh_direct  => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "24",
      "mem"      => "20gb"
    },
  },
  bwa => {
    class              => "Alignment::BWA",
    perform            => 1,
    target_dir         => "${target_dir}/bwa",
    option             => "",
    bwa_index          => $bwa_fasta,
    picard_jar         => $picard_jar,
    source_ref         => [ "cutadapt", ".fastq.gz" ],
    sort_by_coordinate => 1,
    sh_direct          => 0,
    pbs                => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  bwa_cleanbam => {
    class                   => "ATACseq::CleanBam",
    perform                 => 1,
    target_dir              => "${target_dir}/bwa_cleanbam",
    option                  => "",
    source_ref              => "bwa",
    picard_jar              => $picard_jar,
    remove_chromosome       => "M",
    keep_chromosome         => "chr",
    is_sorted_by_coordinate => 1,
    sh_direct               => 0,
    pbs                     => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "240",
      "mem"      => "40gb"
    },
  },
  bwa_bam2bed => {
    class                   => "Format::Bam2Bed",
    perform                 => 1,
    target_dir              => "${target_dir}/bwa_bam2bed",
    option                  => "",
    source_ref              => "bwa_cleanbam",
    blacklist_file          => "/scratch/cqs/shengq1/references/mappable_region/hg19/wgEncodeDacMapabilityConsensusExcludable.bed",
    is_sorted_by_name       => 0,
    is_paired_end           => 1,
    maximum_fragment_length => 1000,
    minimum_fragment_length => 30,
    sh_direct               => 1,
    pbs                     => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  bwa_macs2callpeak => {
    class      => "Chipseq::MACS2Callpeak",
    perform    => 1,
    target_dir => "${target_dir}/bwa_macs2callpeak",
    option     => $macs2call_option_qvalue,
    source_ref => "bwa_bam2bed",
    sh_direct  => 0,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  bwa_macs2callpeak_bradner_rose => {
    class                => "Chipseq::BradnerRose2",
    perform              => 1,
    target_dir           => "${target_dir}/bwa_macs2callpeak_bradner_rose",
    option               => "",
    source_ref           => "bwa_cleanbam",
    groups_ref           => "treatments",
    pipeline_dir         => "/scratch/cqs/shengq1/local/bin/bradnerlab",
    binding_site_bed_ref => [ "bwa_macs2callpeak", ".bed\$" ],
    genome               => "hg19",
    sh_direct            => 1,
    pbs                  => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  bwa_macs2callpeak_bradner_rose_coltron => {
    class              => "Chipseq::Coltron",
    perform            => 1,
    target_dir         => "${target_dir}/bwa_macs2callpeak_bradner_rose_coltron",
    option             => "",
    source_ref         => "bwa_cleanbam",
    groups_ref         => "treatments",
    enhancer_files_ref => [ "bwa_macs2callpeak_bradner_rose", "_AllEnhancers.table.txt" ],
    genome             => "HG19",
    pipeline_dir       => "/scratch/cqs/shengq1/local/bin/bradnerlab",
    sh_direct          => 1,
    pbs                => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  sequencetask => {
    class      => "CQS::SequenceTask",
    perform    => 1,
    target_dir => "${target_dir}/sequencetask",
    option     => "",
    source     => {
      T1 => [ "fastqc_raw",         "cutadapt",               "fastqc_trimmed",                 "fastqlen", "bwa", "bwa_cleanbam", "bwa_bam2bed", "bwa_macs2callpeak" ],
      T2 => [ "fastqc_raw_summary", "fastqc_trimmed_summary", "bwa_macs2callpeak_bradner_rose", "bwa_macs2callpeak_bradner_rose_coltron" ],

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

#performTask( $config, "bowtie2_pretrim" );

1;

