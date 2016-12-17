#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;
use CQS::FileUtils;
use Data::Dumper;

my $target_dir = create_directory_or_die("/scratch/cqs/shengq1/brown/20161215_brown_chipseq_3593_9to10");
my $task       = "chipseq_3593_9to12";

my $fasta_file   = "/scratch/cqs/shengq1/references/gencode/hg19/bowtie_index_1.1.2/GRCh37.p13.genome.fa";
my $bowtie_index = "/scratch/cqs/shengq1/references/gencode/hg19/bowtie_index_1.1.2/GRCh37.p13.genome";
my $cqstools     = "/home/shengq1/cqstools/cqstools.exe";
my $plot_gff     = "/scratch/cqs/shengq1/chipseq/20160823_janathan_chipseq_195R3_gse53999_bamplot/config/H3K27ac.gff";

my $email = "quanhu.sheng\@vanderbilt.edu";

my $config = {
  general => { task_name => $task },
  files   => {
    "H3K27ac_Input" => ["/gpfs21/scratch/cqs/shengq1/brown/data/3593/3593-JDB-9_S1_R1_001.fastq.gz"],
    "H3K27ac"       => ["/gpfs21/scratch/cqs/shengq1/brown/data/3593/3593-JDB-10_S2_R1_001.fastq.gz"],
  },
  treatments => {
    "H3K27ac" => ["H3K27ac"],
  },
  controls => {
    "H3K27ac" => ["H3K27ac_Input"],
  },
  plotgroups => {
    "H3K27ac" => [ "H3K27ac", "H3K27ac_Input" ],
  },
  fastqc_raw => {
    class      => "QC::FastQC",
    perform    => 1,
    target_dir => "${target_dir}/fastqc_raw",
    option     => "",
    source_ref => "files",
    sh_direct  => 0,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "2",
      "mem"      => "10gb"
    },
  },
  fastqc_raw_summary => {
    class      => "QC::FastQCSummary",
    perform    => 1,
    sh_direct  => 1,
    target_dir => "${target_dir}/fastqc_raw",
    cqstools   => $cqstools,
    option     => "",
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
    source_ref => ["files"],
    adapter    => "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",    #trueseq adapter
    extension  => "_clipped.fastq",
    sh_direct  => 0,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "24",
      "mem"      => "20gb"
    },
  },
  fastqc_post_trim => {
    class      => "QC::FastQC",
    perform    => 1,
    target_dir => "${target_dir}/fastqc_post_trim",
    option     => "",
    sh_direct  => 1,
    source_ref => [ "cutadapt", ".fastq.gz" ],
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "2",
      "mem"      => "10gb"
    },
  },
  fastqc_post_trim_summary => {
    class      => "QC::FastQCSummary",
    perform    => 1,
    sh_direct  => 1,
    target_dir => "${target_dir}/fastqc_post_trim",
    cqstools   => $cqstools,
    option     => "",
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "2",
      "mem"      => "10gb"
    },
  },
  fastq_len => {
    class      => "CQS::FastqLen",
    perform    => 1,
    target_dir => "$target_dir/fastq_len",
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
  bowtie1 => {
    class                   => "Alignment::Bowtie1",
    perform                 => 1,
    target_dir              => "${target_dir}/bowtie1",
    option                  => "-v 1 -m 1 --best --strata",
    fasta_file              => $fasta_file,
    source_ref              => [ "cutadapt", ".fastq.gz" ],
    bowtie1_index           => $bowtie_index,
    chromosome_grep_pattern => "\"^chr\"",
    sh_direct               => 0,
    pbs                     => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  macs1callpeak => {
    class        => "Chipseq::MACS",
    perform      => 1,
    target_dir   => "${target_dir}/macs1callpeak",
    option       => "-p 1e-9 -w -S --space=50",
    source_ref   => [ "bowtie1", ".bam\$" ],
    groups_ref   => "treatments",
    controls_ref => "controls",
    sh_direct    => 0,
    pbs          => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  macs1callpeak_bradner_rose2 => {
    class                => "Chipseq::BradnerRose2",
    perform              => 1,
    target_dir           => "${target_dir}/macs1callpeak_bradner_rose2",
    option               => "",
    source_ref           => "bowtie1",
    groups_ref           => "treatments",
    controls_ref         => "controls",
    pipeline_dir         => "/scratch/cqs/shengq1/local/bin/bradnerlab",
    binding_site_bed_ref => [ "macs1callpeak", ".bed\$" ],
    genome               => "hg19",
    sh_direct            => 1,
    pbs                  => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  bamplot => {
    class      => "Visualization::Bamplot",
    perform    => 1,
    target_dir => "${target_dir}/bamplot",
    option     => "-g HG19 -y uniform -r",

    #option        => "-g HG19 -y uniform -r --save-temp",
    source_ref         => "bowtie1",
    groups_ref         => "plotgroups",
    gff_file           => $plot_gff,
    is_rainbow_color   => 0,
    is_draw_individual => 0,
    is_single_pdf      => 1,
    sh_direct          => 1,
    pbs                => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "1",
      "mem"      => "10gb"
    },
  },
  sequencetask => {
    class      => "CQS::SequenceTask",
    perform    => 1,
    target_dir => "${target_dir}/sequencetask",
    option     => "",
    source     => {
      step_1 => [ "fastqc_raw",         "cutadapt",                 "fastqc_post_trim", "fastq_len",                   "bowtie1" ],
      step_2 => [ "fastqc_raw_summary", "fastqc_post_trim_summary", "macs1callpeak",    "macs1callpeak_bradner_rose2", "bamplot" ],
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
