#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;
use CQS::FileUtils;
use Data::Dumper;

my $target_dir = create_directory_or_die("/scratch/cqs/shengq1/brown/20160901_chipseq_ips");
my $task       = "chipseq_ips";

my $fasta_file     = "/scratch/cqs/shengq1/references/gencode/hg19/bowtie_index_1.1.2/GRCh37.p13.genome.fa";
my $bowtie_index   = "/scratch/cqs/shengq1/references/gencode/hg19/bowtie_index_1.1.2/GRCh37.p13.genome";
my $bwa_index      = "/scratch/cqs/shengq1/references/gencode/hg19/bwa_index_0.7.12/GRCh37.p13.genome.fa";
my $cqstools       = "/home/shengq1/cqstools/cqstools.exe";
my $picard_jar     = "/scratch/cqs/shengq1/local/bin/picard/picard.jar";
my $transcript_gtf = "/scratch/cqs/shengq1/references/gencode/hg19/gencode.v19.chr_patch_hapl_scaff.annotation.gtf";
my $qc3_perl       = "/scratch/cqs/shengq1/local/bin/qc3/qc3.pl";
my $plot_gff       = "/scratch/cqs/shengq1/chipseq/20160823_janathan_chipseq_195R3_gse53999_bamplot/config/H3K27ac.gff";

my $email = "quanhu.sheng\@vanderbilt.edu";

my $config = {
  general => { task_name => $task },
  files   => {
    "IPS_EC1_WCE"     => ["/gpfs21/scratch/cqs/shengq1/brown/20160901_chipseq_ips/data/CAGATC-s_2_1_sequence.txt.gz"],
    "IPS_EC1_H3K27ac" => ["/gpfs21/scratch/cqs/shengq1/brown/20160901_chipseq_ips/data/ACAGTG-s_2_1_sequence.txt.gz"],
    "IPS_EC2_WCE"     => ["/gpfs21/scratch/cqs/shengq1/brown/20160901_chipseq_ips/data/ACTTGA-s_2_1_sequence.txt.gz"],
    "IPS_EC2_H3K27ac" => ["/gpfs21/scratch/cqs/shengq1/brown/20160901_chipseq_ips/data/GCCAAT-s_2_1_sequence.txt.gz"],
  },
  treatments => {
    "IPS_EC1_H3K27ac" => ["IPS_EC1_H3K27ac"],
    "IPS_EC2_H3K27ac" => ["IPS_EC2_H3K27ac"],
  },
  controls => {
    "IPS_EC1_H3K27ac" => ["IPS_EC1_WCE"],
    "IPS_EC2_H3K27ac" => ["IPS_EC2_WCE"],
  },
  plotgroups => {
    "IPS_H3K27ac" => [ "IPS_EC1_H3K27ac", "IPS_EC1_WCE", "IPS_EC2_H3K27ac", "IPS_EC2_WCE" ],
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
  bowtie1 => {
    class                   => "Alignment::Bowtie1",
    perform                 => 1,
    target_dir              => "${target_dir}/bowtie1",
    option                  => "-v 1 -m 1 --best --strata",
    fasta_file              => $fasta_file,
    source_ref              => "files",
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
    source_ref   => "bowtie1",
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
    source_ref    => "bowtie1",
    groups_ref    => "plotgroups",
    gff_file      => $plot_gff,
    rainbow_color => 0,
    sh_direct     => 1,
    pbs           => {
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
      step_1 => [ "fastqc_raw",         "bowtie1" ],
      step_2 => [ "fastqc_raw_summary", "macs1callpeak", "macs1callpeak_bradner_rose2", "bamplot" ],
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
