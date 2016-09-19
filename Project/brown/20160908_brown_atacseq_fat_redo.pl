#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;
use CQS::FileUtils;
use Data::Dumper;

my $target_dir = create_directory_or_die("/scratch/cqs/shengq1/brown/20160908_atacseq_fat_redo");

my $email = "quanhu.sheng\@vanderbilt.edu";
my $task  = "fat";

my $macs1call_option        = "-p 1e-9 -w -S --space=50";
my $macs2call_option_qvalue = "-f BEDPE --broad -g mm -B -q 0.01 --broad-cutoff 0.01 --nomodel --slocal 20000 --llocal 20000 --keep-dup all";
my $macs2call_option_pvalue = "-f BEDPE --broad -g mm -B -p 1e-9 --broad-cutoff 1e-9 --nomodel --slocal 20000 --llocal 20000 --keep-dup all";

my $config = {
  general => { task_name => $task },
  files   => {
    "SQ1_CHOW"   => ["/gpfs21/scratch/cqs/shengq1/brown/data/fat/20150911_4895_mm9.noChrM.fix.rmdup.sorted.bam"],
    "SQ2_CHOW"   => ["/gpfs21/scratch/cqs/shengq1/brown/data/fat/20150911_4896_mm9.noChrM.fix.rmdup.sorted.bam"],
    "Visc1_CHOW" => ["/gpfs21/scratch/cqs/shengq1/brown/data/fat/20150911_4897_mm9.noChrM.fix.rmdup.sorted.bam"],
    "Visc2_CHOW" => ["/gpfs21/scratch/cqs/shengq1/brown/data/fat/20150911_4898_mm9.noChrM.fix.rmdup.sorted.bam"],
    "SQ1_HFD"    => ["/gpfs21/scratch/cqs/shengq1/brown/data/fat/20150911_4899_mm9.noChrM.fix.rmdup.sorted.bam"],
    "SQ2_HFD"    => ["/gpfs21/scratch/cqs/shengq1/brown/data/fat/20150911_4900_mm9.noChrM.fix.rmdup.sorted.bam"],
    "Visc1_HFD"  => ["/gpfs21/scratch/cqs/shengq1/brown/data/fat/20150911_4901_mm9.noChrM.fix.rmdup.sorted.bam"],
    "Visc2_HFD"  => ["/gpfs21/scratch/cqs/shengq1/brown/data/fat/20150911_4902_mm9.noChrM.fix.rmdup.sorted.bam"],
  },
  treatments => {
    "SQ1_CHOW"   => "SQ1_CHOW",
    "SQ2_CHOW"   => "SQ2_CHOW",
    "Visc1_CHOW" => "Visc1_CHOW",
    "Visc2_CHOW" => "Visc2_CHOW",
    "SQ1_HFD"    => "SQ1_HFD",
    "SQ2_HFD"    => "SQ2_HFD",
    "Visc1_HFD"  => "Visc1_HFD",
    "Visc2_HFD"  => "Visc2_HFD",
  },
  replicates => {
    "SQ_CHOW"   => [ "SQ1_CHOW",   "SQ2_CHOW" ],
    "Visc_CHOW" => [ "Visc1_CHOW", "Visc2_CHOW" ],
    "SQ_HFD"    => [ "SQ1_HFD",    "SQ2_HFD" ],
    "Visc_HFD"  => [ "Visc1_HFD",  "Visc2_HFD" ],
  },
  replicates_comparison => {
    "SQ_Visc_CHOW" => [ "SQ_CHOW", "Visc_CHOW" ],
    "SQ_Visc_HFD"  => [ "SQ_HFD",  "Visc_HFD" ],
  },
  individual_comparison => {
    "SQ1_Visc1_CHOW" => [ "SQ1_CHOW", "Visc1_CHOW" ],
    "SQ2_Visc2_CHOW" => [ "SQ2_CHOW", "Visc2_CHOW" ],
    "SQ1_Visc1_HFD"  => [ "SQ1_HFD",  "Visc1_HFD" ],
    "SQ2_Visc2_HFD"  => [ "SQ2_HFD",  "Visc2_HFD" ],
    "SQ1_Visc2_CHOW" => [ "SQ1_CHOW", "Visc2_CHOW" ],
    "SQ2_Visc1_CHOW" => [ "SQ2_CHOW", "Visc1_CHOW" ],
    "SQ1_Visc2_HFD"  => [ "SQ1_HFD",  "Visc2_HFD" ],
    "SQ2_Visc1_HFD"  => [ "SQ2_HFD",  "Visc1_HFD" ],
  },
  replicates_group => {
    "SQ_Visc_CHOW" => [ "SQ1_CHOW", "SQ2_CHOW", "Visc1_CHOW", "Visc2_CHOW" ],
    "SQ_Visc_HFD"  => [ "SQ1_HFD",  "SQ2_HFD",  "Visc1_HFD",  "Visc2_HFD" ],
  },
  bam2bed => {
    class                   => "Format::Bam2Bed",
    perform                 => 1,
    target_dir              => "${target_dir}/bam2bed",
    option                  => "",
    source_ref              => "files",
    blacklist_file          => "/scratch/cqs/shengq1/references/mappable_region/mm9/mm9-blacklist.bed",
    is_sorted_by_name       => 0,
    is_paired_end           => 1,
    maximum_fragment_length => 1000,
    minimum_fragment_length => 30,
    sh_direct               => 1,
    pbs                     => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  macs2callpeak_individual_nomodel => {
    class      => "Chipseq::MACS2Callpeak",
    perform    => 1,
    target_dir => "${target_dir}/macs2callpeak_individual_nomodel",
    option     => $macs2call_option_pvalue,
    source_ref => "bam2bed",
    sh_direct  => 0,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  macs2bdgdiff_individual_nomodel => {
    class      => "Chipseq::MACS2Bdgdiff",
    perform    => 1,
    target_dir => "${target_dir}/macs2bdgdiff_individual_nomodel",
    option     => "",
    source_ref => "macs2callpeak_individual_nomodel",
    groups_ref => "individual_comparison",
    sh_direct  => 0,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  macs2callpeak_replicates_nomodel => {
    class      => "Chipseq::MACS2Callpeak",
    perform    => 1,
    target_dir => "${target_dir}/macs2callpeak_replicates_nomodel",
    option     => $macs2call_option_pvalue,
    source_ref => "bam2bed",
    groups_ref => "replicates",
    sh_direct  => 0,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  macs2bdgdiff_replicates_nomodel => {
    class      => "Chipseq::MACS2Bdgdiff",
    perform    => 1,
    target_dir => "${target_dir}/macs2bdgdiff_replicates_nomodel",
    option     => "",
    source_ref => "macs2callpeak_replicates_nomodel",
    groups_ref => "replicates_comparison",
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
    source_ref           => "files",
    groups_ref           => "treatments",
    pipeline_dir         => "/scratch/cqs/shengq1/local/bin/bradnerlab",
    binding_site_bed_ref => [ "macs2callpeak_individual_nomodel", ".bed\$" ],
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
    source_ref         => "files",
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
      "step1" => [ "bam2bed",                         "macs2callpeak_individual_nomodel", "macs2callpeak_replicates_nomodel" ],
      "step2" => [ "macs2bdgdiff_individual_nomodel", "macs2bdgdiff_replicates_nomodel", "bwa_macs2callpeak_bradner_rose", "bwa_macs2callpeak_bradner_rose_coltron" ],
    },
    sh_direct => 0,
    pbs       => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  }
};

performConfig($config);

1;
