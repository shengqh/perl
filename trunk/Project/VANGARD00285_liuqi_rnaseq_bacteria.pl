#!/usr/bin/perl
use strict;
use warnings;

use CQS::QC;
use CQS::DNASeq;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::SomaticMutation;
use CQS::ClassFactory;

my $vangard = "VANGARD00285";

my $target_dir = create_directory_or_die("/scratch/cqs/shengq1/vangard/${vangard}_liuqi_rnaseq_bacteria");

my $email    = "quanhu.sheng\@vanderbilt.edu";
my $cqstools = "/home/shengq1/cqstools/CQS.Tools.exe";

my $config = {
  general    => { task_name => "${vangard}" },
  fastqfiles => {
    "2763-EPS-01" => ["/autofs/blue_sequencer/Runs/projects/2763-EPS/2014-02-06/2763-EPS-1_1.fastq.gz"],
    "2763-EPS-02" => ["/autofs/blue_sequencer/Runs/projects/2763-EPS/2014-02-06/2763-EPS-2_1.fastq.gz"],
    "2763-EPS-03" => ["/autofs/blue_sequencer/Runs/projects/2763-EPS/2014-02-06/2763-EPS-3_1.fastq.gz"],
    "2763-EPS-04" => ["/autofs/blue_sequencer/Runs/projects/2763-EPS/2014-02-06/2763-EPS-4_1.fastq.gz"],
    "2763-EPS-05" => ["/autofs/blue_sequencer/Runs/projects/2763-EPS/2014-02-06/2763-EPS-5_1.fastq.gz"],
    "2763-EPS-06" => ["/autofs/blue_sequencer/Runs/projects/2763-EPS/2014-02-06/2763-EPS-6_1.fastq.gz"],
    "2763-EPS-07" => ["/autofs/blue_sequencer/Runs/projects/2763-EPS/2014-02-06/2763-EPS-7_1.fastq.gz"],
    "2763-EPS-08" => ["/autofs/blue_sequencer/Runs/projects/2763-EPS/2014-02-06/2763-EPS-8_1.fastq.gz"],
    "2763-EPS-09" => ["/autofs/blue_sequencer/Runs/projects/2763-EPS/2014-02-06/2763-EPS-9_1.fastq.gz"],
    "2763-EPS-11" => ["/autofs/blue_sequencer/Runs/projects/2763-EPS/2014-02-06/2763-EPS-11_1.fastq.gz"],
    "2763-EPS-12" => ["/autofs/blue_sequencer/Runs/projects/2763-EPS/2014-02-06/2763-EPS-12_1.fastq.gz"],
    "2763-EPS-14" => ["/autofs/blue_sequencer/Runs/projects/2763-EPS/2014-02-06/2763-EPS-14_1.fastq.gz"],
  },
  groups => {
    "isdGisdI_plgt"      => [ "2763-EPS-01", "2763-EPS-02", "2763-EPS-03" ],
    "isdGisdI_plgt_IsdI" => [ "2763-EPS-04", "2763-EPS-05", "2763-EPS-06" ],
    "isdGisdI_plgt_hmuO" => [ "2763-EPS-07", "2763-EPS-08", "2763-EPS-09" ],
    "isdGisdI_plgt_mhuD" => [ "2763-EPS-11", "2763-EPS-12", "2763-EPS-14" ],
  },
  pairs => {
    "IsdI" => [ "isdGisdI_plgt", "isdGisdI_plgt_IsdI" ],
    "hmuO" => [ "isdGisdI_plgt", "isdGisdI_plgt_hmuO" ],
    "mhuD" => [ "isdGisdI_plgt", "isdGisdI_plgt_mhuD" ],
  },
  fastqc => {
    class      => "FastQC",
    perform    => 1,
    target_dir => "${target_dir}/fastqc",
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
  trimmer => {
    class      => "CQS::FastqTrimmer",
    perform    => 1,
    target_dir => "${target_dir}/trimmer",
    option     => "-n",
    extension  => "_trim.fastq",
    source_ref => "fastqfiles",
    cqstools   => $cqstools,
    sh_direct  => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "2",
      "mem"      => "10gb"
    },
  },
  fastqlen => {
    class      => "FastqLen",
    perform    => 1,
    target_dir => "${target_dir}/trimlen",
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
  rockhopper => {
    class          => "Bacteria::RockHopper",
    perform        => 1,
    target_dir     => "${target_dir}/rockhopper",
    source_ref     => "trimmer",
    groups_ref     => "groups",
    pairs_ref      => "pairs",
    java_option    => "-Xmx10g",
    rockhopper_jar => "/scratch/cqs/shengq1/local/bin/Rockhopper.jar",
    genome_dir     => "/data/cqs/shengq1/reference/bacteria",
    option         => "-p 8 -TIME -v true -c true",
    sh_direct      => 0,
    pbs            => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "10gb"
    },
  },
  bowtie2 => {
    class         => "Bowtie2",
    perform       => 1,
    target_dir    => "${target_dir}/bowtie2",
    source_ref    => "trimmer",
    bowtie2_index => "${target_dir}/genome/bowtie2-index/NC_005945",
    option        => "-p 8",
    sh_direct     => 0,
    pbs           => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "10gb"
    },
  },
  overall => {
    class      => "CQS::SequenceTask",
    perform    => 1,
    target_dir => "${target_dir}/overall",
    option     => "",
    source     => { individual => [ "fastqc", "trimmer", "fastqlen", ], },
    sh_direct  => 0,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },

};

#performConfig($config);
performTask($config, "rockhopper");

1;
