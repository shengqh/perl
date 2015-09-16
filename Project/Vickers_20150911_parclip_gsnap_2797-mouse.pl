#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;
use CQS::PerformSmallRNA;
use CQS::FileUtils;

my $def = {

  #General options
  task_name  => "2797_mouse",
  email      => "quanhu.sheng\@vanderbilt.edu",
  target_dir => create_directory_or_die("/scratch/cqs/shengq1/vickers/Vickers_20150911_parclip_gsnap_2797-mouse/"),
  max_thread => 8,
  cluster    => "slurm",

  cqstools => "/home/shengq1/cqstools/CQS.Tools.exe",

  gsnap_index_directory => "/scratch/cqs/shengq1/references/mm10/gsnap_index_2015-06-23/",
  gsnap_index_name      => "mm10",

  #Data
  files => {
    "RPI47" => [ "/gpfs21/scratch/cqs/shengq1/vickers/20150401_parclip_2797_rat_mouse_human/mm10/identical_NTA/result/RPI47_clipped_identical_NTA.fastq.gz", ],
    "RPI48" => [ "/gpfs21/scratch/cqs/shengq1/vickers/20150401_parclip_2797_rat_mouse_human/mm10/identical_NTA/result/RPI48_clipped_identical_NTA.fastq.gz", ],
  },
  counts => {
    "RPI47" => ["/gpfs21/scratch/cqs/shengq1/vickers/20150401_parclip_2797_rat_mouse_human/mm10/identical_NTA/result/RPI47_clipped_identical_NTA.fastq.gz.dupcount"],
    "RPI48" => ["/gpfs21/scratch/cqs/shengq1/vickers/20150401_parclip_2797_rat_mouse_human/mm10/identical_NTA/result/RPI48_clipped_identical_NTA.fastq.gz.dupcount"],
  },
};

my $config = {
  general => { "task_name" => $def->{task_name}, },
  gsnap   => {
    class                 => "Alignment::Gsnap",
    perform               => 1,
    target_dir            => $def->{target_dir} . "/gsnap",
    option                => "-y 0 -z 0 -Y 0 -Z 0 -m 1 -Q --trim-mismatch-score 0 --trim-indel-score 0 --mode ttoc-nonstranded --gunzip",
    gsnap_index_directory => $def->{gsnap_index_directory},
    gsnap_index_name      => $def->{gsnap_index_name},
    source                => $def->{files},
    sh_direct             => 1,
    pbs                   => {
      "email"    => $def->{email},
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    }
  },
  gsnap_smallRNA_count => {
    class           => 'CQS::SmallRNACount',
    perform         => 1,
    target_dir      => $def->{target_dir} . "/gsnap_smallRNA_count",
    option          => '-s -e 4',
    cluster         => 'slurm',
    sh_direct       => 1,
    fastq_files     => $def->{files},
    coordinate_file => "/scratch/cqs/shengq1/references/smallrna/mm10_miRBase21_ucsc-tRNA_ensembl78.bed",
    fasta_file      => "/scratch/cqs/shengq1/references/smallrna/mm10_miRBase21_ucsc-tRNA_ensembl78.bed.fa",
    source_ref      => 'gsnap',
    cqs_tools       => $def->{cqstools},
    seqcount        => $def->{counts},
    samtools        => '/scratch/cqs/shengq1/local/bin/samtools',
    pbs             => {
      'email'    => 'quanhu.sheng@vanderbilt.edu',
      'walltime' => '72',
      'mem'      => '40gb',
      'nodes'    => '1:ppn=1'
    },
  },
};

performConfig($config);

1;

