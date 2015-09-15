#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;
use CQS::PerformSmallRNA;
use CQS::FileUtils;

my $def = {

  #General options
  task_name  => "2797_rat",
  email      => "quanhu.sheng\@vanderbilt.edu",
  target_dir => create_directory_or_die("/scratch/cqs/shengq1/vickers/Vickers_20150911_parclip_gsnap_2797-rat/"),
  max_thread => 8,
  cluster    => "slurm",

  cqstools => "/home/shengq1/cqstools/CQS.Tools.exe",

  gsnap_index_directory => "/scratch/cqs/shengq1/references/rn5/STAR_index_v38.81_2.4.2a_sjdb49/",
  gsnap_index_name      => "rn5",

  #Data
  files => {
    "RPI40" => [ "/gpfs21/scratch/cqs/shengq1/vickers/20150401_parclip_2797_rat_mouse_human/rn5/identical_NTA/result/RPI40_clipped_identical_NTA.fastq.gz", ],
    "RPI41" => [ "/gpfs21/scratch/cqs/shengq1/vickers/20150401_parclip_2797_rat_mouse_human/rn5/identical_NTA/result/RPI41_clipped_identical_NTA.fastq.gz", ],
    "RPI42" => [ "/gpfs21/scratch/cqs/shengq1/vickers/20150401_parclip_2797_rat_mouse_human/rn5/identical_NTA/result/RPI42_clipped_identical_NTA.fastq.gz", ],
    "RPI43" => [ "/gpfs21/scratch/cqs/shengq1/vickers/20150401_parclip_2797_rat_mouse_human/rn5/identical_NTA/result/RPI43_clipped_identical_NTA.fastq.gz", ],
  },
  counts => {
    "RPI40" => ["/gpfs21/scratch/cqs/shengq1/vickers/20150401_parclip_2797_rat_mouse_human/rn5/identical_NTA/result/RPI40_clipped_identical_NTA.fastq.gz.dupcount"],
    "RPI41" => ["/gpfs21/scratch/cqs/shengq1/vickers/20150401_parclip_2797_rat_mouse_human/rn5/identical_NTA/result/RPI41_clipped_identical_NTA.fastq.gz.dupcount"],
    "RPI42" => ["/gpfs21/scratch/cqs/shengq1/vickers/20150401_parclip_2797_rat_mouse_human/rn5/identical_NTA/result/RPI42_clipped_identical_NTA.fastq.gz.dupcount"],
    "RPI43" => ["/gpfs21/scratch/cqs/shengq1/vickers/20150401_parclip_2797_rat_mouse_human/rn5/identical_NTA/result/RPI43_clipped_identical_NTA.fastq.gz.dupcount"],
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
  'gsnap_smallRNA_count' => {
    'cluster'         => 'slurm',
    'fasta_file'      => "/scratch/cqs/shengq1/references/smallrna/rn5_miRBase21_ucsc-tRNA_ensembl78.bed.fa",
    'sh_direct'       => 1,
    'perform'         => 1,
    'target_dir'      => $def->{target_dir} . "/gsnap_smallRNA_count",
    'fastq_files'     => $def->{files},
    'coordinate_file' => "/scratch/cqs/shengq1/references/smallrna/rn5_miRBase21_ucsc-tRNA_ensembl78.bed",
    'source_ref'      => 'gsnap',
    'cqs_tools'       => $def->{cqstools},
    'seqcount'        => $def->{counts},
    'samtools'        => '/scratch/cqs/shengq1/local/bin/samtools',
    'class'           => 'CQS::SmallRNACount',
    'option'          => '-s -e 4',
    'pbs'             => {
      'email'    => 'quanhu.sheng@vanderbilt.edu',
      'walltime' => '72',
      'mem'      => '40gb',
      'nodes'    => '1:ppn=1'
    },
  },
};

performConfig($config);

1;

