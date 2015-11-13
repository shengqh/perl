#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;
use CQS::PerformSmallRNA;
use CQS::FileUtils;
use Hash::Merge qw( merge );

my $userdef = merge(
  {

    #General options
    task_name  => "2797_mouse",
    email      => "quanhu.sheng\@vanderbilt.edu",
    target_dir => create_directory_or_die("/scratch/cqs/shengq1/vickers/20150911_parclip_gsnap_2797-mouse/"),
    max_thread => 8,
    cluster    => "slurm",

    #Default software parameter (don't change it except you really know it)
    fastq_remove_N => 0,
    adapter        => "TGGAATTCTCGGGTGCCAAGG",
    cqstools       => "/home/shengq1/cqstools/CQS.Tools.exe",

    #Data
    files => {
      "RPI47" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/2797_demultiplexing/2797-KCV-1_RPI47_Ago2MIN6Huh7.fastq.gz"],
      "RPI48" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/2797_demultiplexing/2797-KCV-1_RPI48_Ago3MIN6Huh7.fastq.gz"],
    },
  },
  hg19_3utr()
);

my $def = getSmallRNADefinition( $userdef, mm10_genome() );

#print Dumper($def);

performParclipSmallRNA($def);

1;

