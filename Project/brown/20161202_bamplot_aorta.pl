#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;
use CQS::FileUtils;
use Data::Dumper;

my $target_dir = create_directory_or_die("/scratch/cqs/shengq1/brown/20161202_bamplot_aorta");

my $cqstools = "/home/shengq1/cqstools/cqstools.exe";

my $email    = "quanhu.sheng\@vanderbilt.edu";
my $task     = "chipseq";
my $plot_gff = "/scratch/cqs/shengq1/brown/20161202_bamplot_aorta/document/locus.gff";

my $config = {
  general   => { task_name => $task },
  bam_files => {
    "Aorta_9097"  => ["/gpfs21/scratch/cqs/shengq1/brown/20161202_bamplot_aorta/data/Alignment_Post_Processing_9097.bam"],
    "Aorta_14713" => ["/gpfs21/scratch/cqs/shengq1/brown/20161202_bamplot_aorta/data/Alignment_Post_Processing_14713.bam"],
    "Aorta_14814" => ["/gpfs21/scratch/cqs/shengq1/brown/20161202_bamplot_aorta/data/Alignment_Post_Processing_14814.bam"],
    "Aorta_15725" => ["/gpfs21/scratch/cqs/shengq1/brown/20161202_bamplot_aorta/data/Alignment_Post_Processing_15725.bam"],
  },
  plotgroups => {
    "Aorta" => [ "Aorta_9097", "Aorta_14713", "Aorta_14814", "Aorta_15725" ],
  },
  colormaps => {
    "Aorta_9097"  => "0,0,0",
    "Aorta_14713" => "0,0,0",
    "Aorta_14814" => "0,0,0",
    "Aorta_15725" => "0,0,0",
  },
  bamplot => {
    class      => "Visualization::Bamplot",
    perform    => 1,
    target_dir => "${target_dir}/bamplot",

    #option     => "-g HG19 -y uniform -r",

    option             => "-g HG19 -y uniform -r --save-temp",
    source_ref         => "bam_files",
    groups_ref         => "plotgroups",
    gff_file           => $plot_gff,
    is_rainbow_color   => 0,
    is_single_pdf      => 0,
    is_draw_individual => 0,
    colors_ref         => "colormaps",
    sh_direct          => 1,
    pbs                => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "1",
      "mem"      => "10gb"
    },
  },
};

performConfig($config);
1;
