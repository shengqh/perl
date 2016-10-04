#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;
use CQS::FileUtils;
use Data::Dumper;

my $target_dir = create_directory_or_die("/scratch/cqs/shengq1/brown/20161003_bamplot_tissues");

my $cqstools = "/home/shengq1/cqstools/cqstools.exe";

my $email    = "quanhu.sheng\@vanderbilt.edu";
my $task     = "chipseq";
my $plot_gff = "/scratch/cqs/shengq1/brown/20161003_bamplot_tissues/document/locus.gff";

my $config = {
  general   => { task_name => $task },
  bam_files => {
    "Adipose nuclei" => ["/gpfs21/scratch/cqs/shengq1/brown/20161003_bamplot_tissues/H3K27Ac_multiple_tissues/Alignment_Post_Processing_1500.bam"],
    "Adult CD14"     => ["/gpfs21/scratch/cqs/shengq1/brown/20161003_bamplot_tissues/H3K27Ac_multiple_tissues/Alignment_Post_Processing_543.bam"],
    "Angular gyrus"  => ["/gpfs21/scratch/cqs/shengq1/brown/20161003_bamplot_tissues/H3K27Ac_multiple_tissues/Alignment_Post_Processing_702.bam"],
    "Aorta"          => ["/gpfs21/scratch/cqs/shengq1/brown/20161003_bamplot_tissues/H3K27Ac_multiple_tissues/Alignment_Post_Processing_9097.bam"],

    #"Aorta 2"           => ["/gpfs21/scratch/cqs/shengq1/brown/20161003_bamplot_tissues/H3K27Ac_multiple_tissues/Alignment_Post_Processing_14713.bam"],
    #"Aorta 3"           => ["/gpfs21/scratch/cqs/shengq1/brown/20161003_bamplot_tissues/H3K27Ac_multiple_tissues/Alignment_Post_Processing_14814.bam"],
    "Fibroblast of arm" => ["/gpfs21/scratch/cqs/shengq1/brown/20161003_bamplot_tissues/H3K27Ac_multiple_tissues/Alignment_Post_Processing_6850.bam"],
    "Kidney"            => ["/gpfs21/scratch/cqs/shengq1/brown/20161003_bamplot_tissues/H3K27Ac_multiple_tissues/Alignment_Post_Processing_883.bam"],
    "Lung"              => ["/gpfs21/scratch/cqs/shengq1/brown/20161003_bamplot_tissues/H3K27Ac_multiple_tissues/Alignment_Post_Processing_14834.bam"],
    "Liver"             => ["/gpfs21/scratch/cqs/shengq1/brown/20161003_bamplot_tissues/H3K27Ac_multiple_tissues/Alignment_Post_Processing_15413.bam"],
    "Pancreas"          => ["/gpfs21/scratch/cqs/shengq1/brown/20161003_bamplot_tissues/H3K27Ac_multiple_tissues/Alignment_Post_Processing_9687.bam"],
    "Thyroid gland"     => ["/gpfs21/scratch/cqs/shengq1/brown/20161003_bamplot_tissues/H3K27Ac_multiple_tissues/Alignment_Post_Processing_9685.bam"],
    "Uterus"            => ["/gpfs21/scratch/cqs/shengq1/brown/20161003_bamplot_tissues/H3K27Ac_multiple_tissues/Alignment_Post_Processing_14833.bam"],
  },
  plotgroups => {
    "H3K27ac" => [
      "Adipose nuclei", "Adult CD14", "Angular gyrus", "Aorta",

      #"Aorta 2", "Aorta 3",
      "Fibroblast of arm", "Kidney", "Lung", "Liver", "Pancreas", "Thyroid gland", "Uterus"
    ]
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
    is_single_pdf      => 1,
    is_draw_individual => 0,
    default_color      => "0,0,255",
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
