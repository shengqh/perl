#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;
use CQS::FileUtils;
use Data::Dumper;

my $target_dir = create_directory_or_die("/scratch/cqs/shengq1/brown/20161003_bamplot_aorta");

my $cqstools = "/home/shengq1/cqstools/cqstools.exe";

my $email    = "quanhu.sheng\@vanderbilt.edu";
my $task     = "chipseq";
my $plot_gff = "/scratch/cqs/shengq1/brown/20161003_bamplot_aorta/document/locus.gff";

my $config = {
  general     => { task_name => $task },
  bam_files => {
    "H3K4me3"  => ["/gpfs21/scratch/cqs/shengq1/brown/20161003_bamplot_aorta/All_Histone_Marks_Aorta/Alignment_Post_Processing_9091.bam"],
    "H3K9me3"  => ["/gpfs21/scratch/cqs/shengq1/brown/20161003_bamplot_aorta/All_Histone_Marks_Aorta/Alignment_Post_Processing_9095.bam"],
    "H3K36me3" => ["/gpfs21/scratch/cqs/shengq1/brown/20161003_bamplot_aorta/All_Histone_Marks_Aorta/Alignment_Post_Processing_9096.bam"],
    "H3K27ac"  => ["/gpfs21/scratch/cqs/shengq1/brown/20161003_bamplot_aorta/All_Histone_Marks_Aorta/Alignment_Post_Processing_9097.bam"],
    "Unknown"  => ["/gpfs21/scratch/cqs/shengq1/brown/20161003_bamplot_aorta/All_Histone_Marks_Aorta/Alignment_Post_Processing_9492.bam"],
    "H3K27me3" => ["/gpfs21/scratch/cqs/shengq1/brown/20161003_bamplot_aorta/All_Histone_Marks_Aorta/Alignment_Post_Processing_9499.bam"],
  },
  plotgroups => {
    "Aorta"   => [ "H3K4me3",        "H3K9me3",    "H3K36me3",      "H3K27ac", "H3K27me3" ],
  },
  colormaps => {
    "H3K4me3"  => "0,0,0",
    "H3K9me3"  => "0,0,0",
    "H3K36me3" => "0,0,0",
    "H3K27ac"  => "0,0,255",
    "Unknown"  => "0,0,0",
    "H3K27me3" => "0,0,0",
  },
  bamplot => {
    class      => "Visualization::Bamplot",
    perform    => 1,
    target_dir => "${target_dir}/bamplot",
    #option     => "-g HG19 -y uniform -r",

    option        => "-g HG19 -y uniform -r --save-temp",
    source_ref    => "bam_files",
    groups_ref    => "plotgroups",
    gff_file      => $plot_gff,
    is_rainbow_color => 0,
    is_single_pdf => 0,
    is_draw_individual => 0,
    colors_ref => "colormaps",
    sh_direct     => 1,
    pbs           => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "1",
      "mem"      => "10gb"
    },
  },
};

performConfig($config);
1;
