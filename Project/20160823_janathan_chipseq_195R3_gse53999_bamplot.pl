#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;
use CQS::FileUtils;
use Data::Dumper;

my $target_dir = create_directory_or_die("/scratch/cqs/shengq1/chipseq/20160823_janathan_chipseq_195R3_gse53999_bamplot");

my $cqstools = "/home/shengq1/cqstools/cqstools.exe";

my $email    = "quanhu.sheng\@vanderbilt.edu";
my $task     = "chipseq";
my $plot_gff = "/scratch/cqs/shengq1/chipseq/20160823_janathan_chipseq_195R3_gse53999_bamplot/config/H3K27ac.gff";

my $config = {
  general   => { task_name => $task },
  bam_files => {
    "EC_H3K27AC"                  => ["/gpfs21/scratch/cqs/shengq1/chipseq/20160419_janathan_chipseq_gse53999_hg19/bowtie1/result/EC_H3K27AC_CON/EC_H3K27AC_CON.bam"],
    "EC_WCE"                      => ["/gpfs21/scratch/cqs/shengq1/chipseq/20160419_janathan_chipseq_gse53999_hg19/bowtie1/result/EC_WCE_CON/EC_WCE_CON.bam"],
    "d17_static_ESctrl1_H3K27ac"  => ["/gpfs21/scratch/cqs/shengq1/chipseq/20160302_janathan_chipseq_195R3/bowtie1/result/d17_static_ESctrl1_H3K27ac/d17_static_ESctrl1_H3K27ac.bam"],
    "d17_static_ESctrl2_H3K27ac"  => ["/gpfs21/scratch/cqs/shengq1/chipseq/20160302_janathan_chipseq_195R3/bowtie1/result/d17_static_ESctrl2_H3K27ac/d17_static_ESctrl2_H3K27ac.bam"],
    "d17_static_iPSctrl2_H3K27ac" => ["/scratch/cqs/shengq1/chipseq/20160302_janathan_chipseq_195R3/bowtie1/result/d17_static_iPSctrl2_H3K27ac/d17_static_iPSctrl2_H3K27ac.bam"],
  },
  plotgroups => {
    "H3K27ac" => [ "EC_H3K27AC", "EC_WCE", "d17_static_ESctrl1_H3K27ac", "d17_static_ESctrl2_H3K27ac", "d17_static_iPSctrl2_H3K27ac" ]
  },
  bamplot => {
    class         => "Visualization::Bamplot",
    perform       => 1,
    target_dir    => "${target_dir}/bamplot",
    option        => "-g HG19 -y uniform -r --save-temp",
    source_ref    => "bam_files",
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
};

performConfig($config);
1;
