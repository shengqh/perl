#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;
use CQS::FileUtils;
use Data::Dumper;
use Pipeline::SmallRNAUtils;
use Pipeline::ChipSeqBowtie;

my $def = {
  task_name      => "20161216_chipseq_gse",
  fasta_file     => "/scratch/cqs/shengq1/references/gencode/hg19/bowtie_index_1.1.2/GRCh37.p13.genome.fa",
  bowtie_index   => "/scratch/cqs/shengq1/references/gencode/hg19/bowtie_index_1.1.2/GRCh37.p13.genome",
  cqstools       => "/home/shengq1/cqstools/cqstools.exe",
  plot_gff       => "/scratch/cqs/shengq1/chipseq/20160823_janathan_chipseq_195R3_gse53999_bamplot/config/H3K27ac.gff",
  bamplot_option => "-g HG19 -y uniform -r",
  sra2fastq      => 1,
  email          => "quanhu.sheng\@vanderbilt.edu",
  target_dir     => create_directory_or_die("/scratch/cqs/shengq1/brown/20161216_chipseq_gse"),

  files => {
    "CM1_8_H3K27ac" => ["/gpfs21/scratch/cqs/shengq1/brown/data/20161216_chipseq_gse/GSM2279957.sra"],
    "CM1_8_Input"   => ["/gpfs21/scratch/cqs/shengq1/brown/data/20161216_chipseq_gse/GSM2279961.sra"],
    "CM1_8_Tbx5"    => ["/gpfs21/scratch/cqs/shengq1/brown/data/20161216_chipseq_gse/GSM2279963.sra"],
    "Gata4"         => ["/gpfs21/scratch/cqs/shengq1/brown/data/20161216_chipseq_gse/GSM2279964.sra"],
    "H3K27ac"       => ["/gpfs21/scratch/cqs/shengq1/brown/data/20161216_chipseq_gse/GSM2279965.sra"],
    "H3K27me3"      => ["/gpfs21/scratch/cqs/shengq1/brown/data/20161216_chipseq_gse/GSM2279966.sra"],
    "H3K36me3"      => ["/gpfs21/scratch/cqs/shengq1/brown/data/20161216_chipseq_gse/GSM2279967.sra"],
    "input"         => ["/gpfs21/scratch/cqs/shengq1/brown/data/20161216_chipseq_gse/GSM2279969.sra"],
    "Med1"          => ["/gpfs21/scratch/cqs/shengq1/brown/data/20161216_chipseq_gse/GSM2279970.sra"],
    "Tbx5"          => ["/gpfs21/scratch/cqs/shengq1/brown/data/20161216_chipseq_gse/GSM2279971.sra"],
  },
  treatments => {
    "CM1_8_H3K27ac" => ["CM1_8_H3K27ac"],
    "CM1_8_Tbx5"    => ["CM1_8_Tbx5"],
    "Gata4"         => ["Gata4"],
    "H3K27ac"       => ["H3K27ac"],
    "H3K27me3"      => ["H3K27me3"],
    "H3K36me3"      => ["H3K36me3"],
    "Med1"          => ["Med1"],
    "Tbx5"          => ["Tbx5"],
  },
  controls => {
    "CM1_8_H3K27ac" => ["CM1_8_Input"],
    "CM1_8_Tbx5"    => ["CM1_8_Input"],
    "Gata4"         => ["input"],
    "H3K27ac"       => ["input"],
    "H3K27me3"      => ["input"],
    "H3K36me3"      => ["input"],
    "Med1"          => ["input"],
    "Tbx5"          => ["input"],
  },
  #plotgroups => {
  #  "chipseq" => [ "CM1_8_Input", "CM1_8_H3K27ac", "CM1_8_Tbx5", "input", "Gata4", "H3K27ac", "H3K27me3", "H3K36me3", "Med1", "Tbx5" ],
  #}
};

performChipSeqBowtie($def);
#performChipSeqBowtieTask($def, "bamplot");

1;
