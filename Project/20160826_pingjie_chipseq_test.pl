#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;
use Pipeline::ChIPSeq;

my $input = "/scratch/pingj1/chipseq/tagAlign/InputDNAAll.tagAlign.gz";

my $config = {

  #general options
  task_name       => "chipseq_pipeline",
  email           => "quanhu.sheng\@vanderbilt.edu",
  target_dir      => "/scratch/cqs/shengq1/chipseq/20160826_pingjie_chipseq_test/",
  max_thread      => 8,
  min_read_length => 16,
  cluster         => "slurm",

  #softwares
  cqstools   => "/home/shengq1/cqstools/cqstools.exe",
  picard_jar => "/scratch/cqs/shengq1/local/bin/picard/picard.jar",
  spp_r      => "/home/pingj1/soft/phantompeakqualtools/run_spp.R",

  #options
  bwa_index    => "/scratch/cqs/shengq1/references/gencode/hg19/bwa_index_0.7.12/GRCh37.p13.genome.fa",
  macs2_option => "-f BED -g hs -B -q 0.01",

  #data
  files => {

    #repeat 1
    "REP1_IgG"      => [ "/workspace/pingj1/chipseq/rawdata/REPLICATE_1/3364-WPT-14_1.fq.gz", "/workspace/pingj1/chipseq/rawdata/REPLICATE_1/3364-WPT-14_2.fq.gz" ],
    "REP1_N262"     => [ "/workspace/pingj1/chipseq/rawdata/REPLICATE_1/3364-WPT-15_1.fq.gz", "/workspace/pingj1/chipseq/rawdata/REPLICATE_1/3364-WPT-15_2.fq.gz" ],
    "REP1_Y69"      => [ "/workspace/pingj1/chipseq/rawdata/REPLICATE_1/3364-WPT-16_1.fq.gz", "/workspace/pingj1/chipseq/rawdata/REPLICATE_1/3364-WPT-16_2.fq.gz" ],
    "REP1_WDR5"     => [ "/workspace/pingj1/chipseq/rawdata/REPLICATE_1/3364-WPT-17_1.fq.gz", "/workspace/pingj1/chipseq/rawdata/REPLICATE_1/3364-WPT-17_2.fq.gz" ],
    "REP1_H3K4me3"  => [ "/workspace/pingj1/chipseq/rawdata/REPLICATE_1/3364-WPT-18_1.fq.gz", "/workspace/pingj1/chipseq/rawdata/REPLICATE_1/3364-WPT-18_2.fq.gz" ],
    "REP1_InputDNA" => [ "/workspace/pingj1/chipseq/rawdata/REPLICATE_1/3364-WPT-19_1.fq.gz", "/workspace/pingj1/chipseq/rawdata/REPLICATE_1/3364-WPT-19_2.fq.gz" ],

    #repeat 2
    "REP2_IgG"      => ["/workspace/pingj1/chipseq/rawdata/REPLICATE_2/LTIGG_S47_L008_R1_001.fastq.gz"],
    "REP2_N262"     => ["/workspace/pingj1/chipseq/rawdata/REPLICATE_2/LTN262_S48_L008_R1_001.fastq.gz"],
    "REP2_Y69"      => ["/workspace/pingj1/chipseq/rawdata/REPLICATE_2/LTY69_S49_L008_R1_001.fastq.gz"],
    "REP2_WDR5"     => ["/workspace/pingj1/chipseq/rawdata/REPLICATE_2/LTWDR5_S50_L008_R1_001.fastq.gz"],
    "REP2_H3K4me3"  => ["/workspace/pingj1/chipseq/rawdata/REPLICATE_2/LTH3K4_S51_L008_R1_001.fastq.gz"],
    "REP2_InputDNA" => ["/workspace/pingj1/chipseq/rawdata/REPLICATE_2/LTIN_S52_L008_R1_001.fastq.gz"],

    #repeat 3
    "REP3_IgG"      => ["/workspace/pingj1/chipseq/rawdata/REPLICATE_3/3444-WPT-13_1.fq.gz"],
    "REP3_N262"     => ["/workspace/pingj1/chipseq/rawdata/REPLICATE_3/3444-WPT-14_1.fq.gz"],
    "REP3_Y69"      => ["/workspace/pingj1/chipseq/rawdata/REPLICATE_3/3444-WPT-15_1.fq.gz"],
    "REP3_WDR5"     => ["/workspace/pingj1/chipseq/rawdata/REPLICATE_3/3444-WPT-16_1.fq.gz"],
    "REP3_H3K4me3"  => ["/workspace/pingj1/chipseq/rawdata/REPLICATE_3/3444-WPT-17_1.fq.gz"],
    "REP3_InputDNA" => ["/workspace/pingj1/chipseq/rawdata/REPLICATE_3/3444-WPT-18_1.fq.gz"],
  },
  groups => {

    #repeat 1 DNA
    "REP1_N262_DNA"    => ["REP1_N262"],
    "REP1_Y69_DNA"     => ["REP1_Y69"],
    "REP1_WDR5_DNA"    => ["REP1_WDR5"],
    "REP1_H3K4me3_DNA" => ["REP1_H3K4me3"],

    #repeat 2 DNA
    "REP2_N262_DNA"    => ["REP2_N262"],
    "REP2_Y69_DNA"     => ["REP2_Y69"],
    "REP2_WDR5_DNA"    => ["REP2_WDR5"],
    "REP2_H3K4me3_DNA" => ["REP2_H3K4me3"],

    #repeat 3 DNA
    "REP3_N262_DNA"    => ["REP3_N262"],
    "REP3_Y69_DNA"     => ["REP3_Y69"],
    "REP3_WDR5_DNA"    => ["REP3_WDR5"],
    "REP3_H3K4me3_DNA" => ["REP3_H3K4me3"],

    #repeat 1 IgG
    "REP1_N262_IgG"    => ["REP1_N262"],
    "REP1_Y69_IgG"     => ["REP1_Y69"],
    "REP1_WDR5_IgG"    => ["REP1_WDR5"],
    "REP1_H3K4me3_IgG" => ["REP1_H3K4me3"],

    #repeat 2 IgG
    "REP2_N262_IgG"    => ["REP2_N262"],
    "REP2_Y69_IgG"     => ["REP2_Y69"],
    "REP2_WDR5_IgG"    => ["REP2_WDR5"],
    "REP2_H3K4me3_IgG" => ["REP2_H3K4me3"],

    #repeat 3 IgG
    "REP3_N262_IgG"    => ["REP3_N262"],
    "REP3_Y69_IgG"     => ["REP3_Y69"],
    "REP3_WDR5_IgG"    => ["REP3_WDR5"],
    "REP3_H3K4me3_IgG" => ["REP3_H3K4me3"],
  },

  #merge_tagaligns_files is optional, the key will be used in any task require tagalign files
  merge_tagaligns_files => {
    "InputDNA" => [ "REP1_InputDNA", "REP2_InputDNA", "REP3_InputDNA" ],
    "InputIgG" => [ "REP1_IgG",      "REP2_IgG",      "REP3_IgG" ]
  },

  #spp_inputs is optional, inputs will be used if spp_inputs is not defined.
  spp_inputs => {

    #repeat 1 DNA
    "REP1_N262_DNA"    => ["InputDNA"],
    "REP1_Y69_DNA"     => ["InputDNA"],
    "REP1_WDR5_DNA"    => ["InputDNA"],
    "REP1_H3K4me3_DNA" => ["InputDNA"],

    #repeat 2 DNA
    "REP2_N262_DNA"    => ["InputDNA"],
    "REP2_Y69_DNA"     => ["InputDNA"],
    "REP2_WDR5_DNA"    => ["InputDNA"],
    "REP2_H3K4me3_DNA" => ["InputDNA"],

    #repeat 3 DNA
    "REP3_N262_DNA"    => ["InputDNA"],
    "REP3_Y69_DNA"     => ["InputDNA"],
    "REP3_WDR5_DNA"    => ["InputDNA"],
    "REP3_H3K4me3_DNA" => ["InputDNA"],

    #repeat 1 IgG
    "REP1_N262_IgG"    => ["InputIgG"],
    "REP1_Y69_IgG"     => ["InputIgG"],
    "REP1_WDR5_IgG"    => ["InputIgG"],
    "REP1_H3K4me3_IgG" => ["InputIgG"],

    #repeat 2 IgG
    "REP2_N262_IgG"    => ["InputIgG"],
    "REP2_Y69_IgG"     => ["InputIgG"],
    "REP2_WDR5_IgG"    => ["InputIgG"],
    "REP2_H3K4me3_IgG" => ["InputIgG"],

    #repeat 3 IgG
    "REP3_N262_IgG"    => ["InputIgG"],
    "REP3_Y69_IgG"     => ["InputIgG"],
    "REP3_WDR5_IgG"    => ["InputIgG"],
    "REP3_H3K4me3_IgG" => ["InputIgG"],
  },

  inputs => {

    #repeat 1 DNA
    "REP1_N262_DNA"    => ["REP1_InputDNA"],
    "REP1_Y69_DNA"     => ["REP1_InputDNA"],
    "REP1_WDR5_DNA"    => ["REP1_InputDNA"],
    "REP1_H3K4me3_DNA" => ["REP1_InputDNA"],

    #repeat 2 DNA
    "REP2_N262_DNA"    => ["REP2_InputDNA"],
    "REP2_Y69_DNA"     => ["REP2_InputDNA"],
    "REP2_WDR5_DNA"    => ["REP2_InputDNA"],
    "REP2_H3K4me3_DNA" => ["REP2_InputDNA"],

    #repeat 3 DNA
    "REP3_N262_DNA"    => ["REP3_InputDNA"],
    "REP3_Y69_DNA"     => ["REP3_InputDNA"],
    "REP3_WDR5_DNA"    => ["REP3_InputDNA"],
    "REP3_H3K4me3_DNA" => ["REP3_InputDNA"],

    #repeat 1 IgG
    "REP1_N262_IgG"    => ["REP1_IgG"],
    "REP1_Y69_IgG"     => ["REP1_IgG"],
    "REP1_WDR5_IgG"    => ["REP1_IgG"],
    "REP1_H3K4me3_IgG" => ["REP1_IgG"],

    #repeat 2 IgG
    "REP2_N262_IgG"    => ["REP2_IgG"],
    "REP2_Y69_IgG"     => ["REP2_IgG"],
    "REP2_WDR5_IgG"    => ["REP2_IgG"],
    "REP2_H3K4me3_IgG" => ["REP2_IgG"],

    #repeat 3 IgG
    "REP3_N262_IgG"    => ["REP3_IgG"],
    "REP3_Y69_IgG"     => ["REP3_IgG"],
    "REP3_WDR5_IgG"    => ["REP3_IgG"],
    "REP3_H3K4me3_IgG" => ["REP3_IgG"],
  },
  chipseq_qc_table => {
    "InputDNA" => {
      "REP1_N262" => {
        Sample   => "REP1_N262",
        Input    => "REP1_InputDNA",
        SampleId => "N262",
        RepeatId => "1"
      },
      "REP2_N262" => {
        Sample   => "REP2_N262",
        Input    => "REP2_InputDNA",
        SampleId => "N262",
        RepeatId => "2"
      },
      "REP3_N262" => {
        Sample   => "REP3_N262",
        Input    => "REP3_InputDNA",
        SampleId => "N262",
        RepeatId => "3"
      },
      "REP1_Y69" => {
        Sample   => "REP1_Y69",
        Input    => "REP1_InputDNA",
        SampleId => "Y69",
        RepeatId => "1"
      },
      "REP2_Y69" => {
        Sample   => "REP2_Y69",
        Input    => "REP2_InputDNA",
        SampleId => "Y69",
        RepeatId => "2"
      },
      "REP3_Y69" => {
        Sample   => "REP3_Y69",
        Input    => "REP3_InputDNA",
        SampleId => "Y69",
        RepeatId => "3"
      },
      "REP1_WDR5" => {
        Sample   => "REP1_WDR5",
        Input    => "REP1_InputDNA",
        SampleId => "WDR5",
        RepeatId => "1"
      },
      "REP2_WDR5" => {
        Sample   => "REP2_WDR5",
        Input    => "REP2_InputDNA",
        SampleId => "WDR5",
        RepeatId => "2"
      },
      "REP3_WDR5" => {
        Sample   => "REP3_WDR5",
        Input    => "REP3_InputDNA",
        SampleId => "WDR5",
        RepeatId => "3"
      },
      "REP1_H3K4me3" => {
        Sample   => "REP1_H3K4me3",
        Input    => "REP1_InputDNA",
        SampleId => "WDR5",
        RepeatId => "1"
      },
      "REP2_H3K4me3" => {
        Sample   => "REP2_H3K4me3",
        Input    => "REP2_InputDNA",
        SampleId => "WDR5",
        RepeatId => "2"
      },
      "REP3_H3K4me3" => {
        Sample   => "REP3_H3K4me3",
        Input    => "REP3_InputDNA",
        SampleId => "WDR5",
        RepeatId => "3"
      },
    },
    "InputIgG" => {
      "REP1_N262" => {
        Sample   => "REP1_N262",
        Input    => "REP1_IgG",
        SampleId => "N262",
        RepeatId => "1"
      },
      "REP2_N262" => {
        Sample   => "REP2_N262",
        Input    => "REP2_IgG",
        SampleId => "N262",
        RepeatId => "2"
      },
      "REP3_N262" => {
        Sample   => "REP3_N262",
        Input    => "REP3_IgG",
        SampleId => "N262",
        RepeatId => "3"
      },
      "REP1_Y69" => {
        Sample   => "REP1_Y69",
        Input    => "REP1_IgG",
        SampleId => "Y69",
        RepeatId => "1"
      },
      "REP2_Y69" => {
        Sample   => "REP2_Y69",
        Input    => "REP2_IgG",
        SampleId => "Y69",
        RepeatId => "2"
      },
      "REP3_Y69" => {
        Sample   => "REP3_Y69",
        Input    => "REP3_IgG",
        SampleId => "Y69",
        RepeatId => "3"
      },
      "REP1_WDR5" => {
        Sample   => "REP1_WDR5",
        Input    => "REP1_IgG",
        SampleId => "WDR5",
        RepeatId => "1"
      },
      "REP2_WDR5" => {
        Sample   => "REP2_WDR5",
        Input    => "REP2_IgG",
        SampleId => "WDR5",
        RepeatId => "2"
      },
      "REP3_WDR5" => {
        Sample   => "REP3_WDR5",
        Input    => "REP3_IgG",
        SampleId => "WDR5",
        RepeatId => "3"
      },
      "REP1_H3K4me3" => {
        Sample   => "REP1_H3K4me3",
        Input    => "REP1_IgG",
        SampleId => "WDR5",
        RepeatId => "1"
      },
      "REP2_H3K4me3" => {
        Sample   => "REP2_H3K4me3",
        Input    => "REP2_IgG",
        SampleId => "WDR5",
        RepeatId => "2"
      },
      "REP3_H3K4me3" => {
        Sample   => "REP3_H3K4me3",
        Input    => "REP3_IgG",
        SampleId => "WDR5",
        RepeatId => "3"
      },
    },
  }
};

performChIPSeq($config);

#performChIPSeqTask($config, "spp");

1;

