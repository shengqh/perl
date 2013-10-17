#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;
use CQS::FileUtils;

my $transcript_gtf = "/data/cqs/guoy1/reference/annotation2/mm10/Mus_musculus.GRCm38.68.gtf";
my $smad_gff       = "/data/cqs/shengq1/reference/mm10/smad4/smad4.gff";
my $cqstools       = "/home/shengq1/cqstools/CQS.Tools.exe";
my $bed_file       = "/data/cqs/shengq1/reference/mm10/mm10.bed";

my $email = "quanhu.sheng\@vanderbilt.edu";
my $task  = "smad4";

my $data1 = {
  target_dir => "/scratch/cqs/shengq1/rnaseq/20130828_chenxi_rnaseq_smad4",
  fastqfiles => {
    "2288-RDB-54" => ["/autofs/blue_sequencer/Runs/projects/2288-RDB/2013-08-23/2288-RDB-54_1.fastq.gz"],
    "2288-RDB-55" => ["/autofs/blue_sequencer/Runs/projects/2288-RDB/2013-08-23/2288-RDB-55_1.fastq.gz"],
    "2288-RDB-56" => ["/autofs/blue_sequencer/Runs/projects/2288-RDB/2013-08-23/2288-RDB-56_1.fastq.gz"],
    "2288-RDB-57" => ["/autofs/blue_sequencer/Runs/projects/2288-RDB/2013-08-23/2288-RDB-57_1.fastq.gz"],
    "2288-RDB-58" => ["/autofs/blue_sequencer/Runs/projects/2288-RDB/2013-08-23/2288-RDB-58_1.fastq.gz"],
    "2288-RDB-59" => ["/autofs/blue_sequencer/Runs/projects/2288-RDB/2013-08-23/2288-RDB-59_1.fastq.gz"],
    "2288-RDB-60" => ["/autofs/blue_sequencer/Runs/projects/2288-RDB/2013-08-23/2288-RDB-60_1.fastq.gz"],
    "2288-RDB-61" => ["/autofs/blue_sequencer/Runs/projects/2288-RDB/2013-08-23/2288-RDB-61_1.fastq.gz"],
  },
  groups => {
    "WT" => [ "2288-RDB-54", "2288-RDB-58", "2288-RDB-59", "2288-RDB-61" ],
    "KO" => [ "2288-RDB-55", "2288-RDB-56", "2288-RDB-57", "2288-RDB-60" ],
  },
  pairs => {
    "KO_vs_WT"   => [ "WT",  "KO" ],
    "KO2_vs_WT2" => [ "WT2", "KO2" ]
  },
};

my $data2 = {
  target_dir => "/scratch/cqs/shengq1/rnaseq/20131017_chenxi_rnaseq_smad4",
  fastqfiles => {
    "2288-RDB-81" => ["/autofs/blue_sequencer/Runs/projects/2288-RDB/2013-10-11/2288-RDB-81_1.fastq.gz"],
    "2288-RDB-82" => ["/autofs/blue_sequencer/Runs/projects/2288-RDB/2013-10-11/2288-RDB-82_1.fastq.gz"],
    "2288-RDB-83" => ["/autofs/blue_sequencer/Runs/projects/2288-RDB/2013-10-11/2288-RDB-83_1.fastq.gz"],
    "2288-RDB-84" => ["/autofs/blue_sequencer/Runs/projects/2288-RDB/2013-10-11/2288-RDB-84_1.fastq.gz"],
    "2288-RDB-85" => ["/autofs/blue_sequencer/Runs/projects/2288-RDB/2013-10-11/2288-RDB-85_1.fastq.gz"],
    "2288-RDB-86" => ["/autofs/blue_sequencer/Runs/projects/2288-RDB/2013-10-11/2288-RDB-86_1.fastq.gz"],
  },
  groups => {
    "WT" => [ "2288-RDB-81", "2288-RDB-83", "2288-RDB-85" ],
    "KO" => [ "2288-RDB-82", "2288-RDB-84", "2288-RDB-86" ]
  },
  pairs => { "KO_vs_WT" => [ "WT", "KO" ] },
};

my @defs = (

  #$data1,
  $data2
);

foreach my $def (@defs) {
  my $target_dir = create_directory_or_die( $def->{target_dir} );

  my $config = {
    general    => { task_name => $task },
    fastqfiles => $def->{fastqfiles},
    groups     => $def->{groups},
    pairs      => $def->{pairs},
    fastqc     => {
      class      => "FastQC",
      perform    => 1,
      target_dir => "${target_dir}/fastqc",
      option     => "",
      source_ref => "fastqfiles",
      sh_direct  => 1,
      pbs        => {
        "email"    => $email,
        "nodes"    => "1:ppn=2",
        "walltime" => "2",
        "mem"      => "10gb"
      },
    },
    tophat2 => {
      class         => "Tophat2",
      perform       => 1,
      target_dir    => "${target_dir}/tophat2",
      option        => "--segment-length 25 -r 0 -p 6",
      source_ref    => "fastqfiles",
      bowtie2_index => "/data/cqs/guoy1/reference/mm10/bowtie2_index/mm10",
      sh_direct     => 0,
      pbs           => {
        "email"    => $email,
        "nodes"    => "1:ppn=6",
        "walltime" => "72",
        "mem"      => "30gb"
      },
    },
    sortbam => {
      class         => "Sortbam",
      perform       => 1,
      target_dir    => "${target_dir}/sortname",
      option        => "",
      source_ref    => "tophat2",
      sort_by_query => 1,
      sh_direct     => 1,
      pbs           => {
        "email"    => $email,
        "nodes"    => "1:ppn=8",
        "walltime" => "72",
        "mem"      => "20gb"
      },
    },
    htseqcount => {
      class      => "HTSeqCount",
      perform    => 1,
      target_dir => "${target_dir}/htseqcount",
      option     => "",
      source_ref => "sortbam",
      gff_file   => $transcript_gtf,
      sh_direct  => 1,
      pbs        => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "72",
        "mem"      => "40gb"
      },
    },
    genetable => {
      class         => "CQSDatatable",
      perform       => 1,
      target_dir    => "${target_dir}/genetable",
      option        => "-p ENS --noheader -o ${task}_gene.count",
      source_ref    => "htseqcount",
      name_map_file => "/data/cqs/shengq1/reference/mm10/mm10.gene.map",
      cqs_tools     => $cqstools,
      sh_direct     => 1,
      pbs           => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "10",
        "mem"      => "10gb"
      },
    },
    bedtoolscount => {
      class      => "BedtoolsCount",
      perform    => 0,
      target_dir => "${target_dir}/bedtoolscount",
      option     => "",
      source_ref => "tophat2",
      bed_file   => $bed_file,
      sh_direct  => 1,
      pbs        => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "10",
        "mem"      => "10gb"
      },
    },
    bedtoolstable => {
      class         => "CQSDatatable",
      perform       => 0,
      target_dir    => "${target_dir}/bedtable",
      option        => "-i 1 -v 2 -p ENS --noheader -o ${task}_gene_bedtools.count",
      source_ref    => "bedtoolscount",
      name_map_file => "/data/cqs/shengq1/reference/mm10/mm10.gene.map",
      cqs_tools     => $cqstools,
      sh_direct     => 1,
      pbs           => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "10",
        "mem"      => "10gb"
      },
    },
    deseq2 => {
      class         => "DESeq2",
      perform       => 1,
      target_dir    => "${target_dir}/deseq2",
      option        => "",
      source_ref    => "pairs",
      groups_ref    => "groups",
      countfile_ref => "genetable",
      sh_direct     => 1,
      pbs           => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "10",
        "mem"      => "10gb"
      },
    },
    dexseqcount => {
      class        => "DexseqCount",
      perform      => 1,
      target_dir   => "${target_dir}/dexseqcount",
      option       => "",
      source_ref   => "tophat2",
      gff_file     => $smad_gff,
      dexseq_count => "/home/shengq1/pylibs/bin/dexseq_count.py",
      sh_direct    => 1,
      pbs          => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "72",
        "mem"      => "40gb"
      },
    },
    exontable => {
      class      => "CQSDatatable",
      perform    => 1,
      target_dir => "${target_dir}/exontable",
      option     => "-p ENS --noheader -o ${task}_exon.count",
      source_ref => "dexseqcount",
      cqs_tools  => $cqstools,
      sh_direct  => 1,
      pbs        => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "10",
        "mem"      => "10gb"
      },
    },
  };
  performConfig($config);
}

1;
