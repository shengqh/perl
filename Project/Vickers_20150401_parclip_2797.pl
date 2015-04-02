#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;
use CQS::PerformSmallRNA;

my $def = {

  #General options
  task_name  => "Parclip2797",
  email      => "quanhu.sheng\@vanderbilt.edu",
  target_dir => "/scratch/cqs/shengq1/vickers/20150401_parclip_2797_rat_mouse_human/",
  max_thread => 8,
  cluster    => "slurm",

  #Default software parameter (don't change it except you really know it)
  fastq_remove_N => 0,
  adapter        => "TGGAATTCTCGGGTGCCAAGG",

  samtools => "/scratch/cqs/shengq1/local/bin/samtools",
  cqstools => "/home/shengq1/cqstools/CQS.Tools.exe",

  #Data
  files => {
    "RPI40" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/2797_demultiplexing/2797-KCV-1_RPI40_Ago2INS1Huh7.fastq.gz"],
    "RPI41" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/2797_demultiplexing/2797-KCV-1_RPI41_Ago3INS1Huh7.fastq.gz"],
    "RPI42" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/2797_demultiplexing/2797-KCV-1_RPI42_Ago2INS1HCEAC.fastq.gz"],
    "RPI43" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/2797_demultiplexing/2797-KCV-1_RPI43_Ago3INS1HCEAC.fastq.gz"],
    "RPI47" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/2797_demultiplexing/2797-KCV-1_RPI47_Ago2MIN6Huh7.fastq.gz"],
    "RPI48" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/2797_demultiplexing/2797-KCV-1_RPI48_Ago3MIN6Huh7.fastq.gz"],
  },
};

my $rn5Config = performSmallRNA_rn5( $def, 1 );

my $parclip_def = {
  binding_db => "/data/cqs/shengq1/reference/targetscan/targetscan_v61_hg19.bed",
  utr3_db    => "/data/cqs/shengq1/reference/utr3/20140612_ucsc_hg19_3UTR.txt",

  #for parclip target
  fasta_file   => "/gpfs21/scratch/cqs/shengq1/references/hg19_16569_M/hg19_16569_M.fa",
  refgene_file => "/gpfs21/scratch/cqs/shengq1/references/hg19_16569_M/hg19_refgene.tsv",
  genome_2bit  => "/data/cqs/guoy1/reference/hg19/hg19_rCRS.2bit",
  mirna_db     => "/data/cqs/shengq1/reference/miRBase20/hsa.mature.dna.db",
};

#
#my $config = {
#  general => { "task_name" => $def->{task_name}, },
#  t2c     => {
#    class             => "CQS::ParclipT2CFinder",
#    perform           => 1,
#    target_dir        => $def->{target_dir} . "/t2c_finder",
#    option            => "-p 0.05 -e 0.013",
#    source_config_ref => [ $hg19Config, "bowtie1_genome_1mm_NTA_smallRNA_count", ".mapped.xml\$" ],
#    cqs_tools         => $def->{cqstools},
#    sh_direct         => 1,
#    pbs               => {
#      "email"    => $def->{email},
#      "nodes"    => "1:ppn=1",
#      "walltime" => "72",
#      "mem"      => "20gb"
#    },
#  },
#  utr3_count => {
#    class                  => "CQS::SmallRNACount",
#    perform                => 1,
#    target_dir             => $def->{target_dir} . "/count_3utr",
#    option                 => "-m 0",
#    source_config_ref      => [ $hg19Config, "bowtie1_genome_1mm_NTA", ".bam\$" ],
#    fastq_files_config_ref => [ $hg19Config, "identical_NTA", "fastq.gz\$" ],
#    seqcount_config_ref    => [ $hg19Config, "identical_NTA", ".dupcount\$" ],
#    cqs_tools              => $def->{cqstools},
#    coordinate_file        => $def->{utr3_db},
#    samtools               => $def->{samtools},
#    sh_direct              => 1,
#    pbs                    => {
#      "email"    => $def->{email},
#      "nodes"    => "1:ppn=1",
#      "walltime" => "72",
#      "mem"      => "20gb"
#    },
#  },
#  mirna_target => {
#    class        => "CQS::ParclipMirnaTarget",
#    perform      => 1,
#    target_dir   => $def->{target_dir} . "/mirna_3utr_target",
#    option       => "",
#    source_ref   => [ "t2c", ".xml\$" ],
#    target_ref   => [ "utr3_count", ".xml\$" ],
#    fasta_file   => $def->{fasta_file},
#    refgene_file => $def->{refgene_file},
#    cqs_tools    => $def->{cqstools},
#    sh_direct    => 1,
#    pbs          => {
#      "email"    => $def->{email},
#      "nodes"    => "1:ppn=1",
#      "walltime" => "72",
#      "mem"      => "20gb"
#    },
#  },
#  PARalyzer => {
#    class             => "ParClip::PARalyzer",
#    perform           => 1,
#    target_dir        => $def->{target_dir} . "/paralyzer",
#    option            => "",
#    source_config_ref => [ $hg19Config, "bowtie1_genome_1mm_notidentical", ".bam\$" ],
#    genome2bit        => $def->{genome_2bit},
#    mirna_db          => $def->{mirna_db},
#    sh_direct         => 1,
#    pbs               => {
#      "email"    => $def->{email},
#      "nodes"    => "1:ppn=1",
#      "walltime" => "72",
#      "mem"      => "20gb"
#    },
#  },
#
#  binding_count => {
#    class                  => "CQS::SmallRNACount",
#    perform                => 1,
#    target_dir             => $def->{target_dir} . "/count_binding",
#    option                 => "-m 0",
#    source_config_ref      => [ $hg19Config, "bowtie1_genome_1mm_NTA", ".bam\$" ],
#    fastq_files_config_ref => [ $hg19Config, "identical_NTA", "fastq.gz\$" ],
#    seqcount_config_ref    => [ $hg19Config, "identical_NTA", ".dupcount\$" ],
#    cqs_tools              => $def->{cqstools},
#    coordinate_file        => $def->{binding_db},
#    samtools               => $def->{samtools},
#    sh_direct              => 1,
#    pbs                    => {
#      "email"    => $def->{email},
#      "nodes"    => "1:ppn=1",
#      "walltime" => "72",
#      "mem"      => "20gb"
#    },
#  },
#  sequencetask => {
#    class      => "CQS::SequenceTask",
#    perform    => 1,
#    target_dir => $def->{target_dir} . "/sequencetask",
#    option     => "",
#    source     => { step3 => [ "t2c", "utr3_count", "mirna_target", "PARalyzer", "binding_count" ], },
#    sh_direct  => 0,
#    pbs        => {
#      "email"    => $def->{email},
#      "nodes"    => "1:ppn=1",
#      "walltime" => "72",
#      "mem"      => "20gb"
#    },
#  },
#};
#
#performConfig($config);

1;
