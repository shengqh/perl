#!/usr/bin/perl
use strict;
use warnings;

use CQS::PerformSmallRNA;
use CQS::ClassFactory;

my $def_human = {

  #General options
  task_name  => "3018-KCV-10",
  email      => "quanhu.sheng\@vanderbilt.edu",
  target_dir => "/scratch/cqs/shengq1/vickers/20150224_smallRNA_3018-KCV-10_human",
  max_thread => 8,

  #Default software parameter (don't change it except you really know it)
  fastq_remove_N       => 1,
  smallrnacount_option => "-s --unmapped_fastq",    #export unmapped fastq for following up analysis

  #Data
  files => {
    "3018-KCV-10-37" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-10_Human_Cow/3018-KCV-10-37_CGGAAT_L001_R1_001.fastq.gz"],
    "3018-KCV-10-38" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-10_Human_Cow/3018-KCV-10-38_CTAGCT_L001_R1_001.fastq.gz"],
    "3018-KCV-10-39" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-10_Human_Cow/3018-KCV-10-39_CTATAC_L001_R1_001.fastq.gz"],
    "3018-KCV-10-40" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-10_Human_Cow/3018-KCV-10-40_CTCAGA_L001_R1_001.fastq.gz"],
    "3018-KCV-10-41" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-10_Human_Cow/3018-KCV-10-41_GACGAC_L001_R1_001.fastq.gz"],
    "3018-KCV-10-42" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-10_Human_Cow/3018-KCV-10-42_TAATCG_L001_R1_001.fastq.gz"],
    "3018-KCV-10-43" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-10_Human_Cow/3018-KCV-10-43_TACAGC_L001_R1_001.fastq.gz"],
    "3018-KCV-10-44" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-10_Human_Cow/3018-KCV-10-44_TATAAT_L001_R1_001.fastq.gz"],
    "3018-KCV-10-45" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-10_Human_Cow/3018-KCV-10-45_TCATTC_L001_R1_001.fastq.gz"],
    "3018-KCV-10-46" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-10_Human_Cow/3018-KCV-10-46_TCCCGA_L001_R1_001.fastq.gz"],
    "3018-KCV-10-47" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-10_Human_Cow/3018-KCV-10-47_TCGAAG_L001_R1_001.fastq.gz"],
    "3018-KCV-10-48" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-10_Human_Cow/3018-KCV-10-48_TCGGCA_L001_R1_001.fastq.gz"],
  }
};

my $config_human = performSmallRNA_hg19( $def_human, 0 );

#print Dumper($config_human);

my $def_cow = {

  #General options
  task_name  => "3018-KCV-10",
  email      => "quanhu.sheng\@vanderbilt.edu",
  target_dir => "/scratch/cqs/shengq1/vickers/20150227_smallRNA_3018-KCV-10_cow/",
  max_thread => 8,
  cluster    => "slurm",

  #Default software parameter (don't change it except you really know it)
  bowtie1_option_pm => "-a -m 100 --best --strata -v 0 -p 8",

  #Software and miRBase database options
  samtools => "/scratch/cqs/shengq1/local/bin/samtools",
  cqstools => "/home/shengq1/cqstools/CQS.Tools.exe",

  coordinate    => "/scratch/cqs/shengq1/references/cow/Bos_taurus.UMD3.1.78.gtf",
  bowtie1_index => "/scratch/cqs/shengq1/references/cow/bowtie_index_1.1.1/Bos_taurus_UMD_3.1",
};

my $config_cow = {
  general => { task_name => $def_cow->{task_name}, },

  bowtie1_genome_pm => {
    class             => "Alignment::Bowtie1",
    perform           => 1,
    target_dir        => $def_cow->{target_dir} . "/bowtie1_genome_pm",
    option            => $def_cow->{bowtie1_option_pm},
    source_config_ref => [ $config_human, "bowtie1_genome_1mm_NTA_smallRNA_count", "unmapped.fastq.gz\$" ],
    bowtie1_index     => $def_cow->{bowtie1_index},
    samonly           => 0,
    sh_direct         => 1,
    cluster           => $def_cow->{cluster},
    pbs               => {
      "email"    => $def_cow->{email},
      "nodes"    => "1:ppn=" . $def_cow->{max_thread},
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  bowtie1_genome_pm_count => {
    class               => "CQS::CQSMappedCount",
    perform             => 1,
    target_dir          => $def_cow->{target_dir} . "/bowtie1_genome_pm_count",
    option              => "-s",
    source_ref          => "bowtie1_genome_pm",
    seqcount_config_ref => [ $config_human, "identical", ".dupcount\$" ],
    cqs_tools           => $def_cow->{cqstools},
    gff_file            => $def_cow->{coordinate},
    samtools            => $def_cow->{samtools},
    sh_direct           => 1,
    cluster             => $def_cow->{cluster},
    pbs                 => {
      "email"    => $def_cow->{email},
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  bowtie1_genome_pm_table => {
    class      => "CQS::CQSMappedTable",
    perform    => 1,
    target_dir => $def_cow->{target_dir} . "/bowtie1_genome_pm_table",
    option     => "",
    source_ref => [ "bowtie1_genome_pm_count", ".mapped.xml" ],
    cqs_tools  => $def_cow->{cqstools},
    prefix     => "cow_pm_",
    sh_direct  => 1,
    cluster    => $def_cow->{cluster},
    pbs        => {
      "email"    => $def_cow->{email},
      "nodes"    => "1:ppn=1",
      "walltime" => "10",
      "mem"      => "10gb"
    },
  },
  sequencetask => {
    class      => "CQS::SequenceTask",
    perform    => 0,
    target_dir => $def_cow->{target_dir} . "/sequencetask",
    option     => "",
    source     => {
      one => [ "bowtie1_genome_pm", "bowtie1_genome_pm_count" ],
      all => [

        #NTA table
        "bowtie1_genome_pm_table",
      ],
    },
    sh_direct => 0,
    cluster   => $def_cow->{cluster},
    pbs       => {
      "email"    => $def_cow->{email},
      "nodes"    => "1:ppn=" . $def_cow->{max_thread},
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
};

performConfig($config_cow);

1;

