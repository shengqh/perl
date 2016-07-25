#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;
use CQS::FileUtils;
use Data::Dumper;

my $target_dir = create_directory_or_die("/scratch/cqs/shengq1/brown/20160719_chipseq_3512");

my $fasta_file     = "/scratch/cqs/shengq1/references/gencode/hg19/bowtie_index_1.1.2/GRCh37.p13.genome.fa";
my $bowtie_index   = "/scratch/cqs/shengq1/references/gencode/hg19/bowtie_index_1.1.2/GRCh37.p13.genome";
my $cqstools       = "/home/shengq1/cqstools/cqstools.exe";
my $qc3_perl       = "/scratch/cqs/shengq1/local/bin/qc3/qc3.pl";
my $transcript_gtf = "/scratch/cqs/shengq1/references/gencode/hg19/gencode.v19.chr_patch_hapl_scaff.annotation.gtf";
my $picard_jar     = "/scratch/cqs/shengq1/local/bin/picard/picard.jar";
my $bwa_fasta      = "/scratch/cqs/shengq1/references/gencode/hg19/bwa_index_0.7.12/GRCh37.p13.genome.fa";

my $email = "quanhu.sheng\@vanderbilt.edu";
my $task  = "3512";

my $config = {
  general => { task_name => $task },
  files   => {
    "HUVEC_LSS"               => ["/gpfs21/scratch/cqs/shengq1/brown/data/3512/3512-JDB-1_1_sequence.txt.gz"],
    "HUVEC_Oscillatory"       => ["/gpfs21/scratch/cqs/shengq1/brown/data/3512/3512-JDB-2_1_sequence.txt.gz"],
    "HUVEC_LSS_Input"         => ["/gpfs21/scratch/cqs/shengq1/brown/data/3512/3512-JDB-3_1_sequence.txt.gz"],
    "HUVEC_Oscillatory_Input" => ["/gpfs21/scratch/cqs/shengq1/brown/data/3512/3512-JDB-4_1_sequence.txt.gz"],
  },
  treatments => {
    "HUVEC_LSS"         => ["HUVEC_LSS"],
    "HUVEC_Oscillatory" => ["HUVEC_Oscillatory"],
  },
  controls => {
    "HUVEC_LSS"         => ["HUVEC_LSS_Input"],
    "HUVEC_Oscillatory" => ["HUVEC_Oscillatory_Input"],
  },
  depthgroups => {
    "HUVEC_LSS"         => [ "HUVEC_LSS",         "HUVEC_LSS_Input" ],
    "HUVEC_Oscillatory" => [ "HUVEC_Oscillatory", "HUVEC_Oscillatory_Input" ],
  },
  fastqc_pre_trim => {
    class      => "QC::FastQC",
    perform    => 1,
    target_dir => "${target_dir}/fastqc_pre_trim",
    option     => "",
    source_ref => ["files"],
    sh_direct  => 0,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=2",
      "walltime" => "2",
      "mem"      => "10gb"
    },
  },
  fastqc_pre_trim_summary => {
    class      => "QC::FastQCSummary",
    perform    => 1,
    sh_direct  => 1,
    target_dir => "${target_dir}/fastqc_pre_trim",
    cqstools   => $cqstools,
    option     => "",
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "2",
      "mem"      => "10gb"
    },
  },
  cutadapt => {
    class      => "Trimmer::Cutadapt",
    perform    => 1,
    target_dir => "${target_dir}/cutadapt",
    option     => "-m 30 --trim-n",
    source_ref => ["files"],
    adapter    => "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",    #trueseq adapter
    extension  => "_clipped.fastq",
    sh_direct  => 0,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "24",
      "mem"      => "20gb"
    },
  },
  fastqc_post_trim => {
    class      => "QC::FastQC",
    perform    => 1,
    target_dir => "${target_dir}/fastqc_post_trim",
    option     => "",
    sh_direct  => 1,
    source_ref => [ "cutadapt", ".fastq.gz" ],
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "2",
      "mem"      => "10gb"
    },
  },
  fastqc_post_trim_summary => {
    class      => "QC::FastQCSummary",
    perform    => 1,
    sh_direct  => 1,
    target_dir => "${target_dir}/fastqc_post_trim",
    cqstools   => $cqstools,
    option     => "",
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "2",
      "mem"      => "10gb"
    },
  },
  fastq_len => {
    class      => "CQS::FastqLen",
    perform    => 1,
    target_dir => "$target_dir/fastq_len",
    option     => "",
    source_ref => "cutadapt",
    cqstools   => $cqstools,
    sh_direct  => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "24",
      "mem"      => "20gb"
    },
  },
  bwa => {
    class      => "Alignment::BWA",
    perform    => 1,
    target_dir => "${target_dir}/bwa",
    option     => "",
    bwa_index  => $bwa_fasta,
    source_ref => [ "cutadapt", ".fastq.gz\$" ],
    picard_jar => $picard_jar,
    sh_direct  => 0,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  bowtie1 => {
    class                   => "Alignment::Bowtie1",
    perform                 => 1,
    target_dir              => "${target_dir}/bowtie1",
    option                  => "-v 1 -m 1 --best --strata",
    fasta_file              => $fasta_file,
    source_ref              => [ "cutadapt", ".fastq.gz\$" ],
    bowtie1_index           => $bowtie_index,
    chromosome_grep_pattern => "\"^chr\"",
    sh_direct               => 0,
    pbs                     => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  bowtie1_miss2 => {
    class                   => "Alignment::Bowtie1",
    perform                 => 1,
    target_dir              => "${target_dir}/bowtie1_miss2",
    option                  => "-v 2 -m 1 --best --strata",
    fasta_file              => $fasta_file,
    source_ref              => [ "cutadapt", ".fastq.gz\$" ],
    bowtie1_index           => $bowtie_index,
    chromosome_grep_pattern => "\"^chr\"",
    sh_direct               => 0,
    pbs                     => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  qc3bam => {
    class          => "QC::QC3bam",
    perform        => 1,
    target_dir     => "${target_dir}/qc3bam",
    option         => "",
    transcript_gtf => $transcript_gtf,
    qc3_perl       => $qc3_perl,
    source_ref     => "bowtie1",
    pbs            => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  macs1callpeak => {
    class        => "Chipseq::MACS",
    perform      => 1,
    target_dir   => "${target_dir}/macs1callpeak",
    option       => "-p 1e-9 -w -S --space=50",
    source_ref   => "bowtie1",
    groups_ref   => "treatments",
    controls_ref => "controls",
    sh_direct    => 0,
    pbs          => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  macs1callpeak_bradner_rose2 => {
    class                => "Chipseq::BradnerRose2",
    perform              => 1,
    target_dir           => "${target_dir}/macs1callpeak_bradner_rose2",
    option               => "",
    source_ref           => "bowtie1",
    groups_ref           => "treatments",
    controls_ref         => "controls",
    pipeline_dir         => "/scratch/cqs/shengq1/local/bin/bradnerlab",
    binding_site_bed_ref => [ "macs1callpeak", ".bed\$" ],
    genome               => "hg19",
    sh_direct            => 1,
    pbs                  => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  depth => {
    class         => "Visualization::Depth",
    perform       => 1,
    target_dir    => "${target_dir}/macs1callpeak_depth",
    option        => "",
    source_ref    => [ "macs1callpeak", ".name.bed" ],
    groups_ref    => "depthgroups",
    bam_files_ref => "bowtie1",
    sh_direct     => 1,
    pbs           => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "24",
      "mem"      => "10gb"
    },
  },
  sequencetask => {
    class      => "CQS::SequenceTask",
    perform    => 1,
    target_dir => "${target_dir}/sequencetask",
    option     => "",
    source     => {
      step_1 => [ "fastqc_pre_trim",         "cutadapt",                 "fastqc_post_trim", "fastq_len",     "bowtie1", "bwa" ],
      step_2 => [ "fastqc_pre_trim_summary", "fastqc_post_trim_summary", "qc3bam",           "macs1callpeak", "macs1callpeak_bradner_rose2", "depth" ],
    },
    sh_direct => 0,
    pbs       => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
};

#performConfig($config);
performTask( $config, "depth" );

1;
