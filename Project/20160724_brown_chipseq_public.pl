#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;
use CQS::FileUtils;
use Data::Dumper;

my $target_dir = create_directory_or_die("/scratch/cqs/shengq1/brown/20160724_chipseq_public");
my $task       = "public";

my $fasta_file   = "/scratch/cqs/shengq1/references/gencode/hg19/bowtie_index_1.1.2/GRCh37.p13.genome.fa";
my $bowtie_index = "/scratch/cqs/shengq1/references/gencode/hg19/bowtie_index_1.1.2/GRCh37.p13.genome";
my $cqstools     = "/home/shengq1/cqstools/cqstools.exe";

my $transcript_gtf = "/scratch/cqs/shengq1/references/gencode/hg19/gencode.v19.chr_patch_hapl_scaff.annotation.gtf";
my $qc3_perl       = "/scratch/cqs/shengq1/local/bin/qc3/qc3.pl";

my $email = "quanhu.sheng\@vanderbilt.edu";

my $config = {
  general => { task_name => $task },
  files   => {
    "input"                 => ["/gpfs21/scratch/cqs/shengq1/brown/20160724_chipseq_public/data/GSE77267/input.sra"],
    "iPS_H3K27ac"           => ["/gpfs21/scratch/cqs/shengq1/brown/20160724_chipseq_public/data/GSE77267/H3K27ac.sra"],
    "iPS_H3K27me3"          => ["/gpfs21/scratch/cqs/shengq1/brown/20160724_chipseq_public/data/GSE77267/H3K27me3.sra"],
    "iPS_H3K4me1"           => ["/gpfs21/scratch/cqs/shengq1/brown/20160724_chipseq_public/data/GSE77267/H3K4me1.sra"],
    "iPS_H3K4me3"           => ["/gpfs21/scratch/cqs/shengq1/brown/20160724_chipseq_public/data/GSE77267/H3K4me3.sra"],
    "polII"                 => ["/gpfs21/scratch/cqs/shengq1/brown/20160724_chipseq_public/data/GSE77267/polII.sra"],
    "AdultHeart_acCBP_p300" => ["/gpfs21/scratch/cqs/shengq1/brown/20160724_chipseq_public/data/GSE32587/HumanAdultHeart_acCBP-p300.sra"],
    "FetalHeart_acCBP_p300" => ["/gpfs21/scratch/cqs/shengq1/brown/20160724_chipseq_public/data/GSE32587/HumanFetalHeart_acCBP-p300.sra"],
    "FetalHeart_polII"      => ["/gpfs21/scratch/cqs/shengq1/brown/20160724_chipseq_public/data/GSE32587/HumanFetalHeart_polII.sra"],
  },
  treatments => {
    "iPS_H3K27ac"           => ["iPS_H3K27ac"],
    "iPS_H3K27me3"          => ["iPS_H3K27me3"],
    "iPS_H3K4me1"           => ["iPS_H3K4me1"],
    "iPS_H3K4me3"           => ["iPS_H3K4me3"],
    "polII"                 => ["polII"],
    "AdultHeart_acCBP_p300" => ["AdultHeart_acCBP_p300"],
    "FetalHeart_acCBP_p300" => ["FetalHeart_acCBP_p300"],
    "FetalHeart_polII"      => ["FetalHeart_polII"],
  },
  controls => {
    "iPS_H3K27ac"           => ["input"],
    "iPS_H3K27me3"          => ["input"],
    "iPS_H3K4me1"           => ["input"],
    "iPS_H3K4me3"           => ["input"],
    "polII"                 => ["input"],
    "AdultHeart_acCBP_p300" => ["input"],
    "FetalHeart_acCBP_p300" => ["input"],
    "FetalHeart_polII"      => ["input"],
  },
  depthgroups => {
    "iPS_H3K27ac"           => [ "iPS_H3K27ac",           "input" ],
    "iPS_H3K27me3"          => [ "iPS_H3K27me3",          "input" ],
    "iPS_H3K4me1"           => [ "iPS_H3K4me1",           "input" ],
    "iPS_H3K4me3"           => [ "iPS_H3K4me3",           "input" ],
    "polII"                 => [ "polII",                 "input" ],
    "AdultHeart_acCBP_p300" => [ "AdultHeart_acCBP_p300", "input" ],
    "FetalHeart_acCBP_p300" => [ "FetalHeart_acCBP_p300", "input" ],
    "FetalHeart_polII"      => [ "FetalHeart_polII",      "input" ],
  },
  sra2fastq => {
    class      => "SRA::FastqDump",
    perform    => 1,
    ispaired   => 0,
    target_dir => "${target_dir}/sra2fastq",
    option     => "",
    source_ref => "files",
    sh_direct  => 0,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "10",
      "mem"      => "10gb"
    },
  },
  fastqc_raw => {
    class      => "QC::FastQC",
    perform    => 1,
    target_dir => "${target_dir}/fastqc_raw",
    option     => "",
    source_ref => ["sra2fastq"],
    sh_direct  => 0,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "2",
      "mem"      => "10gb"
    },
  },
  fastqc_raw_summary => {
    class      => "QC::FastQCSummary",
    perform    => 1,
    sh_direct  => 1,
    target_dir => "${target_dir}/fastqc_raw",
    cqstools   => $cqstools,
    option     => "",
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "2",
      "mem"      => "10gb"
    },
  },
  bowtie1 => {
    class                   => "Alignment::Bowtie1",
    perform                 => 1,
    target_dir              => "${target_dir}/bowtie1",
    option                  => "-v 1 -m 1 --best --strata",
    fasta_file              => $fasta_file,
    source_ref              => "sra2fastq",
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
  bowtie1_qc3 => {
    class          => "QC::QC3bam",
    perform        => 1,
    target_dir     => "${target_dir}/bowtie1_qc3",
    option         => "",
    transcript_gtf => $transcript_gtf,
    qc3_perl       => $qc3_perl,
    source_ref     =>  ["bowtie1", ".bam\$" ],
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
  macs1callpeak_depth => {
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
  macs2callpeak => {
    class        => "Chipseq::MACS2Callpeak",
    perform      => 1,
    target_dir   => "${target_dir}/macs2callpeak",
    option       => "-g hs --broad -B -p 1e-9",
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
  macs2callpeak_depth => {
    class         => "Visualization::Depth",
    perform       => 1,
    target_dir    => "${target_dir}/macs2callpeak_depth",
    option        => "",
    source_ref    => [ "macs2callpeak", ".bed\$" ],
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
      step_1 => [ "sra2fastq", "fastqc_raw", "bowtie1" ],
      step_2 => [ "fastqc_raw_summary", "macs1callpeak", "macs1callpeak_depth", "macs2callpeak", "macs2callpeak_depth" ],
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

performTask( $config, "bowtie1_qc3" );

1;
