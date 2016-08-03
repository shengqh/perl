#!/usr/bin/perl
use strict;
use warnings;

use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;
use Hash::Merge qw( merge );

my $target_dir = create_directory_or_die("/scratch/cqs/shengq1/brown/20160803_atacseq_3307_12_13_human");
my $email      = "quanhu.sheng\@vanderbilt.edu";

my $cqstools   = "/home/shengq1/cqstools/cqstools.exe";
my $picard_jar = "/scratch/cqs/shengq1/local/bin/picard/picard.jar";
my $gatk_jar   = "/home/shengq1/local/bin/GATK/GenomeAnalysisTK.jar";

my $macs2call_option = "-f BED -g hg -B -q 0.01 --nomodel";

my $bwa_fasta     = "/scratch/cqs/shengq1/references/gencode/hg19/bwa_index_0.7.12/GRCh37.p13.genome.fa";

my $config = {
  general => { task_name => "HAEC" },
  files   => {
    "HAEC_5percent"  => [ "/scratch/cqs/shengq1/brown/data/3307/3307-JDB-12_1_sequence.txt.gz", "/scratch/cqs/shengq1/brown/data/3307/3307-JDB-12_2_sequence.txt.gz" ],
    "HAEC_10percent" => [ "/scratch/cqs/shengq1/brown/data/3307/3307-JDB-13_1_sequence.txt.gz", "/scratch/cqs/shengq1/brown/data/3307/3307-JDB-13_2_sequence.txt.gz" ],
  },
  comparison => {
    "DIFF_HAEC" => ["HAEC_5percent", "HAEC_10percent"]
  },
  fastqc => {
    class      => "QC::FastQC",
    perform    => 1,
    target_dir => "${target_dir}/fastqc",
    option     => "",
    source_ref => "files",
    sh_direct  => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=2",
      "walltime" => "2",
      "mem"      => "40gb"
    },
  },
  fastqc_summary => {
    class      => "QC::FastQCSummary",
    perform    => 1,
    target_dir => "${target_dir}/fastqc",
    option     => "",
    cqstools   => $cqstools,
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
    option     => "-O 10 -m 50",
    source_ref => "files",
    adapter    => "CTGTCTCTTATA",
    extension  => "_clipped.fastq.gz",
    sh_direct  => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "24",
      "mem"      => "20gb"
    },
  },
  fastqlen => {
    class      => "CQS::FastqLen",
    perform    => 1,
    target_dir => "${target_dir}/fastqlen",
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
    picard_jar => $picard_jar,
    source_ref => ["cutadapt", ".fastq.gz"],
    sh_direct  => 0,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  bwa_cleanbam => {
    class             => "ATACseq::CleanBam",
    perform           => 1,
    target_dir        => "${target_dir}/bwa_cleanbam",
    option            => "",
    source_ref        => "bwa",
    picard_jar        => $picard_jar,
    remove_chromosome => "M",
    sh_direct         => 0,
    pbs               => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "240",
      "mem"      => "40gb"
    },
  },
  bwa_bam2bed => {
    class          => "ATACseq::BamToBed",
    perform        => 1,
    target_dir     => "${target_dir}/bwa_bam2bed",
    option         => "",
    source_ref     => "bwa_cleanbam",
    blacklist_file => "/scratch/cqs/shengq1/references/mappable_region/hg19/wgEncodeDacMapabilityConsensusExcludable.bed",
    sh_direct      => 0,
    pbs            => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  bwa_macs2callpeak => {
    class      => "Chipseq::MACS2Callpeak",
    perform    => 1,
    target_dir => "${target_dir}/bwa_macs2callpeak",
    option     => $macs2call_option,
    source_ref => "bwa_bam2bed",
    sh_direct  => 0,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  bwa_macs2bdgdiff => {
    class      => "Chipseq::MACS2Bdgdiff",
    perform    => 1,
    target_dir => "${target_dir}/bwa_macs2bdgdiff",
    option     => "",
    source_ref => "bwa_macs2callpeak",
    groups_ref => "comparison",
    sh_direct  => 0,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  bwa_macs2bdgdiff_depth => {
    class         => "Visualization::Depth",
    perform       => 1,
    target_dir    => "${target_dir}/bwa_macs2bdgdiff_depth",
    option        => "",
    source_ref    => "bwa_macs2bdgdiff",
    bam_files_ref => "bwa_cleanbam",
    groups_ref    => "comparison",
    sh_direct     => 0,
    pbs           => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },  
  sequencetask => {
    class      => "CQS::SequenceTask",
    perform    => 1,
    target_dir => "${target_dir}/sequencetask",
    option     => "",
    source     => {
      T1 => [ "fastqc", "cutadapt", "fastqlen", "bwa",      "bwa_cleanbam",      "bwa_bam2bed",      "bwa_macs2callpeak" ],
      T2 => [ "fastqc_summary", "bwa_macs2bdgdiff", "bwa_macs2bdgdiff_depth" ],
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

performConfig($config);

#performTask( $config, "bowtie2_pretrim" );

1;

