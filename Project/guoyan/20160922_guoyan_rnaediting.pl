#!/usr/bin/perl
use strict;
use warnings;

use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;
use Hash::Merge qw( merge );

my $target_dir = create_directory_or_die("/workspace/shengq1/guoyan/20160922_rnaediting");
my $email      = "quanhu.sheng\@vanderbilt.edu";

my $cqstools     = "/home/shengq1/cqstools/cqstools.exe";
my $gatk_jar     = "/home/shengq1/local/bin/GATK/GenomeAnalysisTK.jar";
my $picard_jar   = "/scratch/cqs/shengq1/local/bin/picard/picard.jar";
my $bwa_fasta    = "/workspace/shengq1/guoyan/20160922_rnaediting/database/bwa_index_0.7.12/rat_GRM4.fasta";
my $bowtie_fasta = "/workspace/shengq1/guoyan/20160922_rnaediting/database/bowtie_index_1.1.2/rat_GRM4.fasta";
my $bowtie_index = "/workspace/shengq1/guoyan/20160922_rnaediting/database/bowtie_index_1.1.2/rat_GRM4";
my $config       = {
  general => { task_name => "rnaediting" },
  files   => {
    "Cerebellum-Rat1_S01" => [
      "/workspace/guoy1/emeson_lab/RNA-editing/Hussain/2016Apr21/FASTQ/3442-RBE-1-29526771/Cerebellum Rat1-34710992/Data/Intensities/BaseCalls/Cerebellum-Rat1_S1_L001_R1_001.fastq.gz",
      "/workspace/guoy1/emeson_lab/RNA-editing/Hussain/2016Apr21/FASTQ/3442-RBE-1-29526771/Cerebellum Rat1-34710992/Data/Intensities/BaseCalls/Cerebellum-Rat1_S1_L001_R2_001.fastq.gz"
    ],
    "Cerebellum-Rat2_S11" => [
      "/workspace/guoy1/emeson_lab/RNA-editing/Hussain/2016Apr21/FASTQ/3442-RBE-1-29526771/Cerebellum Rat2-34711002/Data/Intensities/BaseCalls/Cerebellum-Rat2_S11_L001_R1_001.fastq.gz",
      "/workspace/guoy1/emeson_lab/RNA-editing/Hussain/2016Apr21/FASTQ/3442-RBE-1-29526771/Cerebellum Rat2-34711002/Data/Intensities/BaseCalls/Cerebellum-Rat2_S11_L001_R2_001.fastq.gz"
    ],
    "Cerebellum-Rat3_S21" => [
      "/workspace/guoy1/emeson_lab/RNA-editing/Hussain/2016Apr21/FASTQ/3442-RBE-1-29526771/Cerebellum Rat3-34711012/Data/Intensities/BaseCalls/Cerebellum-Rat3_S21_L001_R1_001.fastq.gz",
      "/workspace/guoy1/emeson_lab/RNA-editing/Hussain/2016Apr21/FASTQ/3442-RBE-1-29526771/Cerebellum Rat3-34711012/Data/Intensities/BaseCalls/Cerebellum-Rat3_S21_L001_R2_001.fastq.gz"
    ],
    "Colon-Rat1_S07" => [
      "/workspace/guoy1/emeson_lab/RNA-editing/Hussain/2016Apr21/FASTQ/3442-RBE-1-29526771/Colon Rat1-34710998/Data/Intensities/BaseCalls/Colon-Rat1_S7_L001_R1_001.fastq.gz",
      "/workspace/guoy1/emeson_lab/RNA-editing/Hussain/2016Apr21/FASTQ/3442-RBE-1-29526771/Colon Rat1-34710998/Data/Intensities/BaseCalls/Colon-Rat1_S7_L001_R2_001.fastq.gz"
    ],
    "Colon-Rat2_S17" => [
      "/workspace/guoy1/emeson_lab/RNA-editing/Hussain/2016Apr21/FASTQ/3442-RBE-1-29526771/Colon Rat2-34711008/Data/Intensities/BaseCalls/Colon-Rat2_S17_L001_R1_001.fastq.gz",
      "/workspace/guoy1/emeson_lab/RNA-editing/Hussain/2016Apr21/FASTQ/3442-RBE-1-29526771/Colon Rat2-34711008/Data/Intensities/BaseCalls/Colon-Rat2_S17_L001_R2_001.fastq.gz"
    ],
    "Colon-Rat3_S27" => [
      "/workspace/guoy1/emeson_lab/RNA-editing/Hussain/2016Apr21/FASTQ/3442-RBE-1-29526771/Colon Rat3-34711018/Data/Intensities/BaseCalls/Colon-Rat3_S27_L001_R1_001.fastq.gz",
      "/workspace/guoy1/emeson_lab/RNA-editing/Hussain/2016Apr21/FASTQ/3442-RBE-1-29526771/Colon Rat3-34711018/Data/Intensities/BaseCalls/Colon-Rat3_S27_L001_R2_001.fastq.gz"
    ],
    "Cortex-Rat1_S02" => [
      "/workspace/guoy1/emeson_lab/RNA-editing/Hussain/2016Apr21/FASTQ/3442-RBE-1-29526771/Cortex Rat1-34710993/Data/Intensities/BaseCalls/Cortex-Rat1_S2_L001_R1_001.fastq.gz",
      "/workspace/guoy1/emeson_lab/RNA-editing/Hussain/2016Apr21/FASTQ/3442-RBE-1-29526771/Cortex Rat1-34710993/Data/Intensities/BaseCalls/Cortex-Rat1_S2_L001_R2_001.fastq.gz"
    ],
    "Cortex-Rat2_S12" => [
      "/workspace/guoy1/emeson_lab/RNA-editing/Hussain/2016Apr21/FASTQ/3442-RBE-1-29526771/Cortex Rat2-34711003/Data/Intensities/BaseCalls/Cortex-Rat2_S12_L001_R1_001.fastq.gz",
      "/workspace/guoy1/emeson_lab/RNA-editing/Hussain/2016Apr21/FASTQ/3442-RBE-1-29526771/Cortex Rat2-34711003/Data/Intensities/BaseCalls/Cortex-Rat2_S12_L001_R2_001.fastq.gz"
    ],
    "Cortex-Rat3_S22" => [
      "/workspace/guoy1/emeson_lab/RNA-editing/Hussain/2016Apr21/FASTQ/3442-RBE-1-29526771/Cortex Rat3-34711013/Data/Intensities/BaseCalls/Cortex-Rat3_S22_L001_R1_001.fastq.gz",
      "/workspace/guoy1/emeson_lab/RNA-editing/Hussain/2016Apr21/FASTQ/3442-RBE-1-29526771/Cortex Rat3-34711013/Data/Intensities/BaseCalls/Cortex-Rat3_S22_L001_R2_001.fastq.gz"
    ],
    "Hippocampus-Rat1_S03" => [
      "/workspace/guoy1/emeson_lab/RNA-editing/Hussain/2016Apr21/FASTQ/3442-RBE-1-29526771/Hippocampus Rat1-34710994/Data/Intensities/BaseCalls/Hippocampus-Rat1_S3_L001_R1_001.fastq.gz",
      "/workspace/guoy1/emeson_lab/RNA-editing/Hussain/2016Apr21/FASTQ/3442-RBE-1-29526771/Hippocampus Rat1-34710994/Data/Intensities/BaseCalls/Hippocampus-Rat1_S3_L001_R2_001.fastq.gz"
    ],
    "Hippocampus-Rat2_S13" => [
      "/workspace/guoy1/emeson_lab/RNA-editing/Hussain/2016Apr21/FASTQ/3442-RBE-1-29526771/Hippocampus Rat2-34711004/Data/Intensities/BaseCalls/Hippocampus-Rat2_S13_L001_R1_001.fastq.gz",
      "/workspace/guoy1/emeson_lab/RNA-editing/Hussain/2016Apr21/FASTQ/3442-RBE-1-29526771/Hippocampus Rat2-34711004/Data/Intensities/BaseCalls/Hippocampus-Rat2_S13_L001_R2_001.fastq.gz"
    ],
    "Hippocampus-Rat3_S23" => [
      "/workspace/guoy1/emeson_lab/RNA-editing/Hussain/2016Apr21/FASTQ/3442-RBE-1-29526771/Hippocampus Rat3-34711014/Data/Intensities/BaseCalls/Hippocampus-Rat3_S23_L001_R1_001.fastq.gz",
      "/workspace/guoy1/emeson_lab/RNA-editing/Hussain/2016Apr21/FASTQ/3442-RBE-1-29526771/Hippocampus Rat3-34711014/Data/Intensities/BaseCalls/Hippocampus-Rat3_S23_L001_R2_001.fastq.gz"
    ],
    "Hypothalamus-Rat1_S04" => [
      "/workspace/guoy1/emeson_lab/RNA-editing/Hussain/2016Apr21/FASTQ/3442-RBE-1-29526771/Hypothalamus Rat1-34710995/Data/Intensities/BaseCalls/Hypothalamus-Rat1_S4_L001_R1_001.fastq.gz",
      "/workspace/guoy1/emeson_lab/RNA-editing/Hussain/2016Apr21/FASTQ/3442-RBE-1-29526771/Hypothalamus Rat1-34710995/Data/Intensities/BaseCalls/Hypothalamus-Rat1_S4_L001_R2_001.fastq.gz"
    ],
    "Hypothalamus-Rat2_S14" => [
      "/workspace/guoy1/emeson_lab/RNA-editing/Hussain/2016Apr21/FASTQ/3442-RBE-1-29526771/Hypothalamus Rat2-34711005/Data/Intensities/BaseCalls/Hypothalamus-Rat2_S14_L001_R1_001.fastq.gz",
      "/workspace/guoy1/emeson_lab/RNA-editing/Hussain/2016Apr21/FASTQ/3442-RBE-1-29526771/Hypothalamus Rat2-34711005/Data/Intensities/BaseCalls/Hypothalamus-Rat2_S14_L001_R2_001.fastq.gz"
    ],
    "Hypothalamus-Rat3_S24" => [
      "/workspace/guoy1/emeson_lab/RNA-editing/Hussain/2016Apr21/FASTQ/3442-RBE-1-29526771/Hypothalamus Rat3-34711015/Data/Intensities/BaseCalls/Hypothalamus-Rat3_S24_L001_R1_001.fastq.gz",
      "/workspace/guoy1/emeson_lab/RNA-editing/Hussain/2016Apr21/FASTQ/3442-RBE-1-29526771/Hypothalamus Rat3-34711015/Data/Intensities/BaseCalls/Hypothalamus-Rat3_S24_L001_R2_001.fastq.gz"
    ],
    "Kidney-Rat1_S09" => [
      "/workspace/guoy1/emeson_lab/RNA-editing/Hussain/2016Apr21/FASTQ/3442-RBE-1-29526771/Kidney Rat1-34711000/Data/Intensities/BaseCalls/Kidney-Rat1_S9_L001_R1_001.fastq.gz",
      "/workspace/guoy1/emeson_lab/RNA-editing/Hussain/2016Apr21/FASTQ/3442-RBE-1-29526771/Kidney Rat1-34711000/Data/Intensities/BaseCalls/Kidney-Rat1_S9_L001_R2_001.fastq.gz"
    ],
    "Kidney-Rat2_S19" => [
      "/workspace/guoy1/emeson_lab/RNA-editing/Hussain/2016Apr21/FASTQ/3442-RBE-1-29526771/Kidney Rat2-34711010/Data/Intensities/BaseCalls/Kidney-Rat2_S19_L001_R1_001.fastq.gz",
      "/workspace/guoy1/emeson_lab/RNA-editing/Hussain/2016Apr21/FASTQ/3442-RBE-1-29526771/Kidney Rat2-34711010/Data/Intensities/BaseCalls/Kidney-Rat2_S19_L001_R2_001.fastq.gz"
    ],
    "Kidney-Rat3_S29" => [
      "/workspace/guoy1/emeson_lab/RNA-editing/Hussain/2016Apr21/FASTQ/3442-RBE-1-29526771/Kidney Rat3-34711020/Data/Intensities/BaseCalls/Kidney-Rat3_S29_L001_R1_001.fastq.gz",
      "/workspace/guoy1/emeson_lab/RNA-editing/Hussain/2016Apr21/FASTQ/3442-RBE-1-29526771/Kidney Rat3-34711020/Data/Intensities/BaseCalls/Kidney-Rat3_S29_L001_R2_001.fastq.gz"
    ],
    "Lung-Rat1_S06" => [
      "/workspace/guoy1/emeson_lab/RNA-editing/Hussain/2016Apr21/FASTQ/3442-RBE-1-29526771/Lung Rat1-34710997/Data/Intensities/BaseCalls/Lung-Rat1_S6_L001_R1_001.fastq.gz",
      "/workspace/guoy1/emeson_lab/RNA-editing/Hussain/2016Apr21/FASTQ/3442-RBE-1-29526771/Lung Rat1-34710997/Data/Intensities/BaseCalls/Lung-Rat1_S6_L001_R2_001.fastq.gz"
    ],
    "Lung-Rat2_S16" => [
      "/workspace/guoy1/emeson_lab/RNA-editing/Hussain/2016Apr21/FASTQ/3442-RBE-1-29526771/Lung Rat2-34711007/Data/Intensities/BaseCalls/Lung-Rat2_S16_L001_R1_001.fastq.gz",
      "/workspace/guoy1/emeson_lab/RNA-editing/Hussain/2016Apr21/FASTQ/3442-RBE-1-29526771/Lung Rat2-34711007/Data/Intensities/BaseCalls/Lung-Rat2_S16_L001_R2_001.fastq.gz"
    ],
    "Lung-Rat3_S26" => [
      "/workspace/guoy1/emeson_lab/RNA-editing/Hussain/2016Apr21/FASTQ/3442-RBE-1-29526771/Lung Rat3-34711017/Data/Intensities/BaseCalls/Lung-Rat3_S26_L001_R1_001.fastq.gz",
      "/workspace/guoy1/emeson_lab/RNA-editing/Hussain/2016Apr21/FASTQ/3442-RBE-1-29526771/Lung Rat3-34711017/Data/Intensities/BaseCalls/Lung-Rat3_S26_L001_R2_001.fastq.gz"
    ],
    "Pancreas-Rat1_S10" => [
      "/workspace/guoy1/emeson_lab/RNA-editing/Hussain/2016Apr21/FASTQ/3442-RBE-1-29526771/Pancreas Rat1-34711001/Data/Intensities/BaseCalls/Pancreas-Rat1_S10_L001_R1_001.fastq.gz",
      "/workspace/guoy1/emeson_lab/RNA-editing/Hussain/2016Apr21/FASTQ/3442-RBE-1-29526771/Pancreas Rat1-34711001/Data/Intensities/BaseCalls/Pancreas-Rat1_S10_L001_R2_001.fastq.gz"
    ],
    "Pancreas-Rat2_S20" => [
      "/workspace/guoy1/emeson_lab/RNA-editing/Hussain/2016Apr21/FASTQ/3442-RBE-1-29526771/Pancreas Rat2-34711011/Data/Intensities/BaseCalls/Pancreas-Rat2_S20_L001_R1_001.fastq.gz",
      "/workspace/guoy1/emeson_lab/RNA-editing/Hussain/2016Apr21/FASTQ/3442-RBE-1-29526771/Pancreas Rat2-34711011/Data/Intensities/BaseCalls/Pancreas-Rat2_S20_L001_R2_001.fastq.gz"
    ],
    "Pancreas-Rat3_S30" => [
      "/workspace/guoy1/emeson_lab/RNA-editing/Hussain/2016Apr21/FASTQ/3442-RBE-1-29526771/Pancreas Rat3-34711021/Data/Intensities/BaseCalls/Pancreas-Rat3_S30_L001_R1_001.fastq.gz",
      "/workspace/guoy1/emeson_lab/RNA-editing/Hussain/2016Apr21/FASTQ/3442-RBE-1-29526771/Pancreas Rat3-34711021/Data/Intensities/BaseCalls/Pancreas-Rat3_S30_L001_R2_001.fastq.gz"
    ],
    "Stomach-Rat1_S08" => [
      "/workspace/guoy1/emeson_lab/RNA-editing/Hussain/2016Apr21/FASTQ/3442-RBE-1-29526771/Stomach Rat1-34710999/Data/Intensities/BaseCalls/Stomach-Rat1_S8_L001_R1_001.fastq.gz",
      "/workspace/guoy1/emeson_lab/RNA-editing/Hussain/2016Apr21/FASTQ/3442-RBE-1-29526771/Stomach Rat1-34710999/Data/Intensities/BaseCalls/Stomach-Rat1_S8_L001_R2_001.fastq.gz"
    ],
    "Stomach-Rat2_S18" => [
      "/workspace/guoy1/emeson_lab/RNA-editing/Hussain/2016Apr21/FASTQ/3442-RBE-1-29526771/Stomach Rat2-34711009/Data/Intensities/BaseCalls/Stomach-Rat2_S18_L001_R1_001.fastq.gz",
      "/workspace/guoy1/emeson_lab/RNA-editing/Hussain/2016Apr21/FASTQ/3442-RBE-1-29526771/Stomach Rat2-34711009/Data/Intensities/BaseCalls/Stomach-Rat2_S18_L001_R2_001.fastq.gz"
    ],
    "Stomach-Rat3_S28" => [
      "/workspace/guoy1/emeson_lab/RNA-editing/Hussain/2016Apr21/FASTQ/3442-RBE-1-29526771/Stomach Rat3-34711019/Data/Intensities/BaseCalls/Stomach-Rat3_S28_L001_R1_001.fastq.gz",
      "/workspace/guoy1/emeson_lab/RNA-editing/Hussain/2016Apr21/FASTQ/3442-RBE-1-29526771/Stomach Rat3-34711019/Data/Intensities/BaseCalls/Stomach-Rat3_S28_L001_R2_001.fastq.gz"
    ],
    "Striatum-Rat1_S05" => [
      "/workspace/guoy1/emeson_lab/RNA-editing/Hussain/2016Apr21/FASTQ/3442-RBE-1-29526771/Striatum Rat1-34710996/Data/Intensities/BaseCalls/Striatum-Rat1_S5_L001_R1_001.fastq.gz",
      "/workspace/guoy1/emeson_lab/RNA-editing/Hussain/2016Apr21/FASTQ/3442-RBE-1-29526771/Striatum Rat1-34710996/Data/Intensities/BaseCalls/Striatum-Rat1_S5_L001_R2_001.fastq.gz"
    ],
    "Striatum-Rat2_S15" => [
      "/workspace/guoy1/emeson_lab/RNA-editing/Hussain/2016Apr21/FASTQ/3442-RBE-1-29526771/Striatum Rat2-34711006/Data/Intensities/BaseCalls/Striatum-Rat2_S15_L001_R1_001.fastq.gz",
      "/workspace/guoy1/emeson_lab/RNA-editing/Hussain/2016Apr21/FASTQ/3442-RBE-1-29526771/Striatum Rat2-34711006/Data/Intensities/BaseCalls/Striatum-Rat2_S15_L001_R2_001.fastq.gz"
    ],
    "Striatum-Rat3_S25" => [
      "/workspace/guoy1/emeson_lab/RNA-editing/Hussain/2016Apr21/FASTQ/3442-RBE-1-29526771/Striatum Rat3-34711016/Data/Intensities/BaseCalls/Striatum-Rat3_S25_L001_R1_001.fastq.gz",
      "/workspace/guoy1/emeson_lab/RNA-editing/Hussain/2016Apr21/FASTQ/3442-RBE-1-29526771/Striatum Rat3-34711016/Data/Intensities/BaseCalls/Striatum-Rat3_S25_L001_R2_001.fastq.gz"
    ],
  },
  fastqc_raw => {
    class      => "QC::FastQC",
    perform    => 1,
    target_dir => "${target_dir}/fastqc_raw",
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
  fastqc_raw_summary => {
    class      => "QC::FastQCSummary",
    perform    => 1,
    target_dir => "${target_dir}/fastqc_raw",
    option     => "",
    cqstools   => $cqstools,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "2",
      "mem"      => "10gb"
    },
  },
  fastq_join => {
    class      => "Format::FastqJoin",
    perform    => 1,
    target_dir => "${target_dir}/fastq_join",
    option     => "-p 4 -m 100",
    source_ref => "files",
    sh_direct  => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  fastq_join_fastqc => {
    class      => "QC::FastQC",
    perform    => 1,
    target_dir => "${target_dir}/fastq_join_fastqc",
    option     => "",
    source_ref => "fastq_join",
    sh_direct  => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "2",
      "mem"      => "40gb"
    },
  },
  fastq_join_fastqc_summary => {
    class      => "QC::FastQCSummary",
    perform    => 1,
    target_dir => "${target_dir}/fastq_join_fastqc",
    option     => "",
    cqstools   => $cqstools,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "2",
      "mem"      => "10gb"
    },
  },
  fastq_join_bowtie => {
    class                 => "Alignment::Bowtie1",
    perform               => 1,
    target_dir            => "${target_dir}/fastq_join_bowtie",
    option                => "",
    fasta_file            => $bowtie_fasta,
    source_ref            => "fastq_join",
    bowtie1_index         => $bowtie_index,
    add_RG_to_read        => 1,
    picard_jar            => $picard_jar,
    output_to_same_folder => 1,
    sh_direct             => 1,
    pbs                   => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  fastq_join_bowtie_snv => {
    class      => "GATK::UnifiedGenotyperCombine",
    perform    => 1,
    target_dir => "${target_dir}/fastq_join_bowtie_snv",
    option     => "-rf ReassignOneMappingQuality -RMQF 255 -RMQT 60",
    fasta_file => $bowtie_fasta,
    source_ref => "fastq_join_bowtie",
    gatk_jar   => $gatk_jar,
    by_file    => 0,
    sh_direct  => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  fastq_join_bowtie_stat => {
    class                    => "CQS::UniqueR",
    perform                  => 1,
    target_dir               => "${target_dir}/fastq_join_bowtie",
    option                   => "",
    rtemplate                => "samtoolsStatTable.R",
    output_file              => ".MappedStat",
    output_file_ext          => ".Reads.tsv",
    parameterSampleFile1_ref => [ "fastq_join_bowtie", ".stat\$" ],
    sh_direct                => 1,
    pbs                      => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  reads_summary => {
    class              => "CQS::UniqueR",
    perform            => 1,
    target_dir         => "${target_dir}/reads_summary",
    rtemplate          => "rnaeditingReads.R",
    output_file        => ".Summary",
    output_file_ext    => ".Reads.png;.Reads.tsv",
    parameterFile1_ref => [ "fastqc_raw_summary", ".FastQC.summary.reads.tsv\$" ],
    parameterFile2_ref => [ "fastq_join_fastqc_summary", ".FastQC.summary.reads.tsv\$" ],
    parameterFile3_ref => [ "fastq_join_bowtie_stat", ".Reads.tsv\$" ],
    sh_direct          => 1,
    pbs                => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "1",
      "mem"      => "10gb"
    },
  },
  fastq_join_bowtie_mismatch => {
    class        => "QC::BamMismatch",
    perform      => 1,
    target_dir   => "${target_dir}/fastq_join_bowtie_mismatch",
    option       => "",
    source_ref   => "fastq_join_bowtie",
    max_mismatch => 10,
    height_width => "4000 4000",
    sh_direct    => 1,
    pbs          => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  rnaediting_result => {
    class                 => "RNAediting::ParseMutation",
    perform               => 1,
    target_dir            => "${target_dir}/rnaediting_result",
    option                => "",
    source_ref            => "fastq_join_bowtie",
    sequence_fasta        => $bowtie_fasta,
    positions             => "72, 75, 231, 246, 331",
    picture_width         => 6000,
    picture_height        => 3000,
    percentage_thresholds => "1 2 5",
    sh_direct             => 1,
    pbs                   => {
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
      T1 => [ "fastqc_raw",         "fastq_join",                "fastq_join_fastqc",      "fastq_join_bowtie", "fastq_join_bowtie_snv" ],
      T2 => [ "fastqc_raw_summary", "fastq_join_fastqc_summary", "fastq_join_bowtie_stat", "reads_summary",     "rnaediting_result", "fastq_join_bowtie_mismatch" ],
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

1;

