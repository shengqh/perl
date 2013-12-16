#!/usr/bin/perl
use strict;
use warnings;

use CQS::QC;
use CQS::DNASeq;
use CQS::RNASeq;
use CQS::CQSTools;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;
use CQS::Cutadapt;

my $root          = "/scratch/cqs/shengq1/mirna/20131216_liuqi";
my $cqstools      = "/home/shengq1/cqstools/CQS.Tools.exe";
my $hsa_gffs      = "/data/cqs/shengq1/reference/miRBase20/hsa.gff3";
my $hsa_trna_gffs = "/data/cqs/shengq1/reference/trna/hg19_tRNA.bed";

my $target_dir = create_directory_or_die($root);

my $email     = "quanhu.sheng\@vanderbilt.edu";
my $task_name = "mirna";

my $samtools = "/home/shengq1/local/bin/samtools/samtools";

my $bowtie1_option_pm       = "-a -m 100 --best --strata -v 0 -l 12 -p 8";
my $bowtie1_option_1mm      = "-a -m 100 --best --strata -v 1 -l 12 -p 8";
my $bowtie1_option_1mm_trim = "-a -m 100 --best --strata -v 1 -l 12 -p 8 --trim5 2 --trim3 3";
my $bowtie1_option_3mm      = "-a -m 100 --best --strata -v 3 -l 12 -p 8";

my $bowtie1_rat_index   = "/data/cqs/shengq1/reference/rn4/bowtie1_index/rn4";
my $bowtie1_human_index = "/data/cqs/guoy1/reference/hg19/bowtie_index/hg19";
my $bowtie1_mouse_index = "/data/cqs/shengq1/reference/mm10/bowtie_index/mm10";

my $mirnacount_option          = "-s";                                                    #ignore score
my $trnacount_option           = "--length --sequence";
my $mirna_overlap_count_option = "-s --gtf_key miRNA";
my $fasta_file                 = "/data/cqs/shengq1/reference/miRBase20/mature.dna.fa";

my $trna_hg19_fasta = "/data/cqs/shengq1/reference/trna/hg19_tRNA.bed.fa";
my $trna_mm10_fasta = "/data/cqs/shengq1/reference/trna/mm10_tRNA.bed.fa";
my $trna_rn4_fasta  = "/data/cqs/shengq1/reference/trna/rn4_tRNA.bed.fa";

my $shrimp2_option              = "-Q -N 8 -o 1 --qv-offset 33";
my $shrimp2_rat_miRBase_index   = "/data/cqs/shengq1/reference/miRBase19/shrimp_index_2.2.3/rno.mature.dna-ls";
my $shrimp2_human_miRBase_index = "/data/cqs/shengq1/reference/miRBase19/shrimp_index_2.2.3/hsa.mature.dna-ls";
my $shrimp2_mouse_miRBase_index = "/data/cqs/shengq1/reference/miRBase19/shrimp_index_2.2.3/mmu.mature.dna-ls";

my $bwa_option       = "-o 0 -l 8 -n 3 -t 8";
my $bwa_hsammu_fasta = "/data/cqs/shengq1/reference/hg19mm9/bwa_0.7.4_index/hg19mm9.fa";
my $hsammu_gffs      = "/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/smrnapipeline/hsa_mmu_tableL.bed";
my $bwa_clip_option  = "-o 2 -e 3 -l 8 -n 3 -t 8";

my $exrna = {
  mirna => {
    "SL35457" => ["/gpfs21/scratch/liuq6/exrna/cutadpt/result/H148TADXX.s2.0.illumina12index.1.SL35457.clipped.fastq"],
    "SL35458" => ["/gpfs21/scratch/liuq6/exrna/cutadpt/result/H148TADXX.s2.0.illumina12index.2.SL35458.clipped.fastq"],
    "SL35459" => ["/gpfs21/scratch/liuq6/exrna/cutadpt/result/H148TADXX.s2.0.illumina12index.3.SL35459.clipped.fastq"],
    "SL35460" => ["/gpfs21/scratch/liuq6/exrna/cutadpt/result/H148TADXX.s2.0.illumina12index.4.SL35460.clipped.fastq"],
    "SL35461" => ["/gpfs21/scratch/liuq6/exrna/cutadpt/result/H148TADXX.s2.0.illumina12index.5.SL35461.clipped.fastq"],
    "SL35462" => ["/gpfs21/scratch/liuq6/exrna/cutadpt/result/H148TADXX.s2.0.illumina12index.6.SL35462.clipped.fastq"],
    "SL41462" => ["/gpfs21/scratch/liuq6/exrna/cutadpt/result/C315CACXX.s8.0.illumina12index.1.SL41462.clipped.fastq"],
    "SL41463" => ["/gpfs21/scratch/liuq6/exrna/cutadpt/result/C315CACXX.s8.0.illumina12index.2.SL41463.clipped.fastq"],
    "SL41464" => ["/gpfs21/scratch/liuq6/exrna/cutadpt/result/C315CACXX.s8.0.illumina12index.3.SL41464.clipped.fastq"],
    "SL41465" => ["/gpfs21/scratch/liuq6/exrna/cutadpt/result/C315CACXX.s8.0.illumina12index.4.SL41465.clipped.fastq"],
    "SL41466" => ["/gpfs21/scratch/liuq6/exrna/cutadpt/result/C315CACXX.s8.0.illumina12index.5.SL41466.clipped.fastq"],
    "SL41467" => ["/gpfs21/scratch/liuq6/exrna/cutadpt/result/C315CACXX.s8.0.illumina12index.6.SL41467.clipped.fastq"],
    "SL41468" => ["/gpfs21/scratch/liuq6/exrna/cutadpt/result/C315CACXX.s8.0.illumina12index.7.SL41468.clipped.fastq"],
    "SL41469" => ["/gpfs21/scratch/liuq6/exrna/cutadpt/result/C315CACXX.s8.0.illumina12index.8.SL41469.clipped.fastq"],
    "SL41470" => ["/gpfs21/scratch/liuq6/exrna/cutadpt/result/C315CACXX.s8.0.illumina12index.9.SL41470.clipped.fastq"],
    "SL41471" => ["/gpfs21/scratch/liuq6/exrna/cutadpt/result/C315CACXX.s8.0.illumina12index.10.SL41471.clipped.fastq"],
    "SL41473" => ["/gpfs21/scratch/liuq6/exrna/cutadpt/result/C315CACXX.s8.0.illumina12index.12.SL41473.clipped.fastq"],
  },
  coordinate          => $hsa_gffs,
  trna_coordinate     => $hsa_trna_gffs,
  trna_fasta          => $trna_hg19_fasta,
  smallrna_coordinate => "/gpfs21/scratch/cqs/shengq1/references/smallrna/Homo_sapiens.GRCh37.73.smallRNA.bed",
  bowtie1_index       => $bowtie1_human_index,
  shrimp2_index       => $shrimp2_human_miRBase_index,
  target_dir          => $target_dir,
  task_name           => $task_name,
};

my @defs = ($exrna);

foreach my $def (@defs) {
  my $cur_target_dir = $def->{target_dir};
  my $config         = {
    general  => { "task_name" => $def->{task_name}, },
    fastqlen => {
      class      => "FastqLen",
      perform    => 0,
      target_dir => "${cur_target_dir}/fastqlen",
      option     => "",
      source     => $def->{mirna},
      cqstools   => $cqstools,
      sh_direct  => 1,
      pbs        => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "24",
        "mem"      => "20gb"
      },
    },
    identical => {
      class      => "FastqIdentical",
      perform    => 1,
      target_dir => "${cur_target_dir}/identical",
      option     => "-n",
      source     => $def->{mirna},
      cqstools   => $cqstools,
      extension  => "_clipped_identical.fastq",
      sh_direct  => 1,
      pbs        => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "24",
        "mem"      => "20gb"
      },
    },

    #1 mismatch search
    bowtie1_genome_cutadapt_topN_1mm => {
      class         => "Bowtie1",
      perform       => 0,
      target_dir    => "${cur_target_dir}/topN_bowtie1_genome_cutadapt_1mm",
      option        => $bowtie1_option_1mm,
      source_ref    => [ "identical", ".fastq\$" ],
      bowtie1_index => $def->{bowtie1_index},
      samonly       => 0,
      sh_direct     => 1,
      pbs           => {
        "email"    => $email,
        "nodes"    => "1:ppn=8",
        "walltime" => "72",
        "mem"      => "40gb"
      },
    },
    mirna_1mm_count => {
      class           => "MirnaCount",
      perform         => 0,
      target_dir      => "${cur_target_dir}/topN_bowtie1_genome_cutadapt_1mm_count_miRNA",
      option          => $mirnacount_option,
      source_ref      => "bowtie1_genome_cutadapt_topN_1mm",
      fastq_files_ref => "identical",
      seqcount_ref    => [ "identical", ".dupcount\$" ],
      cqs_tools       => $cqstools,
      gff_file        => $def->{coordinate},
      fasta_file      => $fasta_file,
      samtools        => $samtools,
      sh_direct       => 1,
      pbs             => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "72",
        "mem"      => "40gb"
      },
    },
    miRNA_1mm_table => {
      class      => "CQSMirnaTable",
      perform    => 0,
      target_dir => "${cur_target_dir}/topN_bowtie1_genome_cutadapt_1mm_count_miRNA_table",
      option     => "",
      source_ref => "mirna_1mm_count",
      cqs_tools  => $cqstools,
      prefix     => "miRNA_1mm_",
      sh_direct  => 1,
      pbs        => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "10",
        "mem"      => "10gb"
      },
    },
    miRNA_1mm_count_overlap => {
      class           => "CQSMappedCount",
      perform         => 0,
      target_dir      => "${cur_target_dir}/topN_bowtie1_genome_cutadapt_1mm_count_miRNA_overlap",
      option          => $mirna_overlap_count_option,
      source_ref      => "bowtie1_genome_cutadapt_topN_1mm",
      fastq_files_ref => "identical",
      seqcount_ref    => [ "identical", ".dupcount\$" ],
      cqs_tools       => $cqstools,
      gff_file        => $def->{coordinate},
      fasta_file      => $fasta_file,
      samtools        => $samtools,
      sh_direct       => 1,
      pbs             => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "72",
        "mem"      => "40gb"
      },
    },
    miRNA_1mm_overlap_position => {
      class      => "CQSMappedPosition",
      perform    => 0,
      target_dir => "${cur_target_dir}/topN_bowtie1_genome_cutadapt_1mm_count_miRNA_overlap_position",
      option     => "-o " . $def->{task_name} . "_miRNA.position",
      source_ref => "miRNA_1mm_count_overlap",
      cqs_tools  => $cqstools,
      sh_direct  => 1,
      pbs        => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "10",
        "mem"      => "10gb"
      },
    },
    tRNA_1mm_count => {
      class           => "CQSMappedCount",
      perform         => 0,
      target_dir      => "${cur_target_dir}/topN_bowtie1_genome_cutadapt_1mm_count_tRNA",
      option          => $trnacount_option,
      source_ref      => "bowtie1_genome_cutadapt_topN_1mm",
      fastq_files_ref => "identical",
      seqcount_ref    => [ "identical", ".dupcount\$" ],
      cqs_tools       => $cqstools,
      gff_file        => $def->{trna_coordinate},
      fasta_file      => $def->{trna_fasta},
      samtools        => $samtools,
      sh_direct       => 1,
      pbs             => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "72",
        "mem"      => "40gb"
      },
    },
    tRNA_1mm_table => {
      class      => "CQSMappedTable",
      perform    => 0,
      target_dir => "${cur_target_dir}/topN_bowtie1_genome_cutadapt_1mm_count_tRNA_table",
      option     => "",
      source_ref => [ "tRNA_1mm_count", ".xml" ],
      cqs_tools  => $cqstools,
      prefix     => "tRNA_1mm_",
      sh_direct  => 1,
      pbs        => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "10",
        "mem"      => "10gb"
      },
    },
    tRNA_1mm_position => {
      class      => "CQSMappedPosition",
      perform    => 0,
      target_dir => "${cur_target_dir}/topN_bowtie1_genome_cutadapt_1mm_count_tRNA_position",
      option     => "-o " . $def->{task_name} . "_tRNA.position",
      source_ref => "tRNA_1mm_count",
      cqs_tools  => $cqstools,
      sh_direct  => 1,
      pbs        => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "10",
        "mem"      => "10gb"
      },
    },
    smallRNA_1mm_count => {
      class           => "CQSMappedCount",
      perform         => 0,
      target_dir      => "${cur_target_dir}/topN_bowtie1_genome_cutadapt_1mm_count_smallRNA",
      option          => $trnacount_option,
      source_ref      => "bowtie1_genome_cutadapt_topN_1mm",
      fastq_files_ref => "identical",
      seqcount_ref    => [ "identical", ".dupcount\$" ],
      cqs_tools       => $cqstools,
      gff_file        => $def->{smallrna_coordinate},
      samtools        => $samtools,
      sh_direct       => 1,
      pbs             => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "72",
        "mem"      => "40gb"
      },
    },
    smallRNA_1mm_category => {
      class           => "CQSSmallRNACategory",
      perform         => 0,
      target_dir      => "${cur_target_dir}/topN_bowtie1_genome_cutadapt_1mm_count_smallRNA_category",
      option          => "",
      source_ref      => [ "smallRNA_1mm_count", ".mapped.xml\$" ],
      mirna_count_ref => [ "mirna_1mm_count", ".mapped.xml\$" ],
      groups          => $def->{groups},
      cqs_tools       => $cqstools,
      sh_direct       => 1,
      pbs             => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "72",
        "mem"      => "40gb"
      },
    },
  };

  performConfig($config);
}

1;
