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

my $root          = "/scratch/cqs/shengq1/mirna/20131127_guoyan_mirna";
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

my $human = {
  mrna => {
    "2163-CRF-01" => ["/autofs/blue_sequencer/Runs/projects/2163-CRF/2012-06-26/2163-CRF-1_1.fastq.gz"],
    "2163-CRF-02" => ["/autofs/blue_sequencer/Runs/projects/2163-CRF/2012-06-26/2163-CRF-2_1.fastq.gz"],
    "2163-CRF-03" => ["/autofs/blue_sequencer/Runs/projects/2163-CRF/2012-06-26/2163-CRF-3_1.fastq.gz"],
    "2163-CRF-04" => ["/autofs/blue_sequencer/Runs/projects/2163-CRF/2012-06-26/2163-CRF-4_1.fastq.gz"],
    "2163-CRF-05" => ["/autofs/blue_sequencer/Runs/projects/2163-CRF/2012-06-26/2163-CRF-5_1.fastq.gz"],
    "2163-CRF-06" => ["/autofs/blue_sequencer/Runs/projects/2163-CRF/2012-06-26/2163-CRF-6_1.fastq.gz"],
    "2163-CRF-07" => ["/autofs/blue_sequencer/Runs/projects/2163-CRF/2012-06-26/2163-CRF-7_1.fastq.gz"],
    "2163-CRF-08" => ["/autofs/blue_sequencer/Runs/projects/2163-CRF/2012-06-26/2163-CRF-8_1.fastq.gz"],
    "2163-CRF-09" => ["/autofs/blue_sequencer/Runs/projects/2163-CRF/2012-08-03/2163-CRF-9_1.fastq.gz"],
    "2163-CRF-10" => ["/autofs/blue_sequencer/Runs/projects/2163-CRF/2012-08-03/2163-CRF-10_1.fastq.gz"],
    "2163-CRF-11" => ["/autofs/blue_sequencer/Runs/projects/2163-CRF/2012-08-03/2163-CRF-11_1.fastq.gz"],
    "2163-CRF-12" => ["/autofs/blue_sequencer/Runs/projects/2163-CRF/2012-08-03/2163-CRF-12_1.fastq.gz"]
  },
  mirna => {
    "2245-CRF-00" => ["/autofs/blue_sequencer/Runs/projects/2245-CRF/2012-08-03/2245-CRF-0_1.fastq.gz"],
    "2245-CRF-01" => ["/autofs/blue_sequencer/Runs/projects/2245-CRF/2012-08-03/2245-CRF-1_1.fastq.gz"],
    "2245-CRF-02" => ["/autofs/blue_sequencer/Runs/projects/2245-CRF/2012-08-03/2245-CRF-2_1.fastq.gz"],
    "2245-CRF-03" => ["/autofs/blue_sequencer/Runs/projects/2245-CRF/2012-08-03/2245-CRF-3_1.fastq.gz"],
    "2245-CRF-04" => ["/autofs/blue_sequencer/Runs/projects/2245-CRF/2012-08-03/2245-CRF-4_1.fastq.gz"],
    "2245-CRF-05" => ["/autofs/blue_sequencer/Runs/projects/2245-CRF/2012-08-03/2245-CRF-5_1.fastq.gz"],
    "2245-CRF-06" => ["/autofs/blue_sequencer/Runs/projects/2245-CRF/2012-08-03/2245-CRF-6_1.fastq.gz"],
    "2245-CRF-07" => ["/autofs/blue_sequencer/Runs/projects/2245-CRF/2012-08-03/2245-CRF-7_1.fastq.gz"],
    "2245-CRF-08" => ["/autofs/blue_sequencer/Runs/projects/2245-CRF/2012-08-03/2245-CRF-8_1.fastq.gz"],
    "2245-CRF-09" => ["/autofs/blue_sequencer/Runs/projects/2245-CRF/2012-08-03/2245-CRF-9_1.fastq.gz"],
    "2245-CRF-10" => ["/autofs/blue_sequencer/Runs/projects/2245-CRF/2012-08-03/2245-CRF-10_1.fastq.gz"],
    "2245-CRF-11" => ["/autofs/blue_sequencer/Runs/projects/2245-CRF/2012-08-03/2245-CRF-11_1.fastq.gz"],
    "2403-CRF-01" => ["/autofs/blue_sequencer/Runs/projects/2403-CRF/2013-01-22/2403-CRF-1_1.fastq.gz"],
    "2403-CRF-02" => ["/autofs/blue_sequencer/Runs/projects/2403-CRF/2013-01-22/2403-CRF-2_1.fastq.gz"],
    "2403-CRF-03" => ["/autofs/blue_sequencer/Runs/projects/2403-CRF/2013-01-22/2403-CRF-3_1.fastq.gz"],
    "2403-CRF-04" => ["/autofs/blue_sequencer/Runs/projects/2403-CRF/2013-01-22/2403-CRF-4_1.fastq.gz"],
    "2403-CRF-05" => ["/autofs/blue_sequencer/Runs/projects/2403-CRF/2013-01-22/2403-CRF-5_1.fastq.gz"],
    "2403-CRF-06" => ["/autofs/blue_sequencer/Runs/projects/2403-CRF/2013-01-22/2403-CRF-6_1.fastq.gz"],
    "2403-CRF-07" => ["/autofs/blue_sequencer/Runs/projects/2403-CRF/2013-01-22/2403-CRF-7_1.fastq.gz"],
    "2403-CRF-08" => ["/autofs/blue_sequencer/Runs/projects/2403-CRF/2013-01-22/2403-CRF-8_1.fastq.gz"],
    "2403-CRF-09" => ["/autofs/blue_sequencer/Runs/projects/2403-CRF/2013-01-22/2403-CRF-9_1.fastq.gz"],
    "2403-CRF-10" => ["/autofs/blue_sequencer/Runs/projects/2403-CRF/2013-01-22/2403-CRF-10_1.fastq.gz"],
    "2403-CRF-11" => ["/autofs/blue_sequencer/Runs/projects/2403-CRF/2013-01-22/2403-CRF-11_1.fastq.gz"],
    "2403-CRF-12" => ["/autofs/blue_sequencer/Runs/projects/2403-CRF/2013-01-22/2403-CRF-12_1.fastq.gz"],
    "2403-CRF-13" => ["/autofs/blue_sequencer/Runs/projects/2403-CRF/2013-01-22/2403-CRF-13_1.fastq.gz"],
    "2403-CRF-14" => ["/autofs/blue_sequencer/Runs/projects/2403-CRF/2013-01-22/2403-CRF-14_1.fastq.gz"],
    "2403-CRF-15" => ["/autofs/blue_sequencer/Runs/projects/2403-CRF/2013-01-22/2403-CRF-15_1.fastq.gz"],
    "2403-CRF-18" => ["/autofs/blue_sequencer/Runs/projects/2403-CRF/2013-01-22/2403-CRF-18_1.fastq.gz"],
    "2403-CRF-19" => ["/autofs/blue_sequencer/Runs/projects/2403-CRF/2013-01-22/2403-CRF-19_1.fastq.gz"],
    "2403-CRF-20" => ["/autofs/blue_sequencer/Runs/projects/2403-CRF/2013-01-22/2403-CRF-20_1.fastq.gz"],
    "2403-CRF-21" => ["/autofs/blue_sequencer/Runs/projects/2403-CRF/2013-01-22/2403-CRF-21_1.fastq.gz"],
    "2403-CRF-22" => ["/autofs/blue_sequencer/Runs/projects/2403-CRF/2013-01-22/2403-CRF-22_1.fastq.gz"],
    "2403-CRF-23" => ["/autofs/blue_sequencer/Runs/projects/2403-CRF/2013-01-22/2403-CRF-23_1.fastq.gz"],
    "2403-CRF-24" => ["/autofs/blue_sequencer/Runs/projects/2403-CRF/2013-01-22/2403-CRF-24_1.fastq.gz"],
    "2403-CRF-26" => ["/autofs/blue_sequencer/Runs/projects/2403-CRF/2013-01-22/2403-CRF-26_1.fastq.gz"],
    "2403-CRF-28" => ["/autofs/blue_sequencer/Runs/projects/2403-CRF/2013-01-22/2403-CRF-28_1.fastq.gz"],
  },
  coordinate          => $hsa_gffs,
  trna_coordinate     => $hsa_trna_gffs,
  trna_fasta          => $trna_hg19_fasta,
  smallrna_coordinate => "/gpfs21/scratch/cqs/shengq1/references/smallrna/Homo_sapiens.GRCh37.73.smallRNA.bed",
  bowtie1_index       => $bowtie1_human_index,
  shrimp2_index       => $shrimp2_human_miRBase_index,
  target_dir          => $target_dir,
  task_name           => $task_name . "_human",
  mrna_groups              => {
    "2163-CRF" =>
      [ "2163-CRF-01", "2163-CRF-10", "2163-CRF-11", "2163-CRF-12", "2163-CRF-02", "2163-CRF-03", "2163-CRF-04", "2163-CRF-05", "2163-CRF-06", "2163-CRF-07", "2163-CRF-08", "2163-CRF-09" ],
  },
  groups              => {
    "2245-CRF" =>
      [ "2245-CRF-00", "2245-CRF-01", "2245-CRF-10", "2245-CRF-11", "2245-CRF-02", "2245-CRF-03", "2245-CRF-04", "2245-CRF-05", "2245-CRF-06", "2245-CRF-07", "2245-CRF-08", "2245-CRF-09" ],
    "2403-CRF" => [
      "2403-CRF-01", "2403-CRF-10", "2403-CRF-11", "2403-CRF-12", "2403-CRF-13", "2403-CRF-14", "2403-CRF-15", "2403-CRF-18", "2403-CRF-19", "2403-CRF-02", "2403-CRF-20", "2403-CRF-21",
      "2403-CRF-22", "2403-CRF-23", "2403-CRF-24", "2403-CRF-26", "2403-CRF-28", "2403-CRF-03", "2403-CRF-04", "2403-CRF-05", "2403-CRF-06", "2403-CRF-07", "2403-CRF-08", "2403-CRF-09"
    ],
  },
};

my @defs = ($human);

foreach my $def (@defs) {
  my $cur_target_dir = $def->{target_dir};
  my $config         = {
    general  => { "task_name" => $def->{task_name}, },
    cutadapt_len => {
      class      => "Cutadapt",
      perform    => 0,
      target_dir => "${cur_target_dir}/cutadapt_len",
      option     => "-O 10 -m 12",
      source     => $def->{mirna},
      adaptor    => "TGGAATTCTCGGGTGCCAAGG",
      extension  => "_clipped.fastq",
      sh_direct  => 1,
      pbs        => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "24",
        "mem"      => "20gb"
      },
    },
    fastqlen => {
      class      => "FastqLen",
      perform    => 0,
      target_dir => "${cur_target_dir}/fastqlen",
      option     => "",
      source_ref => "cutadapt_len",
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
      perform    => 0,
      target_dir => "${cur_target_dir}/identical",
      option     => "",
      source_ref => ["cutadapt_len", "fastq.gz\$"],
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

    #1 mismatch notidentical search
    bowtie1_genome_cutadapt_topN_1mm_notidentical => {
      class         => "Bowtie1",
      perform       => 1,
      target_dir    => "${cur_target_dir}/topN_bowtie1_genome_cutadapt_1mm_notidentical",
      option        => $bowtie1_option_1mm,
      source_ref    => "cutadapt_len",
      bowtie1_index => $def->{bowtie1_index},
      samonly       => 0,
      sh_direct     => 0,
      pbs           => {
        "email"    => $email,
        "nodes"    => "1:ppn=8",
        "walltime" => "72",
        "mem"      => "40gb"
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
      groups     => $def->{groups},
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
    miRNA_1mm_overlap_table => {
      class      => "CQSMappedTable",
      perform    => 0,
      target_dir => "${cur_target_dir}/topN_bowtie1_genome_cutadapt_1mm_count_miRNA_overlap_table",
      option     => "",
      source_ref => "miRNA_1mm_count_overlap",
      groups     => $def->{groups},
      cqs_tools  => $cqstools,
      prefix     => "miRNA_overlap_1mm_",
      sh_direct  => 1,
      pbs        => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "10",
        "mem"      => "10gb"
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
      groups     => $def->{groups},
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
