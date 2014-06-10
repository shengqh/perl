#!/usr/bin/perl
use strict;
use warnings;

use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;

my $root     = create_directory_or_die("/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201403_smallRNA_2795");
my $cqstools = "/home/shengq1/cqstools/CQS.Tools.exe";

my $rn4_mrna_gff      = "/data/cqs/shengq1/reference/miRBase20/rno.gff3";
my $rn4_trna_bed      = "/data/cqs/guoy1/reference/smallrna/rn4_tRNA_ucsc_ensembl.bed";
my $rn4_trna_fasta    = "/data/cqs/guoy1/reference/smallrna/rn4_tRNA_ucsc_ensembl.bed.fa";
my $rn4_smallrna_bed  = "/data/cqs/guoy1/reference/smallrna/rn4_smallRNA_ucsc_ensembl.bed";
my $rn4_bowtie1_index = "/data/cqs/shengq1/reference/rn4/bowtie1_index/rn4";

my $email     = "quanhu.sheng\@vanderbilt.edu";
my $task_name = "Vicky2795";

my $samtools           = "/home/shengq1/local/bin/samtools/samtools";
my $bowtie1_option_1mm = "-a -m 100 --best --strata -v 1 -l 12 -p 8";

my $mirnacount_option          = "-s";                                                    #ignore score
my $trnacount_option           = "--length --sequence";
my $mirna_overlap_count_option = "-s --gtf_key miRNA";
my $mirna_fasta                = "/data/cqs/shengq1/reference/miRBase20/mature.dna.fa";

my $kcv2795 = {
  files => {
    "2795-KCV-1" => ["/gpfs21/scratch/vantage_repo/Vickers/2795/2795-KCV-1_1_sequence.txt.gz"],
    "2795-KCV-2" => ["/gpfs21/scratch/vantage_repo/Vickers/2795/2795-KCV-2_1_sequence.txt.gz"],
    "2795-KCV-3" => ["/gpfs21/scratch/vantage_repo/Vickers/2795/2795-KCV-3_1_sequence.txt.gz"],
    "2795-KCV-4" => ["/gpfs21/scratch/vantage_repo/Vickers/2795/2795-KCV-4_1_sequence.txt.gz"],
  },
  maps => {
    "2795-KCV-1" => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201403_smallRNA_2795/demultiplexing/pbs/2795-KCV-1.map"],
    "2795-KCV-2" => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201403_smallRNA_2795/demultiplexing/pbs/2795-KCV-2.map"],
    "2795-KCV-3" => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201403_smallRNA_2795/demultiplexing/pbs/2795-KCV-3.map"],
    "2795-KCV-4" => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201403_smallRNA_2795/demultiplexing/pbs/2795-KCV-4.map"],
  },
  task_name           => "2795-KCV",
  mirna_coordinate    => $rn4_mrna_gff,
  trna_coordinate     => $rn4_trna_bed,
  trna_fasta          => $rn4_trna_fasta,
  smallrna_coordinate => $rn4_smallrna_bed,
  bowtie1_index       => $rn4_bowtie1_index,
  target_dir          => $root,
};

my @defs = ($kcv2795);

foreach my $def (@defs) {
  my $target_dir = $def->{target_dir};
  my $config     = {
    general        => { "task_name" => $def->{task_name}, },
    demultiplexing => {
      class      => "Format::Demultiplexing",
      perform    => 0,
      target_dir => "${target_dir}/demultiplexing",
      option     => "",
      source     => $def->{files},
      maps       => $def->{maps},
      cqstools   => $cqstools,
      sh_direct  => 1,
      pbs        => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "24",
        "mem"      => "20gb"
      },
    },
    cutadapt => {
      class      => "Cutadapt",
      perform    => 0,
      target_dir => "${target_dir}/cutadapt",
      option     => "-O 10 -m 12",
      source_ref => "demultiplexing",
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
    identical => {
      class      => "FastqIdentical",
      perform    => 0,
      target_dir => "${target_dir}/identical",
      option     => "",
      source_ref => [ "cutadapt", ".fastq.gz" ],
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

    #not identical, for IGV
    bowtie1_genome_cutadapt_topN_1mm_notidentical => {
      class         => "Bowtie1",
      perform       => 0,
      target_dir    => "${target_dir}/topN_bowtie1_genome_cutadapt_1mm_notidentical",
      option        => $bowtie1_option_1mm,
      source_ref    => [ "cutadapt", ".fastq.gz" ],
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
      target_dir    => "${target_dir}/topN_bowtie1_genome_cutadapt_1mm",
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
      target_dir      => "${target_dir}/topN_bowtie1_genome_cutadapt_1mm_count_miRNA",
      option          => $mirnacount_option,
      source_ref      => "bowtie1_genome_cutadapt_topN_1mm",
      fastq_files_ref => "identical",
      seqcount_ref    => [ "identical", ".dupcount\$" ],
      cqs_tools       => $cqstools,
      gff_file        => $def->{mirna_coordinate},
      fasta_file      => $mirna_fasta,
      samtools        => $samtools,
      sh_direct       => 1,
      pbs             => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "72",
        "mem"      => "40gb"
      },
    },
    less_50_mirna_1mm_count => {
      class           => "MirnaCount",
      perform         => 1,
      target_dir      => "${target_dir}/less50_topN_bowtie1_genome_cutadapt_1mm_count_miRNA",
      option          => $mirnacount_option . " --maxlen 46",
      source_ref      => "bowtie1_genome_cutadapt_topN_1mm",
      fastq_files_ref => "identical",
      seqcount_ref    => [ "identical", ".dupcount\$" ],
      cqs_tools       => $cqstools,
      gff_file        => $def->{mirna_coordinate},
      fasta_file      => $mirna_fasta,
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
      target_dir => "${target_dir}/topN_bowtie1_genome_cutadapt_1mm_count_miRNA_table",
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
      target_dir      => "${target_dir}/topN_bowtie1_genome_cutadapt_1mm_count_miRNA_overlap",
      option          => $mirna_overlap_count_option,
      source_ref      => "bowtie1_genome_cutadapt_topN_1mm",
      fastq_files_ref => "identical",
      seqcount_ref    => [ "identical", ".dupcount\$" ],
      cqs_tools       => $cqstools,
      gff_file        => $def->{mirna_coordinate},
      fasta_file      => $mirna_fasta,
      samtools        => $samtools,
      sh_direct       => 1,
      pbs             => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "72",
        "mem"      => "20gb"
      },
    },
    miRNA_1mm_overlap_position => {
      class      => "CQSMappedPosition",
      perform    => 0,
      target_dir => "${target_dir}/topN_bowtie1_genome_cutadapt_1mm_count_miRNA_overlap_position",
      option     => "-o " . $task_name . "_miRNA.position",
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
      target_dir      => "${target_dir}/topN_bowtie1_genome_cutadapt_1mm_count_tRNA",
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
        "mem"      => "20gb"
      },
    },
    tRNA_1mm_table => {
      class      => "CQSMappedTable",
      perform    => 0,
      target_dir => "${target_dir}/topN_bowtie1_genome_cutadapt_1mm_count_tRNA_table",
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
      target_dir => "${target_dir}/topN_bowtie1_genome_cutadapt_1mm_count_tRNA_position",
      option     => "-o " . $task_name . "_tRNA.position",
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
      target_dir      => "${target_dir}/topN_bowtie1_genome_cutadapt_1mm_count_smallRNA",
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
        "mem"      => "20gb"
      },
    },
    smallRNA_1mm_table => {
      class      => "CQSMappedTable",
      perform    => 1,
      target_dir => "${target_dir}/topN_bowtie1_genome_cutadapt_1mm_count_smallRNA_table",
      option     => "",
      source_ref => [ "smallRNA_1mm_count", ".xml" ],
      cqs_tools  => $cqstools,
      prefix     => "smallRNA_1mm_",
      sh_direct  => 1,
      pbs        => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "10",
        "mem"      => "10gb"
      },
    },
    smallRNA_1mm_category => {
      class           => "CQSSmallRNACategory",
      perform         => 0,
      target_dir      => "${target_dir}/topN_bowtie1_genome_cutadapt_1mm_count_smallRNA_category",
      option          => "",
      source_ref      => [ "smallRNA_1mm_count", ".mapped.xml\$" ],
      mirna_count_ref => [ "mirna_1mm_count", ".mapped.xml\$" ],
      cqs_tools       => $cqstools,
      sh_direct       => 1,
      pbs             => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "72",
        "mem"      => "40gb"
      },
    },
    sequencetask => {
      class      => "CQS::SequenceTask",
      perform    => 0,
      target_dir => "${target_dir}/sequencetask",
      option     => "",
      source     => {
        prepare    => ["demultiplexing"],
        individual => [
          "cutadapt", "fastqlen", "identical",
          "bowtie1_genome_cutadapt_topN_1mm_notidentical",
          "bowtie1_genome_cutadapt_topN_1mm",
          "mirna_1mm_count", "miRNA_1mm_count_overlap", "tRNA_1mm_count", "smallRNA_1mm_count",
        ],
        summary => [ "miRNA_1mm_table", "tRNA_1mm_table", "smallRNA_1mm_table", "smallRNA_1mm_category", "miRNA_1mm_overlap_position", "tRNA_1mm_position" ],
      },
      sh_direct => 1,
      pbs       => {
        "email"    => $email,
        "nodes"    => "1:ppn=8",
        "walltime" => "72",
        "mem"      => "40gb"
      },
    },
  };

  performConfig($config);

}

1;
