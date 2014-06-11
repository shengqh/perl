#!/usr/bin/perl
use strict;
use warnings;

use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;

my $root =
  create_directory_or_die(
  "/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201406_smallRNA_2829/");
my $cqstools = "/home/shengq1/cqstools/CQS.Tools.exe";

my $hg19_mrna_gff = "/data/cqs/shengq1/reference/miRBase20/hsa.gff3";
my $hg19_trna_bed =
  "/data/cqs/guoy1/reference/smallrna/hg19_tRNA_ucsc_ensembl.bed";
my $hg19_trna_fasta =
  "/data/cqs/guoy1/reference/smallrna/hg19_tRNA_ucsc_ensembl.bed.fa";
my $hg19_smallrna_bed =
  "/data/cqs/guoy1/reference/smallrna/hg19_smallRNA_ucsc_ensembl.bed";
my $hg19_bowtie1_index =
  "/data/cqs/guoy1/reference/hg19/bowtie_index_hg19_rCRS_1.0.0/hg19_rCRS";

my $email     = "quanhu.sheng\@vanderbilt.edu";
my $task_name = "Vicky2829";

my $samtools           = "/home/shengq1/local/bin/samtools/samtools";
my $bowtie1_option_1mm = "-a -m 100 --best --strata -v 1 -l 12 -p 8";

my $mirnacount_option = "-s";                  #ignore score
my $trnacount_option  = "--length --sequence";
my $mirna_overlap_count_option = "-s --gtf_key miRNA";
my $mirna_fasta = "/data/cqs/shengq1/reference/miRBase20/mature.dna.fa";

my $kcv2829 = {
  task_name           => "2829-KCV",
  mirna_coordinate    => $hg19_mrna_gff,
  trna_coordinate     => $hg19_trna_bed,
  trna_fasta          => $hg19_trna_fasta,
  smallrna_coordinate => $hg19_smallrna_bed,
  bowtie1_index       => $hg19_bowtie1_index,
  target_dir          => $root,
};

my @defs = ($kcv2829);

foreach my $def (@defs) {
  my $target_dir = $def->{target_dir};
  my $config     = {
    general => { "task_name" => $def->{task_name}, },
    files   => {
      "2829-KCV-1A" => [
"/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201405_smallRNA_2829/raw/2829-KCV-1A_ATCACG_L008_R1_001.fastq.gz"
      ],
      "2829-KCV-1B" => [
"/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201405_smallRNA_2829/raw/2829-KCV-1B_CGATGT_L008_R1_001.fastq.gz"
      ],
      "2829-KCV-1C" => [
"/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201405_smallRNA_2829/raw/2829-KCV-1C_TTAGGC_L008_R1_001.fastq.gz"
      ],
      "2829-KCV-1D" => [
"/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201405_smallRNA_2829/raw/2829-KCV-1D_TGACCA_L008_R1_001.fastq.gz"
      ],
      "2829-KCV-1E" => [
"/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201405_smallRNA_2829/raw/2829-KCV-1E_ACAGTG_L008_R1_001.fastq.gz"
      ],
      "2829-KCV-1F" => [
"/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201405_smallRNA_2829/raw/2829-KCV-1F_GCCAAT_L008_R1_001.fastq.gz"
      ],
      "2829-KCV-1G" => [
"/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201405_smallRNA_2829/raw/2829-KCV-1G_CAGATC_L008_R1_001.fastq.gz"
      ],
      "2829-KCV-1H" => [
"/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201405_smallRNA_2829/raw/2829-KCV-1H_ACTTGA_L008_R1_001.fastq.gz"
      ],
      "2829-KCV-1I" => [
"/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201405_smallRNA_2829/raw/2829-KCV-1I_GATCAG_L008_R1_001.fastq.gz"
      ],
      "2829-KCV-1J" => [
"/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201405_smallRNA_2829/raw/2829-KCV-1J_TAGCTT_L008_R1_001.fastq.gz"
      ],
    },
    trimmer => {
      class      => "CQS::FastqTrimmer",
      perform    => 1,
      target_dir => "${target_dir}/FastqTrimmer",
      option     => "-n -z",
      extension  => "_trim.fastq.gz",
      source_ref => "files",
      cqstools   => $cqstools,
      sh_direct  => 1,
      pbs        => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "2",
        "mem"      => "10gb"
      },
    },
    cutadapt => {
      class      => "Cutadapt",
      perform    => 1,
      target_dir => "${target_dir}/cutadapt",
      option     => "-O 10 -m 12",
      source_ref => "trimmer",
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
    identical => {
      class      => "FastqIdentical",
      perform    => 1,
      target_dir => "${target_dir}/identical",
      option     => "",
      source_ref => [ "cutadapt", ".fastq.gz" ],
      cqstools   => $cqstools,
      extension  => "_clipped_identical.fastq.gz",
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
      class      => "Bowtie1",
      perform    => 1,
      target_dir =>
        "${target_dir}/topN_bowtie1_genome_cutadapt_1mm_notidentical",
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
      perform       => 1,
      target_dir    => "${target_dir}/topN_bowtie1_genome_cutadapt_1mm",
      option        => $bowtie1_option_1mm,
      source_ref    => [ "identical", ".fastq.gz\$" ],
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
      class      => "MirnaCount",
      perform    => 1,
      target_dir =>
        "${target_dir}/topN_bowtie1_genome_cutadapt_1mm_count_miRNA",
      option          => $mirnacount_option,
      source_ref      => "bowtie1_genome_cutadapt_topN_1mm",
      fastq_files_ref => [ "identical", ".fastq.gz\$" ],
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
      perform    => 1,
      target_dir =>
        "${target_dir}/topN_bowtie1_genome_cutadapt_1mm_count_miRNA_table",
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
      class      => "CQSMappedCount",
      perform    => 1,
      target_dir =>
        "${target_dir}/topN_bowtie1_genome_cutadapt_1mm_count_miRNA_overlap",
      option          => $mirna_overlap_count_option,
      source_ref      => "bowtie1_genome_cutadapt_topN_1mm",
      fastq_files_ref => [ "identical", ".fastq.gz\$" ],
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
      perform    => 1,
      target_dir =>
"${target_dir}/topN_bowtie1_genome_cutadapt_1mm_count_miRNA_overlap_position",
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
      class      => "CQSMappedCount",
      perform    => 1,
      target_dir => "${target_dir}/topN_bowtie1_genome_cutadapt_1mm_count_tRNA",
      option     => $trnacount_option,
      source_ref => "bowtie1_genome_cutadapt_topN_1mm",
      fastq_files_ref => [ "identical", ".fastq.gz\$" ],
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
      perform    => 1,
      target_dir =>
        "${target_dir}/topN_bowtie1_genome_cutadapt_1mm_count_tRNA_table",
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
      perform    => 1,
      target_dir =>
        "${target_dir}/topN_bowtie1_genome_cutadapt_1mm_count_tRNA_position",
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
      class      => "CQSMappedCount",
      perform    => 1,
      target_dir =>
        "${target_dir}/topN_bowtie1_genome_cutadapt_1mm_count_smallRNA",
      option          => $trnacount_option,
      source_ref      => "bowtie1_genome_cutadapt_topN_1mm",
      fastq_files_ref => [ "identical", ".fastq.gz\$" ],
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
      target_dir =>
        "${target_dir}/topN_bowtie1_genome_cutadapt_1mm_count_smallRNA_table",
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
      class      => "CQSSmallRNACategory",
      perform    => 1,
      target_dir =>
"${target_dir}/topN_bowtie1_genome_cutadapt_1mm_count_smallRNA_category",
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
      perform    => 1,
      target_dir => "${target_dir}/sequencetask",
      option     => "",
      source     => {
        individual => [
          "cutadapt",
          "fastqlen",
          "identical",
          #"bowtie1_genome_cutadapt_topN_1mm_notidentical",
          "bowtie1_genome_cutadapt_topN_1mm",
          "mirna_1mm_count",
          "miRNA_1mm_count_overlap",
          "tRNA_1mm_count",
          "smallRNA_1mm_count",
        ],
        summary => [
          "miRNA_1mm_table",            "tRNA_1mm_table",
          "smallRNA_1mm_table",         "smallRNA_1mm_category",
          "miRNA_1mm_overlap_position", "tRNA_1mm_position"
        ],
      },
      sh_direct => 1,
      pbs       => {
        "email"    => $email,
        "nodes"    => "1:ppn=8",
        "walltime" => "72",
        "mem"      => "40gb"
      },
    },
    
    less50_cutadapt => {
      class      => "Cutadapt",
      perform    => 1,
      target_dir => "${target_dir}/less50_cutadapt",
      option     => "-O 10 -m 12 -M 48",
      source_ref => "trimmer",
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
    less50_fastqlen => {
      class      => "FastqLen",
      perform    => 1,
      target_dir => "${target_dir}/less50_fastqlen",
      option     => "",
      source_ref => "less50_cutadapt",
      cqstools   => $cqstools,
      sh_direct  => 1,
      pbs        => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "24",
        "mem"      => "20gb"
      },
    },
    less50_identical => {
      class      => "FastqIdentical",
      perform    => 1,
      target_dir => "${target_dir}/less50_identical",
      option     => "",
      source_ref => [ "less50_cutadapt", ".fastq.gz" ],
      cqstools   => $cqstools,
      extension  => "_clipped_identical.fastq.gz",
      sh_direct  => 1,
      pbs        => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "24",
        "mem"      => "20gb"
      },
    },

    #not identical, for IGV
#    bowtie1_genome_cutadapt_topN_1mm_notidentical => {
#      class      => "Bowtie1",
#      perform    => 1,
#      target_dir =>
#        "${target_dir}/topN_bowtie1_genome_cutadapt_1mm_notidentical",
#      option        => $bowtie1_option_1mm,
#      source_ref    => [ "cutadapt", ".fastq.gz" ],
#      bowtie1_index => $def->{bowtie1_index},
#      samonly       => 0,
#      sh_direct     => 0,
#      pbs           => {
#        "email"    => $email,
#        "nodes"    => "1:ppn=8",
#        "walltime" => "72",
#        "mem"      => "40gb"
#      },
#    },

    #1 mismatch search
    less50_bowtie1_genome_cutadapt_topN_1mm => {
      class         => "Bowtie1",
      perform       => 1,
      target_dir    => "${target_dir}/less50_topN_bowtie1_genome_cutadapt_1mm",
      option        => $bowtie1_option_1mm,
      source_ref    => [ "less50_identical", ".fastq.gz\$" ],
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
    less50_mirna_1mm_count => {
      class      => "MirnaCount",
      perform    => 1,
      target_dir =>
        "${target_dir}/less50_topN_bowtie1_genome_cutadapt_1mm_count_miRNA",
      option          => $mirnacount_option,
      source_ref      => "bowtie1_genome_cutadapt_topN_1mm",
      fastq_files_ref => [ "less50_identical", ".fastq.gz\$" ],
      seqcount_ref    => [ "less50_identical", ".dupcount\$" ],
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
    less50_miRNA_1mm_table => {
      class      => "CQSMirnaTable",
      perform    => 1,
      target_dir =>
        "${target_dir}/less50_topN_bowtie1_genome_cutadapt_1mm_count_miRNA_table",
      option     => "",
      source_ref => "less50_mirna_1mm_count",
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
    less50_miRNA_1mm_count_overlap => {
      class      => "CQSMappedCount",
      perform    => 1,
      target_dir =>
        "${target_dir}/less50_topN_bowtie1_genome_cutadapt_1mm_count_miRNA_overlap",
      option          => $mirna_overlap_count_option,
      source_ref      => "less50_bowtie1_genome_cutadapt_topN_1mm",
      fastq_files_ref => [ "less50_identical", ".fastq.gz\$" ],
      seqcount_ref    => [ "less50_identical", ".dupcount\$" ],
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
    less50_miRNA_1mm_overlap_position => {
      class      => "CQSMappedPosition",
      perform    => 1,
      target_dir =>
"${target_dir}/less50_topN_bowtie1_genome_cutadapt_1mm_count_miRNA_overlap_position",
      option     => "-o " . $task_name . "_miRNA.position",
      source_ref => "less50_miRNA_1mm_count_overlap",
      cqs_tools  => $cqstools,
      sh_direct  => 1,
      pbs        => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "10",
        "mem"      => "10gb"
      },
    },
    less50_tRNA_1mm_count => {
      class      => "CQSMappedCount",
      perform    => 1,
      target_dir => "${target_dir}/less50_topN_bowtie1_genome_cutadapt_1mm_count_tRNA",
      option     => $trnacount_option,
      source_ref => "less50_bowtie1_genome_cutadapt_topN_1mm",
      fastq_files_ref => [ "less50_identical", ".fastq.gz\$" ],
      seqcount_ref    => [ "less50_identical", ".dupcount\$" ],
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
    less50_tRNA_1mm_table => {
      class      => "CQSMappedTable",
      perform    => 1,
      target_dir =>
        "${target_dir}/less50_topN_bowtie1_genome_cutadapt_1mm_count_tRNA_table",
      option     => "",
      source_ref => [ "less50_tRNA_1mm_count", ".xml" ],
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
    less50_tRNA_1mm_position => {
      class      => "CQSMappedPosition",
      perform    => 1,
      target_dir =>
        "${target_dir}/less50_topN_bowtie1_genome_cutadapt_1mm_count_tRNA_position",
      option     => "-o " . $task_name . "_tRNA.position",
      source_ref => "less50_tRNA_1mm_count",
      cqs_tools  => $cqstools,
      sh_direct  => 1,
      pbs        => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "10",
        "mem"      => "10gb"
      },
    },
    less50_smallRNA_1mm_count => {
      class      => "CQSMappedCount",
      perform    => 1,
      target_dir =>
        "${target_dir}/less50_topN_bowtie1_genome_cutadapt_1mm_count_smallRNA",
      option          => $trnacount_option,
      source_ref      => "less50_bowtie1_genome_cutadapt_topN_1mm",
      fastq_files_ref => [ "less50_identical", ".fastq.gz\$" ],
      seqcount_ref    => [ "less50_identical", ".dupcount\$" ],
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
    less50_smallRNA_1mm_table => {
      class      => "CQSMappedTable",
      perform    => 1,
      target_dir =>
        "${target_dir}/less50_topN_bowtie1_genome_cutadapt_1mm_count_smallRNA_table",
      option     => "",
      source_ref => [ "less50_smallRNA_1mm_count", ".xml" ],
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
    less50_smallRNA_1mm_category => {
      class      => "CQSSmallRNACategory",
      perform    => 1,
      target_dir =>
"${target_dir}/less50_topN_bowtie1_genome_cutadapt_1mm_count_smallRNA_category",
      option          => "",
      source_ref      => [ "less50_smallRNA_1mm_count", ".mapped.xml\$" ],
      mirna_count_ref => [ "less50_mirna_1mm_count", ".mapped.xml\$" ],
      cqs_tools       => $cqstools,
      sh_direct       => 1,
      pbs             => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "72",
        "mem"      => "40gb"
      },
    },
    less50_sequencetask => {
      class      => "CQS::SequenceTask",
      perform    => 1,
      target_dir => "${target_dir}/less50_sequencetask",
      option     => "",
      source     => {
        individual => [
          "trimmer",
          "less50_cutadapt",
          "less50_fastqlen",
          "less50_identical",
          #"less50_bowtie1_genome_cutadapt_topN_1mm_notidentical",
          "less50_bowtie1_genome_cutadapt_topN_1mm",
          "less50_mirna_1mm_count",
          "less50_miRNA_1mm_count_overlap",
          "less50_tRNA_1mm_count",
          "less50_smallRNA_1mm_count",
        ],
        summary => [
          "less50_miRNA_1mm_table",            "less50_tRNA_1mm_table",
          "less50_smallRNA_1mm_table",         "less50_smallRNA_1mm_category",
          "less50_miRNA_1mm_overlap_position", "less50_tRNA_1mm_position"
        ],
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
  #performTask( $config, "smallRNA_1mm_table" );
}

1;

