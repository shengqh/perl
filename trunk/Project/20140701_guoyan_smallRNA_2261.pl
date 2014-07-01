#!/usr/bin/perl
use strict;
use warnings;

use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;

my $root     = create_directory_or_die("/scratch/cqs/shengq1/smallRNA/20140701_guoyan_smallRNA_2261/");
my $cqstools = "/home/shengq1/cqstools/CQS.Tools.exe";

my $hg19_mrna_gff      = "/data/cqs/shengq1/reference/miRBase20/hsa.gff3";
my $hg19_trna_bed      = "/data/cqs/guoy1/reference/smallrna/hg19_tRNA_ucsc_ensembl.bed";
my $hg19_trna_fasta    = "/data/cqs/guoy1/reference/smallrna/hg19_tRNA_ucsc_ensembl.bed.fa";
my $hg19_smallrna_bed  = "/data/cqs/guoy1/reference/smallrna/hg19_smallRNA_ucsc_ensembl.bed";
my $hg19_bowtie1_index = "/data/cqs/guoy1/reference/hg19/bowtie_index_hg19_rCRS_1.0.0/hg19_rCRS";

my $email     = "quanhu.sheng\@vanderbilt.edu";
my $task_name = "2261";

my $samtools           = "/home/shengq1/local/bin/samtools/samtools";
my $bowtie1_option_1mm = "-a -m 100 --best --strata -v 1 -l 12 -p 8";

my $mirnacount_option          = "-s";                                                    #ignore score
my $trnacount_option           = "--length --sequence";
my $mirna_overlap_count_option = "-s --gtf_key miRNA";
my $mirna_fasta                = "/data/cqs/shengq1/reference/miRBase20/mature.dna.fa";

my $kcv2829 = {
  task_name           => "2261",
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
      "2261-ASK-001"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_1/2261-ASK-1_1.fastq"],
      "2261-ASK-002"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_1/2261-ASK-2_1.fastq"],
      "2261-ASK-003"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_1/2261-ASK-3_1.fastq"],
      "2261-ASK-004"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_1/2261-ASK-4_1.fastq"],
      "2261-ASK-005"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_1/2261-ASK-5_1.fastq"],
      "2261-ASK-006"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_1/2261-ASK-6_1.fastq"],
      "2261-ASK-007"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_1/2261-ASK-7_1.fastq"],
      "2261-ASK-008"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_2/2261-ASK-8_1.fastq"],
      "2261-ASK-009"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_2/2261-ASK-9_1.fastq"],
      "2261-ASK-010"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_2/2261-ASK-10_1.fastq"],
      "2261-ASK-011"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_2/2261-ASK-11_1.fastq"],
      "2261-ASK-012"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_2/2261-ASK-12_1.fastq"],
      "2261-ASK-013"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_2/2261-ASK-13_1.fastq"],
      "2261-ASK-014"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_2/2261-ASK-14_1.fastq"],
      "2261-ASK-015"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_2/2261-ASK-15_1.fastq"],
      "2261-ASK-018"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_2/2261-ASK-18_1.fastq"],
      "2261-ASK-019"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_2/2261-ASK-19_1.fastq"],
      "2261-ASK-020"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_2/2261-ASK-20_1.fastq"],
      "2261-ASK-021"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_2/2261-ASK-21_1.fastq"],
      "2261-ASK-022"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_2/2261-ASK-22_1.fastq"],
      "2261-ASK-023"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_2/2261-ASK-23_1.fastq"],
      "2261-ASK-024"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_2/2261-ASK-24_1.fastq"],
      "2261-ASK-025"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_2/2261-ASK-25_1.fastq"],
      "2261-ASK-026"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_2/2261-ASK-26_1.fastq"],
      "2261-ASK-027"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_2/2261-ASK-27_1.fastq"],
      "2261-ASK-029"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_2/2261-ASK-29_1.fastq"],
      "2261-ASK-030"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_2/2261-ASK-30_1.fastq"],
      "2261-ASK-032"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_2/2261-ASK-32_1.fastq"],
      "2261-ASK-033"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_2/2261-ASK-33_1.fastq"],
      "2261-ASK-034"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_2/2261-ASK-34_1.fastq"],
      "2261-ASK-035"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_2/2261-ASK-35_1.fastq"],
      "2261-ASK-036"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_2/2261-ASK-36_1.fastq"],
      "2261-ASK-037"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_2/2261-ASK-37_1.fastq"],
      "2261-ASK-038"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_2/2261-ASK-38_1.fastq"],
      "2261-ASK-039"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_2/2261-ASK-39_1.fastq"],
      "2261-ASK-040"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_2/2261-ASK-40_1.fastq"],
      "2261-ASK-041"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_2/2261-ASK-41_1.fastq"],
      "2261-ASK-042"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_2/2261-ASK-42_1.fastq"],
      "2261-ASK-043"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_2/2261-ASK-43_1.fastq"],
      "2261-ASK-044"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_2/2261-ASK-44_1.fastq"],
      "2261-ASK-045"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_2/2261-ASK-45_1.fastq"],
      "2261-ASK-046"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_2/2261-ASK-46_1.fastq"],
      "2261-ASK-047"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_2/2261-ASK-47_1.fastq"],
      "2261-ASK-049"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_2/2261-ASK-49_1.fastq"],
      "2261-ASK-050"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_2/2261-ASK-50_1.fastq"],
      "2261-ASK-051"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_2/2261-ASK-51_1.fastq"],
      "2261-ASK-052"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_2/2261-ASK-52_1.fastq"],
      "2261-ASK-054"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_2/2261-ASK-54_1.fastq"],
      "2261-ASK-056"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_2/2261-ASK-56_1.fastq"],
      "2261-ASK-057"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_4/2261-ASK-57_1.fastq"],
      "2261-ASK-058a"  => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_2/2261-ASK-58a_1.fastq"],
      "2261-ASK-058ab" => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_3/2261-ASK-58ab_1.fastq"],
      "2261-ASK-058b"  => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_3/2261-ASK-58b_1.fastq"],
      "2261-ASK-059"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_4/2261-ASK-59_1.fastq"],
      "2261-ASK-060a"  => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_2/2261-ASK-60a_1.fastq"],
      "2261-ASK-060ab" => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_3/2261-ASK-60ab_1.fastq"],
      "2261-ASK-060b"  => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_3/2261-ASK-60b_1.fastq"],
      "2261-ASK-061"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_4/2261-ASK-61_1.fastq"],
      "2261-ASK-062"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_4/2261-ASK-62_1.fastq"],
      "2261-ASK-063"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_4/2261-ASK-63_1.fastq"],
      "2261-ASK-064"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_3/2261-ASK-64_1.fastq"],
      "2261-ASK-065"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_3/2261-ASK-65_1.fastq"],
      "2261-ASK-066"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_4/2261-ASK-66_1.fastq"],
      "2261-ASK-067"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_3/2261-ASK-67_1.fastq"],
      "2261-ASK-068"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_4/2261-ASK-68_1.fastq"],
      "2261-ASK-069"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_3/2261-ASK-69_1.fastq"],
      "2261-ASK-070"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_3/2261-ASK-70_1.fastq"],
      "2261-ASK-071"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_4/2261-ASK-71_1.fastq"],
      "2261-ASK-072"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_3/2261-ASK-72_1.fastq"],
      "2261-ASK-073"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_3/2261-ASK-73_1.fastq"],
      "2261-ASK-074"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_3/2261-ASK-74_1.fastq"],
      "2261-ASK-075"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_3/2261-ASK-75_1.fastq"],
      "2261-ASK-077"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_3/2261-ASK-77_1.fastq"],
      "2261-ASK-078"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_3/2261-ASK-78_1.fastq"],
      "2261-ASK-080"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_3/2261-ASK-80_1.fastq"],
      "2261-ASK-082"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_4/2261-ASK-82_1.fastq"],
      "2261-ASK-083"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_3/2261-ASK-83_1.fastq"],
      "2261-ASK-084"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_3/2261-ASK-84_1.fastq"],
      "2261-ASK-086"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_3/2261-ASK-86_1.fastq"],
      "2261-ASK-087"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_3/2261-ASK-87_1.fastq"],
      "2261-ASK-088"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_3/2261-ASK-88_1.fastq"],
      "2261-ASK-089"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_3/2261-ASK-89_1.fastq"],
      "2261-ASK-090"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_3/2261-ASK-90_1.fastq"],
      "2261-ASK-091"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_3/2261-ASK-91_1.fastq"],
      "2261-ASK-092"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_3/2261-ASK-92_1.fastq"],
      "2261-ASK-093"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_4/2261-ASK-93_1.fastq"],
      "2261-ASK-094"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_4/2261-ASK-94_1.fastq"],
      "2261-ASK-095"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_4/2261-ASK-95_1.fastq"],
      "2261-ASK-096"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_4/2261-ASK-96_1.fastq"],
      "2261-ASK-097"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_3/2261-ASK-97_1.fastq"],
      "2261-ASK-098"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_3/2261-ASK-98_1.fastq"],
      "2261-ASK-099"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_3/2261-ASK-99_1.fastq"],
      "2261-ASK-100"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_3/2261-ASK-100_1.fastq"],
      "2261-ASK-101"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_3/2261-ASK-101_1.fastq"],
      "2261-ASK-102"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_4/2261-ASK-102_1.fastq"],
      "2261-ASK-103"   => ["/gpfs21/scratch/cqs/guom1/2261/rawdata_4/2261-ASK-103_1.fastq"],
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
      class         => "Bowtie1",
      perform       => 1,
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
      class           => "MirnaCount",
      perform         => 1,
      target_dir      => "${target_dir}/topN_bowtie1_genome_cutadapt_1mm_count_miRNA",
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
      perform         => 1,
      target_dir      => "${target_dir}/topN_bowtie1_genome_cutadapt_1mm_count_miRNA_overlap",
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
      perform         => 1,
      target_dir      => "${target_dir}/topN_bowtie1_genome_cutadapt_1mm_count_tRNA",
      option          => $trnacount_option,
      source_ref      => "bowtie1_genome_cutadapt_topN_1mm",
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
      perform    => 1,
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
      perform         => 1,
      target_dir      => "${target_dir}/topN_bowtie1_genome_cutadapt_1mm_count_smallRNA",
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
      perform         => 1,
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
      perform    => 1,
      target_dir => "${target_dir}/sequencetask",
      option     => "",
      source     => {
        individual => [
          "cutadapt",
          "fastqlen",
          "identical",
          "bowtie1_genome_cutadapt_topN_1mm",
          "mirna_1mm_count",
          "miRNA_1mm_count_overlap",
          "tRNA_1mm_count",
          "smallRNA_1mm_count",
          "bowtie1_genome_cutadapt_topN_1mm_notidentical",
        ],
        summary => [ "miRNA_1mm_table", "tRNA_1mm_table", "smallRNA_1mm_table", "smallRNA_1mm_category", "miRNA_1mm_overlap_position", "tRNA_1mm_position" ],
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
}

1;

