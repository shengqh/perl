#!/usr/bin/perl
use strict;
use warnings;

use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;

my $root     = create_directory_or_die("/scratch/cqs/shengq1/smallRNA/20140926_bingshan_smallRNA_human/");
my $cqstools = "/home/shengq1/cqstools/CQS.Tools.exe";

my $email = "quanhu.sheng\@vanderbilt.edu";

my $samtools           = "/home/shengq1/local/bin/samtools/samtools";
my $bowtie1_option_1mm = "-a -m 100 --best --strata -v 1 -l 12 -p 8";

my $mirnacount_option          = "-s";                                                    #ignore score
my $trnacount_option           = "--length --sequence";
my $mirna_overlap_count_option = "-s --gtf_key miRNA";
my $mirna_fasta                = "/data/cqs/shengq1/reference/miRBase20/mature.dna.fa";

my $bowtie1_option_pm = "-a -m 100 --best --strata -v 0 -l 12 -p 8";

my $def = {
  task_name           => "3018_b2",
  mirna_coordinate    => "/data/cqs/shengq1/reference/miRBase20/hsa.gff3",
  trna_coordinate     => "/data/cqs/guoy1/reference/smallrna/hg19_tRNA_ucsc_ensembl.bed",
  trna_fasta          => "/data/cqs/guoy1/reference/smallrna/hg19_tRNA_ucsc_ensembl.bed.fa",
  smallrna_coordinate => "/data/cqs/guoy1/reference/smallrna/hg19_smallRNA_ucsc_ensembl.bed",
  bowtie1_index       => "/data/cqs/guoy1/reference/hg19/bowtie_index_hg19_rCRS_1.0.0/hg19_rCRS",
  target_dir          => $root,
  files               => {
    "2812-BL-01"   => ["/gpfs21/scratch/shengq1/bingshan/batch1_12samples/fastq/2812-BL-1_1_sequence.txt.gz"],
    "2812-BL-10"   => ["/gpfs21/scratch/shengq1/bingshan/batch1_12samples/fastq/2812-BL-10_1_sequence.txt.gz"],
    "2812-BL-11"   => ["/gpfs21/scratch/shengq1/bingshan/batch1_12samples/fastq/2812-BL-11_1_sequence.txt.gz"],
    "2812-BL-12"   => ["/gpfs21/scratch/shengq1/bingshan/batch1_12samples/fastq/2812-BL-12_1_sequence.txt.gz"],
    "2812-BL-02"   => ["/gpfs21/scratch/shengq1/bingshan/batch1_12samples/fastq/2812-BL-2_1_sequence.txt.gz"],
    "2812-BL-03"   => ["/gpfs21/scratch/shengq1/bingshan/batch1_12samples/fastq/2812-BL-3_1_sequence.txt.gz"],
    "2812-BL-04"   => ["/gpfs21/scratch/shengq1/bingshan/batch1_12samples/fastq/2812-BL-4_1_sequence.txt.gz"],
    "2812-BL-05"   => ["/gpfs21/scratch/shengq1/bingshan/batch1_12samples/fastq/2812-BL-5_1_sequence.txt.gz"],
    "2812-BL-06"   => ["/gpfs21/scratch/shengq1/bingshan/batch1_12samples/fastq/2812-BL-6_1_sequence.txt.gz"],
    "2812-BL-07"   => ["/gpfs21/scratch/shengq1/bingshan/batch1_12samples/fastq/2812-BL-7_1_sequence.txt.gz"],
    "2812-BL-08"   => ["/gpfs21/scratch/shengq1/bingshan/batch1_12samples/fastq/2812-BL-8_1_sequence.txt.gz"],
    "2812-BL-09"   => ["/gpfs21/scratch/shengq1/bingshan/batch1_12samples/fastq/2812-BL-9_1_sequence.txt.gz"],
    "3005-BL-01"   => ["/gpfs21/scratch/shengq1/bingshan/batch2_30samples/fastq/3005-BL-1_1_sequence.txt.gz"],
    "3005-BL-10"   => ["/gpfs21/scratch/shengq1/bingshan/batch2_30samples/fastq/3005-BL-10_1_sequence.txt.gz"],
    "3005-BL-11"   => ["/gpfs21/scratch/shengq1/bingshan/batch2_30samples/fastq/3005-BL-11_1_sequence.txt.gz"],
    "3005-BL-12"   => ["/gpfs21/scratch/shengq1/bingshan/batch2_30samples/fastq/3005-BL-12_1_sequence.txt.gz"],
    "3005-BL-13"   => ["/gpfs21/scratch/shengq1/bingshan/batch2_30samples/fastq/3005-BL-13_1_sequence.txt.gz"],
    "3005-BL-14"   => ["/gpfs21/scratch/shengq1/bingshan/batch2_30samples/fastq/3005-BL-14_1_sequence.txt.gz"],
    "3005-BL-15"   => ["/gpfs21/scratch/shengq1/bingshan/batch2_30samples/fastq/3005-BL-15_1_sequence.txt.gz"],
    "3005-BL-16"   => ["/gpfs21/scratch/shengq1/bingshan/batch2_30samples/fastq/3005-BL-16_1_sequence.txt.gz"],
    "3005-BL-17"   => ["/gpfs21/scratch/shengq1/bingshan/batch2_30samples/fastq/3005-BL-17_1_sequence.txt.gz"],
    "3005-BL-18"   => ["/gpfs21/scratch/shengq1/bingshan/batch2_30samples/fastq/3005-BL-18_1_sequence.txt.gz"],
    "3005-BL-19"   => ["/gpfs21/scratch/shengq1/bingshan/batch2_30samples/fastq/3005-BL-19_1_sequence.txt.gz"],
    "3005-BL-02"   => ["/gpfs21/scratch/shengq1/bingshan/batch2_30samples/fastq/3005-BL-2_1_sequence.txt.gz"],
    "3005-BL-20"   => ["/gpfs21/scratch/shengq1/bingshan/batch2_30samples/fastq/3005-BL-20_1_sequence.txt.gz"],
    "3005-BL-21"   => ["/gpfs21/scratch/shengq1/bingshan/batch2_30samples/fastq/3005-BL-21_1_sequence.txt.gz"],
    "3005-BL-22"   => ["/gpfs21/scratch/shengq1/bingshan/batch2_30samples/fastq/3005-BL-22_1_sequence.txt.gz"],
    "3005-BL-23"   => ["/gpfs21/scratch/shengq1/bingshan/batch2_30samples/fastq/3005-BL-23_1_sequence.txt.gz"],
    "3005-BL-24"   => ["/gpfs21/scratch/shengq1/bingshan/batch2_30samples/fastq/3005-BL-24_1_sequence.txt.gz"],
    "3005-BL-25"   => ["/gpfs21/scratch/shengq1/bingshan/batch2_30samples/fastq/3005-BL-25_1_sequence.txt.gz"],
    "3005-BL-26"   => ["/gpfs21/scratch/shengq1/bingshan/batch2_30samples/fastq/3005-BL-26_1_sequence.txt.gz"],
    "3005-BL-27"   => ["/gpfs21/scratch/shengq1/bingshan/batch2_30samples/fastq/3005-BL-27_1_sequence.txt.gz"],
    "3005-BL-28"   => ["/gpfs21/scratch/shengq1/bingshan/batch2_30samples/fastq/3005-BL-28_1_sequence.txt.gz"],
    "3005-BL-29"   => ["/gpfs21/scratch/shengq1/bingshan/batch2_30samples/fastq/3005-BL-29_1_sequence.txt.gz"],
    "3005-BL-03"   => ["/gpfs21/scratch/shengq1/bingshan/batch2_30samples/fastq/3005-BL-3_1_sequence.txt.gz"],
    "3005-BL-30"   => ["/gpfs21/scratch/shengq1/bingshan/batch2_30samples/fastq/3005-BL-30_1_sequence.txt.gz"],
    "3005-BL-04"   => ["/gpfs21/scratch/shengq1/bingshan/batch2_30samples/fastq/3005-BL-4_1_sequence.txt.gz"],
    "3005-BL-05"   => ["/gpfs21/scratch/shengq1/bingshan/batch2_30samples/fastq/3005-BL-5_1_sequence.txt.gz"],
    "3005-BL-06"   => ["/gpfs21/scratch/shengq1/bingshan/batch2_30samples/fastq/3005-BL-6_1_sequence.txt.gz"],
    "3005-BL-07"   => ["/gpfs21/scratch/shengq1/bingshan/batch2_30samples/fastq/3005-BL-7_1_sequence.txt.gz"],
    "3005-BL-08"   => ["/gpfs21/scratch/shengq1/bingshan/batch2_30samples/fastq/3005-BL-8_1_sequence.txt.gz"],
    "3005-BL-09"   => ["/gpfs21/scratch/shengq1/bingshan/batch2_30samples/fastq/3005-BL-9_1_sequence.txt.gz"],
    "3034-BL-1-01" => ["/gpfs21/scratch/shengq1/bingshan/batch3_10samples_pooledByHushan/fastq/3034-BL-1-1_ATCACG_L008_R1_001.fastq.gz"],
    "3034-BL-1-10" => ["/gpfs21/scratch/shengq1/bingshan/batch3_10samples_pooledByHushan/fastq/3034-BL-1-10_TAGCTT_L008_R1_001.fastq.gz"],
    "3034-BL-1-02" => ["/gpfs21/scratch/shengq1/bingshan/batch3_10samples_pooledByHushan/fastq/3034-BL-1-2_CGATGT_L008_R1_001.fastq.gz"],
    "3034-BL-1-03" => ["/gpfs21/scratch/shengq1/bingshan/batch3_10samples_pooledByHushan/fastq/3034-BL-1-3_TTAGGC_L008_R1_001.fastq.gz"],
    "3034-BL-1-04" => ["/gpfs21/scratch/shengq1/bingshan/batch3_10samples_pooledByHushan/fastq/3034-BL-1-4_TGACCA_L008_R1_001.fastq.gz"],
    "3034-BL-1-05" => ["/gpfs21/scratch/shengq1/bingshan/batch3_10samples_pooledByHushan/fastq/3034-BL-1-5_ACAGTG_L008_R1_001.fastq.gz"],
    "3034-BL-1-06" => ["/gpfs21/scratch/shengq1/bingshan/batch3_10samples_pooledByHushan/fastq/3034-BL-1-6_GCCAAT_L008_R1_001.fastq.gz"],
    "3034-BL-1-07" => ["/gpfs21/scratch/shengq1/bingshan/batch3_10samples_pooledByHushan/fastq/3034-BL-1-7_CAGATC_L008_R1_001.fastq.gz"],
    "3034-BL-1-08" => ["/gpfs21/scratch/shengq1/bingshan/batch3_10samples_pooledByHushan/fastq/3034-BL-1-8_ACTTGA_L008_R1_001.fastq.gz"],
    "3034-BL-1-09" => ["/gpfs21/scratch/shengq1/bingshan/batch3_10samples_pooledByHushan/fastq/3034-BL-1-9_GATCAG_L008_R1_001.fastq.gz"],
  },
};

my $target_dir = $def->{target_dir};
my $config     = {
  general => { "task_name" => $def->{task_name}, },
  files   => $def->{files},
  fastqc  => {
    class      => "QC::FastQC",
    perform    => 1,
    target_dir => "${target_dir}/fastqc",
    option     => "",
    source_ref => "files",
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
    source_ref => "files",
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
  fastqc_post => {
    class      => "QC::FastQC",
    perform    => 1,
    target_dir => "${target_dir}/fastqc_posttrim",
    option     => "",
    source_ref => [ "cutadapt", ".fastq.gz" ],
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "2",
      "mem"      => "10gb"
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
    source_ref    => [ "cutadapt", ".fastq.gz\$" ],
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
    perform    => 1,
    target_dir => "${target_dir}/topN_bowtie1_genome_cutadapt_1mm_count_miRNA_overlap_position",
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
    perform         => 1,
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
    perform         => 1,
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

  #2 perfect match search to mirbase only
  bowtie1_genome_cutadapt_topN_genome_pmnames => {
    class      => "Samtools::PerfectMappedReadNames",
    perform    => 1,
    target_dir => "${target_dir}/topN_bowtie1_genome_cutadapt_pmnames",
    option     => "",
    source_ref => "bowtie1_genome_cutadapt_topN_1mm",
    sh_direct  => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  bowtie1_genome_cutadapt_topN_miRbase_pm => {
    class         => "Alignment::Bowtie1",
    perform       => 1,
    target_dir    => "${target_dir}/topN_bowtie1_genome_cutadapt_miRbase_pm",
    option        => $bowtie1_option_pm,
    source_ref    => [ "identical", ".fastq.gz\$" ],
    bowtie1_index => "/data/cqs/shengq1/reference/miRBase21/bowtie_index_1.0.1/mature.dna",
    samonly       => 0,
    sh_direct     => 1,
    pbs           => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  chromosome_count => {
    class                   => "CQS::CQSChromosomeCount",
    perform                 => 1,
    target_dir              => "${target_dir}/topN_bowtie1_genome_cutadapt_miRbase_pm_count",
    option                  => "",
    source_ref              => "bowtie1_genome_cutadapt_topN_miRbase_pm",
    seqcount_ref            => [ "identical", ".dupcount\$" ],
    perfect_mapped_name_ref => "bowtie1_genome_cutadapt_topN_genome_pmnames",
    cqs_tools               => $cqstools,
    samtools                => $samtools,
    sh_direct               => 1,
    pbs                     => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  chromosome_count_table => {
    class      => "CQSChromosomeTable",
    perform    => 1,
    target_dir => "${target_dir}/topN_bowtie1_genome_cutadapt_miRbase_pm_table",
    option     => "",
    source_ref => [ "chromosome_count", ".xml" ],
    cqs_tools  => $cqstools,
    prefix     => "miRBase_pm_",
    sh_direct  => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "10",
      "mem"      => "10gb"
    },
  },
  sequencetask => {
    class      => "CQS::SequenceTask",
    perform    => 1,
    target_dir => "${target_dir}/sequencetask",
    option     => "",
    source     => {
      individual => [
        "fastqc", "cutadapt", "fastqc_post", "fastqlen", "identical", "bowtie1_genome_cutadapt_topN_1mm",
        "mirna_1mm_count", "miRNA_1mm_count_overlap", "tRNA_1mm_count", "smallRNA_1mm_count",
        "bowtie1_genome_cutadapt_topN_genome_pmnames",
        "bowtie1_genome_cutadapt_topN_miRbase_pm",
        "chromosome_count", "bowtie1_genome_cutadapt_topN_1mm_notidentical",
      ],
      summary => [ "miRNA_1mm_table", "tRNA_1mm_table", "smallRNA_1mm_table", "smallRNA_1mm_category", "miRNA_1mm_overlap_position", "tRNA_1mm_position", "chromosome_count_table" ],
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

#performTask($config, "fastqc_post");

1;

