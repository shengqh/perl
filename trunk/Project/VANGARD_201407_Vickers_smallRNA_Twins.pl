#!/usr/bin/perl
use strict;
use warnings;

use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;

my $root     = create_directory_or_die("/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201407_smallRNA_Twins/");
my $cqstools = "/home/shengq1/cqstools/CQS.Tools.exe";

my $hg19_mrna_gff      = "/data/cqs/shengq1/reference/miRBase20/hsa.gff3";
my $hg19_trna_bed      = "/data/cqs/guoy1/reference/smallrna/hg19_tRNA_ucsc_ensembl.bed";
my $hg19_trna_fasta    = "/data/cqs/guoy1/reference/smallrna/hg19_tRNA_ucsc_ensembl.bed.fa";
my $hg19_smallrna_bed  = "/data/cqs/guoy1/reference/smallrna/hg19_smallRNA_ucsc_ensembl.bed";
my $hg19_bowtie1_index = "/data/cqs/guoy1/reference/hg19/bowtie_index_hg19_rCRS_1.0.0/hg19_rCRS";

my $email = "quanhu.sheng\@vanderbilt.edu";

my $samtools           = "/home/shengq1/local/bin/samtools/samtools";
my $bowtie1_option_1mm = "-a -m 100 --best --strata -v 1 -l 12 -p 8";

my $mirnacount_option          = "-s";                                                    #ignore score
my $trnacount_option           = "--length --sequence";
my $mirna_overlap_count_option = "-s --gtf_key miRNA";
my $mirna_fasta                = "/data/cqs/shengq1/reference/miRBase20/mature.dna.fa";

my $def = {
  task_name           => "2501",
  mirna_coordinate    => $hg19_mrna_gff,
  trna_coordinate     => $hg19_trna_bed,
  trna_fasta          => $hg19_trna_fasta,
  smallrna_coordinate => $hg19_smallrna_bed,
  bowtie1_index       => $hg19_bowtie1_index,
  target_dir          => $root,
};

my $target_dir = $def->{target_dir};
my $config     = {
  general => { "task_name" => $def->{task_name}, },
  files   => {
    "DB113_5036" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201407_smallRNA_Twins/raw/Combined_DB113_2_CGTGAT_3.fastq.gz"],
    "DB114_5165" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201407_smallRNA_Twins/raw/Combined_DB114_2_ACATCG_3.fastq.gz"],
    "DB115_5166" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201407_smallRNA_Twins/raw/Combined_DB115_2_GCCTAA_3.fastq.gz"],
    "DB116_5181" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201407_smallRNA_Twins/raw/Combined_DB116_2_TGGTCA_3.fastq.gz"],
    "DB117_5182" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201407_smallRNA_Twins/raw/Combined_DB117_2_CACTGT_3.fastq.gz"],
    "DB118_5215" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201407_smallRNA_Twins/raw/Combined_DB118_2_ATTGGC_3.fastq.gz"],
    "DB119_5216" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201407_smallRNA_Twins/raw/Combined_DB119_2_GATCTG_3.fastq.gz"],
    "DB120_5237" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201407_smallRNA_Twins/raw/Combined_DB120_2_TCAAGT_3.fastq.gz"],
    "DB121_5238" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201407_smallRNA_Twins/raw/Combined_DB121_2_CTGATC_3.fastq.gz"],
    "DB123_9999" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201407_smallRNA_Twins/raw/Combined_DB123_2_GTAGCC_3.fastq.gz"],
  },
  groups => {
    "FIRST"  => [ "DB115_5166", "DB117_5182", "DB119_5216", "DB121_5238", "DB113_5036" ],
    "SECOND" => [ "DB114_5165", "DB116_5181", "DB118_5215", "DB120_5237", "DB123_9999" ],
  },
  paired => {
    "DEATH_ORDER" => {
      groups => [ "FIRST", "SECOND" ],
      paired => [ "P2083", "P2091", "P2108", "P2119", "P2018" ],
    },
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
  fastqc => {
    class      => "CQS::FastQC",
    perform    => 1,
    target_dir => "${target_dir}/fastqc",
    option     => "",
    source_ref => "trimmer",
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
    option     =>
"-q 3 -O 10 -m 12 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATGTATCTCGTATGCCGTCTTCTGCTTG -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACTTAGGCATCTCGTATGCCGTCTTCTGCTTG -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACTGACCAATCTCGTATGCCGTCTTCTGCTTG -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACGCCAATATCTCGTATGCCGTCTTCTGCTTG -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACCAGATCATCTCGTATGCCGTCTTCTGCTTG -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACACTTGAATCTCGTATGCCGTCTTCTGCTTG -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACGATCAGATCTCGTATGCCGTCTTCTGCTTG -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACTAGCTTATCTCGTATGCCGTCTTCTGCTTG -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACGGCTACATCTCGTATGCCGTCTTCTGCTTG -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACCTTGTAATCTCGTATGCCGTCTTCTGCTTG",
    source_ref => "trimmer",
    adaptor    => "GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG",
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
  miRNA_deseq2 => {
    class      => "DESeq2",
    perform    => 1,
    target_dir => "${target_dir}/topN_bowtie1_genome_cutadapt_1mm_count_miRNA_table_deseq2",
    option     => "",
    source_ref => "pairs",
    groups_ref => "groups",
    countfile_ref  => "miRNA_1mm_table",
    sh_direct => 1,
    pbs       => {
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
        "trimmer", "fastqc", "cutadapt", "fastqlen", "identical", "bowtie1_genome_cutadapt_topN_1mm",
        "mirna_1mm_count", "miRNA_1mm_count_overlap", "tRNA_1mm_count", "smallRNA_1mm_count", "bowtie1_genome_cutadapt_topN_1mm_notidentical",
      ],
      summary => [ "miRNA_1mm_table", "miRNA_deseq2", "tRNA_1mm_table", "smallRNA_1mm_table", "smallRNA_1mm_category", "miRNA_1mm_overlap_position", "tRNA_1mm_position" ],
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
performTask( $config, "miRNA_deseq2" );

1;

