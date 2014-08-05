#!/usr/bin/perl
use strict;
use warnings;

use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;

my $root     = create_directory_or_die("/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/");
my $cqstools = "/home/shengq1/cqstools/CQS.Tools.exe";

my $hg19_mrna_gff      = "/data/cqs/shengq1/reference/miRBase20/hsa.gff3";
my $hg19_trna_bed      = "/data/cqs/guoy1/reference/smallrna/hg19_tRNA_ucsc_ensembl.bed";
my $hg19_trna_fasta    = "/data/cqs/guoy1/reference/smallrna/hg19_tRNA_ucsc_ensembl.bed.fa";
my $hg19_smallrna_bed  = "/data/cqs/guoy1/reference/smallrna/hg19_smallRNA_ucsc_ensembl.bed";
my $hg19_bowtie1_index = "/data/cqs/guoy1/reference/hg19/bowtie_index_hg19_rCRS_1.0.0/hg19_rCRS";

my $email = "quanhu.sheng\@vanderbilt.edu";

my $samtools           = "/home/shengq1/local/bin/samtools/samtools";
my $bowtie1_option_1mm = "-a -m 100 --best --strata -v 1 -l 12 -p 8";

my $mirnacount_option          = "-s --no_unmapped_fastq";                                #ignore score
my $trnacount_option           = "--length --sequence --no_unmapped_fastq";
my $mirna_overlap_count_option = "-s --gtf_key miRNA --no_unmapped_fastq";
my $mirna_fasta                = "/data/cqs/shengq1/reference/miRBase20/mature.dna.fa";

my $bowtie1_option_pm = "-a -m 100 --best --strata -v 0 -l 12 -p 8";

my $def = {
  task_name           => "Human_HDL_N10",
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
    "01-18-Post" => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/data/01-018-Post_CTTGTA_L005_R1_001.fastq.gz"],
    "01-18-Pre"  => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/data/01-018-Pre_GGCTAC_L005_R1_001.fastq.gz"],
    "01-28-Post" => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/data/01-28-Post_CGATGT_L004_R1_001.fastq.gz"],
    "01-28-Pre"  => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/data/01-28-Pre_ATCACG_L004_R1_001.fastq.gz"],
    "01-29-Post" => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/data/01-29-Post_TGACCA_L005_R1_001.fastq.gz"],
    "01-29-Pre"  => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/data/01-29-Pre_TTAGGC_L005_R1_001.fastq.gz"],
    "01-31-Post" => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/data/01-031-Post_GGTAGC_L005_R1_001.fastq.gz"],
    "01-31-Pre"  => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/data/01-031-Pre_GAGTGG_L005_R1_001.fastq.gz"],
    "01-36-Post" => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/data/01-36-Post_GCCAAT_L004_R1_001.fastq.gz"],
    "01-36-Pre"  => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/data/01-36-Pre_ACAGTG_L004_R1_001.fastq.gz"],
    "01-61-Post" => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/data/01-061-Post_ATGAGC_L004_R1_001.fastq.gz"],
    "01-61-Pre"  => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/data/01-061-Pre_ACTGAT_L004_R1_001.fastq.gz"],
    "03-07-Post" => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/data/03-007-Post_CAAAAG_L005_R1_001.fastq.gz"],
    "03-07-Pre"  => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/data/03-007-Pre_ATTCCT_L005_R1_001.fastq.gz"],
    "03-11-Post" => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/data/03-011-Post_CACCGG_L004_R1_001.fastq.gz"],
    "03-11-Pre"  => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/data/03-011-Pre_CAACTA_L004_R1_001.fastq.gz"],
    "03-15-Post" => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/data/03-015-Post_CACTCA_L005_R1_001.fastq.gz"],
    "03-15-Pre"  => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/data/03-015-Pre_CACGAT_L005_R1_001.fastq.gz"],
    "03-16-Post" => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/data/03-16-Post_ACTTGA_L005_R1_001.fastq.gz"],
    "03-16-Pre"  => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/data/03-16-Pre_CAGATC_L005_R1_001.fastq.gz"],
    "03-17-Post" => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/data/03-17-Post_TAGCTT_L004_R1_001.fastq.gz"],
    "03-17-Pre"  => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/data/03-17-Pre_GATCAG_L004_R1_001.fastq.gz"],
    "03-18-Post" => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/data/03-018-Post_CGTACG_L004_R1_001.fastq.gz"],
    "03-18-Pre"  => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/data/03-018-Pre_GTTTCG_L004_R1_001.fastq.gz"],
    "03-26-Post" => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/data/03-026-Post_AGTTCC_L004_R1_001.fastq.gz"],
    "03-26-Pre"  => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/data/03-026-Pre_AGTCAA_L004_R1_001.fastq.gz"],
    "03-31-Post" => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/data/03-031-Post_CATGGC_L004_R1_001.fastq.gz"],
    "03-31-Pre"  => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/data/03-031-Pre_CAGGCG_L004_R1_001.fastq.gz"],
    "03-33-Post" => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/data/03-033-Post_CCAACA_L005_R1_001.fastq.gz"],
    "03-33-Pre"  => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/data/03-033-Pre_CATTTT_L005_R1_001.fastq.gz"],
    "03-36-Post" => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/data/03-036-Post_CCGTCC_L005_R1_001.fastq.gz"],
    "03-36-Pre"  => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/data/03-036-Pre_ATGTCA_L005_R1_001.fastq.gz"],
    "03-47-Post" => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/data/03-047-Post_CTAGCT_L004_R1_001.fastq.gz"],
    "03-47-Pre"  => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/data/03-047-Pre_CGGAAT_L004_R1_001.fastq.gz"],
    "03-49-Post" => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/data/03-049-Post_CTCAGA_L005_R1_001.fastq.gz"],
    "03-49-Pre"  => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/data/03-049-Pre_CTATAC_L005_R1_001.fastq.gz"],
    "03-63-Post" => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/data/03-063-Post_GTGGCC_L005_R1_001.fastq.gz"],
    "03-63-Pre"  => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/data/03-063-Pre_GTGAAA_L005_R1_001.fastq.gz"],
    "03-65-Post" => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/data/03-065-Post_GTCCGC_L004_R1_001.fastq.gz"],
    "03-65-Pre"  => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/data/03-065-Pre_GTAGAG_L004_R1_001.fastq.gz"],
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

  #2 perfect match search
  bowtie1_genome_cutadapt_topN_pm => {
    class         => "Bowtie1",
    perform       => 1,
    target_dir    => "${target_dir}/topN_bowtie1_genome_cutadapt_miRbase_pm",
    option        => $bowtie1_option_pm,
    source_ref    => [ "identical", ".fastq\$" ],
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

  sequencetask => {
    class      => "CQS::SequenceTask",
    perform    => 1,
    target_dir => "${target_dir}/sequencetask",
    option     => "",
    source     => {
      individual => [
        "cutadapt", "fastqlen", "identical",
        "bowtie1_genome_cutadapt_topN_1mm_notidentical",
        "bowtie1_genome_cutadapt_topN_1mm",
        "mirna_1mm_count", "miRNA_1mm_count_overlap", "tRNA_1mm_count", "smallRNA_1mm_count",
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

#performConfig($config);
performTask( $config, "bowtie1_genome_cutadapt_topN_pm" );

1;
