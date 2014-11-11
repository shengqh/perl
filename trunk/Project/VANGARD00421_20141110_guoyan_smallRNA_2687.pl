#!/usr/bin/perl
use strict;
use warnings;

use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;

my $root     = create_directory_or_die("/scratch/cqs/shengq1/vangard/VANGARD00421_20141110_guoyan_smallRNA_2687/");
my $cqstools = "/home/shengq1/cqstools/CQS.Tools.exe";

my $email = "quanhu.sheng\@vanderbilt.edu";

my $samtools                = "/home/shengq1/local/bin/samtools/samtools";
my $bowtie1_option_1mm      = "-a -m 100 --best --strata -v 1 -p 8";
my $mirnacount_option          = "-s";                                                    #ignore score
my $trnacount_option           = "--length --sequence";
my $mirna_overlap_count_option = "-s --gtf_key miRNA";
my $mirna_fasta                = "/data/cqs/shengq1/reference/miRBase20/mature.dna.fa";
my $bowtie1_option_pm = "-a -m 100 --best --strata -v 0 -p 8";

my $def = {
  task_name           => "2687",
  mirna_coordinate    => "/data/cqs/shengq1/reference/miRBase20/hsa.gff3",
  trna_coordinate     => "/data/cqs/guoy1/reference/smallrna/hg19_tRNA_ucsc_ensembl.bed",
  trna_fasta          => "/data/cqs/guoy1/reference/smallrna/hg19_tRNA_ucsc_ensembl.bed.fa",
  smallrna_coordinate => "/data/cqs/guoy1/reference/smallrna/hg19_smallRNA_ucsc_ensembl.bed",
  bowtie1_index       => "/data/cqs/guoy1/reference/hg19/bowtie_index_hg19_rCRS_1.0.0/hg19_rCRS",
  target_dir          => $root,
  files               => {
    "2687-GG-20" => ["/scratch/cqs/guoy1/2687/2687-GG-20_1_sequence.txt.gz"],
    "2687-GG-23" => ["/scratch/cqs/guoy1/2687/2687-GG-23_1_sequence.txt.gz"],
    "2687-GG-35" => ["/scratch/cqs/guoy1/2687/2687-GG-35_1_sequence.txt.gz"],
    "2687-GG-59" => ["/scratch/cqs/guoy1/2687/2687-GG-59_1_sequence.txt.gz"],
    "2687-GG-60" => ["/scratch/cqs/guoy1/2687/2687-GG-60_1_sequence.txt.gz"],
    "2687-GG-61" => ["/scratch/cqs/guoy1/2687/2687-GG-61_1_sequence.txt.gz"],
    "2687-GG-62" => ["/scratch/cqs/guoy1/2687/2687-GG-62_1_sequence.txt.gz"],
    "2687-GG-63" => ["/scratch/cqs/guoy1/2687/2687-GG-63_1_sequence.txt.gz"],
    "2687-GG-65" => ["/scratch/cqs/guoy1/2687/2687-GG-65_1_sequence.txt.gz"],
    "2687-GG-66" => ["/scratch/cqs/guoy1/2687/2687-GG-66_1_sequence.txt.gz"],
    "2687-GG-67" => ["/scratch/cqs/guoy1/2687/2687-GG-67_1_sequence.txt.gz"],
    "2687-GG-69" => ["/scratch/cqs/guoy1/2687/2687-GG-69_1_sequence.txt.gz"],
    "2687-GG-70" => ["/scratch/cqs/guoy1/2687/2687-GG-70_1_sequence.txt.gz"],
  },
};

my $target_dir = $def->{target_dir};

my $preprocessing     = {
  general => { "task_name" => $def->{task_name}, },
  files   => $def->{files},
  trimmer => {
    class      => "CQS::FastqTrimmer",
    perform    => 1,
    target_dir => "${target_dir}/fastq_trimN",
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
  fastqc_pre => {
    class      => "QC::FastQC",
    perform    => 1,
    target_dir => "${target_dir}/fastqc_pretrim",
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
    option     => "-O 10 -m 15",
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
  fastq_mirna => {
    class        => "CQS::FastqMirna",
    perform      => 1,
    target_dir   => "${target_dir}/fastq_mirna",
    option       => "",
    source_ref   => [ "identical", ".fastq.gz\$" ],
    seqcount_ref => [ "identical", ".dupcount\$" ],
    cqstools     => $cqstools,
    extension    => "_clipped_identical_mirna.fastq.gz",
    sh_direct    => 1,
    pbs          => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "24",
      "mem"      => "20gb"
    },
  },
  sequencetask => {
    class      => "CQS::SequenceTask",
    perform    => 1,
    target_dir => "${target_dir}/preprocessing_sequencetask",
    option     => "",
    source     => {
      individual => [
        "trimmer", "fastqc_pre", "cutadapt", "fastqc_post", "fastqlen", "identical",
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

performConfig($preprocessing);

my $human     = {
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
      "nodes"    => "1:ppn=6",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },

  #not identical, for IGV
  bowtie1_genome_cutadapt_topN_1mm_notidentical_sRNAbench => {
    class         => "Bowtie1",
    perform       => 1,
    target_dir    => "${target_dir}/bowtie1_genome_cutadapt_topN_1mm_notidentical_sRNAbench_2",
    option        => "--chunkmbs 2000 -p 8 -n 1 -l 19  -a -m 40  --best --strata",
    source_ref    => [ "cutadapt", ".fastq.gz\$" ],
    bowtie1_index => $def->{bowtie1_index},
    samonly       => 0,
    sh_direct     => 0,
    pbs           => {
      "email"    => $email,
      "nodes"    => "1:ppn=6",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },

  #1 mismatch search
  bowtie1_genome_cutadapt_topN_1mm_mirna => {
    class         => "Bowtie1",
    perform       => 1,
    target_dir    => "${target_dir}/topN_bowtie1_genome_cutadapt_1mm_mirna",
    option        => $bowtie1_option_1mm,
    source_ref    => [ "fastq_mirna", ".fastq.gz\$" ],
    bowtie1_index => $def->{bowtie1_index},
    samonly       => 0,
    sh_direct     => 1,
    pbs           => {
      "email"    => $email,
      "nodes"    => "1:ppn=6",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  mirna_1mm_count_NTA => {
    class           => "MirnaCount",
    perform         => 1,
    target_dir      => "${target_dir}/topN_bowtie1_genome_cutadapt_1mm_mirna_count_miRNA",
    option          => $mirnacount_option,
    source_ref      => "bowtie1_genome_cutadapt_topN_1mm_mirna",
    fastq_files_ref => "fastq_mirna",
    seqcount_ref    => [ "fastq_mirna", ".dupcount\$" ],
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
  miRNA_1mm_NTA_table => {
    class      => "CQSMirnaNTATable",
    perform    => 1,
    target_dir => "${target_dir}/topN_bowtie1_genome_cutadapt_1mm_mirna_count_miRNA_table",
    option     => "",
    source_ref => [ "mirna_1mm_count_NTA", ".mapped.xml" ],
    cqs_tools  => $cqstools,
    prefix     => "miRNA_1mm_NTA_",
    sh_direct  => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "10",
      "mem"      => "10gb"
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
      "nodes"    => "1:ppn=6",
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

  smallRNA_1mm_NTA_category => {
    class           => "CQSSmallRNACategory",
    perform         => 1,
    target_dir      => "${target_dir}/topN_bowtie1_genome_cutadapt_1mm_count_smallRNA_category_NTA",
    option          => "",
    source_ref      => [ "smallRNA_1mm_count", ".mapped.xml\$" ],
    mirna_count_ref => [ "mirna_1mm_count_NTA", ".mapped.xml\$" ],
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
      "nodes"    => "1:ppn=6",
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
      "nodes"    => "1:ppn=6",
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
      "nodes"    => "1:ppn=6",
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
        "trimmer", "fastqc_pre", "cutadapt", "fastqc_post", "fastqlen", "identical", "bowtie1_genome_cutadapt_topN_1mm",
        "mirna_1mm_count", "miRNA_1mm_count_overlap", "tRNA_1mm_count", "smallRNA_1mm_count",
        "bowtie1_genome_cutadapt_topN_genome_pmnames",
        "bowtie1_genome_cutadapt_topN_miRbase_pm",
        "chromosome_count", "bowtie1_genome_cutadapt_topN_1mm_notidentical",
      ],
      summary => [ "miRNA_1mm_table", "tRNA_1mm_table", "smallRNA_1mm_table", "smallRNA_1mm_category", "miRNA_1mm_overlap_position", "tRNA_1mm_position", "chromosome_count_table" ],
    },
    sh_direct => 1,
    pbs       => {
      "email"    => $email,
      "nodes"    => "1:ppn=6",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
};

1;

