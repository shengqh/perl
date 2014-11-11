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

my $samtools                   = "/home/shengq1/local/bin/samtools/samtools";
my $bowtie1_option_1mm         = "-a -m 100 --best --strata -v 1 -p 8";
my $mirnacount_option          = "-s";                                                    #ignore score
my $trnacount_option           = "--length --sequence";
my $mirna_overlap_count_option = "-s --gtf_key miRNA";
my $mirna_fasta                = "/data/cqs/shengq1/reference/miRBase20/mature.dna.fa";
my $bowtie1_option_pm          = "-a -m 100 --best --strata -v 0 -p 8";

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

my $preprocessing = {
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
    source     => { individual => [ "trimmer", "fastqc_pre", "cutadapt", "fastqc_post", "fastqlen", "identical", ], },
    sh_direct  => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
};

#performConfig($preprocessing);

my $human_def = {
  task_name           => $def->{task_name},
  mirna_coordinate    => "/data/cqs/shengq1/reference/miRBase20/hsa.gff3",
  trna_coordinate     => "/data/cqs/guoy1/reference/smallrna/hg19_tRNA_ucsc_ensembl.bed",
  trna_fasta          => "/data/cqs/guoy1/reference/smallrna/hg19_tRNA_ucsc_ensembl.bed.fa",
  smallrna_coordinate => "/data/cqs/guoy1/reference/smallrna/hg19_smallRNA_ucsc_ensembl.bed",
  bowtie1_index       => "/data/cqs/guoy1/reference/hg19/bowtie_index_hg19_rCRS_1.0.0/hg19_rCRS",
  target_dir          => $root . "human/",
  fastq_files         => {
    "2687-GG-20" => ["/scratch/cqs/shengq1/vangard/VANGARD00421_20141110_guoyan_smallRNA_2687/identical/result/2687-GG-20_clipped_identical.fastq.gz"],
    "2687-GG-23" => ["/scratch/cqs/shengq1/vangard/VANGARD00421_20141110_guoyan_smallRNA_2687/identical/result/2687-GG-23_clipped_identical.fastq.gz"],
    "2687-GG-61" => ["/scratch/cqs/shengq1/vangard/VANGARD00421_20141110_guoyan_smallRNA_2687/identical/result/2687-GG-61_clipped_identical.fastq.gz"],
    "2687-GG-65" => ["/scratch/cqs/shengq1/vangard/VANGARD00421_20141110_guoyan_smallRNA_2687/identical/result/2687-GG-65_clipped_identical.fastq.gz"],
  },
  count_files => {
    "2687-GG-20" => ["/scratch/cqs/shengq1/vangard/VANGARD00421_20141110_guoyan_smallRNA_2687/identical/result/2687-GG-20_clipped_identical.fastq.dupcount"],
    "2687-GG-23" => ["/scratch/cqs/shengq1/vangard/VANGARD00421_20141110_guoyan_smallRNA_2687/identical/result/2687-GG-23_clipped_identical.fastq.dupcount"],
    "2687-GG-61" => ["/scratch/cqs/shengq1/vangard/VANGARD00421_20141110_guoyan_smallRNA_2687/identical/result/2687-GG-61_clipped_identical.fastq.dupcount"],
    "2687-GG-65" => ["/scratch/cqs/shengq1/vangard/VANGARD00421_20141110_guoyan_smallRNA_2687/identical/result/2687-GG-65_clipped_identical.fastq.dupcount"],
  },
};

my $mouse_def = {
  task_name           => $def->{task_name},
  mirna_coordinate    => "/data/cqs/shengq1/reference/miRBase20/mmu.gff3",
  trna_coordinate     => "/data/cqs/guoy1/reference/smallrna/mm10_tRNA_ucsc_ensembl.bed",
  trna_fasta          => "/data/cqs/guoy1/reference/smallrna/mm10_tRNA_ucsc_ensembl.bed.fa",
  smallrna_coordinate => "/data/cqs/guoy1/reference/smallrna/mm10_smallRNA_ucsc_ensembl.bed",
  bowtie1_index       => "/data/cqs/shengq1/reference/mm10/bowtie_index/mm10",
  target_dir          => $root . "mouse/",
  fastq_files         => {
    "2687-GG-35" => ["/scratch/cqs/shengq1/vangard/VANGARD00421_20141110_guoyan_smallRNA_2687/identical/result/2687-GG-35_clipped_identical.fastq.gz"],
    "2687-GG-59" => ["/scratch/cqs/shengq1/vangard/VANGARD00421_20141110_guoyan_smallRNA_2687/identical/result/2687-GG-59_clipped_identical.fastq.gz"],
    "2687-GG-60" => ["/scratch/cqs/shengq1/vangard/VANGARD00421_20141110_guoyan_smallRNA_2687/identical/result/2687-GG-60_clipped_identical.fastq.gz"],
    "2687-GG-61" => ["/scratch/cqs/shengq1/vangard/VANGARD00421_20141110_guoyan_smallRNA_2687/identical/result/2687-GG-61_clipped_identical.fastq.gz"],
    "2687-GG-62" => ["/scratch/cqs/shengq1/vangard/VANGARD00421_20141110_guoyan_smallRNA_2687/identical/result/2687-GG-62_clipped_identical.fastq.gz"],
    "2687-GG-63" => ["/scratch/cqs/shengq1/vangard/VANGARD00421_20141110_guoyan_smallRNA_2687/identical/result/2687-GG-63_clipped_identical.fastq.gz"],
    "2687-GG-65" => ["/scratch/cqs/shengq1/vangard/VANGARD00421_20141110_guoyan_smallRNA_2687/identical/result/2687-GG-65_clipped_identical.fastq.gz"],
    "2687-GG-66" => ["/scratch/cqs/shengq1/vangard/VANGARD00421_20141110_guoyan_smallRNA_2687/identical/result/2687-GG-66_clipped_identical.fastq.gz"],
    "2687-GG-67" => ["/scratch/cqs/shengq1/vangard/VANGARD00421_20141110_guoyan_smallRNA_2687/identical/result/2687-GG-67_clipped_identical.fastq.gz"],
    "2687-GG-69" => ["/scratch/cqs/shengq1/vangard/VANGARD00421_20141110_guoyan_smallRNA_2687/identical/result/2687-GG-69_clipped_identical.fastq.gz"],
    "2687-GG-70" => ["/scratch/cqs/shengq1/vangard/VANGARD00421_20141110_guoyan_smallRNA_2687/identical/result/2687-GG-70_clipped_identical.fastq.gz"],
  },
  count_files => {
    "2687-GG-35" => ["/scratch/cqs/shengq1/vangard/VANGARD00421_20141110_guoyan_smallRNA_2687/identical/result/2687-GG-35_clipped_identical.fastq.dupcount"],
    "2687-GG-59" => ["/scratch/cqs/shengq1/vangard/VANGARD00421_20141110_guoyan_smallRNA_2687/identical/result/2687-GG-59_clipped_identical.fastq.dupcount"],
    "2687-GG-60" => ["/scratch/cqs/shengq1/vangard/VANGARD00421_20141110_guoyan_smallRNA_2687/identical/result/2687-GG-60_clipped_identical.fastq.dupcount"],
    "2687-GG-61" => ["/scratch/cqs/shengq1/vangard/VANGARD00421_20141110_guoyan_smallRNA_2687/identical/result/2687-GG-61_clipped_identical.fastq.dupcount"],
    "2687-GG-62" => ["/scratch/cqs/shengq1/vangard/VANGARD00421_20141110_guoyan_smallRNA_2687/identical/result/2687-GG-62_clipped_identical.fastq.dupcount"],
    "2687-GG-63" => ["/scratch/cqs/shengq1/vangard/VANGARD00421_20141110_guoyan_smallRNA_2687/identical/result/2687-GG-63_clipped_identical.fastq.dupcount"],
    "2687-GG-65" => ["/scratch/cqs/shengq1/vangard/VANGARD00421_20141110_guoyan_smallRNA_2687/identical/result/2687-GG-65_clipped_identical.fastq.dupcount"],
    "2687-GG-66" => ["/scratch/cqs/shengq1/vangard/VANGARD00421_20141110_guoyan_smallRNA_2687/identical/result/2687-GG-66_clipped_identical.fastq.dupcount"],
    "2687-GG-67" => ["/scratch/cqs/shengq1/vangard/VANGARD00421_20141110_guoyan_smallRNA_2687/identical/result/2687-GG-67_clipped_identical.fastq.dupcount"],
    "2687-GG-69" => ["/scratch/cqs/shengq1/vangard/VANGARD00421_20141110_guoyan_smallRNA_2687/identical/result/2687-GG-69_clipped_identical.fastq.dupcount"],
    "2687-GG-70" => ["/scratch/cqs/shengq1/vangard/VANGARD00421_20141110_guoyan_smallRNA_2687/identical/result/2687-GG-70_clipped_identical.fastq.dupcount"],
  },
};

my @defs = ( $human_def, $mouse_def );

foreach my $def (@defs) {
  my $target_dir = $def->{target_dir};
  my $config     = {
    general => { "task_name" => $def->{task_name} },

    #1 mismatch search
    bowtie1_genome_cutadapt_topN_1mm => {
      class         => "Bowtie1",
      perform       => 1,
      target_dir    => "${target_dir}/topN_bowtie1_genome_cutadapt_1mm",
      option        => $bowtie1_option_1mm,
      source        => $def->{fastq_files},
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
      class       => "MirnaCount",
      perform     => 1,
      target_dir  => "${target_dir}/topN_bowtie1_genome_cutadapt_1mm_count_miRNA",
      option      => $mirnacount_option,
      source_ref  => "bowtie1_genome_cutadapt_topN_1mm",
      fastq_files => $def->{fastq_files},
      seqcount    => $def->{count_files},
      cqs_tools   => $cqstools,
      gff_file    => $def->{mirna_coordinate},
      fasta_file  => $mirna_fasta,
      samtools    => $samtools,
      sh_direct   => 1,
      pbs         => {
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
      class       => "CQSMappedCount",
      perform     => 1,
      target_dir  => "${target_dir}/topN_bowtie1_genome_cutadapt_1mm_count_miRNA_overlap",
      option      => $mirna_overlap_count_option,
      source_ref  => "bowtie1_genome_cutadapt_topN_1mm",
      fastq_files => $def->{fastq_files},
      seqcount    => $def->{count_files},
      cqs_tools   => $cqstools,
      gff_file    => $def->{mirna_coordinate},
      fasta_file  => $mirna_fasta,
      samtools    => $samtools,
      sh_direct   => 1,
      pbs         => {
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
      class       => "CQSMappedCount",
      perform     => 1,
      target_dir  => "${target_dir}/topN_bowtie1_genome_cutadapt_1mm_count_tRNA",
      option      => $trnacount_option,
      source_ref  => "bowtie1_genome_cutadapt_topN_1mm",
      fastq_files => $def->{fastq_files},
      seqcount    => $def->{count_files},
      cqs_tools   => $cqstools,
      gff_file    => $def->{trna_coordinate},
      fasta_file  => $def->{trna_fasta},
      samtools    => $samtools,
      sh_direct   => 1,
      pbs         => {
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
      class       => "CQSMappedCount",
      perform     => 1,
      target_dir  => "${target_dir}/topN_bowtie1_genome_cutadapt_1mm_count_smallRNA",
      option      => $trnacount_option,
      source_ref  => "bowtie1_genome_cutadapt_topN_1mm",
      fastq_files => $def->{fastq_files},
      seqcount    => $def->{count_files},
      cqs_tools   => $cqstools,
      gff_file    => $def->{smallrna_coordinate},
      samtools    => $samtools,
      sh_direct   => 1,
      pbs         => {
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
        individual => [ "bowtie1_genome_cutadapt_topN_1mm", "mirna_1mm_count", "miRNA_1mm_count_overlap", "tRNA_1mm_count",        "smallRNA_1mm_count", ],
        summary    => [ "miRNA_1mm_table",                  "tRNA_1mm_table",  "smallRNA_1mm_table",      "smallRNA_1mm_category", "miRNA_1mm_overlap_position", "tRNA_1mm_position" ],
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

