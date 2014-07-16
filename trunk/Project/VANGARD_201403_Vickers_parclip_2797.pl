#!/usr/bin/perl
use strict;
use warnings;

use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;

my $root        = create_directory_or_die("/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201403_parclip_2797");
my $cqstools    = "/home/shengq1/cqstools/CQS.Tools.exe";
my $samtools    = "/home/shengq1/local/bin/samtools/samtools";
my $mirna_fasta = "/data/cqs/shengq1/reference/miRBase20/mature.dna.fa";
my $email       = "quanhu.sheng\@vanderbilt.edu";

my $mirnacount_option = "-s";                                    #ignore score
my $trnacount_option  = "--length --sequence";
my $utr3_option  = "-m 0";


my $bowtie1_option    = "-a -m 100 --best --strata -v 1 -p 8";

my $demultiplexing_config = {
  general        => { "task_name" => "parclip", },
  demultiplexing => {
    class      => "Format::Demultiplexing",
    perform    => 0,
    target_dir => "${root}/demultiplexing",
    option     => "",
    source     => { "2797-KCV-1" => ["/autofs/blue_sequencer/Runs/projects/2797-KCV/2014-02-06/2797-KCV-1_1.fastq.gz"], },
    maps       => { "2797-KCV-1" => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201403_parclip_2797/demultiplexing/pbs/2797-KCV-1.map"], },
    cqstools   => $cqstools,
    sh_direct  => 0,
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
    target_dir => "${root}/cutadapt",
    option     => "-m 12 -O 10 -e 0.083",
    source_ref => "demultiplexing",
    adaptor    => "TGGAATTCTCGGGTGCCAAGG",
    extension  => "_clipped.fastq",
    sh_direct  => 0,
    gzipped    => 1,
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
    target_dir => "${root}/fastqlen",
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
    target_dir => "${root}/identical",
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
};

#performConfig($demultiplexing_config);

my $kcv2797human = {
  files => {
    "2797-KCV-1_RPI40_Ago2INS1Huh7"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201403_parclip_2797/identical/result/2797-KCV-1_RPI40_Ago2INS1Huh7_clipped_identical.fastq.gz"],
    "2797-KCV-1_RPI41_Ago3INS1Huh7"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201403_parclip_2797/identical/result/2797-KCV-1_RPI41_Ago3INS1Huh7_clipped_identical.fastq.gz"],
    "2797-KCV-1_RPI42_Ago2INS1HCEAC" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201403_parclip_2797/identical/result/2797-KCV-1_RPI42_Ago2INS1HCEAC_clipped_identical.fastq.gz"],
    "2797-KCV-1_RPI43_Ago3INS1HCEAC" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201403_parclip_2797/identical/result/2797-KCV-1_RPI43_Ago3INS1HCEAC_clipped_identical.fastq.gz"],
    "2797-KCV-1_RPI47_Ago2MIN6Huh7"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201403_parclip_2797/identical/result/2797-KCV-1_RPI47_Ago2MIN6Huh7_clipped_identical.fastq.gz"],
    "2797-KCV-1_RPI48_Ago3MIN6Huh7"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201403_parclip_2797/identical/result/2797-KCV-1_RPI48_Ago3MIN6Huh7_clipped_identical.fastq.gz"],
  },
  countfiles => {
    "2797-KCV-1_RPI40_Ago2INS1Huh7"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201403_parclip_2797/identical/result/2797-KCV-1_RPI40_Ago2INS1Huh7_clipped_identical.fastq.dupcount"],
    "2797-KCV-1_RPI41_Ago3INS1Huh7"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201403_parclip_2797/identical/result/2797-KCV-1_RPI41_Ago3INS1Huh7_clipped_identical.fastq.dupcount"],
    "2797-KCV-1_RPI42_Ago2INS1HCEAC" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201403_parclip_2797/identical/result/2797-KCV-1_RPI42_Ago2INS1HCEAC_clipped_identical.fastq.dupcount"],
    "2797-KCV-1_RPI43_Ago3INS1HCEAC" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201403_parclip_2797/identical/result/2797-KCV-1_RPI43_Ago3INS1HCEAC_clipped_identical.fastq.dupcount"],
    "2797-KCV-1_RPI47_Ago2MIN6Huh7"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201403_parclip_2797/identical/result/2797-KCV-1_RPI47_Ago2MIN6Huh7_clipped_identical.fastq.dupcount"],
    "2797-KCV-1_RPI48_Ago3MIN6Huh7"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201403_parclip_2797/identical/result/2797-KCV-1_RPI48_Ago3MIN6Huh7_clipped_identical.fastq.dupcount"],
  },
  task_name           => "2797-KCV-hg19",
  mirna_coordinate    => "/data/cqs/shengq1/reference/miRBase20/hsa.gff3",
  trna_coordinate     => "/data/cqs/guoy1/reference/smallrna/hg19_tRNA_ucsc_ensembl.bed",
  bowtie1_index       => "/data/cqs/guoy1/reference/hg19/bowtie_index_hg19_rCRS_1.0.0/hg19_rCRS",
  genome_2bit         => "/data/cqs/guoy1/reference/hg19/hg19_rCRS.2bit",
  mirna_db            => "/data/cqs/shengq1/reference/miRBase20/hsa.mature.dna.db",
  binding_db          => "/data/cqs/shengq1/reference/targetscan/targetscan_v61_hg19.bed",
  utr3_db             => "/data/cqs/shengq1/reference/utr3/20140612_ucsc_hg19_3UTR.txt",
  smallrna_coordinate => "/data/cqs/guoy1/reference/smallrna/hg19_smallRNA_ucsc_ensembl.bed"
};

my $kcv2797mouse = {
  files => {
    "2797-KCV-1_RPI47_Ago2MIN6Huh7" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201403_parclip_2797/identical/result/2797-KCV-1_RPI47_Ago2MIN6Huh7_clipped_identical.fastq.gz"],
    "2797-KCV-1_RPI48_Ago3MIN6Huh7" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201403_parclip_2797/identical/result/2797-KCV-1_RPI48_Ago3MIN6Huh7_clipped_identical.fastq.gz"],
  },
  countfiles => {
    "2797-KCV-1_RPI47_Ago2MIN6Huh7"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201403_parclip_2797/identical/result/2797-KCV-1_RPI47_Ago2MIN6Huh7_clipped_identical.fastq.dupcount"],
    "2797-KCV-1_RPI48_Ago3MIN6Huh7"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201403_parclip_2797/identical/result/2797-KCV-1_RPI48_Ago3MIN6Huh7_clipped_identical.fastq.dupcount"],
  },
  task_name        => "2797-KCV-mm10",
  mirna_coordinate => "/data/cqs/shengq1/reference/miRBase20/mmu.gff3",
  trna_coordinate  => "/data/cqs/guoy1/reference/smallrna/mm10_tRNA_ucsc_ensembl.bed",
  bowtie1_index    => "/data/cqs/shengq1/reference/mm10/bowtie_index/mm10",
  genome_2bit      => "/data/cqs/guoy1/reference/mm10/mm10.2bit",
  mirna_db         => "/data/cqs/shengq1/reference/miRBase20/mmu.mature.dna.db",
};

my @datasets = (
  $kcv2797human

    #, $kcv2797mouse
);

foreach my $dataset (@datasets) {
  my $target_dir     = create_directory_or_die( $root . "/" . $dataset->{task_name} );
  my $parclip_config = {
    general => { "task_name" => "parclip", },

    #    bowtie1out => {
    #      class         => "Bowtie1",
    #      perform       => 1,
    #      target_dir    => "${target_dir}/bowtie1out",
    #      option        => $bowtie1_option,
    #      source        => $dataset->{files},
    #      bowtie1_index => $dataset->{bowtie1_index},
    #      samformat     => 0,
    #      samonly       => 0,
    #      sh_direct     => 0,
    #      pbs           => {
    #        "email"    => $email,
    #        "nodes"    => "1:ppn=1",
    #        "walltime" => "72",
    #        "mem"      => "20gb"
    #      },
    #    },
    #    PARalyzer => {
    #      class      => "ParClip::PARalyzer",
    #      perform    => 1,
    #      target_dir => "${target_dir}/paralyzer",
    #      option     => "",
    #      source_ref => "bowtie1out",
    #      genome2bit => $dataset->{genome_2bit},
    #      mirna_db   => $dataset->{mirna_db},
    #      sh_direct  => 1,
    #      pbs        => {
    #        "email"    => $email,
    #        "nodes"    => "1:ppn=1",
    #        "walltime" => "72",
    #        "mem"      => "20gb"
    #      },
    #    },
    #    annotation => {
    #      class            => "CQS::ParalyzerClusterAnnotator",
    #      perform          => 1,
    #      target_dir       => "${target_dir}/paralyzer",
    #      option           => "-f miRNA",
    #      source_ref       => [ "PARalyzer", ".cluster.csv" ],
    #      cqstools         => $cqstools,
    #      coordinate_files => [ $dataset->{mirna_coordinate}, $dataset->{trna_coordinate} ],
    #      sh_direct        => 1,
    #      pbs              => {
    #        "email"    => $email,
    #        "nodes"    => "1:ppn=1",
    #        "walltime" => "72",
    #        "mem"      => "20gb"
    #      },
    #    },
    bowtie1bam => {
      class         => "Bowtie1",
      perform       => 1,
      target_dir    => "${target_dir}/bowtie1bam",
      option        => $bowtie1_option,
      source        => $dataset->{files},
      bowtie1_index => $dataset->{bowtie1_index},
      samonly       => 0,
      sh_direct     => 0,
      pbs           => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "72",
        "mem"      => "20gb"
      },
    },
    mirna_count => {
      class           => "MirnaCount",
      perform         => 1,
      target_dir      => "${target_dir}/count_miRNA",
      option          => $mirnacount_option,
      source_ref      => "bowtie1bam",
      fastq_files => $dataset->{files},
      seqcount    => $dataset->{countfiles},
      cqs_tools       => $cqstools,
      gff_file        => $dataset->{mirna_coordinate},
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
    smallRNA_1mm_count => {
      class           => "CQSMappedCount",
      perform         => 1,
      target_dir      => "${target_dir}/count_smallRNA",
      option          => $trnacount_option,
      source_ref      => "bowtie1bam",
      fastq_files => $dataset->{files},
      seqcount    => $dataset->{countfiles},
      cqs_tools       => $cqstools,
      gff_file        => $dataset->{smallrna_coordinate},
      samtools        => $samtools,
      sh_direct       => 1,
      pbs             => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "72",
        "mem"      => "20gb"
      },
    },
    utr3_count => {
      class           => "CQSMappedCount",
      perform         => 1,
      target_dir      => "${target_dir}/count_3utr",
      option          => "-m 0",
      source_ref      => "bowtie1bam",
      fastq_files => $dataset->{files},
      seqcount    => $dataset->{countfiles},
      cqs_tools       => $cqstools,
      gff_file        => $dataset->{utr3_db},
      samtools        => $samtools,
      sh_direct       => 1,
      pbs             => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "72",
        "mem"      => "20gb"
      },
    },
    binding_count => {
      class           => "CQSMappedCount",
      perform         => 1,
      target_dir      => "${target_dir}/count_binding",
      option          => "-m 0",
      source_ref      => "bowtie1bam",
      fastq_files => $dataset->{files},
      seqcount    => $dataset->{countfiles},
      cqs_tools       => $cqstools,
      gff_file        => $dataset->{binding_db},
      samtools        => $samtools,
      sh_direct       => 1,
      pbs             => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "72",
        "mem"      => "20gb"
      },
    },
    sequencetask => {
      class      => "CQS::SequenceTask",
      perform    => 1,
      target_dir => "${target_dir}/sequencetask",
      option     => "",
      source     => {
        T1_individual => [
          #"bowtie1out", "PARalyzer",
          "bowtie1bam", "mirna_count", "smallRNA_1mm_count", "utr3_count", "binding_count"
        ],
        #T2_summary => ["annotation"],
      },
      sh_direct => 0,
      pbs       => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "72",
        "mem"      => "20gb"
      },
    }
  };

  performConfig($parclip_config);
}

1;
