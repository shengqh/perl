#!/usr/bin/perl
use strict;
use warnings;

use CQS::QC;
use CQS::DNASeq;
use CQS::FileUtils;
use CQS::SystemUtils;

my $vangard  ="VANGARD00095";

my $target_dir = create_directory_or_die("/scratch/cqs/shengq1/vangard/${vangard}_jennifer_rnaseq");

my $bwa_dir = "${target_dir}/bwa_refine";

my $transcript_gtf       = "/data/cqs/guoy1/reference/annotation2/hg19/Homo_sapiens.GRCh37.68.gtf";
my $transcript_gtf_index = "/scratch/cqs/shengq1/gtfindex/hg19_GRCh37_68";

my $email = "quanhu.sheng\@vanderbilt.edu";

my $config = {
  general => {
    bowtie2_index        => "/data/cqs/guoy1/reference/hg19/bowtie2_index/hg19",
    transcript_gtf       => $transcript_gtf,
    transcript_gtf_index => $transcript_gtf_index,
    path_file            => "/home/shengq1/local/bin/path.txt",
    task_name            => "${vangard}"
  },
  fastqfiles => {
    "2562-JP-1" => [ "/blue/sequencer/Runs/projects/2562-JP/2013-06-20/2562-JP-1_1.fastq.gz", "/blue/sequencer/Runs/projects/2562-JP/2013-06-20/2562-JP-1_2.fastq.gz" ],
    "2562-JP-2" => [ "/blue/sequencer/Runs/projects/2562-JP/2013-06-20/2562-JP-2_1.fastq.gz", "/blue/sequencer/Runs/projects/2562-JP/2013-06-20/2562-JP-2_2.fastq.gz" ]
  },
  fastqc => {
    target_dir => "${target_dir}/fastqc",
    option     => "",
    source_ref => "fastqfiles",
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "2",
      "mem"      => "10gb"
    },
  },
  bwa_refine => {
    target_dir         => $bwa_dir,
    option             => "-q 15 -t 8",
    option_samse       => "",
    option_sampe       => "",
    option_gatk        => "-Xmx40g",
    fasta_file         => "/data/cqs/shengq1/reference/hg19/hg19_chr.fa",
    source_ref         => "fastqfiles",
    thread_count       => 8,
    vcf_files          => ["/data/cqs/shengq1/reference/snp137/human/00-All.vcf"],
    gatk_jar           => "/home/shengq1/local/bin/GATK/GenomeAnalysisTK.jar",
    markDuplicates_jar => "/home/shengq1/local/bin/picard/MarkDuplicates.jar",
    pbs                => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  muTect => {
    target_dir    => "${target_dir}/muTect",
    option        => "-nt 8",
    source_ref    => "bwa_refine",
    fasta_file    => "/data/cqs/guoy1/reference/hg19/hg19_rCRS/hg19_rCRS.fa",
    cosmic_file   => "/data/cqs/shengq1/reference/cosmic/cosmic_v64_02042013.hg19.vcf",
    dbsnp_file    => "/data/cqs/shengq1/reference/snp137/dbsnp_137.b37.vcf",
    annovar_param => "--buildver hg19 --verdbsnp 137 --ver1000g 1000g2012apr --veresp 6500si --genetype refgene --alltranscript --remove",
    annovar_db    => "/scratch/cqs/shengq1/references/annovar/humandb/",
    sh_direct     => 1,
    execute_file  => "/home/shengq1/local/bin/muTect-1.1.4.jar",
    pbs           => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
};

fastqc_by_pbs( $config, "fastqc" );
bwa_refine( $config, "bwa_refine" );

1;
