#!/usr/bin/perl
use strict;
use warnings;

use CQS::QC;
use CQS::DNASeq;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::SomaticMutation;
use CQS::ClassFactory;

my $vangard = "VANGARD00095";

my $target_dir = create_directory_or_die("/scratch/cqs/shengq1/vangard/${vangard}_jennifer_exome");

my $bwa_dir = "${target_dir}/bwa_refine";

my $transcript_gtf       = "/data/cqs/guoy1/reference/annotation2/hg19/Homo_sapiens.GRCh37.68.gtf";
my $transcript_gtf_index = "/scratch/cqs/shengq1/gtfindex/hg19_GRCh37_68";

my $annovar_param =  "-buildver hg19 -protocol refGene,phastConsElements46way,genomicSuperDups,esp6500si_all,1000g2012apr_all,snp137,ljb2_all,cosmic64 -operation g,r,r,f,f,f,f,f --alltranscript --remove --otherinfo";
my $annovar_db = "/scratch/cqs/shengq1/references/annovar/humandb/";

my $email = "quanhu.sheng\@vanderbilt.edu";

my $config = {
  general    => { task_name => "${vangard}" },
  fastqfiles => {
    "2562-JP-1" => [ "/blue/sequencer/Runs/projects/2562-JP/2013-06-20/2562-JP-1_1.fastq.gz", "/blue/sequencer/Runs/projects/2562-JP/2013-06-20/2562-JP-1_2.fastq.gz" ],
    "2562-JP-2" => [ "/blue/sequencer/Runs/projects/2562-JP/2013-06-20/2562-JP-2_1.fastq.gz", "/blue/sequencer/Runs/projects/2562-JP/2013-06-20/2562-JP-2_2.fastq.gz" ]
  },
  groups => { "2562-JP" => [ "2562-JP-1", "2562-JP-2" ] },
  fastqc => {
    class      => "FastQC",
    perform    => 0,
    target_dir => "${target_dir}/fastqc",
    option     => "",
    source_ref => "fastqfiles",
    sh_direct  => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "2",
      "mem"      => "10gb"
    },
  },
  bwa => {
    class      => "BWA",
    perform    => 0,
    target_dir => "${target_dir}/bwa",
    option     => "-q 15 -t 8",
    fasta_file => "/data/cqs/guoy1/reference/hg19/bwa_index_0.7.4/hg19_chr.fa",
    source_ref => "fastqfiles",
    sh_direct  => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  refine => {
    class              => "GATKRefine",
    perform            => 0,
    target_dir         => "${target_dir}/refine",
    option             => "-Xmx40g",
    fasta_file         => "/data/cqs/guoy1/reference/hg19/bwa_index_0.7.4/hg19_chr.fa",
    source_ref         => "bwa",
    thread_count       => 8,
    vcf_files          => ["/data/cqs/shengq1/reference/snp137/human/00-All.vcf"],
    gatk_jar           => "/home/shengq1/local/bin/GATK/GenomeAnalysisTK.jar",
    markDuplicates_jar => "/home/shengq1/local/bin/picard/MarkDuplicates.jar",
    sh_direct          => 1,
    pbs                => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  muTect => {
    class         => "MuTect",
    perform       => 0,
    target_dir    => "${target_dir}/muTect",
    option        => "-nt 8",
    source_ref    => "refine",
    groups_ref    => "groups",
    java_option   => "-Xmx40g",
    fasta_file    => "/data/cqs/guoy1/reference/hg19/bwa_index_0.7.4/hg19_chr.fa",
    cosmic_file   => "/data/cqs/shengq1/reference/cosmic/cosmic_v65_28052013.hg19.16571.vcf",
    dbsnp_file    => "/data/cqs/shengq1/reference/snp137/human/00-All.vcf",
    sh_direct     => 1,
    muTect_jar    => "/home/shengq1/local/bin/muTect-1.1.4.jar",
    pbs           => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  annovar_mutect => {
    class      => "Annovar",
    perform    => 1,
    target_dir => "${target_dir}/muTect",
    option     => $annovar_param,
    source_ref => ["muTect", "\.vcf\$"],
    annovar_db => $annovar_db,
    sh_direct  => 1,
    isvcf      => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "10gb"
    },
  },
  snpindel => {
    class       => "GATKSNPIndel",
    perform     => 0,
    target_dir  => "${target_dir}/SNPindel",
    option      => "-l INFO -G Standard -stand_call_conf 50.0 -stand_emit_conf 10.0 -dcov 200 -nct 8",
    source_ref  => "refine",
    #groups_ref  => "groups",
    java_option => "-Xmx40g",
    fasta_file  => "/data/cqs/guoy1/reference/hg19/hg19_chr.fa",
    vcf_files   => ["/data/cqs/shengq1/reference/snp137/human/00-All.vcf"],
    gatk_jar    => "/home/shengq1/local/bin/GATK/GenomeAnalysisTK.jar",
    pbs         => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  annovar_snpindel => {
    class      => "Annovar",
    perform    => 1,
    target_dir => "${target_dir}/SNPindel",
    source_ref => "snpindel",
    option     => $annovar_param,
    annovar_db => $annovar_db,
    sh_direct  => 1,
    isvcf      => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "10gb"
    },
  },
};

performConfig($config);

1;
