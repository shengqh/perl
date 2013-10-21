#!/usr/bin/perl
use strict;
use warnings;

use CQS::QC;
use CQS::DNASeq;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::SomaticMutation;
use CQS::ClassFactory;

my $vangard = "VANGARD00180";

my $target_dir = create_directory_or_die("/scratch/cqs/shengq1/vangard/${vangard}_liuqi_exome_lungcancer");

my $bwa_dir = "${target_dir}/bwa_refine";

my $fasta_file           = "/data/cqs/guoy1/reference/hg19/bwa_index_0.7.4/hg19_chr.fa";
my $transcript_gtf       = "/scratch/cqs/shengq1/references/hg19/Homo_sapiens.GRCh37.73.gtf";
my $transcript_gtf_index = "/scratch/cqs/shengq1/gtfindex/hg19_GRCh37_73";

#my $annovar_param =  "-buildver hg19 -protocol refGene,phastConsElements46way,genomicSuperDups,esp6500si_all,1000g2012apr_all,snp137,ljb2_all,cosmic64 -operation g,r,r,f,f,f,f,f --remove --otherinfo";
my $annovar_param = "-protocol refGene,snp137,cosmic64,esp6500si_all,1000g2012apr_all -operation g,f,f,f,f --remove --otherinfo";
my $annovar_db    = "/scratch/cqs/shengq1/references/annovar/humandb/";

my $email = "quanhu.sheng\@vanderbilt.edu";

my $config = {
  general    => { task_name => "${vangard}" },
  fastqfiles => {
    "2055-PM-00" => [ "/autofs/blue_sequencer/Runs/projects/2055-PM/2013-09-24/2055-PM-0_1.fastq.gz",  "/autofs/blue_sequencer/Runs/projects/2055-PM/2013-09-24/2055-PM-0_2.fastq.gz" ],
    "2055-PM-01" => [ "/autofs/blue_sequencer/Runs/projects/2055-PM/2013-09-24/2055-PM-1_1.fastq.gz",  "/autofs/blue_sequencer/Runs/projects/2055-PM/2013-09-24/2055-PM-1_2.fastq.gz" ],
    "2055-PM-02" => [ "/autofs/blue_sequencer/Runs/projects/2055-PM/2013-09-24/2055-PM-2_1.fastq.gz",  "/autofs/blue_sequencer/Runs/projects/2055-PM/2013-09-24/2055-PM-2_2.fastq.gz" ],
    "2055-PM-04" => [ "/autofs/blue_sequencer/Runs/projects/2055-PM/2013-09-24/2055-PM-4_1.fastq.gz",  "/autofs/blue_sequencer/Runs/projects/2055-PM/2013-09-24/2055-PM-4_2.fastq.gz" ],
    "2055-PM-05" => [ "/autofs/blue_sequencer/Runs/projects/2055-PM/2013-09-24/2055-PM-5_1.fastq.gz",  "/autofs/blue_sequencer/Runs/projects/2055-PM/2013-09-24/2055-PM-5_2.fastq.gz" ],
    "2055-PM-06" => [ "/autofs/blue_sequencer/Runs/projects/2055-PM/2013-09-24/2055-PM-6_1.fastq.gz",  "/autofs/blue_sequencer/Runs/projects/2055-PM/2013-09-24/2055-PM-6_2.fastq.gz" ],
    "2055-PM-07" => [ "/autofs/blue_sequencer/Runs/projects/2055-PM/2013-09-24/2055-PM-7_1.fastq.gz",  "/autofs/blue_sequencer/Runs/projects/2055-PM/2013-09-24/2055-PM-7_2.fastq.gz" ],
    "2055-PM-08" => [ "/autofs/blue_sequencer/Runs/projects/2055-PM/2013-09-24/2055-PM-8_1.fastq.gz",  "/autofs/blue_sequencer/Runs/projects/2055-PM/2013-09-24/2055-PM-8_2.fastq.gz" ],
    "2055-PM-09" => [ "/autofs/blue_sequencer/Runs/projects/2055-PM/2013-07-12/2055-PM-9_1.fastq.gz",  "/autofs/blue_sequencer/Runs/projects/2055-PM/2013-07-12/2055-PM-9_2.fastq.gz" ],
    "2055-PM-10" => [ "/autofs/blue_sequencer/Runs/projects/2055-PM/2013-07-12/2055-PM-10_1.fastq.gz", "/autofs/blue_sequencer/Runs/projects/2055-PM/2013-07-12/2055-PM-10_2.fastq.gz" ],
    "2055-PM-11" => [ "/autofs/blue_sequencer/Runs/projects/2055-PM/2013-07-12/2055-PM-11_1.fastq.gz", "/autofs/blue_sequencer/Runs/projects/2055-PM/2013-07-12/2055-PM-11_2.fastq.gz" ],
    "2055-PM-12" => [ "/autofs/blue_sequencer/Runs/projects/2055-PM/2013-07-12/2055-PM-12_1.fastq.gz", "/autofs/blue_sequencer/Runs/projects/2055-PM/2013-07-12/2055-PM-12_2.fastq.gz" ],
    "2055-PM-13" => [ "/autofs/blue_sequencer/Runs/projects/2055-PM/2013-07-12/2055-PM-13_1.fastq.gz", "/autofs/blue_sequencer/Runs/projects/2055-PM/2013-07-12/2055-PM-13_2.fastq.gz" ],
  },
  groups => {
    "P4413_1_NL"       => [ "2055-PM-00", "2055-PM-01" ],
    "P4413_2_SEV_D"    => [ "2055-PM-00", "2055-PM-02" ],
    "P4413_3_CLS"      => [ "2055-PM-00", "2055-PM-13" ],
    "P4413_4_INVASIVE" => [ "2055-PM-00", "2055-PM-04" ],
    "P6708"            => [ "2055-PM-11", "2055-PM-05" ],
    "P6578"            => [ "2055-PM-09", "2055-PM-06" ],
    "P6186"            => [ "2055-PM-10", "2055-PM-07" ],
    "P6741"            => [ "2055-PM-12", "2055-PM-08" ],
  },
  fastqc => {
    class      => "FastQC",
    perform    => 1,
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
    perform    => 1,
    target_dir => "${target_dir}/bwa",
    option     => "-q 15 -t 8",
    fasta_file => $fasta_file,
    source_ref => "fastqfiles",
    sh_direct  => 0,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  refine => {
    class              => "GATKRefine",
    perform            => 1,
    target_dir         => "${target_dir}/refine",
    option             => "-Xmx40g",
    fasta_file         => $fasta_file,
    source_ref         => "bwa",
    thread_count       => 8,
    vcf_files          => ["/data/cqs/shengq1/reference/snp137/human/00-All.vcf"],
    gatk_jar           => "/home/shengq1/local/bin/GATK/GenomeAnalysisTK.jar",
    markDuplicates_jar => "/home/shengq1/local/bin/picard/MarkDuplicates.jar",
    sh_direct          => 0,
    pbs                => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  muTect => {
    class       => "MuTect",
    perform     => 1,
    target_dir  => "${target_dir}/muTect",
    option      => "",
    source_ref  => "refine",
    groups_ref  => "groups",
    java_option => "-Xmx40g",
    fasta_file  => $fasta_file,
    cosmic_file => "/data/cqs/shengq1/reference/cosmic/cosmic_v65_28052013.hg19.16571.vcf",
    dbsnp_file  => "/data/cqs/shengq1/reference/snp137/human/00-All.vcf",
    sh_direct   => 0,
    muTect_jar  => "/home/shengq1/local/bin/muTect-1.1.4.jar",
    pbs         => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  annovar_mutect => {
    class      => "Annovar",
    perform    => 0,
    target_dir => "${target_dir}/muTect",
    option     => $annovar_param,
    source_ref => [ "muTect", "\.vcf\$" ],
    annovar_db => $annovar_db,
    buildver   => "hg19",
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
    class      => "GATKSNPIndel",
    perform    => 0,
    target_dir => "${target_dir}/SNPindel",
    option     => "-l INFO -G Standard -stand_call_conf 50.0 -stand_emit_conf 10.0 -dcov 200 -nct 8",
    source_ref => "refine",

    #groups_ref  => "groups",
    java_option => "-Xmx40g",
    fasta_file  => $fasta_file,
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
    perform    => 0,
    target_dir => "${target_dir}/SNPindel",
    source_ref => "snpindel",
    option     => $annovar_param,
    annovar_db => $annovar_db,
    buildver   => "hg19",
    sh_direct  => 1,
    isvcf      => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "10gb"
    },
  },
  rsmc => {
    class            => "RSMC",
    perform          => 0,
    target_dir       => "${target_dir}/rsmc",
    option           => "-c 8",                                              #thread mode
    source_ref       => "refine",
    groups_ref       => "groups",
    source_type      => "bam",                                               #source_type can be bam/mpileup
    fasta_file       => $fasta_file,
    annovar_buildver => "hg19",
    rnaediting_db    => "/data/cqs/shengq1/reference/rnaediting/hg19.txt",
    sh_direct        => 0,
    execute_file     => "/home/shengq1/rsmc/rsmc.exe",
    pbs              => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "10gb"
    },
  },
};

performConfig($config);

1;
