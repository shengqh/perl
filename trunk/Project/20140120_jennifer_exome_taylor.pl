#!/usr/bin/perl
use strict;
use warnings;

use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ClassFactory;
use CQS::CNV;

my $target_dir = create_directory_or_die("/scratch/cqs/shengq1/dnaseq/20140120_jennifer_exome_taylor");

my $bwa_dir = "${target_dir}/bwa_refine";

my $fasta_file           = "/data/cqs/guoy1/reference/hg19/bwa_index_0.7.4/hg19_chr.fa";
my $transcript_gtf       = "/data/cqs/guoy1/reference/annotation2/hg19/Homo_sapiens.GRCh37.68.gtf";
my $transcript_gtf_index = "/scratch/cqs/shengq1/gtfindex/hg19_GRCh37_68";
my $cosmic               = "/data/cqs/shengq1/reference/cosmic/cosmic_v67_20131024.hg19.16571.vcf";
my $dbsnp                = "/data/cqs/shengq1/reference/snp137/human/00-All.vcf";
my $gatk                 = "/home/shengq1/local/bin/GATK/GenomeAnalysisTK.jar";
my $picard_dir           = "/home/shengq1/local/bin/picard/";
my $mutect               = "/home/shengq1/local/bin/muTect-1.1.4.jar";
my $bedfile              = "/scratch/cqs/lij17/cnv/SureSelect_XT_Human_All_Exon_V4_withoutchr_withoutY_lite.bed";
my $conifer              = "/home/shengq1/pylibs/bin/conifer.py";

#my $annovar_param =  "-buildver hg19 -protocol refGene,phastConsElements46way,genomicSuperDups,esp6500si_all,1000g2012apr_all,snp137,ljb2_all,cosmic64 -operation g,r,r,f,f,f,f,f --remove --otherinfo";
my $annovar_param = "-protocol refGene,snp137,cosmic64,esp6500si_all,1000g2012apr_all -operation g,f,f,f,f --remove --otherinfo";
my $annovar_db    = "/scratch/cqs/shengq1/references/annovar/humandb/";

my $email = "quanhu.sheng\@vanderbilt.edu";

my $config = {
  general    => { task_name => "exome_taylor" },
  fastqfiles => {
    "2476-JP-01" => [ "/autofs/blue_sequencer/Runs/projects/2476-JP/2013-09-24/2476-JP-1_1.fastq.gz",  "/autofs/blue_sequencer/Runs/projects/2476-JP/2013-09-24/2476-JP-1_2.fastq.gz" ],
    "2476-JP-02" => [ "/autofs/blue_sequencer/Runs/projects/2476-JP/2013-09-24/2476-JP-2_1.fastq.gz",  "/autofs/blue_sequencer/Runs/projects/2476-JP/2013-09-24/2476-JP-2_2.fastq.gz" ],
    "2476-JP-04" => [ "/autofs/blue_sequencer/Runs/projects/2476-JP/2013-09-24/2476-JP-4_1.fastq.gz",  "/autofs/blue_sequencer/Runs/projects/2476-JP/2013-09-24/2476-JP-4_2.fastq.gz" ],
    "2476-JP-05" => [ "/autofs/blue_sequencer/Runs/projects/2476-JP/2013-09-24/2476-JP-5_1.fastq.gz",  "/autofs/blue_sequencer/Runs/projects/2476-JP/2013-09-24/2476-JP-5_2.fastq.gz" ],
    "2476-JP-06" => [ "/autofs/blue_sequencer/Runs/projects/2476-JP/2013-09-24/2476-JP-6_1.fastq.gz",  "/autofs/blue_sequencer/Runs/projects/2476-JP/2013-09-24/2476-JP-6_2.fastq.gz" ],
    "2476-JP-07" => [ "/autofs/blue_sequencer/Runs/projects/2476-JP/2013-09-24/2476-JP-7_1.fastq.gz",  "/autofs/blue_sequencer/Runs/projects/2476-JP/2013-09-24/2476-JP-7_2.fastq.gz" ],
    "2476-JP-08" => [ "/autofs/blue_sequencer/Runs/projects/2476-JP/2013-09-24/2476-JP-8_1.fastq.gz",  "/autofs/blue_sequencer/Runs/projects/2476-JP/2013-09-24/2476-JP-8_2.fastq.gz" ],
    "2476-JP-09" => [ "/autofs/blue_sequencer/Runs/projects/2476-JP/2013-09-24/2476-JP-9_1.fastq.gz",  "/autofs/blue_sequencer/Runs/projects/2476-JP/2013-09-24/2476-JP-9_2.fastq.gz" ],
    "2476-JP-10" => [ "/autofs/blue_sequencer/Runs/projects/2476-JP/2013-09-24/2476-JP-10_1.fastq.gz", "/autofs/blue_sequencer/Runs/projects/2476-JP/2013-09-24/2476-JP-10_2.fastq.gz" ],
    "2476-JP-11" => [ "/autofs/blue_sequencer/Runs/projects/2476-JP/2013-09-24/2476-JP-11_1.fastq.gz", "/autofs/blue_sequencer/Runs/projects/2476-JP/2013-09-24/2476-JP-11_2.fastq.gz" ],
    "2476-JP-12" => [ "/autofs/blue_sequencer/Runs/projects/2476-JP/2013-09-24/2476-JP-12_1.fastq.gz", "/autofs/blue_sequencer/Runs/projects/2476-JP/2013-09-24/2476-JP-12_2.fastq.gz" ],
    "2476-JP-13" => [ "/autofs/blue_sequencer/Runs/projects/2476-JP/2013-09-24/2476-JP-13_1.fastq.gz", "/autofs/blue_sequencer/Runs/projects/2476-JP/2013-09-24/2476-JP-13_2.fastq.gz" ],
    "2476-JP-14" => [ "/autofs/blue_sequencer/Runs/projects/2476-JP/2013-09-24/2476-JP-14_1.fastq.gz", "/autofs/blue_sequencer/Runs/projects/2476-JP/2013-09-24/2476-JP-14_2.fastq.gz" ],
    "2476-JP-15" => [ "/autofs/blue_sequencer/Runs/projects/2476-JP/2013-09-24/2476-JP-15_1.fastq.gz", "/autofs/blue_sequencer/Runs/projects/2476-JP/2013-09-24/2476-JP-15_2.fastq.gz" ],
    "2476-JP-16" => [ "/autofs/blue_sequencer/Runs/projects/2476-JP/2013-09-24/2476-JP-16_1.fastq.gz", "/autofs/blue_sequencer/Runs/projects/2476-JP/2013-09-24/2476-JP-16_2.fastq.gz" ],
    "2476-JP-17" => [ "/autofs/blue_sequencer/Runs/projects/2476-JP/2013-09-24/2476-JP-17_1.fastq.gz", "/autofs/blue_sequencer/Runs/projects/2476-JP/2013-09-24/2476-JP-17_2.fastq.gz" ],
    "2476-JP-18" => [ "/autofs/blue_sequencer/Runs/projects/2476-JP/2013-11-12/2476-JP-18_1.fastq.gz", "/autofs/blue_sequencer/Runs/projects/2476-JP/2013-11-12/2476-JP-18_2.fastq.gz" ],
    "2476-JP-19" => [ "/autofs/blue_sequencer/Runs/projects/2476-JP/2013-12-07/2476-JP-19_1.fastq.gz", "/autofs/blue_sequencer/Runs/projects/2476-JP/2013-12-07/2476-JP-19_2.fastq.gz" ]
  },
  groups => {
    "2476-JP-0102" => [ "2476-JP-01", "2476-JP-02" ],
    "2476-JP-0506" => [ "2476-JP-05", "2476-JP-06" ],
    "2476-JP-0708" => [ "2476-JP-07", "2476-JP-08" ],
    "2476-JP-0910" => [ "2476-JP-09", "2476-JP-10" ],
    "2476-JP-1112" => [ "2476-JP-11", "2476-JP-12" ],
    "2476-JP-1314" => [ "2476-JP-13", "2476-JP-14" ],
    "2476-JP-1516" => [ "2476-JP-15", "2476-JP-16" ]
  },
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
    option     => "-t 8",
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
    perform            => 0,
    target_dir         => "${target_dir}/refine",
    option             => "-Xmx40g",
    fasta_file         => $fasta_file,
    source_ref         => "bwa",
    thread_count       => 8,
    vcf_files          => [$dbsnp],
    gatk_jar           => $gatk,
    markDuplicates_jar => "${picard_dir}/MarkDuplicates.jar",
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
    java_option => "-Xmx20g",
    fasta_file  => $fasta_file,
    cosmic_file => $cosmic,
    dbsnp_file  => $dbsnp,
    sh_direct   => 0,
    muTect_jar  => $mutect,
    pbs         => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "20gb"
    },
  },
  annovar_mutect => {
    class      => "Annovar",
    perform    => 1,
    target_dir => "${target_dir}/muTect",
    option     => $annovar_param,
    source_ref => [ "muTect", "\.vcf\$" ],
    annovar_db => $annovar_db,
    buildver   => "hg19",
    sh_direct  => 0,
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
    java_option => "-Xmx40g",
    fasta_file  => $fasta_file,
    vcf_files   => [$dbsnp],
    gatk_jar    => $gatk,
    sh_direct   => 0,
    pbs         => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  cnmops => {
    class       => "CNV::cnMops",
    perform     => 0,
    target_dir  => "${target_dir}/cnmops",
    option      => "",
    source_ref  => "refine",
    bedfile     => $bedfile,
    pairmode    => "paired",
    isbamsorted => 1,
    pbs         => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "720",
      "mem"      => "40gb"
    },
  },
  conifer => {
    class       => "CNV::Conifer",
    perform     => 0,
    target_dir  => "${target_dir}/conifer",
    option      => "",
    source_ref  => "refine",
    conifer     => $conifer,
    bedfile     => $bedfile,
    isbamsorted => 1,
    pbs         => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "720",
      "mem"      => "10gb"
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
};

performConfig($config);

1;
