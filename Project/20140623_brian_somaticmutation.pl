#!/usr/bin/perl
use strict;
use warnings;

use CQS::FileUtils;
use CQS::ClassFactory;

my $target_dir = create_directory_or_die("/scratch/cqs/shengq1/rnaseq/20140623_brian_rnaseq_somaticmutation");

my $email    = "quanhu.sheng\@vanderbilt.edu";
my $cqstools = "/home/shengq1/cqstools/CQS.Tools.exe";
my $samtools = "/home/shengq1/local/bin/samtools/samtools";

##hg19.16569.MT###
my $fasta_file_16569_MT    = "/data/cqs/shengq1/reference/hg19_16569_MT/bwa_index_0.7.8/hg19_16569_MT.fa";
my $cosmic_file_16569_MT   = "/data/cqs/shengq1/reference/cosmic/cosmic_v69_hg19_16569_MT.vcf";
my $snp_file_16569_MT      = "/data/cqs/shengq1/reference/dbsnp/human_GRCh37_v141_16569_MT.vcf";

my $annovar_param = "-protocol refGene,snp138,cosmic68,esp6500si_all,1000g2012feb -operation g,f,f,f,f --remove";
my $annovar_db    = "/scratch/cqs/shengq1/references/annovar/humandb/";

my $config = {
  general => { task_name => "somaticmutation" },

  files => {

    #cqstools file_def -i . -r -f t.sorted.bam$ -d
    "2944-JP-01" => ["/gpfs21/scratch/cqs/shengq1/rnaseq/20140623_brian_rnaseq_somaticmutation/star/2944-JP-1/Aligned.out.sorted.bam"],
    "2944-JP-02" => ["/gpfs21/scratch/cqs/shengq1/rnaseq/20140623_brian_rnaseq_somaticmutation/star/2944-JP-2/Aligned.out.sorted.bam"],
    "2944-JP-03" => ["/gpfs21/scratch/cqs/shengq1/rnaseq/20140623_brian_rnaseq_somaticmutation/star/2944-JP-3/Aligned.out.sorted.bam"],
    "2944-JP-04" => ["/gpfs21/scratch/cqs/shengq1/rnaseq/20140623_brian_rnaseq_somaticmutation/star/2944-JP-4/Aligned.out.sorted.bam"],
    "2944-JP-05" => ["/gpfs21/scratch/cqs/shengq1/rnaseq/20140623_brian_rnaseq_somaticmutation/star/2944-JP-5/Aligned.out.sorted.bam"],
    "2944-JP-06" => ["/gpfs21/scratch/cqs/shengq1/rnaseq/20140623_brian_rnaseq_somaticmutation/star/2944-JP-6/Aligned.out.sorted.bam"],
    "2944-JP-07" => ["/gpfs21/scratch/cqs/shengq1/rnaseq/20140623_brian_rnaseq_somaticmutation/star/2944-JP-7/Aligned.out.sorted.bam"],
    "2944-JP-08" => ["/gpfs21/scratch/cqs/shengq1/rnaseq/20140623_brian_rnaseq_somaticmutation/star/2944-JP-8/Aligned.out.sorted.bam"],
    "2944-JP-09" => ["/gpfs21/scratch/cqs/shengq1/rnaseq/20140623_brian_rnaseq_somaticmutation/star/2944-JP-9/Aligned.out.sorted.bam"],
    "2944-JP-10" => ["/gpfs21/scratch/cqs/shengq1/rnaseq/20140623_brian_rnaseq_somaticmutation/star/2944-JP-10/Aligned.out.sorted.bam"],
  },
  groups => {
    "2944-JP-0102" => [ "2944-JP-01", "2944-JP-02" ],
    "2944-JP-0304" => [ "2944-JP-03", "2944-JP-04" ],
    "2944-JP-0506" => [ "2944-JP-05", "2944-JP-06" ],
    "2944-JP-0708" => [ "2944-JP-07", "2944-JP-08" ],
    "2944-JP-0910" => [ "2944-JP-09", "2944-JP-10" ],
  },
  varscan2 => {
    class           => "VarScan2::Somatic",
    perform         => 0,
    target_dir      => "${target_dir}/varscan2",
    option          => "--min-coverage 10",
    mpileup_options => "-A -q 20 -Q 20",
    java_option     => "-Xmx40g",
    source_ref      =>  "files" ,
    groups_ref      =>  "groups",
    fasta_file      => $fasta_file_16569_MT,
    somatic_p_value => 0.05,
    sh_direct       => 0,
    VarScan2_jar    => "/home/shengq1/local/bin/VarScan.v2.3.5.jar",
    pbs             => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  annovar_varscan2 => {
    class      => "Annovar",
    perform    => 1,
    target_dir => "${target_dir}/varscan2",
    source_ref => "varscan2",
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

#performTask($config, "rsmc_16569");

1;

