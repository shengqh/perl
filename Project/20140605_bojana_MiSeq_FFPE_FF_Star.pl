#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;

my $task = "20140605_bojana_FFPE_FF";

my $target_dir = "/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF";

my $fasta_file_16569_MT  = "/data/cqs/shengq1/reference/hg19_16569_MT/bwa_index_0.7.8/hg19_16569_MT.fa";
my $transcript_gtf       = "/scratch/cqs/shengq1/references/hg19/Homo_sapiens.GRCh37.73.gtf";
my $transcript_gtf_index = "/scratch/cqs/shengq1/gtfindex/hg19_GRCh37_73";
my $bowtie2_index        = "/data/cqs/shengq1/reference/hg19_16569_MT/bowtie2_index_2.1.0/hg19_16569_MT";

my $annovar_param = "-protocol refGene,snp138,cosmic68,esp6500si_all -operation g,f,f,f --remove";
my $annovar_db    = "/scratch/cqs/shengq1/references/annovar/humandb/";

my $email = "quanhu.sheng\@vanderbilt.edu";

my $config = {
  general    => { task_name => $task },
  fastqfiles => {
    "IG-01" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-1_1.fastq.gz",  "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-1_2.fastq.gz" ],
    "IG-02" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-2_1.fastq.gz",  "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-2_2.fastq.gz" ],
    "IG-03" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-3_1.fastq.gz",  "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-3_2.fastq.gz" ],
    "IG-04" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-4_1.fastq.gz",  "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-4_2.fastq.gz" ],
    "IG-05" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-5_1.fastq.gz",  "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-5_2.fastq.gz" ],
    "IG-06" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-6_1.fastq.gz",  "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-6_2.fastq.gz" ],
    "IG-07" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-7_1.fastq.gz",  "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-7_2.fastq.gz" ],
    "IG-08" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-8_1.fastq.gz",  "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-8_2.fastq.gz" ],
    "IG-09" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-9_1.fastq.gz",  "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-9_2.fastq.gz" ],
    "IG-10" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-10_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-10_2.fastq.gz" ],
    "IG-11" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-11_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-11_2.fastq.gz" ],
    "IG-12" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-12_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-12_2.fastq.gz" ],
    "IG-13" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-13_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-13_2.fastq.gz" ],
    "IG-14" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-14_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-14_2.fastq.gz" ],
    "IG-15" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-15_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-15_2.fastq.gz" ],
    "IG-16" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-16_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-16_2.fastq.gz" ],
    "IG-17" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-17_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-17_2.fastq.gz" ],
    "IG-18" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-18_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-18_2.fastq.gz" ],
    "IG-19" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-19_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-19_2.fastq.gz" ],
    "IG-20" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-20_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-20_2.fastq.gz" ],
    "IG-21" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-21_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-21_2.fastq.gz" ],
    "IG-22" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-22_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-22_2.fastq.gz" ],
    "IG-23" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-23_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-23_2.fastq.gz" ],
    "IG-24" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-24_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-24_2.fastq.gz" ],
    "IG-25" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-25_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-25_2.fastq.gz" ],
    "IG-26" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-26_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-26_2.fastq.gz" ],
    "IG-27" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-27_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-27_2.fastq.gz" ],
    "IG-28" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-28_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-28_2.fastq.gz" ],
    "IG-29" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-29_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-29_2.fastq.gz" ],
    "IG-33" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-33_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-33_2.fastq.gz" ],
    "IG-34" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-34_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-34_2.fastq.gz" ],
    "IG-39" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-39_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-39_2.fastq.gz" ],
    "IG-40" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-40_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-40_2.fastq.gz" ],
    "IG-41" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-41_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-41_2.fastq.gz" ],
    "IG-42" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-42_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-42_2.fastq.gz" ],
    "IG-43" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-43_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-43_2.fastq.gz" ],
    "IG-44" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-44_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-44_2.fastq.gz" ],
    "IG-45" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-45_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-45_2.fastq.gz" ],
    "IG-46" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-46_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-46_2.fastq.gz" ],
    "IG-47" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-47_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-47_2.fastq.gz" ],
    "IG-49" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-49_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-49_2.fastq.gz" ],
    "IG-50" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-50_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-50_2.fastq.gz" ],
    "IG-51" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-51_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-51_2.fastq.gz" ],
    "IG-52" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-52_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-52_2.fastq.gz" ],
    "IG-53" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-53_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-53_2.fastq.gz" ],
    "IG-54" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-54_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-54_2.fastq.gz" ],
    "IG-55" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-55_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-55_2.fastq.gz" ],
    "IG-56" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-56_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-56_2.fastq.gz" ],
    "IG-57" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-57_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-57_2.fastq.gz" ],
    "IG-58" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-58_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-58_2.fastq.gz" ],
    "IG-59" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-59_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-59_2.fastq.gz" ],
    "IG-60" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-60_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-60_2.fastq.gz" ],
    "IG-61" => [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-61_1.fastq.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/fastq/IG-61_2.fastq.gz" ],
  },
  bamfiles => {
    "IG-01" => ["/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-1/Aligned.out.sorted.bam"],
    "IG-02" => ["/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-2/Aligned.out.sorted.bam"],
    "IG-03" => ["/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-3/Aligned.out.sorted.bam"],
    "IG-04" => ["/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-4/Aligned.out.sorted.bam"],
    "IG-05" => ["/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-5/Aligned.out.sorted.bam"],
    "IG-06" => ["/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-6/Aligned.out.sorted.bam"],
    "IG-07" => ["/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-7/Aligned.out.sorted.bam"],
    "IG-08" => ["/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-8/Aligned.out.sorted.bam"],
    "IG-09" => ["/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-9/Aligned.out.sorted.bam"],
    "IG-10" => ["/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-10/Aligned.out.sorted.bam"],
    "IG-11" => ["/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-11/Aligned.out.sorted.bam"],
    "IG-12" => ["/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-12/Aligned.out.sorted.bam"],
    "IG-13" => ["/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-13/Aligned.out.sorted.bam"],
    "IG-14" => ["/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-14/Aligned.out.sorted.bam"],
    "IG-15" => ["/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-15/Aligned.out.sorted.bam"],
    "IG-16" => ["/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-16/Aligned.out.sorted.bam"],
    "IG-17" => ["/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-17/Aligned.out.sorted.bam"],
    "IG-18" => ["/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-18/Aligned.out.sorted.bam"],
    "IG-19" => ["/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-19/Aligned.out.sorted.bam"],
    "IG-20" => ["/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-20/Aligned.out.sorted.bam"],
    "IG-21" => ["/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-21/Aligned.out.sorted.bam"],
    "IG-22" => ["/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-22/Aligned.out.sorted.bam"],
    "IG-23" => ["/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-23/Aligned.out.sorted.bam"],
    "IG-24" => ["/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-24/Aligned.out.sorted.bam"],
    "IG-25" => ["/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-25/Aligned.out.sorted.bam"],
    "IG-26" => ["/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-26/Aligned.out.sorted.bam"],
    "IG-27" => ["/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-27/Aligned.out.sorted.bam"],
    "IG-28" => ["/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-28/Aligned.out.sorted.bam"],
    "IG-29" => ["/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-29/Aligned.out.sorted.bam"],
    "IG-33" => ["/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-33/Aligned.out.sorted.bam"],
    "IG-34" => ["/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-34/Aligned.out.sorted.bam"],
    "IG-39" => ["/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-39/Aligned.out.sorted.bam"],
    "IG-40" => ["/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-40/Aligned.out.sorted.bam"],
    "IG-41" => ["/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-41/Aligned.out.sorted.bam"],
    "IG-42" => ["/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-42/Aligned.out.sorted.bam"],
    "IG-43" => ["/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-43/Aligned.out.sorted.bam"],
    "IG-44" => ["/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-44/Aligned.out.sorted.bam"],
    "IG-45" => ["/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-45/Aligned.out.sorted.bam"],
    "IG-46" => ["/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-46/Aligned.out.sorted.bam"],
    "IG-47" => ["/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-47/Aligned.out.sorted.bam"],
    "IG-49" => ["/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-49/Aligned.out.sorted.bam"],
    "IG-50" => ["/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-50/Aligned.out.sorted.bam"],
    "IG-51" => ["/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-51/Aligned.out.sorted.bam"],
    "IG-52" => ["/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-52/Aligned.out.sorted.bam"],
    "IG-53" => ["/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-53/Aligned.out.sorted.bam"],
    "IG-54" => ["/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-54/Aligned.out.sorted.bam"],
    "IG-55" => ["/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-55/Aligned.out.sorted.bam"],
    "IG-56" => ["/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-56/Aligned.out.sorted.bam"],
    "IG-57" => ["/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-57/Aligned.out.sorted.bam"],
    "IG-58" => ["/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-58/Aligned.out.sorted.bam"],
    "IG-59" => ["/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-59/Aligned.out.sorted.bam"],
    "IG-60" => ["/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-60/Aligned.out.sorted.bam"],
    "IG-61" => ["/scratch/cqs/shengq1/rnaseq/20140605_bojana_MiSeq_FFPE_FF/star/bams/IG-61/Aligned.out.sorted.bam"],
  },
  star_varscan2 => {
    class           => "VarScan2::Mpileup2snp",
    perform         => 1,
    target_dir      => "${target_dir}/star_varscan2",
    option          => "",
    mpileup_options => "",
    java_option     => "-Xmx40g",
    source_ref      => "bamfiles",
    fasta_file      => $fasta_file_16569_MT,
    sh_direct       => 0,
    VarScan2_jar    => "/home/shengq1/local/bin/VarScan.v2.3.5.jar",
    pbs             => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  star_annovar_varscan2 => {
    class      => "Annovar",
    perform    => 1,
    target_dir => "${target_dir}/star_varscan2",
    source_ref => "star_varscan2",
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
  fastqc => {
    class      => "FastQC",
    perform    => 1,
    target_dir => "${target_dir}/fastqc",
    option     => "",
    source_ref => "fastqfiles",
    sh_direct  => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=2",
      "walltime" => "2",
      "mem"      => "10gb"
    },
  },
  tophat2 => {
    class                => "Tophat2",
    perform              => 1,
    target_dir           => "${target_dir}/tophat2",
    option               => "--segment-length 25 -r 0 -p 8",
    source_ref           => "fastqfiles",
    bowtie2_index        => $bowtie2_index,
    transcript_gtf       => $transcript_gtf,
    transcript_gtf_index => $transcript_gtf_index,
    rename_bam           => 1,
    sh_direct            => 0,
    pbs                  => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "30gb"
    },
  },
  tophat2_varscan2 => {
    class           => "VarScan2::Mpileup2snp",
    perform         => 1,
    target_dir      => "${target_dir}/tophat2_varscan2",
    option          => "",
    mpileup_options => "",
    java_option     => "-Xmx40g",
    source_ref      => "tophat2",
    fasta_file      => $fasta_file_16569_MT,
    sh_direct       => 0,
    VarScan2_jar    => "/home/shengq1/local/bin/VarScan.v2.3.5.jar",
    pbs             => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  tophat2_annovar_varscan2 => {
    class      => "Annovar",
    perform    => 1,
    target_dir => "${target_dir}/tophat2_varscan2",
    source_ref => "tophat2_varscan2",
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
  sequencetask => {
    class      => "CQS::SequenceTask",
    perform    => 1,
    target_dir => "${target_dir}/sequencetask",
    option     => "",
    source     => { individual => [ "star_varscan2", "star_annovar_varscan2", "tophat2", "tophat2_varscan2", "tophat2_annovar_varscan2", "fastqc" ], },
    sh_direct  => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
};

performConfig($config);

1;
