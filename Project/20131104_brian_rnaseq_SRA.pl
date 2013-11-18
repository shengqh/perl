#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;
use CQS::FileUtils;

my $target_dir = create_directory_or_die("/scratch/cqs/shengq1/rnaseq/20131104_brian_rnaseq_SRA");

my $fasta_file           = "/data/cqs/guoy1/reference/hg19/bwa_index_0.7.4/hg19_chr.fa";
my $transcript_gtf       = "/scratch/cqs/shengq1/references/hg19/Homo_sapiens.GRCh37.73.gtf";
my $transcript_gtf_index = "/scratch/cqs/shengq1/gtfindex/hg19_GRCh37_73";

my $hg19_gff = "/scratch/cqs/shengq1/references/hg19/dexseq_gff/Homo_sapiens.GRCh37.73.dexseq.gff";
my $hg19_map = "/scratch/cqs/shengq1/references/hg19/Homo_sapiens.GRCh37.73.map";

my $bowtie2_index = "/data/cqs/guoy1/reference/hg19/bowtie2_index/hg19";

my $annovar_param = "-protocol refGene,snp137,cosmic64,esp6500si_all,1000g2012apr_all -operation g,f,f,f,f --remove --otherinfo";
my $annovar_db    = "/scratch/cqs/shengq1/references/annovar/humandb/";

my $cqstools = "/home/shengq1/cqstools/CQS.Tools.exe";

my $email = "quanhu.sheng\@vanderbilt.edu";
my $task  = "SRA";

my $config = {
  general => { task_name => $task },
  files   => {
    "SRR925687" => ["/gpfs20/data/lehmanbd/RNAseq_Gray/SRR925687.sra"],
    "SRR925688" => ["/gpfs20/data/lehmanbd/RNAseq_Gray/SRR925688.sra"],
    "SRR925689" => ["/gpfs20/data/lehmanbd/RNAseq_Gray/SRR925689.sra"],
    "SRR925690" => ["/gpfs20/data/lehmanbd/RNAseq_Gray/SRR925690.sra"],
    "SRR925691" => ["/gpfs20/data/lehmanbd/RNAseq_Gray/SRR925691.sra"],
    "SRR925692" => ["/gpfs20/data/lehmanbd/RNAseq_Gray/SRR925692.sra"],
    "SRR925693" => ["/gpfs20/data/lehmanbd/RNAseq_Gray/SRR925693.sra"],
    "SRR925694" => ["/gpfs20/data/lehmanbd/RNAseq_Gray/SRR925694.sra"],
    "SRR925695" => ["/gpfs20/data/lehmanbd/RNAseq_Gray/SRR925695.sra"],
    "SRR925696" => ["/gpfs20/data/lehmanbd/RNAseq_Gray/SRR925696.sra"],
    "SRR925697" => ["/gpfs20/data/lehmanbd/RNAseq_Gray/SRR925697.sra"],
    "SRR925698" => ["/gpfs20/data/lehmanbd/RNAseq_Gray/SRR925698.sra"],
    "SRR925699" => ["/gpfs20/data/lehmanbd/RNAseq_Gray/SRR925699.sra"],
    "SRR925700" => ["/gpfs20/data/lehmanbd/RNAseq_Gray/SRR925700.sra"],
    "SRR925701" => ["/gpfs20/data/lehmanbd/RNAseq_Gray/SRR925701.sra"],
    "SRR925702" => ["/gpfs20/data/lehmanbd/RNAseq_Gray/SRR925702.sra"],
    "SRR925703" => ["/gpfs20/data/lehmanbd/RNAseq_Gray/SRR925703.sra"],
    "SRR925704" => ["/gpfs20/data/lehmanbd/RNAseq_Gray/SRR925704.sra"],
    "SRR925705" => ["/gpfs20/data/lehmanbd/RNAseq_Gray/SRR925705.sra"],
    "SRR925706" => ["/gpfs20/data/lehmanbd/RNAseq_Gray/SRR925706.sra"],
    "SRR925707" => ["/gpfs20/data/lehmanbd/RNAseq_Gray/SRR925707.sra"],
    "SRR925708" => ["/gpfs20/data/lehmanbd/RNAseq_Gray/SRR925708.sra"],
    "SRR925710" => ["/gpfs20/data/lehmanbd/RNAseq_Gray/SRR925710.sra"],
    "SRR925711" => ["/gpfs20/data/lehmanbd/RNAseq_Gray/SRR925711.sra"],
    "SRR925713" => ["/gpfs20/data/lehmanbd/RNAseq_Gray/SRR925713.sra"],
    "SRR925714" => ["/gpfs20/data/lehmanbd/RNAseq_Gray/SRR925714.sra"],
    "SRR925715" => ["/gpfs20/data/lehmanbd/RNAseq_Gray/SRR925715.sra"],
    "SRR925716" => ["/gpfs20/data/lehmanbd/RNAseq_Gray/SRR925716.sra"],
    "SRR925717" => ["/gpfs20/data/lehmanbd/RNAseq_Gray/SRR925717.sra"],
    "SRR925718" => ["/gpfs20/data/lehmanbd/RNAseq_Gray/SRR925718.sra"],
    "SRR925719" => ["/gpfs20/data/lehmanbd/RNAseq_Gray/SRR925719.sra"],
    "SRR925720" => ["/gpfs20/data/lehmanbd/RNAseq_Gray/SRR925720.sra"],
    "SRR925721" => ["/gpfs20/data/lehmanbd/RNAseq_Gray/SRR925721.sra"],
    "SRR925723" => ["/gpfs20/data/lehmanbd/RNAseq_Gray/SRR925723.sra"],
    "SRR925724" => ["/gpfs20/data/lehmanbd/RNAseq_Gray/SRR925724.sra"],
    "SRR925725" => ["/gpfs20/data/lehmanbd/RNAseq_Gray/SRR925725.sra"],
    "SRR925726" => ["/gpfs20/data/lehmanbd/RNAseq_Gray/SRR925726.sra"],
    "SRR925727" => ["/gpfs20/data/lehmanbd/RNAseq_Gray/SRR925727.sra"],
    "SRR925728" => ["/gpfs20/data/lehmanbd/RNAseq_Gray/SRR925728.sra"],
    "SRR925729" => ["/gpfs20/data/lehmanbd/RNAseq_Gray/SRR925729.sra"],
    "SRR925732" => ["/gpfs20/data/lehmanbd/RNAseq_Gray/SRR925732.sra"],
    "SRR925734" => ["/gpfs20/data/lehmanbd/RNAseq_Gray/SRR925734.sra"],
    "SRR925735" => ["/scratch/cqs/shengq1/rnaseq/20131104_brian_rnaseq_SRA/raw/SRR925735.sra"],
    "SRR925736" => ["/gpfs20/data/lehmanbd/RNAseq_Gray/SRR925736.sra"],
    "SRR925737" => ["/gpfs20/data/lehmanbd/RNAseq_Gray/SRR925737.sra"],
    "SRR925740" => ["/gpfs20/data/lehmanbd/RNAseq_Gray/SRR925740.sra"],
    "SRR925741" => ["/gpfs20/data/lehmanbd/RNAseq_Gray/SRR925741.sra"],
    "SRR925742" => ["/gpfs20/data/lehmanbd/RNAseq_Gray/SRR925742.sra"],
    "SRR934631" => ["/gpfs20/data/lehmanbd/RNAseq_Gray/SRR934631.sra"],
    "SRR934632" => ["/gpfs20/data/lehmanbd/RNAseq_Gray/SRR934632.sra"],
    "SRR934633" => ["/gpfs20/data/lehmanbd/RNAseq_Gray/SRR934633.sra"],
    "SRR934634" => ["/gpfs20/data/lehmanbd/RNAseq_Gray/SRR934634.sra"],
    "SRR934635" => ["/gpfs20/data/lehmanbd/RNAseq_Gray/SRR934635.sra"],
    "SRR934636" => ["/gpfs20/data/lehmanbd/RNAseq_Gray/SRR934636.sra"],
    "SRR934637" => ["/gpfs20/data/lehmanbd/RNAseq_Gray/SRR934637.sra"],
    "SRR934639" => ["/gpfs20/data/lehmanbd/RNAseq_Gray/SRR934639.sra"],
    "SRR934641" => ["/gpfs20/data/lehmanbd/RNAseq_Gray/SRR934641.sra"],
    "SRR934643" => ["/gpfs20/data/lehmanbd/RNAseq_Gray/SRR934643.sra"],
    "SRR934644" => ["/gpfs20/data/lehmanbd/RNAseq_Gray/SRR934644.sra"],
    "SRR934645" => ["/gpfs20/data/lehmanbd/RNAseq_Gray/SRR934645.sra"]
  },
  sra2fastq => {
    class      => "SRA::FastqDump",
    perform    => 0,
    ispaired   => 1,
    target_dir => "${target_dir}/FastqDump",
    option     => "",
    source_ref => "files",
    sh_direct  => 0,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "10",
      "mem"      => "10gb"
    },
  },
  fastqc => {
    class      => "FastQC",
    perform    => 0,
    target_dir => "${target_dir}/fastqc",
    option     => "",
    source_ref => "sra2fastq",
    sh_direct  => 0,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=2",
      "walltime" => "2",
      "mem"      => "10gb"
    },
  },
  tophat2 => {
    class                => "Tophat2",
    perform              => 0,
    target_dir           => "${target_dir}/tophat2",
    option               => "--segment-length 25 -r 0 -p 6",
    source_ref           => "sra2fastq",
    bowtie2_index        => $bowtie2_index,
    transcript_gtf       => $transcript_gtf,
    transcript_gtf_index => $transcript_gtf_index,
    rename_bam           => 1,
    sh_direct            => 0,
    pbs                  => {
      "email"    => $email,
      "nodes"    => "1:ppn=6",
      "walltime" => "72",
      "mem"      => "30gb"
    },
  },
  sortbam => {
    class         => "Sortbam",
    perform       => 0,
    target_dir    => "${target_dir}/sortname",
    option        => "",
    source_ref    => "tophat2",
    sort_by_query => 1,
    sh_direct     => 0,
    pbs           => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "20gb"
    },
  },
  htseqcount => {
    class      => "HTSeqCount",
    perform    => 0,
    target_dir => "${target_dir}/htseqcount",
    option     => "",
    source_ref => "sortbam",
    gff_file   => $transcript_gtf,
    sh_direct  => 0,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  genetable => {
    class         => "CQSDatatable",
    perform       => 0,
    target_dir    => "${target_dir}/genetable",
    option        => "-p ENS --noheader -o ${task}_gene.count",
    source_ref    => "htseqcount",
    name_map_file => $hg19_map,
    cqs_tools     => $cqstools,
    sh_direct     => 1,
    pbs           => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "10",
      "mem"      => "10gb"
    },
  },
  dexseqcount => {
    class        => "DexseqCount",
    perform      => 0,
    target_dir   => "${target_dir}/dexseqcount",
    option       => "",
    source_ref   => "tophat2",
    gff_file     => $hg19_gff,
    dexseq_count => "/home/shengq1/pylibs/bin/dexseq_count.py",
    sh_direct    => 0,
    pbs          => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  exontable => {
    class         => "CQSDatatable",
    perform       => 0,
    target_dir    => "${target_dir}/exontable",
    option        => "-p ENS --noheader -o ${task}_exon.count",
    name_map_file => $hg19_map,
    source_ref    => "dexseqcount",
    cqs_tools     => $cqstools,
    sh_direct     => 1,
    pbs           => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "10",
      "mem"      => "10gb"
    },
  },
  varscan2 => {
    class           => "VarScan2::Mpileup2snp",
    perform         => 0,
    target_dir      => "${target_dir}/varscan2",
    option          => "--min-coverage 10",
    mpileup_options => "-q 20",
    java_option     => "-Xmx40g",
    source_ref      => "tophat2",
    fasta_file      => $fasta_file,
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
    option     => $annovar_param,
    source_ref => [ "varscan2", "\.vcf\$" ],
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
