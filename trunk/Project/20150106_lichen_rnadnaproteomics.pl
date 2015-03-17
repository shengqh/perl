#!/usr/bin/perl
use strict;
use warnings;

use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;

my $target_dir = "/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics";
my $cqstools   = "/home/shengq1/cqstools/CQS.Tools.exe";
my $email      = "quanhu.sheng\@vanderbilt.edu";

my $bwa_fasta      = "/scratch/cqs/shengq1/references/hg19_16569_M/bwa_index_0.7.12/hg19_16569_M.fa";
my $star_index     = "/scratch/cqs/shengq1/references/hg19_16569_M/STAR_index_v37.75_2.4.0j";
my $transcript_gtf = "/scratch/cqs/shengq1/references/ensembl_gtf/v75/Homo_sapiens.GRCh37.75.M.gtf";
my $name_map_file  = "/scratch/cqs/shengq1/references/ensembl_gtf/v75/Homo_sapiens.GRCh37.75.M.map";
my $cosmic         = "/scratch/cqs/shengq1/references/cosmic/cosmic_v71_hg19_16569_MT.vcf";
my $dbsnp          = "/scratch/cqs/shengq1/references/dbsnp/human_GRCh37_v142_16569_M.vcf";
my $annovar_param  = "-protocol refGene,snp138,cosmic70 -operation g,f,f --remove";
my $annovar_db     = "/scratch/cqs/shengq1/references/annovar/humandb/";
my $gatk_jar       = "/home/shengq1/local/bin/GATK/GenomeAnalysisTK.jar";
my $picard_jar     = "/scratch/cqs/shengq1/local/bin/picard/picar.jar";
my $mutect_jar     = "/home/shengq1/local/bin/muTect-1.1.4.jar";

my $cluster = "slurm";

my $config = {
  general => { task_name => "lichen" },

  files => {
    "DNA_273_SRR926701" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/DNA_273_SRR926701.sra"],
    "DNA_273_SRR926707" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/DNA_273_SRR926707.sra"],
    "DNA_273_SRR926711" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/DNA_273_SRR926711.sra"],
    "DNA_273_SRR926713" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/DNA_273_SRR926713.sra"],
    "DNA_273_SRR926719" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/DNA_273_SRR926719.sra"],
    "DNA_273_SRR926721" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/DNA_273_SRR926721.sra"],
    "DNA_273_SRR926727" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/DNA_273_SRR926727.sra"],
    "DNA_273_SRR926729" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/DNA_273_SRR926729.sra"],
    "DNA_283_SRR926702" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/DNA_283_SRR926702.sra"],
    "DNA_283_SRR926706" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/DNA_283_SRR926706.sra"],
    "DNA_283_SRR926710" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/DNA_283_SRR926710.sra"],
    "DNA_283_SRR926714" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/DNA_283_SRR926714.sra"],
    "DNA_283_SRR926718" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/DNA_283_SRR926718.sra"],
    "DNA_283_SRR926722" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/DNA_283_SRR926722.sra"],
    "DNA_283_SRR926726" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/DNA_283_SRR926726.sra"],
    "DNA_283_SRR926730" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/DNA_283_SRR926730.sra"],
    "DNA_311_SRR926703" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/DNA_311_SRR926703.sra"],
    "DNA_311_SRR926705" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/DNA_311_SRR926705.sra"],
    "DNA_311_SRR926709" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/DNA_311_SRR926709.sra"],
    "DNA_311_SRR926715" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/DNA_311_SRR926715.sra"],
    "DNA_311_SRR926717" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/DNA_311_SRR926717.sra"],
    "DNA_311_SRR926723" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/DNA_311_SRR926723.sra"],
    "DNA_311_SRR926725" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/DNA_311_SRR926725.sra"],
    "DNA_311_SRR926731" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/DNA_311_SRR926731.sra"],
    "DNA_321_SRR926704" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/DNA_321_SRR926704.sra"],
    "DNA_321_SRR926708" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/DNA_321_SRR926708.sra"],
    "DNA_321_SRR926712" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/DNA_321_SRR926712.sra"],
    "DNA_321_SRR926716" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/DNA_321_SRR926716.sra"],
    "DNA_321_SRR926720" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/DNA_321_SRR926720.sra"],
    "DNA_321_SRR926724" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/DNA_321_SRR926724.sra"],
    "DNA_321_SRR926728" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/DNA_321_SRR926728.sra"],
    "DNA_321_SRR926732" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/DNA_321_SRR926732.sra"],
    "RNA_273_SRR926693" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/RNA_273_SRR926693.sra"],
    "RNA_273_SRR926694" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/RNA_273_SRR926694.sra"],
    "RNA_273_SRR926695" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/RNA_273_SRR926695.sra"],
    "RNA_273_SRR926696" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/RNA_273_SRR926696.sra"],
    "RNA_283_SRR926685" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/RNA_283_SRR926685.sra"],
    "RNA_283_SRR926686" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/RNA_283_SRR926686.sra"],
    "RNA_283_SRR926687" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/RNA_283_SRR926687.sra"],
    "RNA_283_SRR926688" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/RNA_283_SRR926688.sra"],
    "RNA_311_SRR926697" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/RNA_311_SRR926697.sra"],
    "RNA_311_SRR926698" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/RNA_311_SRR926698.sra"],
    "RNA_311_SRR926699" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/RNA_311_SRR926699.sra"],
    "RNA_311_SRR926700" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/RNA_311_SRR926700.sra"],
    "RNA_321_SRR926689" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/RNA_321_SRR926689.sra"],
    "RNA_321_SRR926690" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/RNA_321_SRR926690.sra"],
    "RNA_321_SRR926691" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/RNA_321_SRR926691.sra"],
    "RNA_321_SRR926692" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/raw/RNA_321_SRR926692.sra"],
  },
  sra2fastq => {
    class      => "SRA::FastqDump",
    perform    => 0,
    ispaired   => 1,
    target_dir => "${target_dir}/FastqDump",
    option     => "",
    source_ref => "files",
    sh_direct  => 0,

    #cluster => "slurm",
    pbs => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "10",
      "mem"      => "10gb"
    },
  }
};

performConfig($config);

my $rna_config = {
  general => { task_name => "lichen" },

  files => {
    "RNA_273" => [
      "/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/FastqDump/result/RNA_273_1.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/FastqDump/result/RNA_273_2.fastq.gz"
    ],
    "RNA_283" => [
      "/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/FastqDump/result/RNA_283_1.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/FastqDump/result/RNA_283_2.fastq.gz"
    ],
    "RNA_311" => [
      "/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/FastqDump/result/RNA_311_1.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/FastqDump/result/RNA_311_2.fastq.gz"
    ],
    "RNA_321" => [
      "/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/FastqDump/result/RNA_321_1.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/FastqDump/result/RNA_321_2.fastq.gz"
    ],
  },
  groups => {
    "RNA_GLC" => [ "RNA_283", "RNA_273" ],
    "RNA_DZC" => [ "RNA_321", "RNA_311" ],
  },
  fastqc => {
    class      => "QC::FastQC",
    perform    => 0,
    target_dir => "${target_dir}/fastqc",
    option     => "",
    source_ref => "files",
    sh_direct  => 1,
    cluster    => $cluster,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=2",
      "walltime" => "2",
      "mem"      => "40gb"
    },
  },
  fastqc_summary => {
    class      => "QC::FastQCSummary",
    perform    => 0,
    target_dir => "${target_dir}/fastqc",
    option     => "",
    cluster    => $cluster,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "2",
      "mem"      => "10gb"
    },
  },
  star => {
    class      => "Alignment::STAR",
    perform    => 1,
    target_dir => "${target_dir}/rna_star",
    option     => "",
    source_ref => "files",
    genome_dir => $star_index,
    sh_direct  => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "30gb"
    },
  },
};

performConfig($rna_config);

my $dna_config = {
  general => { task_name => "lichen" },

  files => {
    "DNA_273" => [
      "/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/FastqDump/result/DNA_273_1.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/FastqDump/result/DNA_273_2.fastq.gz"
    ],
    "DNA_283" => [
      "/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/FastqDump/result/DNA_283_1.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/FastqDump/result/DNA_283_2.fastq.gz"
    ],
    "DNA_311" => [
      "/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/FastqDump/result/DNA_311_1.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/FastqDump/result/DNA_311_2.fastq.gz"
    ],
    "DNA_321" => [
      "/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/FastqDump/result/DNA_321_1.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20150106-lichen-rnadnaproteomics/FastqDump/result/DNA_321_2.fastq.gz"
    ],
  },
  groups => {
    "DNA_GLC" => [ "DNA_283", "DNA_273" ],
    "DNA_DZC" => [ "DNA_321", "DNA_311" ],
  },
  fastqc => {
    class      => "QC::FastQC",
    perform    => 0,
    target_dir => "${target_dir}/dna_fastqc",
    option     => "",
    source_ref => "files",
    sh_direct  => 1,
    cluster    => $cluster,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=2",
      "walltime" => "2",
      "mem"      => "40gb"
    },
  },
  fastqc_summary => {
    class      => "QC::FastQCSummary",
    perform    => 0,
    target_dir => "${target_dir}/dna_fastqc",
    option     => "",
    cluster    => $cluster,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "2",
      "mem"      => "10gb"
    },
  },
  bwa => {
    class      => "Alignment::BWA",
    perform    => 1,
    target_dir => "${target_dir}/dna_bwa",
    option     => "-T 15",
    fasta_file => $bwa_fasta,
    source_ref => "files",
    picard_jar => $picard_jar,
    sh_direct  => 0,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  bwa_refine => {
    class        => "GATK::Refine",
    perform      => 1,
    target_dir   => "${target_dir}/dna_bwa_refine",
    option       => "-Xmx40g",
    gatk_option  => "--fix_misencoded_quality_scores",
    fasta_file   => $bwa_fasta,
    source_ref   => "bwa",
    thread_count => 8,
    vcf_files    => [$dbsnp],
    gatk_jar     => $gatk_jar,
    picard_jar   => $picard_jar,
    sh_direct    => 1,
    pbs          => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  dna_muTect => {
    class       => "GATK::MuTect",
    perform     => 1,
    target_dir  => "${target_dir}/dna_muTect",
    option      => "-nt 8",
    source_ref  => "bwa_refine",
    groups_ref  => "groups",
    java_option => "-Xmx40g",
    fasta_file  => $bwa_fasta,
    cosmic_file => $cosmic,
    dbsnp_file  => $dbsnp,
    sh_direct   => 1,
    muTect_jar  => $mutect_jar,
    pbs         => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  dna_muTect_annovar => {
    class      => "Annotation::Annovar",
    perform    => 1,
    target_dir => "${target_dir}/dna_muTect",
    option     => $annovar_param,
    source_ref => [ "dna_muTect", "\.vcf\$" ],
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
  dna_snpindel => {
    class       => "GATK::SNPIndel",
    perform     => 1,
    target_dir  => "${target_dir}/dna_SNPindel",
    option      => "-l INFO -G Standard -stand_call_conf 50.0 -stand_emit_conf 10.0 -dcov 200 -nct 8",
    source_ref  => "bwa_refine",
    java_option => "-Xmx40g",
    fasta_file  => $bwa_fasta,
    vcf_files   => [$dbsnp],
    gatk_jar    => $gatk_jar,
    pbs         => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  dna_snpindel_annovar => {
    class      => "Annotation::Annovar",
    perform    => 1,
    target_dir => "${target_dir}/dna_SNPindel",
    source_ref => "dna_snpindel",
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

performConfig($dna_config);

1
