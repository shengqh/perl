#!/usr/bin/perl
use strict;
use warnings;

use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;
use Hash::Merge qw( merge );

my $dataset_dir = create_directory_or_die("/scratch/cqs/shengq1/pipelines/exomeseq_somatic_mutation");

my $email = "quanhu.sheng\@vanderbilt.edu";

my $cqstools   = "/home/shengq1/cqstools/cqstools.exe";
my $mutect     = "/home/shengq1/local/bin/mutect-1.1.7.jar";
my $gatk_jar   = "/home/shengq1/local/bin/GATK/GenomeAnalysisTK.jar";
my $picard_jar = "/scratch/cqs/shengq1/local/bin/picard/picard.jar";

my $bwa_fasta     = "/scratch/cqs/shengq1/references/mm10_sorted_M/bwa_index_0.7.12/mm10.fa";
my $dbsnp         = "/scratch/cqs/shengq1/references/dbsnp/mm10/mouse_GRCm38_v142_M.vcf";
my $indel_vcf     = "/scratch/cqs/shengq1/references/mm10/mpg.v5.indels.pass.reordered.vcf";
my $capture_bed   = "/scratch/cqs/shengq1/references/sureselect/S0276129_Mouse_All_Exon_V1/S0276129_mm10_All_Exon_V1_M.bed";
my $name_map_file = "/scratch/cqs/shengq1/references/ensembl_gtf/v75/Mus_musculus.GRCm38.75.M.map";

my $annovar_protocol  = "refGene";
my $annovar_operation = "g";
my $annovar_param     = "-protocol ${annovar_protocol} -operation ${annovar_operation} --remove";
my $annovar_db        = "/scratch/cqs/shengq1/references/annovar/mm10db/";

#minimum quality score 10, minimum overlap 4 bases, remove reads with length less than 30
my $cutadapt_option = "-q 10 -O 4 -m 30";

my $cluster = "slurm";

my $config = {
  general => {
    task_name => "WES"
  },
  files => {
    "N04" => [
      "/gpfs21/scratch/cqs/shengq1/pipelines/exomeseq_somatic_mutation/data/N04_DUSP4flox_LACZ_19.1.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/pipelines/exomeseq_somatic_mutation/data/N04_DUSP4flox_LACZ_19.2.fastq.gz"
    ],
    "N05" => [
      "/gpfs21/scratch/cqs/shengq1/pipelines/exomeseq_somatic_mutation/data/N05_DUSP4flox_Trp53null1_LACZ_19.1.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/pipelines/exomeseq_somatic_mutation/data/N05_DUSP4flox_Trp53null1_LACZ_19.2.fastq.gz"
    ],
  },
  groups => {
    "N05_vs_N04" => [ "N04", "N05" ],
  },

  fastqc => {
    class      => "QC::FastQC",
    perform    => 1,
    target_dir => "${dataset_dir}/fastqc",
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
    perform    => 1,
    target_dir => "${dataset_dir}/fastqc",
    option     => "",
    cluster    => $cluster,
    cqstools   => $cqstools,
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
    target_dir => "${dataset_dir}/bwa",
    option     => "",
    bwa_index  => $bwa_fasta,
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
    class             => "GATK::Refine",
    perform           => 1,
    init_command      => "setpkgs -a java_1.8",
    target_dir        => "${dataset_dir}/bwa_refine",
    option            => "-Xmx40g",
    gatk_option       => "--fix_misencoded_quality_scores",
    fasta_file        => $bwa_fasta,
    source_ref        => "bwa",
    indel_vcf_files   => [$indel_vcf],
    known_vcf_files   => [$dbsnp],
    gatk_jar          => $gatk_jar,
    picard_jar        => $picard_jar,
    indel_realignment => 1,
    slim_print_reads  => 1,
    sh_direct         => 0,
    sorted            => 1,
    pbs               => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "240",
      "mem"      => "40gb"
    },
  },
  muTect => {
    class        => "GATK::MuTect",
    perform      => 1,
    init_command => "setpkgs -a java",
    target_dir   => "${dataset_dir}/muTect",
    option       => "--min_qscore 20 --filter_reads_with_N_cigar",
    java_option  => "-Xmx40g",
    source_ref   => "bwa_refine",
    groups_ref   => "groups",
    fasta_file   => $bwa_fasta,
    dbsnp_file   => $dbsnp,
    bychromosome => 0,
    sh_direct    => 0,
    muTect_jar   => $mutect,
    pbs          => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "240",
      "mem"      => "40gb"
    },
  },
  muTect_annovar => {
    class      => "Annotation::Annovar",
    perform    => 1,
    target_dir => "${dataset_dir}/muTect_annovar",
    option     => $annovar_param,
    source_ref => [ "muTect", ".pass.vcf\$" ],
    annovar_db => $annovar_db,
    buildver   => "mm10",
    cqstools   => $cqstools,
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
    target_dir => "${dataset_dir}/sequencetask",
    option     => "",
    source     => {
      step1 => [ "fastqc",         "bwa",    "bwa_refine" ],
      step2 => [ "fastqc_summary", "muTect", "muTect_annovar" ],
    },
    sh_direct => 0,
    pbs       => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  }
};

#performConfig($config);
performTask( $config, "bwa_refine" );

1;

