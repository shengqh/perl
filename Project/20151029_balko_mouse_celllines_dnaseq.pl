#!/usr/bin/perl
use strict;
use warnings;

use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;
use Hash::Merge qw( merge );

my $target_dir = create_directory_or_die("/scratch/cqs/shengq1/dnaseq/20151029_balko_mouse_celllines");
my $workspace  = "/workspace/shengq1/dnaseq/";

my $email = "quanhu.sheng\@vanderbilt.edu";

my $cqstools      = "/home/shengq1/cqstools/CQS.Tools.exe";
my $glmvc         = "/home/shengq1/glmvc/glmvc.exe";
my $glmvc_1_3_6   = "/home/shengq1/glmvc_old/glmvc.x64.1.3.6/glmvc.exe";
my $mutect        = "/home/shengq1/local/bin/mutect-1.1.7.jar";
my $varscan2      = "/home/shengq1/local/bin/VarScan.v2.4.1.jar";
my $gatk_jar      = "/home/shengq1/local/bin/GATK/GenomeAnalysisTK.jar";
my $picard_jar    = "/scratch/cqs/shengq1/local/bin/picard/picard.jar";
my $conifer       = "/home/shengq1/pylibs/bin/conifer.py";
my $qc3_perl      = "/scratch/cqs/shengq1/local/bin/qc3/qc3.pl";
my $dexseq_script = "/home/shengq1/R/x86_64-pc-linux-gnu-library/3.2/DEXSeq/python_scripts/dexseq_count.py";

my $bwa_fasta        = "/scratch/cqs/shengq1/references/mm10_sorted_M/bwa_index_0.7.12/mm10.fa";
my $dbsnp            = "/scratch/cqs/shengq1/references/dbsnp/mm10/mouse_GRCm38_v142_M.vcf";
my $capture_bed      = "/scratch/cqs/shengq1/references/sureselect/S0276129_Mouse_All_Exon_V1/S0276129_mm10_All_Exon_V1_M.bed";
my $capture_slim_bed = "/scratch/cqs/shengq1/references/sureselect/S0276129_Mouse_All_Exon_V1/S0276129_mm10_All_Exon_V1_slim.bed";
my $gene_bed         = "/scratch/cqs/shengq1/dnaseq/20151029_balko_mouse_celllines/config/dusp4_tp53_myc.bed";
my $exclude_bed      = "/scratch/cqs/shengq1/dnaseq/20151029_balko_mouse_celllines/config/myc.bed";
my $dexseq_gff       = "/scratch/cqs/shengq1/references/ensembl_gtf/v75/Mus_musculus.GRCm38.75.M.dexseq.gff";
my $name_map_file    = "/scratch/cqs/shengq1/references/ensembl_gtf/v75/Mus_musculus.GRCm38.75.M.map";

my $annovar_protocol  = "refGene";
my $annovar_operation = "g";
my $annovar_param     = "-protocol ${annovar_protocol} -operation ${annovar_operation} --remove";
my $annovar_db        = "/scratch/cqs/shengq1/references/annovar/mm10db/";

#minimum quality score 10, minimum overlap 4 bases, remove reads with length less than 30
my $cutadapt_option = "-q 10 -O 4 -m 30";

my $glmvc_option = "--max_normal_percentage 0.01 --min_tumor_percentage 0.1 --min_tumor_read 5 --glm_pvalue 0.1 --exclude_bed $exclude_bed -r 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,X,M";
my $glmvc_option_pvalue =
  "--max_normal_percentage 0.01 --min_tumor_percentage 0.1 --min_tumor_read 5 --glm_pvalue 0.01 --glm_use_raw_pvalue --exclude_bed $exclude_bed -r 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,X,M";

my $cluster = "slurm";

my $groups = {
  "N07_DUSP4flox_MYC"             => [ "N04_DUSP4flox_LACZ", "N07_DUSP4flox_MYC" ],
  "N05_DUSP4flox_Trp53null1_LACZ" => [ "N04_DUSP4flox_LACZ", "N05_DUSP4flox_Trp53null1_LACZ" ],
  "N09_DUSP4flox_Trp53null3_MYC"  => [ "N04_DUSP4flox_LACZ", "N09_DUSP4flox_Trp53null3_MYC" ],
  "N13_DUSP4null_LACZ"            => [ "N04_DUSP4flox_LACZ", "N13_DUSP4null_LACZ" ],
  "N16_DUSP4null_MYC"             => [ "N04_DUSP4flox_LACZ", "N16_DUSP4null_MYC" ],
  "N15_DUSP4null_Trp53null3_LACZ" => [ "N04_DUSP4flox_LACZ", "N15_DUSP4null_Trp53null3_LACZ" ],
  "N17_DUSP4null_Trp53null1_MYC"  => [ "N04_DUSP4flox_LACZ", "N17_DUSP4null_Trp53null1_MYC" ],
};

my $wes = {
  task_name        => "WES",
  no_adapter_files => {
    "N04_DUSP4flox_LACZ" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20151029_balko_mouse_celllines/WES/data/3162-JMB-1_1_sequence.txt.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20151029_balko_mouse_celllines/WES/data/3162-JMB-1_2_sequence.txt.gz"
    ],
    "N05_DUSP4flox_Trp53null1_LACZ" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20151029_balko_mouse_celllines/WES/data/3162-JMB-2_1_sequence.txt.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20151029_balko_mouse_celllines/WES/data/3162-JMB-2_2_sequence.txt.gz"
    ],
    "N07_DUSP4flox_MYC" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20151029_balko_mouse_celllines/WES/data/3162-JMB-3_1_sequence.txt.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20151029_balko_mouse_celllines/WES/data/3162-JMB-3_2_sequence.txt.gz"
    ],
    "N09_DUSP4flox_Trp53null3_MYC" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20151029_balko_mouse_celllines/WES/data/3162-JMB-4_1_sequence.txt.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20151029_balko_mouse_celllines/WES/data/3162-JMB-4_2_sequence.txt.gz"
    ],
    "N13_DUSP4null_LACZ" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20151029_balko_mouse_celllines/WES/data/3162-JMB-5_1_sequence.txt.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20151029_balko_mouse_celllines/WES/data/3162-JMB-5_2_sequence.txt.gz"
    ],
    "N15_DUSP4null_Trp53null3_LACZ" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20151029_balko_mouse_celllines/WES/data/3162-JMB-6_1_sequence.txt.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20151029_balko_mouse_celllines/WES/data/3162-JMB-6_2_sequence.txt.gz"
    ],
    "N16_DUSP4null_MYC" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20151029_balko_mouse_celllines/WES/data/3162-JMB-7_1_sequence.txt.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20151029_balko_mouse_celllines/WES/data/3162-JMB-7_2_sequence.txt.gz"
    ],
    "N17_DUSP4null_Trp53null1_MYC" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20151029_balko_mouse_celllines/WES/data/3162-JMB-8_1_sequence.txt.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20151029_balko_mouse_celllines/WES/data/3162-JMB-8_2_sequence.txt.gz"
    ],
  },
  adapter_files => {
    "N08_DUSP4flox_Trp53null1_MYC" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20151029_balko_mouse_celllines/WES/data/3162-JMB-9_1_merged_sequence.txt.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20151029_balko_mouse_celllines/WES/data/3162-JMB-9_2_merged_sequence.txt.gz"
    ],
    "N18_DUSP4null_Trp53null3_MYC" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20151029_balko_mouse_celllines/WES/data/3162-JMB-10_1_merged_sequence.txt.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20151029_balko_mouse_celllines/WES/data/3162-JMB-10_2_merged_sequence.txt.gz"
    ],
  },
  groups => merge(
    $groups,
    {
      "N08_DUSP4flox_Trp53null1_MYC" => [ "N04_DUSP4flox_LACZ", "N08_DUSP4flox_Trp53null1_MYC" ],
      "N18_DUSP4null_Trp53null3_MYC" => [ "N04_DUSP4flox_LACZ", "N18_DUSP4null_Trp53null3_MYC" ],
    }
  ),
  depthgroups => {
    "WES" => [
      "N04_DUSP4flox_LACZ",           "N05_DUSP4flox_Trp53null1_LACZ", "N07_DUSP4flox_MYC",             "N08_DUSP4flox_Trp53null1_MYC",
      "N09_DUSP4flox_Trp53null3_MYC", "N13_DUSP4null_LACZ",            "N15_DUSP4null_Trp53null3_LACZ", "N16_DUSP4null_MYC",
      "N17_DUSP4null_Trp53null1_MYC", "N18_DUSP4null_Trp53null3_MYC"
    ]
  },
  cutadapt         => 1,
  capture_bed      => $capture_bed,
  capture_slim_bed => $capture_slim_bed
};

my $wgs = {
  task_name        => "WGS",
  no_adapter_files => {

  },
  adapter_files => {
    "N04_DUSP4flox_LACZ" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20151029_balko_mouse_celllines/WGS/data/3162-JB-1_1_sequence.txt.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20151029_balko_mouse_celllines/WGS/data/3162-JB-1_2_sequence.txt.gz"
    ],
    "N05_DUSP4flox_Trp53null1_LACZ" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20151029_balko_mouse_celllines/WGS/data/3162-JB-2_1_sequence.txt.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20151029_balko_mouse_celllines/WGS/data/3162-JB-2_2_sequence.txt.gz"
    ],
    "N07_DUSP4flox_MYC" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20151029_balko_mouse_celllines/WGS/data/3162-JB-3_1_sequence.txt.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20151029_balko_mouse_celllines/WGS/data/3162-JB-3_2_sequence.txt.gz"
    ],
    "N09_DUSP4flox_Trp53null3_MYC" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20151029_balko_mouse_celllines/WGS/data/3162-JB-4_1_sequence.txt.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20151029_balko_mouse_celllines/WGS/data/3162-JB-4_2_sequence.txt.gz"
    ],
    "N13_DUSP4null_LACZ" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20151029_balko_mouse_celllines/WGS/data/3162-JB-5_1_sequence.txt.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20151029_balko_mouse_celllines/WGS/data/3162-JB-5_2_sequence.txt.gz"
    ],
    "N15_DUSP4null_Trp53null3_LACZ" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20151029_balko_mouse_celllines/WGS/data/3162-JB-6_1_sequence.txt.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20151029_balko_mouse_celllines/WGS/data/3162-JB-6_2_sequence.txt.gz"
    ],
    "N16_DUSP4null_MYC" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20151029_balko_mouse_celllines/WGS/data/3162-JB-7_1_sequence.txt.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20151029_balko_mouse_celllines/WGS/data/3162-JB-7_2_sequence.txt.gz"
    ],
    "N17_DUSP4null_Trp53null1_MYC" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20151029_balko_mouse_celllines/WGS/data/3162-JB-8_1_sequence.txt.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20151029_balko_mouse_celllines/WGS/data/3162-JB-8_2_sequence.txt.gz"
    ],
  },
  groups      => $groups,
  depthgroups => {
    "WGS" => [
      "N04_DUSP4flox_LACZ", "N05_DUSP4flox_Trp53null1_LACZ", "N07_DUSP4flox_MYC", "N09_DUSP4flox_Trp53null3_MYC",
      "N13_DUSP4null_LACZ", "N15_DUSP4null_Trp53null3_LACZ", "N16_DUSP4null_MYC", "N17_DUSP4null_Trp53null1_MYC"
    ]
  },
  cutadapt => 1
};

my @datasets = ( $wes, $wgs );

for my $dataset (@datasets) {
  my $dataset_dir = create_directory_or_die( "${target_dir}/" . $dataset->{task_name} );
  my $config      = {
    general          => { task_name => $dataset->{task_name} },
    adapter_files    => $dataset->{adapter_files},
    no_adapter_files => $dataset->{no_adapter_files},
    groups           => $dataset->{groups},
    fastqc           => {
      class      => "QC::FastQC",
      perform    => 1,
      target_dir => "${dataset_dir}/fastqc",
      option     => "",
      source_ref => [ "adapter_files", "no_adapter_files" ],
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
      target_dir => "${target_dir}/" . $dataset->{task_name} . "/fastqc",
      option     => "",
      cluster    => $cluster,
      cqstools   => $cqstools,
      pbs        => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "2",
        "mem"      => "10gb"
      },
    }
  };

  my @step1  = ("fastqc");
  my @step2  = ("fastqc_summary");
  my @step3  = ();
  my $source = "files";
  if ( defined $dataset->{cutadapt} && $dataset->{cutadapt} ) {
    $config = merge(
      $config,
      {
        cutadapt => {
          class      => "Trimmer::Cutadapt",
          perform    => 1,
          target_dir => "${target_dir}/" . $dataset->{task_name} . "/cutadapt",
          option     => $cutadapt_option,
          source_ref => "adapter_files",
          adapter    => "AGATCGGAAGAGC",
          extension  => "_clipped.fastq.gz",
          sh_direct  => 1,
          cluster    => $cluster,
          pairend    => 1,
          pbs        => {
            "email"    => $email,
            "nodes"    => "1:ppn=1",
            "walltime" => "24",
            "mem"      => "20gb"
          },
        },
        fastqlen => {
          class      => "CQS::FastqLen",
          perform    => 1,
          target_dir => "${target_dir}/" . $dataset->{task_name} . "/fastqlen",
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
        cutadapt_fastqc => {
          class      => "QC::FastQC",
          perform    => 1,
          target_dir => "${target_dir}/" . $dataset->{task_name} . "/cutadapt_fastqc",
          option     => "",
          source_ref => "cutadapt",
          sh_direct  => 1,
          cluster    => $cluster,
          pbs        => {
            "email"    => $email,
            "nodes"    => "1:ppn=2",
            "walltime" => "2",
            "mem"      => "40gb"
          },
        },
        cutadapt_fastqc_summary => {
          class      => "QC::FastQCSummary",
          perform    => 1,
          target_dir => "${target_dir}/" . $dataset->{task_name} . "/cutadapt_fastqc",
          option     => "",
          cluster    => $cluster,
          cqstools   => $cqstools,
          pbs        => {
            "email"    => $email,
            "nodes"    => "1:ppn=1",
            "walltime" => "2",
            "mem"      => "10gb"
          },
        }
      }
    );
    $source = [ "no_adapter_files", "cutadapt" ];
    push @step1, ( "cutadapt", "fastqlen", "cutadapt_fastqc" );
    push @step2, ("cutadapt_fastqc_summary");

    #performTask( $config, "cutadapt" );
  }

  $config = merge(
    $config,
    {
      bwa => {
        class      => "Alignment::BWA",
        perform    => 1,
        target_dir => "${target_dir}/" . $dataset->{task_name} . "/bwa",
        option     => "",
        bwa_index  => $bwa_fasta,
        source_ref => $source,
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
        class            => "GATK::Refine",
        perform          => 1,
        target_dir       => "${target_dir}/" . $dataset->{task_name} . "/bwa_refine",
        option           => "-Xmx40g",
        gatk_option      => "--fix_misencoded_quality_scores",
        fasta_file       => $bwa_fasta,
        source_ref       => "bwa",
        indel_vcf_files  => [$dbsnp],
        known_vcf_files  => [$dbsnp],
        gatk_jar         => $gatk_jar,
        picard_jar       => $picard_jar,
        slim_print_reads => 1,
        sh_direct        => 0,
        sorted           => 1,
        pbs              => {
          "email"    => $email,
          "nodes"    => "1:ppn=8",
          "walltime" => "240",
          "mem"      => "40gb"
        },
      },
      bwa_refine_genes_bam => {
        class      => "Samtools::View",
        perform    => 1,
        target_dir => "${target_dir}/" . $dataset->{task_name} . "/bwa_refine_genes_bam",
        option     => "-h -b -L $gene_bed",
        extension  => ".bam",
        source_ref => "bwa_refine",
        sh_direct  => 1,
        sorted     => 1,
        pbs        => {
          "email"    => $email,
          "nodes"    => "1:ppn=1",
          "walltime" => "24",
          "mem"      => "10gb"
        },
      },
      muTect => {
        class        => "GATK::MuTect",
        perform      => 1,
        target_dir   => "${target_dir}/" . $dataset->{task_name} . "/muTect",
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
        target_dir => "${target_dir}/" . $dataset->{task_name} . "/muTect_annovar",
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
      glmvc_noMYC_fdr => {
        class             => "Variants::GlmvcCall",
        perform           => 1,
        target_dir        => "${target_dir}/" . $dataset->{task_name} . "/glmvc_noMYC_fdr",
        option            => $glmvc_option,
        source_type       => "BAM",
        source_ref        => "bwa_refine",
        groups_ref        => "groups",
        fasta_file        => $bwa_fasta,
        annovar_buildver  => "mm10",
        annovar_protocol  => $annovar_protocol,
        annovar_operation => $annovar_operation,
        annovar_db        => $annovar_db,
        sh_direct         => 0,
        execute_file      => $glmvc,
        pbs               => {
          "email"    => $email,
          "nodes"    => "1:ppn=6",
          "walltime" => "72",
          "mem"      => "40gb"
        },
      },
      glmvc_noMYC_rawpvalue => {
        class   => "Variants::GlmvcCall",
        perform => 1,

        #target_dir        => "${target_dir}/" . $dataset->{task_name} . "/glmvc_noMYC_rawpvalue",
        target_dir        => "/workspace/shengq1/dnaseq/" . $dataset->{task_name} . "/glmvc_noMYC_rawpvalue",
        option            => $glmvc_option_pvalue,
        source_type       => "BAM",
        source_ref        => "bwa_refine",
        groups_ref        => "groups",
        fasta_file        => $bwa_fasta,
        annovar_buildver  => "mm10",
        annovar_protocol  => $annovar_protocol,
        annovar_operation => $annovar_operation,
        annovar_db        => $annovar_db,
        sh_direct         => 1,
        execute_file      => $glmvc_1_3_6,
        pbs               => {
          "email"    => $email,
          "nodes"    => "1:ppn=12",
          "walltime" => "72",
          "mem"      => "40gb"
        },
      },
      glmvc_noMYC_v1_3_6_fdr => {
        class             => "Variants::GlmvcCall",
        perform           => 0,
        target_dir        => "${target_dir}/" . $dataset->{task_name} . "/glmvc_noMYC_v1_3_6_fdr",
        option            => $glmvc_option,
        source_type       => "BAM",
        source_ref        => "bwa_refine",
        groups_ref        => "groups",
        fasta_file        => $bwa_fasta,
        annovar_buildver  => "mm10",
        annovar_protocol  => $annovar_protocol,
        annovar_operation => $annovar_operation,
        annovar_db        => $annovar_db,
        sh_direct         => 0,
        execute_file      => $glmvc_1_3_6,
        pbs               => {
          "email"    => $email,
          "nodes"    => "1:ppn=6",
          "walltime" => "72",
          "mem"      => "40gb"
        },
      },
      glmvc_noMYC_v1_3_6_rawpvalue => {
        class             => "Variants::GlmvcCall",
        perform           => 0,
        target_dir        => "${target_dir}/" . $dataset->{task_name} . "/glmvc_noMYC_v1_3_6_rawpvalue",
        option            => $glmvc_option_pvalue,
        source_type       => "BAM",
        source_ref        => "bwa_refine",
        groups_ref        => "groups",
        fasta_file        => $bwa_fasta,
        annovar_buildver  => "mm10",
        annovar_protocol  => $annovar_protocol,
        annovar_operation => $annovar_operation,
        annovar_db        => $annovar_db,
        sh_direct         => 0,
        execute_file      => $glmvc_1_3_6,
        pbs               => {
          "email"    => $email,
          "nodes"    => "1:ppn=6",
          "walltime" => "72",
          "mem"      => "40gb"
        },
      },
      glmvc_noMYC_table => {
        class   => "Variants::GlmvcTable",
        perform => 1,

        #target_dir   => "${target_dir}/" . $dataset->{task_name} . "/glmvc_noMYC_table",
        target_dir   => "/workspace/shengq1/dnaseq/" . $dataset->{task_name} . "/glmvc_noMYC_table",
        option       => "",
        source_ref   => [ "glmvc_noMYC_rawpvalue", "annotation.tsv" ],
        sh_direct    => 1,
        execute_file => $glmvc,
        pbs          => {
          "email"    => $email,
          "nodes"    => "1:ppn=8",
          "walltime" => "72",
          "mem"      => "40gb"
        },
      }
    }
  );

  push @step1, ( "bwa", "bwa_refine", "bwa_refine_genes_bam" );
  push @step2, ( "muTect", "muTect_annovar", "glmvc_noMYC_rawpvalue", "glmvc_noMYC_table" );

  if ( defined $dataset->{capture_bed} ) {
    $config = merge(
      $config,
      {
        cnmops => {
          class       => "CNV::cnMops",
          perform     => 1,
          target_dir  => $workspace . $dataset->{task_name} . "/cnmops",
          option      => "",
          source_ref  => "bwa_refine",
          bedfile     => $dataset->{capture_bed},
          refnames    => ["N04_DUSP4flox_LACZ"],
          pairmode    => "paired",
          isbamsorted => 1,
          sh_direct   => 1,
          cqstools    => $cqstools,
          pbs         => {
            "email"    => $email,
            "nodes"    => "1:ppn=1",
            "walltime" => "72",
            "mem"      => "40gb"
          },
        },
        cnmops_depth => {
          class          => "Visualization::Depth",
          perform        => 1,
          target_dir     => "${target_dir}/" . $dataset->{task_name} . "/cnmops_depth",
          option         => "-h -b -L " . "${target_dir}/" . $dataset->{task_name} . "/cnmops/result/" . $dataset->{task_name} . ".call.bed",
          source_ref     => [ "cnmops", ".call.tsv.cnvr\$" ],
          groups_ref     => $dataset->{depthgroups},
          bam_files_ref  => "bwa_refine",
          cnvr_files_ref => [ "cnmops", ".call.tsv.cnvr\$" ],
          pbs            => {
            "email"    => $email,
            "nodes"    => "1:ppn=1",
            "walltime" => "24",
            "mem"      => "10gb"
          },
        },
      }
    );
    push @step2, ( "cnmops", "cnmops_depth" );
  }
  else {
    $config = merge(
      $config,
      {
        cnmops => {
          class         => "CNV::cnMops",
          perform       => 1,
          target_dir    => $workspace . $dataset->{task_name} . "/cnmops",
          option        => "",
          source_ref    => "bwa_refine",
          pairmode      => "paired",
          isbamsorted   => 1,
          ref_seq_names => [ "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "X", "Y", "M" ],
          cqstools      => $cqstools,
          sh_direct     => 1,
          pbs           => {
            "email"    => $email,
            "nodes"    => "1:ppn=1",
            "walltime" => "72",
            "mem"      => "40gb"
          },
        },
        glmvc_WES_validation => {
          class            => "Variants::GlmvcValidate",
          perform          => 0,
          target_dir       => "${target_dir}/" . $dataset->{task_name} . "/glmvc_WES_validation",
          option           => "",
          source_type      => "BAM",
          source_ref       => "bwa_refine",
          validation_files => $target_dir . "/WES/glmvc_noMYC_table/result/WES.tsv",
          groups_ref       => "groups",
          fasta_file       => $bwa_fasta,
          sh_direct        => 1,
          execute_file     => $glmvc,
          pbs              => {
            "email"    => $email,
            "nodes"    => "1:ppn=1",
            "walltime" => "72",
            "mem"      => "40gb"
          },
        }
      }
    );
    push @step2, ( "cnmops", "glmvc_WES_validation" );
  }

  $config->{varscan2_copynumber} = {
    class           => "VarScan2::Copynumber",
    perform         => 1,
    target_dir      => "${target_dir}/" . $dataset->{task_name} . "/varscan2_copynumber",
    option          => "",
    java_option     => "-Xmx40g",
    mpileup_options => "-q 1",
    call_options    => "--amp-threshold 0.3 --del-threshold 0.3",
    source_ref      => "bwa_refine",
    VarScan2_jar    => $varscan2,
    groups          => $dataset->{groups},
    fasta_file      => $bwa_fasta,
    sh_direct       => 0,
    pbs             => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  };
  push @step2, ("varscan2_copynumber");

  $config->{sequencetask} = {
    class      => "CQS::SequenceTask",
    perform    => 1,
    target_dir => "${target_dir}/" . $dataset->{task_name} . "/sequencetask",
    option     => "",
    source     => {
      step1 => \@step1,
      step2 => \@step2,
      step3 => \@step3,
    },
    sh_direct => 0,
    pbs       => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  };

  #  performConfig($config);
  #performTask($config, "glmvc_noMYC_fdr");
  #performTask( $config, "glmvc_noMYC_rawpvalue" );

  #performTask($config, "glmvc_noMYC_v1_3_6_rawpvalue");
  #performTask( $config, "glmvc_noMYC_table" );

  performTask( $config, "cnmops" );
  if ( defined $dataset->{capture_bed} ) {
    performTask( $config, "cnmops_depth" );
  }
}

1;

