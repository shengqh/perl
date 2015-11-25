#!/usr/bin/perl
use strict;
use warnings;

use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;
use Hash::Merge qw( merge );

my $target_dir = "/scratch/cqs/shengq1/dnaseq/20151029_balko_mouse_celllines";
my $email      = "quanhu.sheng\@vanderbilt.edu";

my $cqstools   = "/home/shengq1/cqstools/CQS.Tools.exe";
my $glmvc      = "/home/shengq1/glmvc/glmvc.exe";
my $mutect     = "/home/shengq1/local/bin/mutect-1.1.7.jar";
my $gatk_jar   = "/home/shengq1/local/bin/GATK/GenomeAnalysisTK.jar";
my $picard_jar = "/scratch/cqs/shengq1/local/bin/picard/picard.jar";
my $conifer    = "/home/shengq1/pylibs/bin/conifer.py";
my $qc3_perl   = "/scratch/cqs/shengq1/local/bin/qc3/qc3.pl";

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
  task_name => "WES",
  files     => {
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
  groups => $groups,

  capture_bed      => $capture_bed,
  capture_slim_bed => $capture_slim_bed
};

my $wgs = {
  task_name => "WGS",
  files     => {
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
  groups   => $groups,
  cutadapt => 1
};

my @datasets = ( $wes, $wgs );

for my $dataset (@datasets) {
  my $config = {
    general => { task_name => $dataset->{task_name} },
    files   => $dataset->{files},
    groups  => $dataset->{groups},
    fastqc  => {
      class      => "QC::FastQC",
      perform    => 1,
      target_dir => "${target_dir}/" . $dataset->{task_name} . "/fastqc",
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

  my @individuals = ("fastqc");
  my @all         = ("fastqc_summary");
  my $source      = "files";
  if ( defined $dataset->{cutadapt} && $dataset->{cutadapt} ) {
    $config = merge(
      $config,
      {
        cutadapt => {
          class      => "Trimmer::Cutadapt",
          perform    => 1,
          target_dir => "${target_dir}/" . $dataset->{task_name} . "/cutadapt",
          option     => $cutadapt_option,
          source_ref => "files",
          adapter    => "AGATCGGAAGAGC",
          extension  => "_clipped.fastq.gz",
          sh_direct  => 1,
          cluster    => $cluster,
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
    $source = "cutadapt";
    push @individuals, ( "cutadapt", "fastqlen", "cutadapt_fastqc" );
    push @all, ("cutadapt_fastqc_summary");
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
      gene_bam => {
        class      => "Samtools::View",
        perform    => 1,
        target_dir => "${target_dir}/" . $dataset->{task_name} . "/gene_bam",
        option     => "-h -b -L $gene_bed",
        extension  => ".bam",
        source_ref => "bwa",
        sh_direct  => 1,
        sorted     => 1,
        pbs        => {
          "email"    => $email,
          "nodes"    => "1:ppn=1",
          "walltime" => "24",
          "mem"      => "10gb"
        },
      },
      bwa_refine => {
        class       => "GATK::Refine",
        perform     => 1,
        target_dir  => "${target_dir}/" . $dataset->{task_name} . "/bwa_refine",
        option      => "-Xmx40g",
        gatk_option => "--fix_misencoded_quality_scores",
        fasta_file  => $bwa_fasta,
        source_ref  => "bwa",
        vcf_files   => [$dbsnp],
        gatk_jar    => $gatk_jar,
        picard_jar  => $picard_jar,
        sh_direct   => 0,
        sorted      => 1,
        pbs         => {
          "email"    => $email,
          "nodes"    => "1:ppn=8",
          "walltime" => "240",
          "mem"      => "40gb"
        },
      },
      bwa_dexseqcount => {
        class        => "Count::DexseqCount",
        perform      => 1,
        target_dir   => "${target_dir}/" . $dataset->{task_name} . "/bwa_dexseqcount",
        option       => "",
        source_ref   => ["bwa"],
        gff_file     => $dexseq_gff,
        dexseq_count => "/home/shengq1/pylibs/bin/dexseq_count.py",
        sh_direct    => 0,
        pbs          => {
          "email"    => $email,
          "nodes"    => "1:ppn=1",
          "walltime" => "72",
          "mem"      => "40gb"
        },
      },
      bwa_exontable => {
        class         => "CQS::CQSDatatable",
        perform       => 1,
        target_dir    => "${target_dir}/" . $dataset->{task_name} . "/bwa_exontable",
        option        => "-p ENS --noheader -o " . $dataset->{task_name} . "_exon.count",
        source_ref    => "bwa_dexseqcount",
        name_map_file => $name_map_file,
        cqs_tools     => $cqstools,
        sh_direct     => 1,
        pbs           => {
          "email"    => $email,
          "nodes"    => "1:ppn=1",
          "walltime" => "10",
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
      glmvc => {
        class             => "Variants::GlmvcCall",
        perform           => 1,
        target_dir        => "${target_dir}/" . $dataset->{task_name} . "/glmvc",
        option            => "--glm_pvalue 0.1",
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
          "nodes"    => "1:ppn=8",
          "walltime" => "72",
          "mem"      => "40gb"
        },
      },
      glmvc_noMYC => {
        class             => "Variants::GlmvcCall",
        perform           => 1,
        target_dir        => "${target_dir}/" . $dataset->{task_name} . "/glmvc_noMYC",
        option            => "--glm_pvalue 0.1 --exclude_bed $exclude_bed -r 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,X,M --glm_use_raw_pvalue",
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
          "nodes"    => "1:ppn=8",
          "walltime" => "72",
          "mem"      => "40gb"
        },
      },
      glmvc_noMYC_table => {
        class        => "Variants::GlmvcTable",
        perform      => 1,
        target_dir   => "${target_dir}/" . $dataset->{task_name} . "/glmvc_noMYC_table",
        option       => "",
        source_ref   => [ "glmvc_noMYC", "annotation.tsv" ],
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

  push @individuals, ( "bwa", "bwa_refine" );

  if ( defined $dataset->{capture_bed} ) {
    $config = merge(
      $config,
      {
        conifer => {
          class       => "CNV::Conifer",
          perform     => 1,
          target_dir  => "${target_dir}/" . $dataset->{task_name} . "/conifer",
          option      => "",
          source_ref  => "bwa",
          conifer     => $conifer,
          bedfile     => $dataset->{capture_slim_bed},
          isbamsorted => 1,
          sh_direct   => 1,
          pbs         => {
            "email"    => $email,
            "nodes"    => "1:ppn=1",
            "walltime" => "72",
            "mem"      => "10gb"
          },
        },
        cnmops => {
          class       => "CNV::cnMops",
          perform     => 1,
          target_dir  => "${target_dir}/" . $dataset->{task_name} . "/cnmops",
          option      => "",
          source_ref  => "bwa",
          bedfile     => $dataset->{capture_bed},
          refnames    => ["N04_DUSP4flox_LACZ"],
          pairmode    => "paired",
          isbamsorted => 1,
          sh_direct   => 1,
          pbs         => {
            "email"    => $email,
            "nodes"    => "1:ppn=1",
            "walltime" => "72",
            "mem"      => "40gb"
          },
        },
        cnmops_bam => {
          class      => "Samtools::View",
          perform    => 1,
          target_dir => "${target_dir}/" . $dataset->{task_name} . "/cnmops_bam",
          option     => "-h -b -L " . "${target_dir}/" . $dataset->{task_name} . "/cnmops/result/" . $dataset->{task_name} . ".call.bed",
          extension  => ".bam",
          source_ref => "bwa",
          sh_direct  => 1,
          sorted     => 1,
          pbs        => {
            "email"    => $email,
            "nodes"    => "1:ppn=1",
            "walltime" => "24",
            "mem"      => "10gb"
          },
        },
      }
    );
    push @all, ( "conifer", "bwa_dexseqcount" );

    #performTask( $config, "conifer" );
    #performTask( $config, "cnmops_bam" );
  }
  else {
    $config->{glmvc_WES_validation} = {
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
    };

    #performTask( $config, "glmvc_WES_validation" );
  }

  $config->{sequencetask} = {
    class      => "CQS::SequenceTask",
    perform    => 1,
    target_dir => "${target_dir}/" . $dataset->{task_name} . "/sequencetask",
    option     => "",
    source     => {
      prepare => \@individuals,
      sm      => [ "muTect", "muTect_annovar", "glmvc" ],
      all     => \@all
    },
    sh_direct => 0,
    pbs       => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  };

  #performConfig($config);

  if ( $dataset == $wes ) {
    performTask( $config, "bwa_dexseqcount" );
    performTask( $config, "bwa_exontable" );
  }
}

1;

