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

my $bwa_fasta      = "/scratch/cqs/shengq1/references/hg19_16569_MT/bwa_index_0.7.12/hg19_16569_MT.fa";
my $transcript_gtf = "/scratch/cqs/shengq1/references/ensembl_gtf/v75/Homo_sapiens.GRCh37.75.MT.gtf";
my $name_map_file  = "/scratch/cqs/shengq1/references/ensembl_gtf/v75/Homo_sapiens.GRCh37.75.MT.map";

my $dbsnp  = "/scratch/cqs/shengq1/references/dbsnp/human_GRCh37_v142_16569_MT.vcf";
my $hapmap = "/scratch/cqs/shengq1/references/gatk/b37/hapmap_3.3.b37.vcf";
my $omni   = "/scratch/cqs/shengq1/references/gatk/b37/1000G_omni2.5.b37.vcf";
my $g1000  = "/scratch/cqs/shengq1/references/gatk/b37/1000G_phase1.snps.high_confidence.b37.vcf";
my $mills  = "/scratch/cqs/shengq1/references/gatk/b37/Mills_and_1000G_gold_standard.indels.b37.vcf";

my $cosmic    = "/scratch/cqs/shengq1/references/cosmic/cosmic_v71_hg19_16569_MT.vcf";
my $affy_file = "/data/cqs/shengq1/reference/affy/HG-U133_Plus_2.na33.annot.csv";

my $annovar_protocol  = "refGene,snp138,cosmic70";
my $annovar_operation = "g,f,f";
my $annovar_param     = "-protocol ${annovar_protocol} -operation ${annovar_operation} --remove";
my $annovar_db        = "/scratch/cqs/shengq1/references/annovar/humandb/";
my $rnaediting_db     = "/data/cqs/shengq1/reference/rnaediting/hg19.txt";

#minimum quality score 10, minimum overlap 4 bases, remove reads with length less than 30
my $cutadapt_option = "-q 10 -O 4 -m 30";

my $cluster = "slurm";

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

    #    "N16_DUSP4null_MYC" => [
    #      "/gpfs21/scratch/cqs/shengq1/dnaseq/20151029_balko_mouse_celllines/WES/data/3162-JMB-7_1_sequence.txt.gz",
    #      "/gpfs21/scratch/cqs/shengq1/dnaseq/20151029_balko_mouse_celllines/WES/data/3162-JMB-7_2_sequence.txt.gz"
    #    ],
    "N17_DUSP4null_Trp53null1_MYC" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20151029_balko_mouse_celllines/WES/data/3162-JMB-8_1_sequence.txt.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20151029_balko_mouse_celllines/WES/data/3162-JMB-8_2_sequence.txt.gz"
    ],
  },
  groups => {
    "N07_DUSP4flox_MYC"             => [ "N04_DUSP4flox_LACZ", "N07_DUSP4flox_MYC" ],
    "N05_DUSP4flox_Trp53null1_LACZ" => [ "N04_DUSP4flox_LACZ", "N05_DUSP4flox_Trp53null1_LACZ" ],
    "N09_DUSP4flox_Trp53null3_MYC"  => [ "N04_DUSP4flox_LACZ", "N09_DUSP4flox_Trp53null3_MYC" ],
    "N13_DUSP4null_LACZ"            => [ "N04_DUSP4flox_LACZ", "N13_DUSP4null_LACZ" ],

    #"N16_DUSP4null_MYC"             => [ "N04_DUSP4flox_LACZ", "N16_DUSP4null_MYC" ],
    "N15_DUSP4null_Trp53null3_LACZ" => [ "N04_DUSP4flox_LACZ", "N15_DUSP4null_Trp53null3_LACZ" ],
    "N17_DUSP4null_Trp53null1_MYC"  => [ "N04_DUSP4flox_LACZ", "N17_DUSP4null_Trp53null1_MYC" ],
  },

  covered_bed => "/scratch/cqs/shengq1/references/sureselect/S0276129_Mouse_All_Exon_V1/S0276129_Regions.bed"
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
  groups => {
    "N07_DUSP4flox_MYC"             => [ "N04_DUSP4flox_LACZ", "N07_DUSP4flox_MYC" ],
    "N05_DUSP4flox_Trp53null1_LACZ" => [ "N04_DUSP4flox_LACZ", "N05_DUSP4flox_Trp53null1_LACZ" ],
    "N09_DUSP4flox_Trp53null3_MYC"  => [ "N04_DUSP4flox_LACZ", "N09_DUSP4flox_Trp53null3_MYC" ],
    "N13_DUSP4null_LACZ"            => [ "N04_DUSP4flox_LACZ", "N13_DUSP4null_LACZ" ],
    "N16_DUSP4null_MYC"             => [ "N04_DUSP4flox_LACZ", "N16_DUSP4null_MYC" ],
    "N15_DUSP4null_Trp53null3_LACZ" => [ "N04_DUSP4flox_LACZ", "N15_DUSP4null_Trp53null3_LACZ" ],
    "N17_DUSP4null_Trp53null1_MYC"  => [ "N04_DUSP4flox_LACZ", "N17_DUSP4null_Trp53null1_MYC" ],
  },
  cutadapt => 1
};

my @datasets = ( $wes, $wgs );

for my $dataset (@datasets) {
  my $qc = {
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
  my $source      = "files";
  if ( defined $dataset->{cutadapt} && $dataset->{cutadapt} ) {
    $qc->{cutadapt} = {
      class      => "Trimmer::Cutadapt",
      perform    => 1,
      target_dir => "${target_dir}/" . $dataset->{task_name} . "/cutadapt",
      option     => $cutadapt_option,
      source_ref => "files",
      adapter    => "GATCGGAAGAGC",
      extension  => "_clipped.fastq.gz",
      sh_direct  => 1,
      cluster    => $cluster,
      pbs        => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "24",
        "mem"      => "20gb"
      },
    };
    $qc->{fastqlen} = {
      class      => "CQS::FastqLen",
      perform    => 1,
      target_dir => "${target_dir}/fastqlen",
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
    };
    $source = "cutadapt";
    push @individuals, ( "cutadapt", "fastqlen" );
  }

  my $config = merge(
    $qc,
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
      muTect => {
        class        => "GATK::MuTect",
        perform      => 1,
        target_dir   => "${target_dir}/" . $dataset->{task_name} . "/muTect",
        option       => "--min_qscore 20 --filter_reads_with_N_cigar",
        java_option  => "-Xmx40g",
        source_ref   => "bwa_refine",
        groups_ref   => "groups",
        fasta_file   => $bwa_fasta,
        cosmic_file  => $cosmic,
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
        target_dir => "${target_dir}/" . $dataset->{task_name} . "/muTect",
        option     => $annovar_param,
        source_ref => [ "muTect", ".pass.vcf\$" ],
        annovar_db => $annovar_db,
        buildver   => "hg19",
        cqstools   => $cqstools,
        affy_file  => $affy_file,
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
        annovar_buildver  => "hg19",
        annovar_protocol  => $annovar_protocol,
        annovar_operation => $annovar_operation,
        rnaediting_db     => $rnaediting_db,
        distance_exon_gtf => $transcript_gtf,
        sh_direct         => 0,
        execute_file      => $glmvc,
        pbs               => {
          "email"    => $email,
          "nodes"    => "1:ppn=8",
          "walltime" => "72",
          "mem"      => "40gb"
        },
      }
    }
  );

  push @individuals, ( "bwa", "bwa_refine" );
  my @all = ("fastqc_summary");

  if ( defined $dataset->{covered_bed} ) {
    $config->{conifer} = {
      class       => "CNV::Conifer",
      perform     => 1,
      target_dir  => "${target_dir}/" . $dataset->{task_name} . "/conifer",
      option      => "",
      source_ref  => "bwa",
      conifer     => $conifer,
      bedfile     => $dataset->{covered_bed},
      isbamsorted => 1,
      sh_direct   => 1,
      pbs         => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "720",
        "mem"      => "10gb"
      },
    };
    $config->{cnmops} = {
      class       => "CNV::cnMops",
      perform     => 1,
      target_dir  => "${target_dir}/" . $dataset->{task_name} . "/cnmops",
      option      => "",
      source_ref  => "bwa",
      bedfile     => $dataset->{convered_bed},
      pairmode    => "paired",
      isbamsorted => 1,
      sh_direct   => 1,
      pbs         => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "720",
        "mem"      => "40gb"
      },
    };
    push @all, ( "conifer", "cnmops" );
  }

  $config->{sequencetask} = {
    class      => "CQS::SequenceTask",
    perform    => 1,
    target_dir => "${target_dir}/" . $dataset->{task_name} . "/sequencetask",
    option     => "",
    source     => {
      prepare => [ "fastqc", "bwa",            "bwa_refine" ],
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

  performConfig($config);
}

1;

