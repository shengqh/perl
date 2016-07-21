#!/usr/bin/perl

######################
#Hg20 not finished yet
######################

use strict;
use warnings;

use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;
use Hash::Merge qw( merge );
my $target_dir = "/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38";
my $cqstools   = "/home/shengq1/cqstools/CQS.Tools.exe";
my $email      = "shilin.zhao\@vanderbilt.edu";

my $adapter         = "AGATCGGAAGAG";
my $cutadapt_option = "-m 50 ";

my $bwa_fasta_hg38 = "/scratch/cqs/zhaos/reference/hg38/hg38bundle/bwa_index_0.7.12/Homo_sapiens_assembly38.fasta";

#GATK realign
my $dbsnp_hg38              = "/scratch/cqs/zhaos/reference/hg38/hg38bundle/Homo_sapiens_assembly38.dbsnp138.vcf";
my $mills_g1000_indel_hg38  = "/scratch/cqs/zhaos/reference/hg38/hg38bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf";
my $g1000_phase1_indel_hg38 = "/scratch/cqs/zhaos/reference/hg38/hg38bundle/Homo_sapiens_assembly38.known_indels.vcf";

#GATK call SNV
my $hapmap_hg38    = "/scratch/cqs/zhaos/reference/hg38/hg38bundle/hapmap_3.3.hg38.vcf";
my $omni_hg38      = "/scratch/cqs/zhaos/reference/hg38/hg38bundle/1000G_omni2.5.hg38.vcf";
my $g1000_snp_hg38 = "/scratch/cqs/zhaos/reference/hg38/hg38bundle/1000G_phase1.snps.high_confidence.hg38.vcf";

#Other programe
my $annovar_db_hg38 = "/scratch/cqs/zhaos/reference/annovardb/humandb";
my $annovar_param   = "-protocol refGene,avsnp144,cosmic70,1000g2015aug_all,clinvar_20160302,exac03nontcga,esp6500siv2_all, -operation g,f,f,f,f,f,f --remove";
my $gatk_jar        = "/home/zhaos/bin/GenomeAnalysisTK36.jar";
my $picard_jar      = "/scratch/cqs/shengq1/local/bin/picard/picard.jar";
my $cosmic_hg38     = "/scratch/cqs/zhaos/reference/COSMIC/CosmicCodingMutsChrSortGatk.vcf";

my $covered_bed_hg20 = "/scratch/cqs/zhaos/Pierre/3199_bedfile/131212_HG19_cancer_TS_EZ_TiledOnly_95pToHg38.bed";

my $qc3_perl = "/home/zhaos/source/perl_cqs/software_QC/RC2/qc3.pl";

my $cluster = "slurm";

my $config = {
  general => { task_name => "3199" },

  files => {

    "T106"     => [ '/gpfs21/scratch/cqs/zhaos/Pierre/3199/run_1/3199-PM-100_1_sequence.txt.gz', '/gpfs21/scratch/cqs/zhaos/Pierre/3199/run_1/3199-PM-100_2_sequence.txt.gz' ],
    "T107"     => [ '/gpfs21/scratch/cqs/zhaos/Pierre/3199/run_1/3199-PM-101_1_sequence.txt.gz', '/gpfs21/scratch/cqs/zhaos/Pierre/3199/run_1/3199-PM-101_2_sequence.txt.gz' ],
    "T108"     => [ '/gpfs21/scratch/cqs/zhaos/Pierre/3199/run_1/3199-PM-102_1_sequence.txt.gz', '/gpfs21/scratch/cqs/zhaos/Pierre/3199/run_1/3199-PM-102_2_sequence.txt.gz' ],
    "T109"     => [ '/gpfs21/scratch/cqs/zhaos/Pierre/3199/run_1/3199-PM-103_1_sequence.txt.gz', '/gpfs21/scratch/cqs/zhaos/Pierre/3199/run_1/3199-PM-103_2_sequence.txt.gz' ],
  },

  cutadapt => {
    class      => "Trimmer::Cutadapt",
    perform    => 1,
    target_dir => "${target_dir}/cutadapt",
    option     => $cutadapt_option,
    source_ref => "files",
    adapter    => $adapter,
    extension  => "_clipped.fastq",
    sh_direct  => 0,
    cluster    => $cluster,
    pairend    => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "2",
      "mem"      => "20gb"
    },
  },

  fastqcPostTrim => {
    class      => "QC::FastQC",
    perform    => 1,
    target_dir => "${target_dir}/fastqcPostTrim",
    option     => "",
    source_ref => "cutadapt",
    sh_direct  => 1,
    cluster    => $cluster,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=2",
      "walltime" => "1",
      "mem"      => "40gb"
    },
  },
  fastqc_summary_PostTrim => {
    class      => "QC::FastQCSummary",
    perform    => 1,
    target_dir => "${target_dir}/fastqcPostTrim",
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

  bwa_hg38 => {
    class      => "Alignment::BWA",
    perform    => 1,
    target_dir => "${target_dir}/bwa",
    option     => "",
    bwa_index  => $bwa_fasta_hg38,
    source_ref => "cutadapt",
    picard_jar => $picard_jar,
    sh_direct  => 0,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "8",
      "mem"      => "40gb"
    },
  },

  bwa_refine_hg38 => {
    class            => "GATK::Refine",
    perform          => 1,
    target_dir       => "${target_dir}/bwa_refine",
    option           => "-Xmx40g",
    gatk_option      => "--fix_misencoded_quality_scores",
    fasta_file       => $bwa_fasta_hg38,
    source_ref       => "bwa_hg38",
    indel_vcf_files  => [ $mills_g1000_indel_hg38, $g1000_phase1_indel_hg38 ],
    known_vcf_files  => [$dbsnp_hg38],
    bed_file         => $covered_bed_hg20,
    interval_padding => 100,
    gatk_jar         => $gatk_jar,
    picard_jar       => $picard_jar,
    sh_direct        => 0,
    sorted           => 1,
    pbs              => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },

  bwa_refine_KeepDupReads_hg38 => {
    class            => "GATK::Refine",
    perform          => 1,
    target_dir       => "${target_dir}/bwa_refine_KeepDupReads",
    option           => "-Xmx40g",
    gatk_option      => "--fix_misencoded_quality_scores",
    fasta_file       => $bwa_fasta_hg38,
    source_ref       => "bwa_hg38",
    indel_vcf_files  => [ $mills_g1000_indel_hg38, $g1000_phase1_indel_hg38 ],
    known_vcf_files  => [$dbsnp_hg38],
    bed_file         => $covered_bed_hg20,
    interval_padding => 100,
    remove_duplicate => 0,
    gatk_jar         => $gatk_jar,
    picard_jar       => $picard_jar,
    sh_direct        => 0,
    sorted           => 1,
    pbs              => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },

  bwa_refine_hc_gvcf_hg38 => {
    class            => "GATK::HaplotypeCallerGVCF",
    perform          => 1,
    target_dir       => "${target_dir}/bwa_refine_hc_gvcf",
    option           => "",
    source_ref       => "bwa_refine_hg38",
    java_option      => "",
    fasta_file       => $bwa_fasta_hg38,
    dbsnp_vcf        => $dbsnp_hg38,
    bed_file         => $covered_bed_hg20,
    interval_padding => 100,
    gatk_jar         => $gatk_jar,
    extension        => ".gvcf",
    sh_direct        => 0,
    pbs              => {
      "email"    => $email,
      "nodes"    => "1:ppn=4",
      "walltime" => "60",
      "mem"      => "70gb"
    },
  },

  bwa_refine_hc_gvcf_KeepDupReads_hg38 => {
    class            => "GATK::HaplotypeCallerGVCF",
    perform          => 1,
    target_dir       => "${target_dir}/bwa_refine_hc_gvcf_KeepDupReads",
    option           => "",
    source_ref       => "bwa_refine_KeepDupReads_hg38",
    java_option      => "",
    fasta_file       => $bwa_fasta_hg38,
    dbsnp_vcf        => $dbsnp_hg38,
    bed_file         => $covered_bed_hg20,
    interval_padding => 100,
    gatk_jar         => $gatk_jar,
    extension        => ".gvcf",
    sh_direct        => 0,
    pbs              => {
      "email"    => $email,
      "nodes"    => "1:ppn=4",
      "walltime" => "60",
      "mem"      => "70gb"
    },
  },

  bwa_refine_hc_gvcf_vqsr_hg38 => {
    class       => "GATK::VariantFilter",
    perform     => 1,
    target_dir  => "${target_dir}/bwa_refine_hc_gvcf_vqsr",
    option      => "",
    source_ref  => "bwa_refine_hc_gvcf_hg38",
    java_option => "",
    fasta_file  => $bwa_fasta_hg38,
    dbsnp_vcf   => $dbsnp_hg38,
    hapmap_vcf  => $hapmap_hg38,
    omni_vcf    => $omni_hg38,
    g1000_vcf   => $g1000_snp_hg38,
    mills_vcf   => $mills_g1000_indel_hg38,
    bed_file    => $covered_bed_hg20,
    gatk_jar    => $gatk_jar,
    cqstools    => $cqstools,
    sh_direct   => 1,
    pbs         => {
      "email"    => $email,
      "nodes"    => "1:ppn=12",
      "walltime" => "96",
      "mem"      => "70gb"
    },
  },

  bwa_refine_hc_gvcf_hardFilter_hg38 => {
    class       => "GATK::VariantFilter",
    perform     => 1,
    target_dir  => "${target_dir}/bwa_refine_hc_gvcf_hardFilter",
    option      => "",
    source_ref  => "bwa_refine_hc_gvcf_hg38",
    java_option => "",
    fasta_file  => $bwa_fasta_hg38,
    dbsnp_vcf   => $dbsnp_hg38,
    is_rna      => 0,

    #		hapmap_vcf  => $hapmap_hg38,
    #		omni_vcf    => $omni_hg38,
    #		g1000_vcf   => $g1000_snp_hg38,
    #		mills_vcf   => $mills_g1000_indel_hg38,
    bed_file  => $covered_bed_hg20,
    gatk_jar  => $gatk_jar,
    cqstools  => $cqstools,
    sh_direct => 1,
    pbs       => {
      "email"    => $email,
      "nodes"    => "1:ppn=12",
      "walltime" => "4",
      "mem"      => "70gb"
    },
  },

  qc3vcf_hg38_snp => {
    class      => "QC::QC3vcf",
    perform    => 1,
    target_dir => "${target_dir}/bwa_refine_hc_gvcf_vqsr_qc3",
    option     => "",
    qc3_perl   => $qc3_perl,
    source_ref => [ "bwa_refine_hc_gvcf_vqsr_hg38", "snp" ],
    annovar_db => $annovar_db_hg38,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },

  bwa_refine_hc_gvcf_vqsr_annovar_hg38_snp => {
    class      => "Annotation::Annovar",
    perform    => 1,
    target_dir => "${target_dir}/bwa_refine_hc_gvcf_vqsr_annovar",
    source_ref => [ "bwa_refine_hc_gvcf_vqsr_hg38", "snp" ],
    option     => $annovar_param,
    annovar_db => $annovar_db_hg38,
    buildver   => "hg38",
    sh_direct  => 1,
    isvcf      => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "10gb"
    },
  },

  bwa_refine_hc_gvcf_hardfilter_annovar_hg38_snp => {
    class      => "Annotation::Annovar",
    perform    => 1,
    target_dir => "${target_dir}/bwa_refine_hc_gvcf_hardfilter_annovar",
    source_ref => [ "bwa_refine_hc_gvcf_hardFilter_hg38", "snp" ],
    option     => $annovar_param,
    annovar_db => $annovar_db_hg38,
    buildver   => "hg38",
    sh_direct  => 1,
    isvcf      => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "10gb"
    },
  },

  bwa_refine_hc_gvcf_vqsr_annovar_hg38_indel => {
    class      => "Annotation::Annovar",
    perform    => 1,
    target_dir => "${target_dir}/bwa_refine_hc_gvcf_vqsr_annovar_indel",
    source_ref => [ "bwa_refine_hc_gvcf_vqsr_hg38", "indel" ],
    option     => $annovar_param,
    annovar_db => $annovar_db_hg38,
    buildver   => "hg38",
    sh_direct  => 1,
    isvcf      => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "10gb"
    },
  },

  muTect2 => {
    class       => "GATK::MuTect2",
    perform     => 1,
    target_dir  => "${target_dir}/muTect2",
    option      => "",
    java_option => "-Xmx60g",
    source_ref  => "bwa_refine_hg38",

    #		groups_ref   => "groups",
    bed_file    => $covered_bed_hg20,
    fasta_file  => $bwa_fasta_hg38,
    dbsnp_file  => $dbsnp_hg38,
    cosmic_file => $cosmic_hg38,
    sh_direct   => 0,
    gatk_jar    => $gatk_jar,
    pbs         => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "240",
      "mem"      => "60gb"
    },
  },

  sequencetask => {
    class      => "CQS::SequenceTask",
    perform    => 1,
    target_dir => "${target_dir}/sequencetask",
    option     => "",
    source     => {
      step1 => [ "cutadapt", "fastqcPostTrim", "bwa_hg38", "bwa_refine_hg38", "bwa_refine_hc_gvcf_hg38" ],
      step2 => ["fastqc_summary_PostTrim"],
      step4 => ["bwa_refine_hc_gvcf_vqsr_hg38"],
      step5 => [ "qc3vcf_hg38_snp", "bwa_refine_hc_gvcf_vqsr_annovar_hg38_snp", "bwa_refine_hc_gvcf_vqsr_annovar_hg38_indel" ],
    },
    sh_direct => 0,
    pbs       => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "36",
      "mem"      => "70gb"
    },
  },

  #	bam_KRAS_mpileup_uniqueRun => {
  #		class                       => "CQS::UniqueRunProgram",
  #		perform                     => 1,
  #		target_dir                  => "${target_dir}/bam_mpileup_uniqueRun",
  #		parameterSampleFile1Label   => "--bam_files ",
  #		parameterSampleFile1SepEach => ",",
  #		parameterSampleFile1_ref    => "bwa_refine_hg38",
  #		parameterFile2Label   => "--bed_file ",
  #		parameterFile2 =>
  #		  "/scratch/cqs/zhaos/Pierre/3199_bedfile/TargetKrasHg38EachPos.bed",
  #		parameterSampleFile3Label => "--bam_names ",
  #		parameterSampleFile3_ref  => "bwa_hg38",
  #		parameterSampleFile3SepEach => ",",
  #		parameterSampleFile3Regex => "([a-zA-Z]+\\d+).bam\$",
  #		runProgram                => "/scratch/cqs/shengq1/mono/bin/mono",
  #		option                    => "/home/shengq1/glmvc/glmvc.exe extract --fasta /scratch/cqs/zhaos/reference/hg38/hg38bundle/bwa_index_0.7.12/Homo_sapiens_assembly38.fasta",
  #		output_file_label         => "-o ",
  #		output_file               => ".bam.KRAS",
  #		output_file_ext           => ".mpileup",
  #		sh_direct                 => 1,
  #		pbs                       => {
  #			"email"    => $email,
  #			"nodes"    => "1:ppn=1",
  #			"walltime" => "40",
  #			"mem"      => "60gb"
  #		},
  #	},
  bam_TP53_mpileup_uniqueRun => {
    class                       => "CQS::UniqueRunProgram",
    perform                     => 1,
    target_dir                  => "${target_dir}/bam_mpileup_uniqueRun",
    parameterSampleFile1Label   => "--bam_files ",
    parameterSampleFile1SepEach => ",",
    parameterSampleFile1_ref    => "bwa_refine_hg38",
    parameterFile2Label         => "--bed_file ",
    parameterFile2              => "/scratch/cqs/zhaos/Pierre/3199_bedfile/TargetTP53Hg38EachPos.bed",
    parameterSampleFile3Label   => "--bam_names ",
    parameterSampleFile3_ref    => "bwa_hg38",
    parameterSampleFile3SepEach => ",",
    parameterSampleFile3Regex   => "([a-zA-Z]+\\d+).bam\$",
    runProgram                  => "/scratch/cqs/shengq1/mono/bin/mono",
    option                      => "/home/shengq1/glmvc/glmvc.exe extract --fasta /scratch/cqs/zhaos/reference/hg38/hg38bundle/bwa_index_0.7.12/Homo_sapiens_assembly38.fasta",
    output_file_label           => "-o ",
    output_file                 => ".bam.TP53",
    output_file_ext             => ".mpileup",
    sh_direct                   => 1,
    pbs                         => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "40",
      "mem"      => "60gb"
    },
  },

  bam_EGFR_mpileup_uniqueRun => {
    class                       => "CQS::UniqueRunProgram",
    perform                     => 1,
    target_dir                  => "${target_dir}/bam_mpileup_uniqueRun",
    parameterSampleFile1Label   => "--bam_files ",
    parameterSampleFile1SepEach => ",",
    parameterSampleFile1_ref    => "bwa_refine_hg38",
    parameterFile2Label         => "--bed_file ",
    parameterFile2              => "/scratch/cqs/zhaos/Pierre/3199_bedfile/TargetEGFRHg38EachPos.bed",
    parameterSampleFile3Label   => "--bam_names ",
    parameterSampleFile3_ref    => "bwa_hg38",
    parameterSampleFile3SepEach => ",",
    parameterSampleFile3Regex   => "([a-zA-Z]+\\d+).bam\$",
    runProgram                  => "/scratch/cqs/shengq1/mono/bin/mono",
    option                      => "/home/shengq1/glmvc/glmvc.exe extract --fasta /scratch/cqs/zhaos/reference/hg38/hg38bundle/bwa_index_0.7.12/Homo_sapiens_assembly38.fasta",
    output_file_label           => "-o ",
    output_file                 => ".bam.EGFR",
    output_file_ext             => ".mpileup",
    sh_direct                   => 1,
    pbs                         => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "40",
      "mem"      => "60gb"
    },
  },

  bam_KRAS_mpileup_Run => {
    class        => "CQS::RunProgram",
    perform      => 1,
    target_dir   => "${target_dir}/bam_mpileup_Run",
    source1Label => "--bam_files ",
    source1_ref  => "bwa_hg38",

    #		source2Label => "--bed_files ",
    #		source2 => "/scratch/cqs/zhaos/Pierre/3199_bedfile/TargetKrasHg38.bed",
    source3Label => "--bam_names ",
    source3_ref  => "bwa_hg38",
    source3Regex => "([a-zA-Z]+\\d+).bam\$",
    runProgram   => "/scratch/cqs/shengq1/mono/bin/mono",
    option =>
"/home/shengq1/glmvc/glmvc.exe extract --bed_file /scratch/cqs/zhaos/Pierre/3199_bedfile/TargetKrasHg38EachPos.bed --fasta /scratch/cqs/zhaos/reference/hg38/hg38bundle/bwa_index_0.7.12/Homo_sapiens_assembly38.fasta",
    output_file_label => "-o ",
    output_file       => ".bam.KRAS",
    output_file_ext   => ".mpileup",
    sh_direct         => 1,
    pbs               => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "1",
      "mem"      => "10gb"
    },
  },

  #	 bam_KRAS_mpileup_Refined_Run => {
  #        class        => "CQS::RunProgram",
  #        perform      => 1,
  #        target_dir   => "${target_dir}/bam_mpileup_Refined_Run",
  #        source1Label => "--bam_files ",
  #        source1_ref  => "bwa_refine_hg38",
##       source2Label => "--bed_files ",
##       source2 => "/scratch/cqs/zhaos/Pierre/3199_bedfile/TargetKrasHg38.bed",
#        source3Label      => "--bam_names ",
#        source3_ref       => "bwa_hg38",
#        source3Regex      => "([a-zA-Z]+\\d+).bam\$",
#        runProgram        => "/scratch/cqs/shengq1/mono/bin/mono",
#        option            => "/home/shengq1/glmvc/glmvc.exe extract --bed_file /scratch/cqs/zhaos/Pierre/3199_bedfile/TargetKrasHg38EachPos.bed --fasta /scratch/cqs/zhaos/reference/hg38/hg38bundle/bwa_index_0.7.12/Homo_sapiens_assembly38.fasta",
#        output_file_label => "-o ",
#        output_file       => ".refined.bam.KRAS",
#        output_file_ext   => ".mpileup",
#        sh_direct         => 1,
#        pbs               => {
#            "email"    => $email,
#            "nodes"    => "1:ppn=1",
#            "walltime" => "1",
#            "mem"      => "10gb"
#        },
#    },

#	bam_TP53_mpileup_Refined_Run => {
#		class        => "CQS::RunProgram",
#		perform      => 1,
#		target_dir   => "${target_dir}/bam_mpileup_Refined_Run",
#		source1Label => "--bam_files ",
#		source1_ref  => "bwa_refine_hg38",
#		source3Label => "--bam_names ",
#		source3_ref  => "bwa_hg38",
#		source3Regex => "([a-zA-Z]+\\d+).bam\$",
#		runProgram   => "/scratch/cqs/shengq1/mono/bin/mono",
#		option =>
#"/home/shengq1/glmvc/glmvc.exe extract --bed_file /scratch/cqs/zhaos/Pierre/3199_bedfile/TargetTP53Hg38EachPos.bed --fasta /scratch/cqs/zhaos/reference/hg38/hg38bundle/bwa_index_0.7.12/Homo_sapiens_assembly38.fasta",
#		output_file_label => "-o ",
#		output_file       => ".refined.bam.TP53",
#		output_file_ext   => ".mpileup",
#		sh_direct         => 1,
#		pbs               => {
#			"email"    => $email,
#			"nodes"    => "1:ppn=1",
#			"walltime" => "1",
#			"mem"      => "10gb"
#		},
#	},

  bam_EGFR_mpileup_Refined_Run => {
    class        => "CQS::RunProgram",
    perform      => 1,
    target_dir   => "${target_dir}/bam_mpileup_Refined_Run",
    source1Label => "--bam_files ",
    source1_ref  => "bwa_refine_hg38",
    source3Label => "--bam_names ",
    source3_ref  => "bwa_hg38",
    source3Regex => "([a-zA-Z]+\\d+).bam\$",
    runProgram   => "/scratch/cqs/shengq1/mono/bin/mono",
    option =>
"/home/shengq1/glmvc/glmvc.exe extract --bed_file /scratch/cqs/zhaos/Pierre/3199_bedfile/TargetEGFRHg38EachPos.bed --fasta /scratch/cqs/zhaos/reference/hg38/hg38bundle/bwa_index_0.7.12/Homo_sapiens_assembly38.fasta",
    output_file_label => "-o ",
    output_file       => ".refined.bam.EGFR",
    output_file_ext   => ".mpileup",
    sh_direct         => 1,
    pbs               => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "1",
      "mem"      => "10gb"
    },
  },
  bam_KRAS_subBam_Run => {
    class             => "CQS::RunProgram",
    perform           => 1,
    target_dir        => "${target_dir}/bam_subBam_Refined_Run",
    source1Label      => "",
    source1_ref       => "bwa_refine_hg38",
    runProgram        => "/scratch/cqs/shengq1/local/bin/samtools",
    option            => " view -L /scratch/cqs/zhaos/Pierre/3199_bedfile/TargetKrasHg38.bed",
    output_file_label => "> ",
    output_file       => ".KRAS",
    output_file_ext   => ".bam",
    sh_direct         => 1,
    pbs               => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "1",
      "mem"      => "10gb"
    },
  },
};

performConfig($config);

1;
