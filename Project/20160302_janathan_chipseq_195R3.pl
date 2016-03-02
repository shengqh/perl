#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;
use CQS::FileUtils;
use Data::Dumper;

my $target_dir = create_directory_or_die("/scratch/cqs/shengq1/chipseq/20160302_janathan_chipseq_195R3");

my $fasta_file   = "/scratch/cqs/shengq1/references/gencode/hg19/bowtie_index_1.1.2/GRCh37.p13.genome.fa";
my $bowtie_index = "/scratch/cqs/shengq1/references/gencode/hg19/bowtie_index_1.1.2/GRCh37.p13.genome";
my $cqstools     = "/home/shengq1/cqstools/cqstools.exe";

my $email = "quanhu.sheng\@vanderbilt.edu";
my $task  = "195R3";

my $config = {
  general => { task_name => $task },
  files   => {
    "CT480-1-16_bc03_GGACCC_L005_R1_001.fastq"      => ["/gpfs21/scratch/cqs/shengq1/chipseq/195R3/195R/day_0/CT480-1-16_bc03_GGACCC_L005_R1_001.fastq.gz"],
    "CT480-1-16_bc03_GGACCC_L006_R1_001.fastq"      => ["/gpfs21/scratch/cqs/shengq1/chipseq/195R3/195R/day_0/CT480-1-16_bc03_GGACCC_L006_R1_001.fastq.gz"],
    "CT480-1-16_bc03_GGACCC_L007_R1_001.fastq"      => ["/gpfs21/scratch/cqs/shengq1/chipseq/195R3/195R/day_0/CT480-1-16_bc03_GGACCC_L007_R1_001.fastq.gz"],
    "CT480-1-16_bc03_GGACCC_L008_R1_001.fastq"      => ["/gpfs21/scratch/cqs/shengq1/chipseq/195R3/195R/day_0/CT480-1-16_bc03_GGACCC_L008_R1_001.fastq.gz"],
    "CT480-1-16_bc04_TTCAGC_L005_R1_001.fastq"      => ["/gpfs21/scratch/cqs/shengq1/chipseq/195R3/195R/day_0/CT480-1-16_bc04_TTCAGC_L005_R1_001.fastq.gz"],
    "CT480-1-16_bc04_TTCAGC_L006_R1_001.fastq"      => ["/gpfs21/scratch/cqs/shengq1/chipseq/195R3/195R/day_0/CT480-1-16_bc04_TTCAGC_L006_R1_001.fastq.gz"],
    "CT480-1-16_bc04_TTCAGC_L007_R1_001.fastq"      => ["/gpfs21/scratch/cqs/shengq1/chipseq/195R3/195R/day_0/CT480-1-16_bc04_TTCAGC_L007_R1_001.fastq.gz"],
    "CT480-1-16_bc04_TTCAGC_L008_R1_001.fastq"      => ["/gpfs21/scratch/cqs/shengq1/chipseq/195R3/195R/day_0/CT480-1-16_bc04_TTCAGC_L008_R1_001.fastq.gz"],
    "CT480-1-16_bc11_GTGTTA_L005_R1_001.fastq"      => ["/gpfs21/scratch/cqs/shengq1/chipseq/195R3/195R/day_0/CT480-1-16_bc11_GTGTTA_L005_R1_001.fastq.gz"],
    "CT480-1-16_bc11_GTGTTA_L006_R1_001.fastq"      => ["/gpfs21/scratch/cqs/shengq1/chipseq/195R3/195R/day_0/CT480-1-16_bc11_GTGTTA_L006_R1_001.fastq.gz"],
    "CT480-1-16_bc11_GTGTTA_L007_R1_001.fastq"      => ["/gpfs21/scratch/cqs/shengq1/chipseq/195R3/195R/day_0/CT480-1-16_bc11_GTGTTA_L007_R1_001.fastq.gz"],
    "CT480-1-16_bc11_GTGTTA_L008_R1_001.fastq"      => ["/gpfs21/scratch/cqs/shengq1/chipseq/195R3/195R/day_0/CT480-1-16_bc11_GTGTTA_L008_R1_001.fastq.gz"],
    "CT480-1-16_bc12_TGTGAA_L005_R1_001.fastq"      => ["/gpfs21/scratch/cqs/shengq1/chipseq/195R3/195R/day_0/CT480-1-16_bc12_TGTGAA_L005_R1_001.fastq.gz"],
    "CT480-1-16_bc12_TGTGAA_L006_R1_001.fastq"      => ["/gpfs21/scratch/cqs/shengq1/chipseq/195R3/195R/day_0/CT480-1-16_bc12_TGTGAA_L006_R1_001.fastq.gz"],
    "CT480-1-16_bc12_TGTGAA_L007_R1_001.fastq"      => ["/gpfs21/scratch/cqs/shengq1/chipseq/195R3/195R/day_0/CT480-1-16_bc12_TGTGAA_L007_R1_001.fastq.gz"],
    "CT480-1-16_bc12_TGTGAA_L008_R1_001.fastq"      => ["/gpfs21/scratch/cqs/shengq1/chipseq/195R3/195R/day_0/CT480-1-16_bc12_TGTGAA_L008_R1_001.fastq.gz"],
    "CT644-20_ACTTGA_L005_R1_001.fastq"             => ["/gpfs21/scratch/cqs/shengq1/chipseq/195R3/195R/day_17_shear/CT644-20_ACTTGA_L005_R1_001.fastq.gz"],
    "CT644-20_ACTTGA_L005_R1_002.fastq"             => ["/gpfs21/scratch/cqs/shengq1/chipseq/195R3/195R/day_17_shear/CT644-20_ACTTGA_L005_R1_002.fastq.gz"],
    "CT644-20_ACTTGA_L005_R1_003.fastq"             => ["/gpfs21/scratch/cqs/shengq1/chipseq/195R3/195R/day_17_shear/CT644-20_ACTTGA_L005_R1_003.fastq.gz"],
    "CT644-20_ACTTGA_L005_R1_004.fastq"             => ["/gpfs21/scratch/cqs/shengq1/chipseq/195R3/195R/day_17_shear/CT644-20_ACTTGA_L005_R1_004.fastq.gz"],
    "CT644-20_ACTTGA_L005_R1_005.fastq"             => ["/gpfs21/scratch/cqs/shengq1/chipseq/195R3/195R/day_17_shear/CT644-20_ACTTGA_L005_R1_005.fastq.gz"],
    "CT644-20_ACTTGA_L005_R1_006.fastq"             => ["/gpfs21/scratch/cqs/shengq1/chipseq/195R3/195R/day_17_shear/CT644-20_ACTTGA_L005_R1_006.fastq.gz"],
    "CT644-20_ACTTGA_L005_R1_007.fastq"             => ["/gpfs21/scratch/cqs/shengq1/chipseq/195R3/195R/day_17_shear/CT644-20_ACTTGA_L005_R1_007.fastq.gz"],
    "CT644-20_ACTTGA_L005_R1_008.fastq"             => ["/gpfs21/scratch/cqs/shengq1/chipseq/195R3/195R/day_17_shear/CT644-20_ACTTGA_L005_R1_008.fastq.gz"],
    "CT644-23_GGCTAC_L006_R1_001.fastq"             => ["/gpfs21/scratch/cqs/shengq1/chipseq/195R3/195R/day_17_static/CT644-23_GGCTAC_L006_R1_001.fastq.gz"],
    "CT644-23_GGCTAC_L006_R1_002.fastq"             => ["/gpfs21/scratch/cqs/shengq1/chipseq/195R3/195R/day_17_static/CT644-23_GGCTAC_L006_R1_002.fastq.gz"],
    "CT644-23_GGCTAC_L006_R1_003.fastq"             => ["/gpfs21/scratch/cqs/shengq1/chipseq/195R3/195R/day_17_static/CT644-23_GGCTAC_L006_R1_003.fastq.gz"],
    "CT644-23_GGCTAC_L006_R1_004.fastq"             => ["/gpfs21/scratch/cqs/shengq1/chipseq/195R3/195R/day_17_static/CT644-23_GGCTAC_L006_R1_004.fastq.gz"],
    "CT644-23_GGCTAC_L006_R1_005.fastq"             => ["/gpfs21/scratch/cqs/shengq1/chipseq/195R3/195R/day_17_static/CT644-23_GGCTAC_L006_R1_005.fastq.gz"],
    "CT644-23_GGCTAC_L006_R1_006.fastq"             => ["/gpfs21/scratch/cqs/shengq1/chipseq/195R3/195R/day_17_static/CT644-23_GGCTAC_L006_R1_006.fastq.gz"],
    "CT644-23_GGCTAC_L006_R1_007.fastq"             => ["/gpfs21/scratch/cqs/shengq1/chipseq/195R3/195R/day_17_static/CT644-23_GGCTAC_L006_R1_007.fastq.gz"],
    "CT644-23_GGCTAC_L006_R1_008.fastq"             => ["/gpfs21/scratch/cqs/shengq1/chipseq/195R3/195R/day_17_static/CT644-23_GGCTAC_L006_R1_008.fastq.gz"],
    "CT644-23_GGCTAC_L006_R1_009.fastq"             => ["/gpfs21/scratch/cqs/shengq1/chipseq/195R3/195R/day_17_static/CT644-23_GGCTAC_L006_R1_009.fastq.gz"],
    "CT644-23_GGCTAC_L006_R1_010.fastq"             => ["/gpfs21/scratch/cqs/shengq1/chipseq/195R3/195R/day_17_static/CT644-23_GGCTAC_L006_R1_010.fastq.gz"],
    "CT644-24_CTTGTA_L006_R1_001.fastq"             => ["/gpfs21/scratch/cqs/shengq1/chipseq/195R3/195R/day_17_static/CT644-24_CTTGTA_L006_R1_001.fastq.gz"],
    "CT644-24_CTTGTA_L006_R1_002.fastq"             => ["/gpfs21/scratch/cqs/shengq1/chipseq/195R3/195R/day_17_static/CT644-24_CTTGTA_L006_R1_002.fastq.gz"],
    "CT644-24_CTTGTA_L006_R1_003.fastq"             => ["/gpfs21/scratch/cqs/shengq1/chipseq/195R3/195R/day_17_static/CT644-24_CTTGTA_L006_R1_003.fastq.gz"],
    "CT644-24_CTTGTA_L006_R1_004.fastq"             => ["/gpfs21/scratch/cqs/shengq1/chipseq/195R3/195R/day_17_static/CT644-24_CTTGTA_L006_R1_004.fastq.gz"],
    "CT644-24_CTTGTA_L006_R1_005.fastq"             => ["/gpfs21/scratch/cqs/shengq1/chipseq/195R3/195R/day_17_static/CT644-24_CTTGTA_L006_R1_005.fastq.gz"],
    "CT644-24_CTTGTA_L006_R1_006.fastq"             => ["/gpfs21/scratch/cqs/shengq1/chipseq/195R3/195R/day_17_static/CT644-24_CTTGTA_L006_R1_006.fastq.gz"],
    "CT644-24_CTTGTA_L006_R1_007.fastq"             => ["/gpfs21/scratch/cqs/shengq1/chipseq/195R3/195R/day_17_static/CT644-24_CTTGTA_L006_R1_007.fastq.gz"],
    "CT644-24_CTTGTA_L006_R1_008.fastq"             => ["/gpfs21/scratch/cqs/shengq1/chipseq/195R3/195R/day_17_static/CT644-24_CTTGTA_L006_R1_008.fastq.gz"],
    "CT644-24_CTTGTA_L006_R1_009.fastq"             => ["/gpfs21/scratch/cqs/shengq1/chipseq/195R3/195R/day_17_static/CT644-24_CTTGTA_L006_R1_009.fastq.gz"],
    "CT644-24_CTTGTA_L006_R1_010.fastq"             => ["/gpfs21/scratch/cqs/shengq1/chipseq/195R3/195R/day_17_static/CT644-24_CTTGTA_L006_R1_010.fastq.gz"],
    "CT644-24_CTTGTA_L006_R1_011.fastq"             => ["/gpfs21/scratch/cqs/shengq1/chipseq/195R3/195R/day_17_static/CT644-24_CTTGTA_L006_R1_011.fastq.gz"],
    "CT663p15-30_index4_TGACCA_L001_R1_001.fastq"   => ["/gpfs21/scratch/cqs/shengq1/chipseq/195R3/195R/day_17_shear/CT663p15-30_index4_TGACCA_L001_R1_001.fastq.gz"],
    "CT663p15-30_index4_TGACCA_L001b_R1_001.fastq"  => ["/gpfs21/scratch/cqs/shengq1/chipseq/195R3/195R/day_17_shear/CT663p15-30_index4_TGACCA_L001b_R1_001.fastq.gz"],
    "CT663p15-30_index7_CAGATC_L001_R1_001.fastq"   => ["/gpfs21/scratch/cqs/shengq1/chipseq/195R3/195R/day_17_static/CT663p15-30_index7_CAGATC_L001_R1_001.fastq.gz"],
    "CT663p15-30_index7_CAGATC_L001b_R1_001.fastq"  => ["/gpfs21/scratch/cqs/shengq1/chipseq/195R3/195R/day_17_static/CT663p15-30_index7_CAGATC_L001b_R1_001.fastq.gz"],
    "CT663p15-30_index8_ACTTGA_L001_R1_001.fastq"   => ["/gpfs21/scratch/cqs/shengq1/chipseq/195R3/195R/day_17_static/CT663p15-30_index8_ACTTGA_L001_R1_001.fastq.gz"],
    "CT663p15-30_index8_ACTTGA_L001b_R1_001.fastq"  => ["/gpfs21/scratch/cqs/shengq1/chipseq/195R3/195R/day_17_static/CT663p15-30_index8_ACTTGA_L001b_R1_001.fastq.gz"],
    "CT663p39-46_index7_CAGATC_L003_R1_001.fastq"   => ["/gpfs21/scratch/cqs/shengq1/chipseq/195R3/195R/day_17_static/CT663p39-46_index7_CAGATC_L003_R1_001.fastq.gz"],
    "CT663p39-46_index7_CAGATC_L003b_R1_001.fastq"  => ["/gpfs21/scratch/cqs/shengq1/chipseq/195R3/195R/day_17_static/CT663p39-46_index7_CAGATC_L003b_R1_001.fastq.gz"],
    "CT663p39-46_index8_ACTTGA_L003_R1_001.fastq"   => ["/gpfs21/scratch/cqs/shengq1/chipseq/195R3/195R/day_17_static/CT663p39-46_index8_ACTTGA_L003_R1_001.fastq.gz"],
    "CT663p39-46_index8_ACTTGA_L003b_R1_001.fastq"  => ["/gpfs21/scratch/cqs/shengq1/chipseq/195R3/195R/day_17_static/CT663p39-46_index8_ACTTGA_L003b_R1_001.fastq.gz"],
    "CT663p47-54_index3_TTAGGC_L004_R1_001.fastq"   => ["/gpfs21/scratch/cqs/shengq1/chipseq/195R3/195R/day_17_static/CT663p47-54_index3_TTAGGC_L004_R1_001.fastq.gz"],
    "CT663p47-54_index3_TTAGGC_L004b_R1_001.fastq"  => ["/gpfs21/scratch/cqs/shengq1/chipseq/195R3/195R/day_17_static/CT663p47-54_index3_TTAGGC_L004b_R1_001.fastq.gz"],
    "CT663p47-54_index4_TGACCA_L004_R1_001.fastq"   => ["/gpfs21/scratch/cqs/shengq1/chipseq/195R3/195R/day_17_static/CT663p47-54_index4_TGACCA_L004_R1_001.fastq.gz"],
    "CT663p47-54_index4_TGACCA_L004b_R1_001.fastq"  => ["/gpfs21/scratch/cqs/shengq1/chipseq/195R3/195R/day_17_static/CT663p47-54_index4_TGACCA_L004b_R1_001.fastq.gz"],
    "CT663p55-62_index11_GGCTAC_L005_R1_001.fastq"  => ["/gpfs21/scratch/cqs/shengq1/chipseq/195R3/195R/day_17_static/CT663p55-62_index11_GGCTAC_L005_R1_001.fastq.gz"],
    "CT663p55-62_index11_GGCTAC_L005b_R1_001.fastq" => ["/gpfs21/scratch/cqs/shengq1/chipseq/195R3/195R/day_17_static/CT663p55-62_index11_GGCTAC_L005b_R1_001.fastq.gz"],
    "CT663p55-62_index12_CTTGTA_L005_R1_001.fastq"  => ["/gpfs21/scratch/cqs/shengq1/chipseq/195R3/195R/day_17_static/CT663p55-62_index12_CTTGTA_L005_R1_001.fastq.gz"],
    "CT663p55-62_index12_CTTGTA_L005b_R1_001.fastq" => ["/gpfs21/scratch/cqs/shengq1/chipseq/195R3/195R/day_17_static/CT663p55-62_index12_CTTGTA_L005b_R1_001.fastq.gz"],
    "CT663p55-62_index7_CAGATC_L005_R1_001.fastq"   => ["/gpfs21/scratch/cqs/shengq1/chipseq/195R3/195R/day_17_shear/CT663p55-62_index7_CAGATC_L005_R1_001.fastq.gz"],
    "CT663p55-62_index7_CAGATC_L005b_R1_001.fastq"  => ["/gpfs21/scratch/cqs/shengq1/chipseq/195R3/195R/day_17_shear/CT663p55-62_index7_CAGATC_L005b_R1_001.fastq.gz"],
    "CT663p55-62_index8_ACTTGA_L005_R1_001.fastq"   => ["/gpfs21/scratch/cqs/shengq1/chipseq/195R3/195R/day_17_shear/CT663p55-62_index8_ACTTGA_L005_R1_001.fastq.gz"],
    "CT663p55-62_index8_ACTTGA_L005b_R1_001.fastq"  => ["/gpfs21/scratch/cqs/shengq1/chipseq/195R3/195R/day_17_shear/CT663p55-62_index8_ACTTGA_L005b_R1_001.fastq.gz"],
    "H3K27ac_ActMo_3s5F_D0_2_4_13CS.fastq"          => ["/gpfs21/scratch/cqs/shengq1/chipseq/195R3/195R/day_0/H3K27ac_ActMo_3s5F_D0_2_4_13CS.fastq.gz"],
    "H3K27ac_H9_Lot016_062513.fastq"                => ["/gpfs21/scratch/cqs/shengq1/chipseq/195R3/195R/day_0/H3K27ac_H9_Lot016_062513.fastq.gz"],
    "H3K27me3_H9_Lot219_062513.fastq"               => ["/gpfs21/scratch/cqs/shengq1/chipseq/195R3/195R/day_0/H3K27me3_H9_Lot219_062513.fastq.gz"],
    "H3K27me3_Mil_3s5F_D0_2_4_13CS.fastq"           => ["/gpfs21/scratch/cqs/shengq1/chipseq/195R3/195R/day_0/H3K27me3_Mil_3s5F_D0_2_4_13CS.fastq.gz"],
  },
  fastqc_pre_trim => {
    class      => "QC::FastQC",
    perform    => 1,
    target_dir => "${target_dir}/fastqc_pre_trim",
    option     => "",
    source_ref => "files",
    sh_direct  => 0,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=2",
      "walltime" => "2",
      "mem"      => "10gb"
    },
  },
  fastqc_pre_trim_summary => {
    class      => "QC::FastQCSummary",
    perform    => 1,
    sh_direct  => 1,
    target_dir => "${target_dir}/fastqc_pre_trim",
    cqstools   => $cqstools,
    option     => "",
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "2",
      "mem"      => "10gb"
    },
  },
  cutadapt => {
    class      => "Trimmer::Cutadapt",
    perform    => 1,
    target_dir => "${target_dir}/cutadapt",
    option     => "-O 10 -m 30",
    source_ref => "files",
    adapter    => "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",    #trueseq adapter
    extension  => "_clipped.fastq",
    sh_direct  => 0,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "24",
      "mem"      => "20gb"
    },
  },
  fastqc_post_trim => {
    class      => "QC::FastQC",
    perform    => 1,
    target_dir => "${target_dir}/fastqc_post_trim",
    option     => "",
    sh_direct  => 1,
    source_ref => [ "cutadapt", ".fastq.gz" ],
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "2",
      "mem"      => "10gb"
    },
  },
  fastqc_post_trim_summary => {
    class      => "QC::FastQCSummary",
    perform    => 1,
    sh_direct  => 1,
    target_dir => "${target_dir}/fastqc_post_trim",
    cqstools   => $cqstools,
    option     => "",
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "2",
      "mem"      => "10gb"
    },
  },
  fastq_len => {
    class      => "CQS::FastqLen",
    perform    => 1,
    target_dir => "$target_dir/fastq_len",
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
  bowtie1 => {
    class         => "Alignment::Bowtie1",
    perform       => 1,
    target_dir    => "${target_dir}/bowtie1",
    option        => "-v 1 -m 1 --best --strata",
    fasta_file    => $fasta_file,
    source_ref    => [ "cutadapt", ".fastq.gz\$" ],
    bowtie1_index => $bowtie_index,
    sh_direct     => 0,
    pbs           => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  MACS => {
    class      => "Chipseq::MACS",
    perform    => 0,
    target_dir => "${target_dir}/MACS",
    option     => "-w",
    source_ref => "bowtie1",
    #groups_ref => "groups",
    #pairs_ref  => "pairs",
    sh_direct  => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  bradner_rose2 => {
    class                 => "Chipseq::BradnerRose2",
    perform               => 0,
    target_dir            => "${target_dir}/BradnerRose2",
    option                => "",
    source_ref            => "bowtie1",
    #groups_ref            => "groups",
    pipeline_dir          => "/scratch/cqs/shengq1/local/bin/bradnerlab",
    binding_site_file_ref => [ "MACS", ".bed\$" ],
    genome                => "hg19",
    sh_direct             => 1,
    pbs                   => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  sequencetask => {
    class      => "CQS::SequenceTask",
    perform    => 0,
    target_dir => "${target_dir}/sequencetask",
    option     => "",
    source     => {
      step_1 => [  "fastqc_pre_trim",          "cutadapt", "fastqc_post_trim", "bowtie1" ],
      step_2 => [ "fastqc_pre_trim_summary", "fastqc_post_trim_summary", "MACS",     "bradner_rose2" ],
    },
    sh_direct => 0,
    pbs       => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
};

performConfig($config);
#performTask( $config, "MACS" );

#print Dumper($config);

#performTask( $config, "bradner_rose2" );

1;
