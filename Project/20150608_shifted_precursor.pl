#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;

my $target_dir     = "/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor";
my $msgf_jar       = "/scratch/cqs/shengq1/local/bin/MSGFPlus/MSGFPlus.jar";
my $mod_file       = "/scratch/cqs/shengq1/local/bin/MSGFPlus/Mods.txt";
my $database_human = "/gpfs21/scratch/cqs/shengq1/proteomics/shifted/rev_Human_uniprot_sprot_v20120613.fasta";
my $database_yeast = "/gpfs21/scratch/cqs/shengq1/proteomics/shifted/rev_Yeast_uniprot_v20120613.fasta";
my $database_ecoli = "/gpfs21/scratch/cqs/shengq1/proteomics/shifted/rev_Ecoli_uniprot_v20120613_P4431.fasta";
my $email          = "quanhu.sheng\@vanderbilt.edu";

my $config = {
  general           => { task_name => "ShiftedTargetDecoy" },
  Elite_CIDIT_Human_MSGFPlus => {
    class      => "Proteomics::Engine::MSGFPlus",
    perform    => 1,
    target_dir => "${target_dir}/Elite_CIDIT_Human",
    option     => "-t 20ppm -ti \"0,1\" -tda 0 -m 1 -inst 0 -e 1 -protocol 5 -ntt 2 -n 1 -addFeatures 1",
    source     => {
      "Elite_CIDIT_Human.plus0.1dalton" => ["/gpfs21/scratch/cqs/shengq1/proteomics/shifted/Elite_CIDIT_Human.plus0.1dalton.mgf"],
      "Elite_CIDIT_Human.plus10dalton"  => ["/gpfs21/scratch/cqs/shengq1/proteomics/shifted/Elite_CIDIT_Human.plus10dalton.mgf"],
    },
    msgf_jar  => $msgf_jar,
    mod_file  => $mod_file,
    database  => $database_human,
    sh_direct => 0,
    pbs       => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  Fusion_CIDIT_Human_MSGFPlus => {
    class      => "Proteomics::Engine::MSGFPlus",
    perform    => 1,
    target_dir => "${target_dir}/Fusion_CIDIT_Human",
    option     => "-t 20ppm -ti \"0,1\" -tda 0 -m 1 -inst 0 -e 1 -protocol 5 -ntt 2 -n 1 -addFeatures 1",
    source     => {
      "Fusion_CIDIT_Human.plus0.1dalton" => ["/gpfs21/scratch/cqs/shengq1/proteomics/shifted/Fusion_CIDIT_Human.plus0.1dalton.mgf"],
      "Fusion_CIDIT_Human.plus10dalton"  => ["/gpfs21/scratch/cqs/shengq1/proteomics/shifted/Fusion_CIDIT_Human.plus10dalton.mgf"],
    },
    msgf_jar  => $msgf_jar,
    mod_file  => $mod_file,
    database  => $database_human,
    sh_direct => 0,
    pbs       => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  Fusion_HCDIT_Yeast_MSGFPlus => {
    class      => "Proteomics::Engine::MSGFPlus",
    perform    => 1,
    target_dir => "${target_dir}/Fusion_HCDIT_Yeast",
    option     => "-t 20ppm -ti \"0,1\" -tda 0 -m 3 -inst 0 -e 1 -protocol 5 -ntt 2 -n 1 -addFeatures 1",
    source     => {
      "Fusion_HCDIT_Yeast.plus0.1dalton" => ["/gpfs21/scratch/cqs/shengq1/proteomics/shifted/Fusion_HCDIT_Yeast.plus0.1dalton.mgf"],
      "Fusion_HCDIT_Yeast.plus10dalton"  => ["/gpfs21/scratch/cqs/shengq1/proteomics/shifted/Fusion_HCDIT_Yeast.plus10dalton.mgf"],
    },
    msgf_jar  => $msgf_jar,
    mod_file  => $mod_file,
    database  => $database_yeast,
    sh_direct => 0,
    pbs       => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  Fusion_HCDOT_Human_MSGFPlus => {
    class      => "Proteomics::Engine::MSGFPlus",
    perform    => 1,
    target_dir => "${target_dir}/Fusion_HCDOT_Human",
    option     => "-t 20ppm -ti \"0,1\" -tda 0 -m 3 -inst 3 -e 1 -protocol 5 -ntt 2 -n 1 -addFeatures 1",
    source     => {
      "Fusion_HCDOT_Human.plus0.1dalton" => ["/gpfs21/scratch/cqs/shengq1/proteomics/shifted/Fusion_HCDOT_Human.plus0.1dalton.mgf"],
      "Fusion_HCDOT_Human.plus10dalton"  => ["/gpfs21/scratch/cqs/shengq1/proteomics/shifted/Fusion_HCDOT_Human.plus10dalton.mgf"],
    },
    msgf_jar  => $msgf_jar,
    mod_file  => $mod_file,
    database  => $database_human,
    sh_direct => 0,
    pbs       => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  QExactive_HCDOT_Human_MSGFPlus => {
    class      => "Proteomics::Engine::MSGFPlus",
    perform    => 1,
    target_dir => "${target_dir}/QExactive_HCDOT_Human",
    option     => "-t 20ppm -ti \"0,1\" -tda 0 -m 3 -inst 3 -e 1 -protocol 5 -ntt 2 -n 1 -addFeatures 1",
    source     => {
      "QExactive_HCDOT_Human.plus0.1dalton" => ["/gpfs21/scratch/cqs/shengq1/proteomics/shifted/QExactive_HCDOT_Human.plus0.1dalton.mgf"],
      "QExactive_HCDOT_Human.plus10dalton"  => ["/gpfs21/scratch/cqs/shengq1/proteomics/shifted/QExactive_HCDOT_Human.plus10dalton.mgf"],
    },
    msgf_jar  => $msgf_jar,
    mod_file  => $mod_file,
    database  => $database_human,
    sh_direct => 0,
    pbs       => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  QTOF_Ecoli_MSGFPlus => {
    class      => "Proteomics::Engine::MSGFPlus",
    perform    => 1,
    target_dir => "${target_dir}/QTOF_Ecoli",
    option     => "-t 20ppm -ti \"0,1\" -tda 0 -m 3 -inst 2 -e 1 -protocol 5 -ntt 2 -n 1 -addFeatures 1",
    source     => {
      "QTOF_Ecoli.plus0.1dalton" => ["/gpfs21/scratch/cqs/shengq1/proteomics/shifted/QTOF_Ecoli.plus0.1dalton.mgf"],
      "QTOF_Ecoli.plus10dalton"  => ["/gpfs21/scratch/cqs/shengq1/proteomics/shifted/QTOF_Ecoli.plus10dalton.mgf"],
    },
    msgf_jar  => $msgf_jar,
    mod_file  => $mod_file,
    database  => $database_ecoli,
    sh_direct => 0,
    pbs       => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  MSGFPlus_Distiller => {
    class      => "Proteomics::Distiller::PSMDistiller",
    perform    => 1,
    target_dir => "${target_dir}/PSMDistiller",
    option     => "-e MSGF -t DTA",
    source_ref     => [ "Elite_CIDIT_Human_MSGFPlus", "Fusion_CIDIT_Human_MSGFPlus", "Fusion_HCDIT_Yeast_MSGFPlus", "Fusion_HCDIT_Yeast_MSGFPlus", "Fusion_HCDOT_Human_MSGFPlus", "QExactive_HCDOT_Human_MSGFPlus", "QTOF_Ecoli_MSGFPlus" ],
    proteomicstools  => $msgf_jar,
    sh_direct => 1,
    pbs       => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  sequencetask => {
    class      => "CQS::SequenceTask",
    perform    => 1,
    target_dir => "${target_dir}/sequencetask",
    option     => "",
    source     => { step_1 => [ "Elite_CIDIT_Human_MSGFPlus", "Fusion_CIDIT_Human_MSGFPlus", "Fusion_HCDIT_Yeast_MSGFPlus", "Fusion_HCDIT_Yeast_MSGFPlus", "Fusion_HCDOT_Human_MSGFPlus", "QExactive_HCDOT_Human_MSGFPlus", "QTOF_Ecoli_MSGFPlus" ], },
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

