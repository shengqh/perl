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
  Elite_CIDIT_Human => {
    class      => "Proteomics::Engine::MSGFPlus",
    perform    => 1,
    target_dir => "${target_dir}/Elite_CIDIT_Human",
    option     => "-t 20ppm -ti \"0,1\" -tda 0 -m 1 -inst 0 -e 1 -protocol 5 -ntt 2 -n 1 -addFeatures 1",
    source     => {
      "Elite_CIDIT_Human.0.1dalton" => ["/gpfs21/scratch/cqs/shengq1/proteomics/shifted/Elite_CIDIT_Human.plus0.1dalton.mgf"],
      "Elite_CIDIT_Human.10dalton"  => ["/gpfs21/scratch/cqs/shengq1/proteomics/shifted/Elite_CIDIT_Human.plus10dalton.mgf"],
    },
    msgf_jar  => $msgf_jar,
    mod_file  => $mod_file,
    database  => $database_human,
    sh_direct => 1,
    pbs       => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  Fusion_CIDIT_Human => {
    class      => "Proteomics::Engine::MSGFPlus",
    perform    => 1,
    target_dir => "${target_dir}/Fusion_CIDIT_Human",
    option     => "-t 20ppm -ti \"0,1\" -tda 0 -m 1 -inst 0 -e 1 -protocol 5 -ntt 2 -n 1 -addFeatures 1",
    source     => {
      "Fusion_CIDIT_Human.0.1dalton" => ["/gpfs21/scratch/cqs/shengq1/proteomics/shifted/Fusion_CIDIT_Human.plus0.1dalton.mgf"],
      "Fusion_CIDIT_Human.10dalton"  => ["/gpfs21/scratch/cqs/shengq1/proteomics/shifted/Fusion_CIDIT_Human.plus10dalton.mgf"],
    },
    msgf_jar  => $msgf_jar,
    mod_file  => $mod_file,
    database  => $database_human,
    sh_direct => 1,
    pbs       => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  Fusion_HCDIT_Yeast => {
    class      => "Proteomics::Engine::MSGFPlus",
    perform    => 1,
    target_dir => "${target_dir}/Fusion_HCDIT_Yeast",
    option     => "-t 20ppm -ti \"0,1\" -tda 0 -m 3 -inst 0 -e 1 -protocol 5 -ntt 2 -n 1 -addFeatures 1",
    source     => {
      "Fusion_HCDIT_Yeast.0.1dalton" => ["/gpfs21/scratch/cqs/shengq1/proteomics/shifted/Fusion_HCDIT_Yeast.plus0.1dalton.mgf"],
      "Fusion_HCDIT_Yeast.10dalton"  => ["/gpfs21/scratch/cqs/shengq1/proteomics/shifted/Fusion_HCDIT_Yeast.plus10dalton.mgf"],
    },
    msgf_jar  => $msgf_jar,
    mod_file  => $mod_file,
    database  => $database_yeast,
    sh_direct => 1,
    pbs       => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  Fusion_HCDOT_Human => {
    class      => "Proteomics::Engine::MSGFPlus",
    perform    => 1,
    target_dir => "${target_dir}/Fusion_HCDOT_Human",
    option     => "-t 20ppm -ti \"0,1\" -tda 0 -m 3 -inst 3 -e 1 -protocol 5 -ntt 2 -n 1 -addFeatures 1",
    source     => {
      "Fusion_HCDOT_Human.0.1dalton" => ["/gpfs21/scratch/cqs/shengq1/proteomics/shifted/Fusion_HCDOT_Human.plus0.1dalton.mgf"],
      "Fusion_HCDOT_Human.10dalton"  => ["/gpfs21/scratch/cqs/shengq1/proteomics/shifted/Fusion_HCDOT_Human.plus10dalton.mgf"],
    },
    msgf_jar  => $msgf_jar,
    mod_file  => $mod_file,
    database  => $database_human,
    sh_direct => 1,
    pbs       => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  QExactive_HCDOT_Human => {
    class      => "Proteomics::Engine::MSGFPlus",
    perform    => 1,
    target_dir => "${target_dir}/QExactive_HCDOT_Human",
    option     => "-t 20ppm -ti \"0,1\" -tda 0 -m 3 -inst 3 -e 1 -protocol 5 -ntt 2 -n 1 -addFeatures 1",
    source     => {
      "QExactive_HCDOT_Human.0.1dalton" => ["/gpfs21/scratch/cqs/shengq1/proteomics/shifted/QExactive_HCDOT_Human.plus0.1dalton.mgf"],
      "QExactive_HCDOT_Human.10dalton"  => ["/gpfs21/scratch/cqs/shengq1/proteomics/shifted/QExactive_HCDOT_Human.plus10dalton.mgf"],
    },
    msgf_jar  => $msgf_jar,
    mod_file  => $mod_file,
    database  => $database_human,
    sh_direct => 1,
    pbs       => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  QTOF_Ecoli => {
    class      => "Proteomics::Engine::MSGFPlus",
    perform    => 1,
    target_dir => "${target_dir}/QTOF_Ecoli",
    option     => "-t 20ppm -ti \"0,1\" -tda 0 -m 3 -inst 2 -e 1 -protocol 5 -ntt 2 -n 1 -addFeatures 1",
    source     => {
      "QTOF_Ecoli.0.1dalton" => ["/gpfs21/scratch/cqs/shengq1/proteomics/shifted/QTOF_Ecoli.plus0.1dalton.mgf"],
      "QTOF_Ecoli.10dalton"  => ["/gpfs21/scratch/cqs/shengq1/proteomics/shifted/QTOF_Ecoli.plus10dalton.mgf"],
    },
    msgf_jar  => $msgf_jar,
    mod_file  => $mod_file,
    database  => $database_ecoli,
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
    source     => { step_1 => [ "Elite_CIDIT_Human", "Fusion_CIDIT_Human", "Fusion_HCDIT_Yeast", "Fusion_HCDIT_Yeast", "Fusion_HCDOT_Human", "QExactive_HCDOT_Human", "QTOF_Ecoli" ], },
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

