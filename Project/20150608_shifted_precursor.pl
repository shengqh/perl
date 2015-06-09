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
  general     => { task_name => "ShiftedTargetDecoy" },
  CIDIT_Human => {
    class      => "Proteomics::Engine::MSGFPlus",
    perform    => 1,
    target_dir => "${target_dir}/CIDIT_Human",
    option     => "-t 20ppm -ti \"0,1\" -tda 0 -m 1 -inst 0 -e 1 -protocol 5 -ntt 2 -n 1 -addFeatures 1",
    source     => {
      "Elite_CIDIT_Human" =>
        [ "/gpfs21/scratch/cqs/shengq1/proteomics/shifted/Elite_CIDIT_Human_plus10dalton.mgf", "/gpfs21/scratch/cqs/shengq1/proteomics/shifted/Elite_CIDIT_Human_plus0.1dalton.mgf" ],
      "Fusion_CIDIT_Human" =>
        [ "/gpfs21/scratch/cqs/shengq1/proteomics/shifted/Fusion_CIDIT_Human.plus10dalton.mgf", "/gpfs21/scratch/cqs/shengq1/proteomics/shifted/Fusion_CIDIT_Human.plus0.1dalton.mgf" ],
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
      "Fusion_HCDIT_Yeast" =>
        [ "/gpfs21/scratch/cqs/shengq1/proteomics/shifted/Fusion_HCDIT_Yeast.plus10dalton.mgf", "/gpfs21/scratch/cqs/shengq1/proteomics/shifted/Fusion_HCDIT_Yeast.plus0.1dalton.mgf" ],
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
  sequencetask => {
    class      => "CQS::SequenceTask",
    perform    => 1,
    target_dir => "${target_dir}/sequencetask",
    option     => "",
    source     => { step_1 => [ "CIDIT_Human", "Fusion_HCDIT_Yeast" ], },
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

