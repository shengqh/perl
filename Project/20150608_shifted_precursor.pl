#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use Data::Dumper;

my $target_dir      = "/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor";
my $msgf_jar        = "/scratch/cqs/shengq1/local/bin/MSGFPlus/MSGFPlus.jar";
my $msamanda        = "/scratch/cqs/shengq1/local/bin/MSAmanda_Standalone_LinuxMac_1.0.0.4485/MSAmanda.exe";
my $proteomicstools = "/home/shengq1/proteomicstools/ProteomicsTools.exe";
my $mod_file        = "/scratch/cqs/shengq1/local/bin/MSGFPlus/Mods.txt";
my $database_human  = "/scratch/cqs/shengq1/proteomics/shifted/rev_Human_uniprot_sprot_v20120613.fasta";
my $database_yeast  = "/scratch/cqs/shengq1/proteomics/shifted/rev_Yeast_uniprot_v20120613.fasta";
my $database_ecoli  = "/scratch/cqs/shengq1/proteomics/shifted/rev_Ecoli_uniprot_v20120613_P4431.fasta";
my $email           = "quanhu.sheng\@vanderbilt.edu";

my $datasets = {
  Elite_CIDIT_Human => {
    source => {
      "Elite_CIDIT_Human_1"               => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/Elite_CIDIT_Human/mgf/20141017Test_DJMa_Cell_04.mgf"],
      "Elite_CIDIT_Human_1.minus10dalton" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/Elite_CIDIT_Human/mgf/20141017Test_DJMa_Cell_04.minus10dalton.mgf"],
      "Elite_CIDIT_Human_1.plus0.1dalton" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/Elite_CIDIT_Human/mgf/20141017Test_DJMa_Cell_04.plus0.1dalton.mgf"],
      "Elite_CIDIT_Human_1.plus10dalton"  => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/Elite_CIDIT_Human/mgf/20141017Test_DJMa_Cell_04.plus10dalton.mgf"],
      "Elite_CIDIT_Human_2"               => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/Elite_CIDIT_Human/mgf/20141017Test_DJMa_Cell_05.mgf"],
      "Elite_CIDIT_Human_2.minus10dalton" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/Elite_CIDIT_Human/mgf/20141017Test_DJMa_Cell_05.minus10dalton.mgf"],
      "Elite_CIDIT_Human_2.plus0.1dalton" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/Elite_CIDIT_Human/mgf/20141017Test_DJMa_Cell_05.plus0.1dalton.mgf"],
      "Elite_CIDIT_Human_2.plus10dalton"  => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/Elite_CIDIT_Human/mgf/20141017Test_DJMa_Cell_05.plus10dalton.mgf"],
      "Elite_CIDIT_Human_3"               => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/Elite_CIDIT_Human/mgf/20141017Test_DJMa_Cell_06.mgf"],
      "Elite_CIDIT_Human_3.minus10dalton" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/Elite_CIDIT_Human/mgf/20141017Test_DJMa_Cell_06.minus10dalton.mgf"],
      "Elite_CIDIT_Human_3.plus0.1dalton" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/Elite_CIDIT_Human/mgf/20141017Test_DJMa_Cell_06.plus0.1dalton.mgf"],
      "Elite_CIDIT_Human_3.plus10dalton"  => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/Elite_CIDIT_Human/mgf/20141017Test_DJMa_Cell_06.plus10dalton.mgf"],
    },
    database              => $database_human,
    MSGF_option           => "-t 20ppm -ti \"0,1\" -tda 0 -m 1 -inst 0 -e 1 -protocol 5 -ntt 2 -n 2 -addFeatures 1",
    comet_config_file     => "/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/parameters/comet-it.params",
    myrimatch_config_file => "/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/parameters/myrimatch-it.config",
    msamanda_config_file  => "/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/parameters/msamanda-it.xml",
  },
  Fusion_CIDIT_Human => {
    source => {
      "Fusion_CIDIT_Human_1"               => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/Fusion_CIDIT_Human/mgf/B01_03_14050113_HCC_Hela_Qu_CID_IT_pepID.mgf"],
      "Fusion_CIDIT_Human_1.minus10dalton" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/Fusion_CIDIT_Human/mgf/B01_03_14050113_HCC_Hela_Qu_CID_IT_pepID.minus10dalton.mgf"],
      "Fusion_CIDIT_Human_1.plus0.1dalton" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/Fusion_CIDIT_Human/mgf/B01_03_14050113_HCC_Hela_Qu_CID_IT_pepID.plus0.1dalton.mgf"],
      "Fusion_CIDIT_Human_1.plus10dalton"  => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/Fusion_CIDIT_Human/mgf/B01_03_14050113_HCC_Hela_Qu_CID_IT_pepID.plus10dalton.mgf"],
      "Fusion_CIDIT_Human_2"               => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/Fusion_CIDIT_Human/mgf/B11_03_140524_HCC_Hela_Qu_CID_IT_pepID.mgf"],
      "Fusion_CIDIT_Human_2.minus10dalton" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/Fusion_CIDIT_Human/mgf/B11_03_140524_HCC_Hela_Qu_CID_IT_pepID.minus10dalton.mgf"],
      "Fusion_CIDIT_Human_2.plus0.1dalton" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/Fusion_CIDIT_Human/mgf/B11_03_140524_HCC_Hela_Qu_CID_IT_pepID.plus0.1dalton.mgf"],
      "Fusion_CIDIT_Human_2.plus10dalton"  => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/Fusion_CIDIT_Human/mgf/B11_03_140524_HCC_Hela_Qu_CID_IT_pepID.plus10dalton.mgf"],
      "Fusion_CIDIT_Human_3"               => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/Fusion_CIDIT_Human/mgf/B13_03_140612_HCC3_Hela_Qu_CID_IT_pepID.mgf"],
      "Fusion_CIDIT_Human_3.minus10dalton" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/Fusion_CIDIT_Human/mgf/B13_03_140612_HCC3_Hela_Qu_CID_IT_pepID.minus10dalton.mgf"],
      "Fusion_CIDIT_Human_3.plus0.1dalton" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/Fusion_CIDIT_Human/mgf/B13_03_140612_HCC3_Hela_Qu_CID_IT_pepID.plus0.1dalton.mgf"],
      "Fusion_CIDIT_Human_3.plus10dalton"  => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/Fusion_CIDIT_Human/mgf/B13_03_140612_HCC3_Hela_Qu_CID_IT_pepID.plus10dalton.mgf"],
    },
    database              => $database_human,
    MSGF_option           => "-t 20ppm -ti \"0,1\" -tda 0 -m 1 -inst 0 -e 1 -protocol 5 -ntt 2 -n 2 -addFeatures 1",
    comet_config_file     => "/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/parameters/comet-it.params",
    myrimatch_config_file => "/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/parameters/myrimatch-it.config",
    msamanda_config_file  => "/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/parameters/msamanda-it.xml",
  },
  Fusion_HCDIT_Yeast => {
    source => {
      "Fusion_HCDIT_Yeast_1"               => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/Fusion_HCDIT_Yeast/mgf/10sep2013_yeast_control_1.mgf"],
      "Fusion_HCDIT_Yeast_1.minus10dalton" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/Fusion_HCDIT_Yeast/mgf/10sep2013_yeast_control_1.minus10dalton.mgf"],
      "Fusion_HCDIT_Yeast_1.plus0.1dalton" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/Fusion_HCDIT_Yeast/mgf/10sep2013_yeast_control_1.plus0.1dalton.mgf"],
      "Fusion_HCDIT_Yeast_1.plus10dalton"  => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/Fusion_HCDIT_Yeast/mgf/10sep2013_yeast_control_1.plus10dalton.mgf"],
      "Fusion_HCDIT_Yeast_2"               => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/Fusion_HCDIT_Yeast/mgf/10sep2013_yeast_control_2.mgf"],
      "Fusion_HCDIT_Yeast_2.minus10dalton" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/Fusion_HCDIT_Yeast/mgf/10sep2013_yeast_control_2.minus10dalton.mgf"],
      "Fusion_HCDIT_Yeast_2.plus0.1dalton" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/Fusion_HCDIT_Yeast/mgf/10sep2013_yeast_control_2.plus0.1dalton.mgf"],
      "Fusion_HCDIT_Yeast_2.plus10dalton"  => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/Fusion_HCDIT_Yeast/mgf/10sep2013_yeast_control_2.plus10dalton.mgf"],
      "Fusion_HCDIT_Yeast_3"               => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/Fusion_HCDIT_Yeast/mgf/10sep2013_yeast_control_4.mgf"],
      "Fusion_HCDIT_Yeast_3.minus10dalton" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/Fusion_HCDIT_Yeast/mgf/10sep2013_yeast_control_4.minus10dalton.mgf"],
      "Fusion_HCDIT_Yeast_3.plus0.1dalton" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/Fusion_HCDIT_Yeast/mgf/10sep2013_yeast_control_4.plus0.1dalton.mgf"],
      "Fusion_HCDIT_Yeast_3.plus10dalton"  => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/Fusion_HCDIT_Yeast/mgf/10sep2013_yeast_control_4.plus10dalton.mgf"],
    },
    database              => $database_yeast,
    MSGF_option           => "-t 20ppm -ti \"0,1\" -tda 0 -m 3 -inst 0 -e 1 -protocol 5 -ntt 2 -n 2 -addFeatures 1",
    comet_config_file     => "/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/parameters/comet-it.params",
    myrimatch_config_file => "/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/parameters/myrimatch-it.config",
    msamanda_config_file  => "/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/parameters/msamanda-it.xml",
  },
  Fusion_HCDOT_Human => {
    source => {
      "Fusion_HCDOT_Human_1"               => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/Fusion_HCDOT_Human/mgf/B00_01_140613_zhuxu_Qu_HCD_OT_2hr.mgf"],
      "Fusion_HCDOT_Human_1.minus10dalton" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/Fusion_HCDOT_Human/mgf/B00_01_140613_zhuxu_Qu_HCD_OT_2hr.minus10dalton.mgf"],
      "Fusion_HCDOT_Human_1.plus0.1dalton" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/Fusion_HCDOT_Human/mgf/B00_01_140613_zhuxu_Qu_HCD_OT_2hr.plus0.1dalton.mgf"],
      "Fusion_HCDOT_Human_1.plus10dalton"  => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/Fusion_HCDOT_Human/mgf/B00_01_140613_zhuxu_Qu_HCD_OT_2hr.plus10dalton.mgf"],
      "Fusion_HCDOT_Human_2"               => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/Fusion_HCDOT_Human/mgf/B00_01_140621_zhuxu_2_Qu_HCD_OT_2hr.mgf"],
      "Fusion_HCDOT_Human_2.minus10dalton" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/Fusion_HCDOT_Human/mgf/B00_01_140621_zhuxu_2_Qu_HCD_OT_2hr.minus10dalton.mgf"],
      "Fusion_HCDOT_Human_2.plus0.1dalton" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/Fusion_HCDOT_Human/mgf/B00_01_140621_zhuxu_2_Qu_HCD_OT_2hr.plus0.1dalton.mgf"],
      "Fusion_HCDOT_Human_2.plus10dalton"  => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/Fusion_HCDOT_Human/mgf/B00_01_140621_zhuxu_2_Qu_HCD_OT_2hr.plus10dalton.mgf"],
      "Fusion_HCDOT_Human_3"               => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/Fusion_HCDOT_Human/mgf/B00_01_140622_zhuxu_2_Qu_HCD_OT_2hr.mgf"],
      "Fusion_HCDOT_Human_3.minus10dalton" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/Fusion_HCDOT_Human/mgf/B00_01_140622_zhuxu_2_Qu_HCD_OT_2hr.minus10dalton.mgf"],
      "Fusion_HCDOT_Human_3.plus0.1dalton" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/Fusion_HCDOT_Human/mgf/B00_01_140622_zhuxu_2_Qu_HCD_OT_2hr.plus0.1dalton.mgf"],
      "Fusion_HCDOT_Human_3.plus10dalton"  => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/Fusion_HCDOT_Human/mgf/B00_01_140622_zhuxu_2_Qu_HCD_OT_2hr.plus10dalton.mgf"],
    },
    database              => $database_human,
    MSGF_option           => "-t 20ppm -ti \"0,1\" -tda 0 -m 3 -inst 3 -e 1 -protocol 5 -ntt 2 -n 2 -addFeatures 1",
    comet_config_file     => "/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/parameters/comet-ot.params",
    myrimatch_config_file => "/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/parameters/myrimatch-ot.config",
    msamanda_config_file  => "/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/parameters/msamanda-ot.xml",
  },
  QExactive_HCDOT_Human => {
    source => {
      "QExactive_HCDOT_Human_1"               => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/QExactive_HCDOT_Human/mgf/QExactive1.mgf"],
      "QExactive_HCDOT_Human_1.minus10dalton" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/QExactive_HCDOT_Human/mgf/QExactive1.minus10dalton.mgf"],
      "QExactive_HCDOT_Human_1.plus0.1dalton" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/QExactive_HCDOT_Human/mgf/QExactive1.plus0.1dalton.mgf"],
      "QExactive_HCDOT_Human_1.plus10dalton"  => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/QExactive_HCDOT_Human/mgf/QExactive1.plus10dalton.mgf"],
      "QExactive_HCDOT_Human_2"               => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/QExactive_HCDOT_Human/mgf/QExactive2.mgf"],
      "QExactive_HCDOT_Human_2.minus10dalton" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/QExactive_HCDOT_Human/mgf/QExactive2.minus10dalton.mgf"],
      "QExactive_HCDOT_Human_2.plus0.1dalton" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/QExactive_HCDOT_Human/mgf/QExactive2.plus0.1dalton.mgf"],
      "QExactive_HCDOT_Human_2.plus10dalton"  => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/QExactive_HCDOT_Human/mgf/QExactive2.plus10dalton.mgf"],
      "QExactive_HCDOT_Human_3"               => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/QExactive_HCDOT_Human/mgf/QExactive3.mgf"],
      "QExactive_HCDOT_Human_3.minus10dalton" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/QExactive_HCDOT_Human/mgf/QExactive3.minus10dalton.mgf"],
      "QExactive_HCDOT_Human_3.plus0.1dalton" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/QExactive_HCDOT_Human/mgf/QExactive3.plus0.1dalton.mgf"],
      "QExactive_HCDOT_Human_3.plus10dalton"  => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/QExactive_HCDOT_Human/mgf/QExactive3.plus10dalton.mgf"],
    },
    database              => $database_human,
    MSGF_option           => "-t 20ppm -ti \"0,1\" -tda 0 -m 3 -inst 3 -e 1 -protocol 5 -ntt 2 -n 2 -addFeatures 1",
    comet_config_file     => "/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/parameters/comet-ot.params",
    myrimatch_config_file => "/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/parameters/myrimatch-ot.config",
    msamanda_config_file  => "/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/parameters/msamanda-ot.xml",
  },
  QTOF_Ecoli => {
    source => {
      "QTOF_Ecoli_1"               => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/QTOF_Ecoli/mgf/ecoli-0500-r001.mgf"],
      "QTOF_Ecoli_1.minus10dalton" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/QTOF_Ecoli/mgf/ecoli-0500-r001.minus10dalton.mgf"],
      "QTOF_Ecoli_1.plus0.1dalton" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/QTOF_Ecoli/mgf/ecoli-0500-r001.plus0.1dalton.mgf"],
      "QTOF_Ecoli_1.plus10dalton"  => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/QTOF_Ecoli/mgf/ecoli-0500-r001.plus10dalton.mgf"],
      "QTOF_Ecoli_2"               => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/QTOF_Ecoli/mgf/ecoli-0500-r002.mgf"],
      "QTOF_Ecoli_2.minus10dalton" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/QTOF_Ecoli/mgf/ecoli-0500-r002.minus10dalton.mgf"],
      "QTOF_Ecoli_2.plus0.1dalton" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/QTOF_Ecoli/mgf/ecoli-0500-r002.plus0.1dalton.mgf"],
      "QTOF_Ecoli_2.plus10dalton"  => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/QTOF_Ecoli/mgf/ecoli-0500-r002.plus10dalton.mgf"],
      "QTOF_Ecoli_3"               => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/QTOF_Ecoli/mgf/ecoli-0500-r003.mgf"],
      "QTOF_Ecoli_3.minus10dalton" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/QTOF_Ecoli/mgf/ecoli-0500-r003.minus10dalton.mgf"],
      "QTOF_Ecoli_3.plus0.1dalton" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/QTOF_Ecoli/mgf/ecoli-0500-r003.plus0.1dalton.mgf"],
      "QTOF_Ecoli_3.plus10dalton"  => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/QTOF_Ecoli/mgf/ecoli-0500-r003.plus10dalton.mgf"],
    },
    database              => $database_ecoli,
    MSGF_option           => "-t 20ppm -ti \"0,1\" -tda 0 -m 3 -inst 2 -e 1 -protocol 5 -ntt 2 -n 2 -addFeatures 1",
    comet_config_file     => "/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/parameters/comet-qtof.params",
    myrimatch_config_file => "/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/parameters/myrimatch-qtof.config",
    msamanda_config_file  => "/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/parameters/msamanda-qtof.xml",
  },
};

my $configall = { general => { task_name => "ShiftedTargetDecoy" }, };

my @sections     = ();
my @psmdistiller = ();
my @msamanda     = ();

for my $dataset ( sort keys %{$datasets} ) {
  $configall->{ $dataset . "_MSGF" } = {
    class      => "Proteomics::Engine::MSGFPlus",
    task_name  => $dataset,
    perform    => 1,
    target_dir => "${target_dir}/$dataset/MSGF",
    option     => $datasets->{$dataset}->{MSGF_option},
    source     => $datasets->{$dataset}->{source},
    msgf_jar   => $msgf_jar,
    mod_file   => $mod_file,
    database   => $datasets->{$dataset}->{database},
    sh_direct  => 0,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  };

  $configall->{ $dataset . "_MSGF_PSM" } = {
    class           => "Proteomics::Distiller::PSMDistiller",
    task_name       => $dataset,
    perform         => 1,
    target_dir      => "${target_dir}/$dataset/MSGF",
    option          => "-e MSGF -t DTA",
    source_ref      => $dataset . "_MSGF",
    proteomicstools => $proteomicstools,
    sh_direct       => 1,
    pbs             => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  };

  $configall->{ $dataset . "_MSGF_PSM_rank2" } = {
    class           => "Proteomics::Distiller::PSMDistiller",
    task_name       => $dataset,
    perform         => 1,
    target_dir      => "${target_dir}/$dataset/MSGF_PSM_rank2",
    option          => "-e MSGF -t DTA --rank2",
    source_ref      => $dataset . "_MSGF",
    proteomicstools => $proteomicstools,
    sh_direct       => 1,
    pbs             => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  };

  push( @sections,     $dataset . "_MSGF" );
  push( @psmdistiller, $dataset . "_MSGF_PSM" );
  push( @psmdistiller, $dataset . "_MSGF_PSM_rank2" );

  $configall->{ $dataset . "_comet" } = {
    class           => "Proteomics::Engine::Comet",
    task_name       => $dataset,
    perform         => 1,
    target_dir      => "${target_dir}/$dataset/comet",
    option          => "",
    source          => $datasets->{$dataset}->{source},
    database        => $datasets->{$dataset}->{database},
    param_file      => $datasets->{$dataset}->{comet_config_file},
    proteomicstools => $proteomicstools,
    titleformat     => "DTA",
    sh_direct       => 0,
    pbs             => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  };

  $configall->{ $dataset . "_comet_PSM" } = {
    class           => "Proteomics::Distiller::PSMDistiller",
    task_name       => $dataset,
    perform         => 1,
    target_dir      => "${target_dir}/$dataset/comet",
    option          => "-e Comet -t DTA",
    source_ref      => $dataset . "_comet",
    proteomicstools => $proteomicstools,
    sh_direct       => 1,
    pbs             => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  };

  $configall->{ $dataset . "_comet_PSM_rank2" } = {
    class           => "Proteomics::Distiller::PSMDistiller",
    task_name       => $dataset . "_rank2",
    perform         => 1,
    target_dir      => "${target_dir}/$dataset/comet_PSM_rank2",
    option          => "-e Comet -t DTA --rank2",
    source_ref      => [ $dataset . "_comet", "" ],
    proteomicstools => $proteomicstools,
    sh_direct       => 1,
    pbs             => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  };

  push( @sections,     $dataset . "_comet" );
  push( @psmdistiller, $dataset . "_comet_PSM" );
  push( @psmdistiller, $dataset . "_comet_PSM_rank2" );

  $configall->{ $dataset . "_myrimatch" } = {
    class      => "Proteomics::Engine::Myrimatch",
    task_name  => $dataset,
    perform    => 1,
    target_dir => "${target_dir}/$dataset/myrimatch",
    option     => "",
    source     => $datasets->{$dataset}->{source},
    cfgfile    => $datasets->{$dataset}->{myrimatch_config_file},
    database   => $datasets->{$dataset}->{database},
    sh_direct  => 0,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  };

  $configall->{ $dataset . "_myrimatch_PSM" } = {
    class           => "Proteomics::Distiller::PSMDistiller",
    task_name       => $dataset,
    perform         => 1,
    target_dir      => "${target_dir}/$dataset/myrimatch",
    option          => "-e MyriMatch -t DTA",
    source_ref      => $dataset . "_myrimatch",
    proteomicstools => $proteomicstools,
    sh_direct       => 1,
    pbs             => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  };

  $configall->{ $dataset . "_myrimatch_PSM_rank2" } = {
    class           => "Proteomics::Distiller::PSMDistiller",
    task_name       => $dataset,
    perform         => 1,
    target_dir      => "${target_dir}/$dataset/myrimatch_PSM_rank2",
    option          => "-e MyriMatch -t DTA --rank2",
    source_ref      => $dataset . "_myrimatch",
    proteomicstools => $proteomicstools,
    sh_direct       => 1,
    pbs             => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  };

  push( @sections,     $dataset . "_myrimatch" );
  push( @psmdistiller, $dataset . "_myrimatch_PSM" );
  push( @psmdistiller, $dataset . "_myrimatch_PSM_rank2" );

  $configall->{ $dataset . "_msamanda" } = {
    class      => "Proteomics::Engine::Msamanda",
    task_name  => $dataset,
    perform    => 1,
    target_dir => "${target_dir}/$dataset/msamanda",
    option     => "",
    source     => $datasets->{$dataset}->{source},
    executable => $msamanda,
    cfgfile    => $datasets->{$dataset}->{msamanda_config_file},
    database   => $datasets->{$dataset}->{database},
    sh_direct  => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  };

  $configall->{ $dataset . "_msamanda_PSM" } = {
    class           => "Proteomics::Distiller::PSMDistiller",
    task_name       => $dataset,
    perform         => 1,
    target_dir      => "${target_dir}/$dataset/msamanda",
    option          => "-e MSAmanda -t DTA",
    source_ref      => $dataset . "_msamanda",
    proteomicstools => $proteomicstools,
    sh_direct       => 1,
    pbs             => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  };

  $configall->{ $dataset . "_msamanda_PSM_rank2" } = {
    class           => "Proteomics::Distiller::PSMDistiller",
    task_name       => $dataset,
    perform         => 1,
    target_dir      => "${target_dir}/$dataset/msamanda_PSM_rank2",
    option          => "-e MSAmanda -t DTA --rank2",
    source_ref      => $dataset . "_msamanda",
    proteomicstools => $proteomicstools,
    sh_direct       => 1,
    pbs             => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  };

  push( @msamanda, $dataset . "_msamanda" );
  push( @msamanda, $dataset . "_msamanda_PSM" );
  push( @msamanda, $dataset . "_msamanda_PSM_rank2" );
}

$configall->{sequencetask} = {
  class      => "CQS::SequenceTask",
  perform    => 1,
  target_dir => "${target_dir}/sequencetask",
  option     => "",
  source     => {
    step_1 => \@sections,
    step_2 => \@psmdistiller
  },
  sh_direct => 0,
  pbs       => {
    "email"    => $email,
    "nodes"    => "1:ppn=1",
    "walltime" => "72",
    "mem"      => "20gb"
  },
};

$configall->{msamanda_sequencetask} = {
  class      => "CQS::SequenceTask",
  perform    => 1,
  target_dir => "${target_dir}/msamanda_sequencetask",
  option     => "",
  source     => { step_1 => \@msamanda },
  sh_direct  => 1,
  pbs        => {
    "email"    => $email,
    "nodes"    => "1:ppn=1",
    "walltime" => "72",
    "mem"      => "20gb"
  },
};

#print Dumper($configall);
performConfig($configall);

1;
