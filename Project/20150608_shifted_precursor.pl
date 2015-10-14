#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use Data::Dumper;
use Hash::Merge qw( merge );

my $target_dir      = "/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor";
my $msgf_jar        = "/scratch/cqs/shengq1/local/bin/MSGFPlus/MSGFPlus.jar";
my $msamanda        = "/scratch/cqs/shengq1/local/bin/MSAmanda_Standalone_LinuxMac_1.0.0.4485/MSAmanda.exe";
my $proteomicstools = "/home/shengq1/proteomicstools/ProteomicsTools.exe";
my $mod_file        = "/scratch/cqs/shengq1/local/bin/MSGFPlus/Mods.txt";
my $database_human  = "${target_dir}/database/rev_Human_uniprot_sprot_v20120613.fasta";
my $database_yeast  = "${target_dir}/database/rev_Yeast_uniprot_v20120613.fasta";
my $database_ecoli  = "${target_dir}/database/rev_Ecoli_uniprot_v20120613_P4431.fasta";
my $email           = "quanhu.sheng\@vanderbilt.edu";

my $buildsummary_msgf_human_file = "/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/parameters/buildsummary_msgf_human.param";
my $buildsummary_msgf_yeast_file = "/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/parameters/buildsummary_msgf_yeast.param";

my $datasets = {
  Elite_CIDIT_Human => {
    source => {
      "Elite_CIDIT_Human_1" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/Elite_CIDIT_Human/mgf/20141017Test_DJMa_Cell_04.mgf"],
      "Elite_CIDIT_Human_2" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/Elite_CIDIT_Human/mgf/20141017Test_DJMa_Cell_05.mgf"],
      "Elite_CIDIT_Human_3" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/Elite_CIDIT_Human/mgf/20141017Test_DJMa_Cell_06.mgf"],
    },
    database              => $database_human,
    MSGF_option           => "-t 20ppm -ti \"0,1\" -tda 0 -m 1 -inst 0 -e 1 -protocol 5 -ntt 2 -n 2 -addFeatures 1",
    comet_config_file     => "/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/parameters/comet-it.params",
    myrimatch_config_file => "/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/parameters/myrimatch-it.config",
    msamanda_config_file  => "/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/parameters/msamanda-it.xml",
  },
  Fusion_CIDIT_Human => {
    source => {
      "Fusion_CIDIT_Human_1" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/Fusion_CIDIT_Human/mgf/B01_03_14050113_HCC_Hela_Qu_CID_IT_pepID.mgf"],
      "Fusion_CIDIT_Human_2" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/Fusion_CIDIT_Human/mgf/B11_03_140524_HCC_Hela_Qu_CID_IT_pepID.mgf"],
      "Fusion_CIDIT_Human_3" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/Fusion_CIDIT_Human/mgf/B13_03_140612_HCC3_Hela_Qu_CID_IT_pepID.mgf"],
    },
    database              => $database_human,
    MSGF_option           => "-t 20ppm -ti \"0,1\" -tda 0 -m 1 -inst 0 -e 1 -protocol 5 -ntt 2 -n 2 -addFeatures 1",
    comet_config_file     => "/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/parameters/comet-it.params",
    myrimatch_config_file => "/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/parameters/myrimatch-it.config",
    msamanda_config_file  => "/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/parameters/msamanda-it.xml",
  },
  Fusion_HCDIT_Yeast => {
    source => {
      "Fusion_HCDIT_Yeast_1" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/Fusion_HCDIT_Yeast/mgf/10sep2013_yeast_control_1.mgf"],
      "Fusion_HCDIT_Yeast_2" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/Fusion_HCDIT_Yeast/mgf/10sep2013_yeast_control_2.mgf"],
      "Fusion_HCDIT_Yeast_3" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/Fusion_HCDIT_Yeast/mgf/10sep2013_yeast_control_4.mgf"],
    },
    database              => $database_yeast,
    MSGF_option           => "-t 20ppm -ti \"0,1\" -tda 0 -m 3 -inst 0 -e 1 -protocol 5 -ntt 2 -n 2 -addFeatures 1",
    comet_config_file     => "/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/parameters/comet-it.params",
    myrimatch_config_file => "/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/parameters/myrimatch-it.config",
    msamanda_config_file  => "/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/parameters/msamanda-it.xml",
  },
  Fusion_HCDOT_Human => {
    source => {
      "Fusion_HCDOT_Human_1" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/Fusion_HCDOT_Human/mgf/B00_01_140613_zhuxu_Qu_HCD_OT_2hr.mgf"],
      "Fusion_HCDOT_Human_2" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/Fusion_HCDOT_Human/mgf/B00_01_140621_zhuxu_2_Qu_HCD_OT_2hr.mgf"],
      "Fusion_HCDOT_Human_3" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/Fusion_HCDOT_Human/mgf/B00_01_140622_zhuxu_2_Qu_HCD_OT_2hr.mgf"],
    },
    database              => $database_human,
    MSGF_option           => "-t 20ppm -ti \"0,1\" -tda 0 -m 3 -inst 3 -e 1 -protocol 5 -ntt 2 -n 2 -addFeatures 1",
    comet_config_file     => "/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/parameters/comet-ot.params",
    myrimatch_config_file => "/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/parameters/myrimatch-ot.config",
    msamanda_config_file  => "/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/parameters/msamanda-ot.xml",
  },
  QExactive_HCDOT_Human => {
    source => {
      "QExactive_HCDOT_Human_1" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/QExactive_HCDOT_Human/mgf/QExactive1.mgf"],
      "QExactive_HCDOT_Human_2" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/QExactive_HCDOT_Human/mgf/QExactive2.mgf"],
      "QExactive_HCDOT_Human_3" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/QExactive_HCDOT_Human/mgf/QExactive3.mgf"],
    },
    database              => $database_human,
    MSGF_option           => "-t 20ppm -ti \"0,1\" -tda 0 -m 3 -inst 3 -e 1 -protocol 5 -ntt 2 -n 2 -addFeatures 1",
    comet_config_file     => "/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/parameters/comet-ot.params",
    myrimatch_config_file => "/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/parameters/myrimatch-ot.config",
    msamanda_config_file  => "/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/parameters/msamanda-ot.xml",
  }
};

my $configall = { general => { task_name => "td" }, };

my @individual = ();
my @summary    = ();

for my $dataset ( sort keys %{$datasets} ) {
  my $def   = $datasets->{$dataset};
  my $shift = {
    "${dataset}_source"         => $def->{source},
    "${dataset}_shift_minus_15" => {
      class           => "Proteomics::Format::PrecursorShiftProcessor",
      perform         => 0,
      target_dir      => "${target_dir}/${dataset}/shift_minus_15",
      option          => "",
      source_ref      => "${dataset}_source",
      proteomicstools => $proteomicstools,
      shiftmass       => -15,
      shiftscan       => 10000000,
      titleformat     => "DTA",
      sh_direct       => 1,
      pbs             => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "72",
        "mem"      => "10gb"
      },
    },
    "${dataset}_shift_minus_10" => {
      class           => "Proteomics::Format::PrecursorShiftProcessor",
      perform         => 0,
      target_dir      => "${target_dir}/${dataset}/shift_minus_10",
      option          => "",
      source_ref      => "${dataset}_source",
      proteomicstools => $proteomicstools,
      shiftmass       => -10,
      shiftscan       => 10000000,
      titleformat     => "DTA",
      sh_direct       => 1,
      pbs             => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "72",
        "mem"      => "10gb"
      },
    },
    "${dataset}_shift_minus_5" => {
      class           => "Proteomics::Format::PrecursorShiftProcessor",
      perform         => 0,
      target_dir      => "${target_dir}/${dataset}/shift_minus_5",
      option          => "",
      source_ref      => "${dataset}_source",
      proteomicstools => $proteomicstools,
      shiftmass       => -5,
      shiftscan       => 10000000,
      titleformat     => "DTA",
      sh_direct       => 1,
      pbs             => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "72",
        "mem"      => "10gb"
      },
    },
    "${dataset}_shift_minus_0.5" => {
      class           => "Proteomics::Format::PrecursorShiftProcessor",
      perform         => 0,
      target_dir      => "${target_dir}/${dataset}/shift_minus_0.5",
      option          => "",
      source_ref      => "${dataset}_source",
      proteomicstools => $proteomicstools,
      shiftmass       => -0.5,
      shiftscan       => 10000000,
      titleformat     => "DTA",
      sh_direct       => 1,
      pbs             => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "72",
        "mem"      => "10gb"
      },
    },
    "${dataset}_shift_plus_0.5" => {
      class           => "Proteomics::Format::PrecursorShiftProcessor",
      perform         => 0,
      target_dir      => "${target_dir}/${dataset}/shift_plus_0.5",
      option          => "",
      source_ref      => "${dataset}_source",
      proteomicstools => $proteomicstools,
      shiftmass       => 0.5,
      shiftscan       => 10000000,
      titleformat     => "DTA",
      sh_direct       => 1,
      pbs             => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "72",
        "mem"      => "10gb"
      },
    },
    "${dataset}_shift_plus_5" => {
      class           => "Proteomics::Format::PrecursorShiftProcessor",
      perform         => 0,
      target_dir      => "${target_dir}/${dataset}/shift_plus_5",
      option          => "",
      source_ref      => "${dataset}_source",
      proteomicstools => $proteomicstools,
      shiftmass       => 5,
      shiftscan       => 10000000,
      titleformat     => "DTA",
      sh_direct       => 1,
      pbs             => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "72",
        "mem"      => "10gb"
      },
    },
    "${dataset}_shift_plus_10" => {
      class           => "Proteomics::Format::PrecursorShiftProcessor",
      perform         => 0,
      target_dir      => "${target_dir}/${dataset}/shift_plus_10",
      option          => "",
      source_ref      => "${dataset}_source",
      proteomicstools => $proteomicstools,
      shiftmass       => 10,
      shiftscan       => 10000000,
      titleformat     => "DTA",
      sh_direct       => 1,
      pbs             => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "72",
        "mem"      => "10gb"
      },
    },
    "${dataset}_shift_plus_15" => {
      class           => "Proteomics::Format::PrecursorShiftProcessor",
      perform         => 0,
      target_dir      => "${target_dir}/${dataset}/shift_plus_15",
      option          => "",
      source_ref      => "${dataset}_source",
      proteomicstools => $proteomicstools,
      shiftmass       => 15,
      shiftscan       => 10000000,
      titleformat     => "DTA",
      sh_direct       => 1,
      pbs             => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "72",
        "mem"      => "10gb"
      },
    },
  };

  $configall = merge( $configall, $shift );
  my @mgf = [
    "${dataset}_source",       "${dataset}_shift_minus_15", "${dataset}_shift_minus_10", "${dataset}_shift_minus_5", "${dataset}_shift_minus_0.5", "${dataset}_shift_plus_0.5",
    "${dataset}_shift_plus_5", "${dataset}_shift_plus_10",  "${dataset}_shift_plus_15"
  ];
  push(
    @individual,
    (
      "${dataset}_shift_minus_15", "${dataset}_shift_minus_10", "${dataset}_shift_minus_5", "${dataset}_shift_minus_0.5",
      "${dataset}_shift_plus_0.5", "${dataset}_shift_plus_5",   "${dataset}_shift_plus_10", "${dataset}_shift_plus_15"
    )
  );
  my $search = {
    "${dataset}_MSGF" => {
      class      => "Proteomics::Engine::MSGFPlus",
      task_name  => $dataset,
      perform    => 0,
      target_dir => "${target_dir}/$dataset/MSGF",
      option     => $datasets->{$dataset}->{MSGF_option},
      source_ref => @mgf,
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
    },
    "${dataset}_MSGF_PSM" => {
      class           => "Proteomics::Distiller::PSMDistiller",
      task_name       => $dataset,
      perform         => 0,
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
    },

    "${dataset}_MSGF_buildsummary_-10daltons" => {
      task_name       => "${dataset}_MSGF_-10daltons",
      class           => "Proteomics::Summary::BuildSummary",
      perform         => 1,
      target_dir      => "${target_dir}/$dataset/MSGF_buildsummary",
      option          => "",
      source_ref      => [ "${dataset}_MSGF", "_??.msgf.mzid", "${dataset}_MSGF", "-10daltons.msgf.mzid", ],
      parameter_file  => $buildsummary_msgf_human_file,
      proteomicstools => $proteomicstools,
      sh_direct       => 0,
      pbs             => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "72",
        "mem"      => "40gb"
      },
    },

    "${dataset}_comet" => {
      class           => "Proteomics::Engine::Comet",
      task_name       => $dataset,
      perform         => 0,
      target_dir      => "${target_dir}/$dataset/comet",
      option          => "",
      source_ref      => @mgf,
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
    },
    "${dataset}_comet_PSM" => {
      class           => "Proteomics::Distiller::PSMDistiller",
      task_name       => $dataset,
      perform         => 0,
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
    },
    "${dataset}_myrimatch" => {
      class      => "Proteomics::Engine::Myrimatch",
      task_name  => $dataset,
      perform    => 0,
      target_dir => "${target_dir}/$dataset/myrimatch",
      option     => "",
      source_ref => @mgf,
      cfgfile    => $datasets->{$dataset}->{myrimatch_config_file},
      database   => $datasets->{$dataset}->{database},
      sh_direct  => 0,
      pbs        => {
        "email"    => $email,
        "nodes"    => "1:ppn=8",
        "walltime" => "72",
        "mem"      => "40gb"
      },
    },
    "${dataset}_myrimatch_PSM" => {
      class           => "Proteomics::Distiller::PSMDistiller",
      task_name       => $dataset,
      perform         => 0,
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
    },
  };

  push( @individual, ( "${dataset}_MSGF", "${dataset}_MSGF_PSM", "${dataset}_comet", "${dataset}_comet_PSM", "${dataset}_myrimatch", "${dataset}_myrimatch_PSM", ) );
  push( @summary, ("${dataset}_MSGF_buildsummary_-10daltons") );

  $configall = merge( $configall, $search );
}

$configall->{sequencetask} = {
  class      => "CQS::SequenceTask",
  perform    => 1,
  target_dir => $target_dir . "/sequencetask",
  option     => "",
  source     => {
    step1 => \@individual,
    step2 => \@summary,
  },
  sh_direct => 0,
  pbs       => {
    "email"    => $email,
    "nodes"    => "1:ppn=8",
    "walltime" => "72",
    "mem"      => "40gb"
  },
};

#print Dumper($configall);
performConfig($configall);

1;

