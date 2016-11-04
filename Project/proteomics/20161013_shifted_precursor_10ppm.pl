#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use Data::Dumper;
use Hash::Merge qw( merge );

#my $target_dir      = "E:/temp";
my $target_dir      = create_directory_or_die("/scratch/cqs/shengq1/proteomics/20161013_shifted_precursor_10ppm");
my $msgf_jar        = "/scratch/cqs/shengq1/local/bin/MSGFPlus/MSGFPlus.jar";
my $proteomicstools = "/home/shengq1/proteomicstools/ProteomicsTools.exe";
my $mod_file        = "/scratch/cqs/shengq1/proteomics/20161013_shifted_precursor_10ppm/parameters/msgf_mod.txt";
my $database_human  = "/scratch/cqs/shengq1/proteomics/database/uniprot_human_201609.REVERSED.fasta";
my $database_yeast  = "/scratch/cqs/shengq1/proteomics/database/uniprot_yeast_201609.REVERSED.fasta";

#my $database_ecoli  = "/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/database/rev_Ecoli_uniprot_v20120613_P4431.fasta";
my $email = "quanhu.sheng\@vanderbilt.edu";

my $buildsummary_msgf_human_file  = "/scratch/cqs/shengq1/proteomics/20161013_shifted_precursor_10ppm/parameters/buildsummary_msgf_human.param";
my $buildsummary_msgf_yeast_file  = "/scratch/cqs/shengq1/proteomics/20161013_shifted_precursor_10ppm/parameters/buildsummary_msgf_yeast.param";
my $buildsummary_comet_human_file = "/scratch/cqs/shengq1/proteomics/20161013_shifted_precursor_10ppm/parameters/buildsummary_comet_human.param";
my $buildsummary_comet_yeast_file = "/scratch/cqs/shengq1/proteomics/20161013_shifted_precursor_10ppm/parameters/buildsummary_comet_yeast.param";

my $datasets = {
  Elite_CIDIT_Human => {
    source => {
      "Elite_CIDIT_Human" => [
        "/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/Elite_CIDIT_Human/mgf/20141017Test_DJMa_Cell_04.mgf",
        "/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/Elite_CIDIT_Human/mgf/20141017Test_DJMa_Cell_05.mgf",
        "/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/Elite_CIDIT_Human/mgf/20141017Test_DJMa_Cell_06.mgf"
      ],
    },
    database                       => $database_human,
    MSGF_option                    => "-t 10ppm -ti \"0,1\" -tda 0 -m 1 -inst 0 -e 1 -protocol 5 -ntt 2 -n 2 -addFeatures 1",
    comet_config_file              => "/scratch/cqs/shengq1/proteomics/20161013_shifted_precursor_10ppm/parameters/comet-it.params",
    myrimatch_config_file          => "/scratch/cqs/shengq1/proteomics/20161013_shifted_precursor_10ppm/parameters/myrimatch-it.config",
    msamanda_config_file           => "/scratch/cqs/shengq1/proteomics/20161013_shifted_precursor_10ppm/parameters/msamanda-it.xml",
    buildsummary_msgf_config_file  => $buildsummary_msgf_human_file,
    buildsummary_comet_config_file => $buildsummary_comet_human_file,
  },
  Fusion_CIDIT_Human => {
    source => {
      "Fusion_CIDIT_Human" => [
        "/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/Fusion_CIDIT_Human/mgf/B01_03_14050113_HCC_Hela_Qu_CID_IT_pepID.mgf",
        "/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/Fusion_CIDIT_Human/mgf/B11_03_140524_HCC_Hela_Qu_CID_IT_pepID.mgf",
        "/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/Fusion_CIDIT_Human/mgf/B13_03_140612_HCC3_Hela_Qu_CID_IT_pepID.mgf"
      ],
    },
    database                       => $database_human,
    MSGF_option                    => "-t 10ppm -ti \"0,1\" -tda 0 -m 1 -inst 0 -e 1 -protocol 5 -ntt 2 -n 2 -addFeatures 1",
    comet_config_file              => "/scratch/cqs/shengq1/proteomics/20161013_shifted_precursor_10ppm/parameters/comet-it.params",
    myrimatch_config_file          => "/scratch/cqs/shengq1/proteomics/20161013_shifted_precursor_10ppm/parameters/myrimatch-it.config",
    msamanda_config_file           => "/scratch/cqs/shengq1/proteomics/20161013_shifted_precursor_10ppm/parameters/msamanda-it.xml",
    buildsummary_msgf_config_file  => $buildsummary_msgf_human_file,
    buildsummary_comet_config_file => $buildsummary_comet_human_file,
  },
  Fusion_HCDIT_Yeast => {
    source => {
      "Fusion_HCDIT_Yeast" => [
        "/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/Fusion_HCDIT_Yeast/mgf/10sep2013_yeast_control_1.mgf",
        "/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/Fusion_HCDIT_Yeast/mgf/10sep2013_yeast_control_2.mgf",
        "/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/Fusion_HCDIT_Yeast/mgf/10sep2013_yeast_control_4.mgf"
      ],
    },
    database                       => $database_yeast,
    MSGF_option                    => "-t 10ppm -ti \"0,1\" -tda 0 -m 3 -inst 0 -e 1 -protocol 5 -ntt 2 -n 2 -addFeatures 1",
    comet_config_file              => "/scratch/cqs/shengq1/proteomics/20161013_shifted_precursor_10ppm/parameters/comet-it.params",
    myrimatch_config_file          => "/scratch/cqs/shengq1/proteomics/20161013_shifted_precursor_10ppm/parameters/myrimatch-it.config",
    msamanda_config_file           => "/scratch/cqs/shengq1/proteomics/20161013_shifted_precursor_10ppm/parameters/msamanda-it.xml",
    buildsummary_msgf_config_file  => $buildsummary_msgf_yeast_file,
    buildsummary_comet_config_file => $buildsummary_comet_yeast_file,
  },
  Fusion_HCDOT_Human => {
    source => {
      "Fusion_HCDOT_Human" => [
        "/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/Fusion_HCDOT_Human/mgf/B00_01_140613_zhuxu_Qu_HCD_OT_2hr.mgf",
        "/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/Fusion_HCDOT_Human/mgf/B00_01_140621_zhuxu_2_Qu_HCD_OT_2hr.mgf",
        "/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/Fusion_HCDOT_Human/mgf/B00_01_140622_zhuxu_2_Qu_HCD_OT_2hr.mgf"
      ],
    },
    database                       => $database_human,
    MSGF_option                    => "-t 10ppm -ti \"0,1\" -tda 0 -m 3 -inst 3 -e 1 -protocol 5 -ntt 2 -n 2 -addFeatures 1",
    comet_config_file              => "/scratch/cqs/shengq1/proteomics/20161013_shifted_precursor_10ppm/parameters/comet-ot.params",
    myrimatch_config_file          => "/scratch/cqs/shengq1/proteomics/20161013_shifted_precursor_10ppm/parameters/myrimatch-ot.config",
    msamanda_config_file           => "/scratch/cqs/shengq1/proteomics/20161013_shifted_precursor_10ppm/parameters/msamanda-ot.xml",
    buildsummary_msgf_config_file  => $buildsummary_msgf_human_file,
    buildsummary_comet_config_file => $buildsummary_comet_human_file,
  },
  QExactive_HCDOT_Human => {
    source => {
      "QExactive_HCDOT_Human" => [
        "/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/QExactive_HCDOT_Human/mgf/QExactive1.mgf",
        "/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/QExactive_HCDOT_Human/mgf/QExactive2.mgf",
        "/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/QExactive_HCDOT_Human/mgf/QExactive3.mgf"
      ],
    },
    database                       => $database_human,
    MSGF_option                    => "-t 10ppm -ti \"0,1\" -tda 0 -m 3 -inst 3 -e 1 -protocol 5 -ntt 2 -n 2 -addFeatures 1",
    comet_config_file              => "/scratch/cqs/shengq1/proteomics/20161013_shifted_precursor_10ppm/parameters/comet-ot.params",
    myrimatch_config_file          => "/scratch/cqs/shengq1/proteomics/20161013_shifted_precursor_10ppm/parameters/myrimatch-ot.config",
    msamanda_config_file           => "/scratch/cqs/shengq1/proteomics/20161013_shifted_precursor_10ppm/parameters/msamanda-ot.xml",
    buildsummary_msgf_config_file  => $buildsummary_msgf_human_file,
    buildsummary_comet_config_file => $buildsummary_comet_human_file,
  }
};

my $configall = { general => { task_name => "td" }, };

my @individual = ();
my @summary    = ();

my @deltas = (7.014395);

for my $dataset ( sort keys %{$datasets} ) {
  my $def   = $datasets->{$dataset};
  my $shift = {
    "${dataset}_source"      => $def->{source},
    "${dataset}_MSGF_source" => {
      class      => "Proteomics::Engine::MSGFPlus",
      perform    => 1,
      target_dir => "${target_dir}/${dataset}_MSGF",
      option     => $datasets->{$dataset}->{MSGF_option},
      source_ref => "${dataset}_source",
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
      perform         => 1,
      target_dir      => "${target_dir}/${dataset}_MSGF",
      option          => "-e MSGF -t DTA",
      source_ref      => "${dataset}_MSGF_source",
      proteomicstools => $proteomicstools,
      sh_direct       => 1,
      pbs             => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "72",
        "mem"      => "40gb"
      },
    },
    "${dataset}_comet_source" => {
      class           => "Proteomics::Engine::Comet",
      perform         => 1,
      target_dir      => "${target_dir}/${dataset}_comet",
      option          => "",
      source_ref      => "${dataset}_source",
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
      perform         => 1,
      target_dir      => "${target_dir}/${dataset}_comet",
      option          => "-e Comet -t DTA",
      source_ref      => "${dataset}_comet_source",
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

  push( @individual, ( "${dataset}_MSGF_source", "${dataset}_MSGF_PSM", "${dataset}_comet_source", "${dataset}_comet_PSM" ) );
  $configall = merge( $configall, $shift );

  for my $deltaValue (@deltas) {
    my $delta       = int($deltaValue);
    my $deltaconfig = {
      "${dataset}_shift_${delta}" => {
        suffix          => "_${delta}",
        class           => "Proteomics::Format::PrecursorShiftProcessor",
        perform         => 1,
        target_dir      => $target_dir . "/${dataset}_shift",
        option          => "",
        source_ref      => "${dataset}_source",
        proteomicstools => $proteomicstools,
        shiftmass       => $deltaValue,
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
      "${dataset}_MSGF_${delta}" => {
        suffix     => "_${delta}",
        class      => "Proteomics::Engine::MSGFPlus",
        perform    => 1,
        target_dir => "${target_dir}/${dataset}_MSGF",
        option     => $datasets->{$dataset}->{MSGF_option},
        source_ref => "${dataset}_shift_${delta}",
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
      "${dataset}_MSGF_${delta}_msgf2pin" => {
        task_name     => "${dataset}_msgf2pin",
        class         => "Proteomics::Engine::Msgf2Pin",
        perform       => 1,
        target_dir    => "${target_dir}/${dataset}_MSGF_${delta}_percolator",
        option        => "",
        source_ref    => [ "${dataset}_MSGF_source", "${dataset}_MSGF_${delta}" ],
        decoy_pattern => "REVERSED_",
        sh_direct     => 0,
        pbs           => {
          "email"    => $email,
          "nodes"    => "1:ppn=1",
          "walltime" => "72",
          "mem"      => "40gb"
        },
      },
      "${dataset}_MSGF_${delta}_percolator" => {
        task_name     => "${dataset}_percolator",
        class         => "Proteomics::Engine::Percolator",
        perform       => 1,
        target_dir    => "${target_dir}/${dataset}_MSGF_${delta}_percolator",
        option        => "-U",
        source_ref    => "${dataset}_MSGF_${delta}_msgf2pin",
        sh_direct     => 1,
        pbs           => {
          "email"    => $email,
          "nodes"    => "1:ppn=1",
          "walltime" => "72",
          "mem"      => "40gb"
        },
      },
      "${dataset}_MSGF_${delta}_PSM" => {
        suffix          => "_${delta}",
        class           => "Proteomics::Distiller::PSMDistiller",
        perform         => 1,
        target_dir      => "${target_dir}/${dataset}_MSGF",
        option          => "-e MSGF -t DTA",
        source_ref      => "${dataset}_MSGF_${delta}",
        proteomicstools => $proteomicstools,
        sh_direct       => 1,
        pbs             => {
          "email"    => $email,
          "nodes"    => "1:ppn=1",
          "walltime" => "72",
          "mem"      => "40gb"
        },
      },

      "${dataset}_MSGF_${delta}_buildsummary" => {
        task_name       => "${dataset}_MSGF_${delta}",
        class           => "Proteomics::Summary::BuildSummary",
        perform         => 1,
        target_dir      => "${target_dir}/buildsummary",
        option          => "",
        source_ref      => [ "${dataset}_MSGF_source", "${dataset}_MSGF_${delta}" ],
        parameter_file  => $datasets->{$dataset}->{buildsummary_msgf_config_file},
        proteomicstools => $proteomicstools,
        sh_direct       => 0,
        pbs             => {
          "email"    => $email,
          "nodes"    => "1:ppn=1",
          "walltime" => "72",
          "mem"      => "40gb"
        },
      },

      "${dataset}_comet_${delta}" => {
        suffix          => "_${delta}",
        class           => "Proteomics::Engine::Comet",
        perform         => 1,
        target_dir      => "${target_dir}/${dataset}_comet",
        option          => "",
        source_ref      => "${dataset}_shift_${delta}",
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
      "${dataset}_comet_${delta}_PSM" => {
        suffix          => "_${delta}",
        class           => "Proteomics::Distiller::PSMDistiller",
        perform         => 1,
        target_dir      => "${target_dir}/${dataset}_comet",
        option          => "-e Comet -t DTA",
        source_ref      => "${dataset}_comet_${delta}",
        proteomicstools => $proteomicstools,
        sh_direct       => 1,
        pbs             => {
          "email"    => $email,
          "nodes"    => "1:ppn=1",
          "walltime" => "72",
          "mem"      => "40gb"
        },
      },
      "${dataset}_comet_${delta}_buildsummary" => {
        task_name       => "${dataset}_comet_${delta}",
        class           => "Proteomics::Summary::BuildSummary",
        perform         => 1,
        target_dir      => "${target_dir}/buildsummary",
        option          => "",
        source_ref      => [ "${dataset}_comet_source", "${dataset}_comet_${delta}" ],
        parameter_file  => $datasets->{$dataset}->{buildsummary_comet_config_file},
        proteomicstools => $proteomicstools,
        sh_direct       => 0,
        pbs             => {
          "email"    => $email,
          "nodes"    => "1:ppn=1",
          "walltime" => "72",
          "mem"      => "40gb"
        },
      },
    };

    push( @individual, ( "${dataset}_shift_${delta}", "${dataset}_MSGF_${delta}", "${dataset}_MSGF_${delta}_PSM", "${dataset}_comet_${delta}", "${dataset}_comet_${delta}_PSM" ) );
    push( @summary, ( "${dataset}_MSGF_${delta}_msgf2pin","${dataset}_MSGF_${delta}_percolator", "${dataset}_MSGF_${delta}_buildsummary", "${dataset}_comet_${delta}_buildsummary" ) );

    $configall = merge( $configall, $deltaconfig );
  }
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
