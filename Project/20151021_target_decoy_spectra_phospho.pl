#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use Data::Dumper;

my $target_dir            = create_directory_or_die("/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho");
my $proteomicstools       = "/home/shengq1/proteomicstools/ProteomicsTools.exe";
my $email                 = "quanhu.sheng\@vanderbilt.edu";
my $msgf_jar              = "/scratch/cqs/shengq1/local/bin/MSGFPlus/MSGFPlus.jar";
my $msgf_mod_file         = "/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/config/msgf_mods.txt";
my $msgf_option           = "-t 20ppm -ti \"0,1\" -tda 0 -m 1 -inst 1 -e 1 -protocol 5 -ntt 2 -n 2 -minLength 7 -addFeatures 1";
my $target_database       = "/scratch/cqs/shengq1/proteomics/20151001_target_decoy_spectra/database/humanRefSeq_Version54_with_tryp.fasta";
my $target_decoy_database = "/scratch/cqs/shengq1/proteomics/20151001_target_decoy_spectra/database/humanRefSeq_Version54_with_tryp_DECOY.fasta";

my $buildsummary_msgf_target_file               = "/scratch/cqs/shengq1/proteomics/20151001_target_decoy_spectra/config/buildsummary_msgf_target.param";
my $buildsummary_msgf_target_decoy_file         = "/scratch/cqs/shengq1/proteomics/20151001_target_decoy_spectra/config/buildsummary_msgf_target_decoy.param";
my $buildsummary_msgf_target_unique2_file       = "/scratch/cqs/shengq1/proteomics/20151001_target_decoy_spectra/config/buildsummary_msgf_target_unique2.param";
my $buildsummary_msgf_target_decoy_unique2_file = "/scratch/cqs/shengq1/proteomics/20151001_target_decoy_spectra/config/buildsummary_msgf_target_decoy_unique2.param";
my $buildsummary_msgf_target_ratio2_file         = "/scratch/cqs/shengq1/proteomics/20151001_target_decoy_spectra/config/buildsummary_msgf_target_ratio2.param";
my $buildsummary_msgf_target_unique2_ratio2_file = "/scratch/cqs/shengq1/proteomics/20151001_target_decoy_spectra/config/buildsummary_msgf_target_unique2_ratio2.param";
my $buildsummary_bins                              = [ 1, 2, 3, 4, 5 ];

my $config = {
  general         => { task_name => "td" },
  files           => {},
  shift_precursor => {
    class           => "Proteomics::Format::PrecursorShiftProcessor",
    perform         => 1,
    target_dir      => "${target_dir}/shift_precursor",
    option          => "",
    source_ref      => "files",
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
  msgf_target => {
    class      => "Proteomics::Engine::MSGFPlus",
    perform    => 1,
    target_dir => "${target_dir}/msgf_target",
    option     => $msgf_option,
    source_ref => [ "files", "shift_precursor" ],
    msgf_jar   => $msgf_jar,
    mod_file   => $msgf_mod_file,
    database   => $target_database,
    sh_direct  => 0,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  msgf_target_psm => {
    class           => "Proteomics::Distiller::PSMDistiller",
    perform         => 1,
    target_dir      => "${target_dir}/msgf_target",
    option          => "-e MSGF -t DTA",
    source_ref      => "msgf_target",
    proteomicstools => $proteomicstools,
    sh_direct       => 1,
    pbs             => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "24",
      "mem"      => "10gb"
    },
  },
  msgf_target_accumulate_buildsummary => {
    class           => "Proteomics::Summary::AccumulateBuildSummary",
    perform         => 1,
    target_dir      => "${target_dir}/msgf_target_accumulate_buildsummary",
    option          => "",
    source_ref      => ["msgf_target"],
    parameter_file  => $buildsummary_msgf_target_file,
    proteomicstools => $proteomicstools,
    sh_direct       => 0,
    bin_size        => 5,
    bins            => $buildsummary_bins,
    pbs             => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  msgf_target_accumulate_buildsummary_ratio2 => {
    class           => "Proteomics::Summary::AccumulateBuildSummary",
    perform         => 1,
    target_dir      => "${target_dir}/msgf_target_accumulate_buildsummary_ratio2",
    option          => "",
    source_ref      => ["msgf_target"],
    parameter_file  => $buildsummary_msgf_target_ratio2_file,
    proteomicstools => $proteomicstools,
    sh_direct       => 0,
    bin_size        => 5,
    bins            => $buildsummary_bins,
    pbs             => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  msgf_target_accumulate_buildsummary_unique2 => {
    class           => "Proteomics::Summary::AccumulateBuildSummary",
    perform         => 1,
    target_dir      => "${target_dir}/msgf_target_accumulate_buildsummary_unique2",
    option          => "",
    source_ref      => ["msgf_target"],
    parameter_file  => $buildsummary_msgf_target_unique2_file,
    proteomicstools => $proteomicstools,
    sh_direct       => 0,
    bin_size        => 5,
    bins            => $buildsummary_bins,
    pbs             => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  msgf_target_accumulate_buildsummary_unique2_ratio2 => {
    class           => "Proteomics::Summary::AccumulateBuildSummary",
    perform         => 1,
    target_dir      => "${target_dir}/msgf_target_accumulate_buildsummary_unique2_ratio2",
    option          => "",
    source_ref      => ["msgf_target"],
    parameter_file  => $buildsummary_msgf_target_unique2_ratio2_file,
    proteomicstools => $proteomicstools,
    sh_direct       => 0,
    bin_size        => 5,
    bins            => $buildsummary_bins,
    pbs             => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  msgf_target_decoy => {
    class      => "Proteomics::Engine::MSGFPlus",
    perform    => 1,
    target_dir => "${target_dir}/msgf_target_decoy",
    option     => $msgf_option,
    source_ref => "files",
    msgf_jar   => $msgf_jar,
    mod_file   => $msgf_mod_file,
    database   => $target_decoy_database,
    sh_direct  => 0,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  msgf_target_decoy_psm => {
    class           => "Proteomics::Distiller::PSMDistiller",
    perform         => 1,
    target_dir      => "${target_dir}/msgf_target_decoy",
    option          => "-e MSGF -t DTA",
    source_ref      => "msgf_target_decoy",
    proteomicstools => $proteomicstools,
    sh_direct       => 1,
    pbs             => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "24",
      "mem"      => "10gb"
    },
  },
  msgf_target_decoy_accumulate_buildsummary => {
    class           => "Proteomics::Summary::AccumulateBuildSummary",
    perform         => 1,
    target_dir      => "${target_dir}/msgf_target_decoy_accumulate_buildsummary",
    option          => "",
    source_ref      => ["msgf_target_decoy"],
    parameter_file  => $buildsummary_msgf_target_decoy_file,
    proteomicstools => $proteomicstools,
    bin_size        => 5,
    bins            => $buildsummary_bins,
    sh_direct       => 0,
    pbs             => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  msgf_target_decoy_accumulate_buildsummary_unique2 => {
    class           => "Proteomics::Summary::AccumulateBuildSummary",
    perform         => 1,
    target_dir      => "${target_dir}/msgf_target_decoy_accumulate_buildsummary_unique2",
    option          => "",
    source_ref      => ["msgf_target_decoy"],
    parameter_file  => $buildsummary_msgf_target_decoy_unique2_file,
    proteomicstools => $proteomicstools,
    bin_size        => 5,
    bins            => $buildsummary_bins,
    sh_direct       => 0,
    pbs             => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  sequencetask => {
    class      => "CQS::SequenceTask",
    perform    => 1,
    target_dir => "${target_dir}/sequencetask",
    option     => "",
    source     => {
      T1_individual => [
        "shift_precursor", "msgf_target", "msgf_target_psm", "msgf_target_decoy", "msgf_target_decoy_psm"
      ],
    },
    sh_direct => 0,
    pbs       => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  }
};

performConfig($config);

1;
