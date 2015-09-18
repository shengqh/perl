#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use Data::Dumper;

my $target_dir      = create_directory_or_die("/scratch/cqs/shengq1/proteomics/20150918_shifted_precursor");
my $proteomicstools = "/home/shengq1/proteomicstools/ProteomicsTools.exe";
my $email           = "quanhu.sheng\@vanderbilt.edu";

my $config = {
  general => { task_name => "ShiftedTargetDecoy" },
  files   => {
    "Elite_CIDIT_Human_1" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/Elite_CIDIT_Human/mgf/20141017Test_DJMa_Cell_04.mgf"],
    "Elite_CIDIT_Human_2" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/Elite_CIDIT_Human/mgf/20141017Test_DJMa_Cell_05.mgf"],
    "Elite_CIDIT_Human_3" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20150608_shifted_precursor/Elite_CIDIT_Human/mgf/20141017Test_DJMa_Cell_06.mgf"],
  },
  shift_precursor => {
    class           => "Proteomics::Format::PrecursorShiftProcessor",
    perform         => 1,
    target_dir      => "${target_dir}/shift_precursor",
    option          => "",
    source_ref      => "files",
    proteomicstools => $proteomicstools,
    sh_direct       => 0,
    pbs             => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  }
};
performConfig($config);

1;
