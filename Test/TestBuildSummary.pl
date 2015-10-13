#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use Data::Dumper;

my $target_dir      = create_directory_or_die("E:/shengquanhu/projects/20151001_target_decoy_spectra");
my $proteomicstools = "E:/sqh/programs/csharp/OmicsLabCSharp/RCPA.Tools/bin/x64/Release/proteomicstools.exe";
my $email           = "quanhu.sheng\@vanderbilt.edu";

my $buildsummary_msgf_target_file = "E:/shengquanhu/projects/20151001_target_decoy_spectra/config/buildsummary_msgf_target.param";

my $config = {
  general => { task_name => "ShiftedTargetDecoy" },
  files   => {
    "TCGA-AA-3534-01A-22" => [
      "/dors/bioinfo/zhanglab/tcga_colon_proteomics/TCGA-AA-3534-01A-22_W_VU_20130208_A0218_10B_R_FR01.mz5",
      "/dors/bioinfo/zhanglab/tcga_colon_proteomics/TCGA-AA-3534-01A-22_W_VU_20130208_A0218_10B_R_FR02.mz5",
      "/dors/bioinfo/zhanglab/tcga_colon_proteomics/TCGA-AA-3534-01A-22_W_VU_20130208_A0218_10B_R_FR03.mz5",
      "/dors/bioinfo/zhanglab/tcga_colon_proteomics/TCGA-AA-3534-01A-22_W_VU_20130208_A0218_10B_R_FR04.mz5",
      "/dors/bioinfo/zhanglab/tcga_colon_proteomics/TCGA-AA-3534-01A-22_W_VU_20130208_A0218_10B_R_FR05.mz5",
      "/dors/bioinfo/zhanglab/tcga_colon_proteomics/TCGA-AA-3534-01A-22_W_VU_20130208_A0218_10B_R_FR06.mz5",
      "/dors/bioinfo/zhanglab/tcga_colon_proteomics/TCGA-AA-3534-01A-22_W_VU_20130208_A0218_10B_R_FR07.mz5",
      "/dors/bioinfo/zhanglab/tcga_colon_proteomics/TCGA-AA-3534-01A-22_W_VU_20130208_A0218_10B_R_FR08.mz5",
      "/dors/bioinfo/zhanglab/tcga_colon_proteomics/TCGA-AA-3534-01A-22_W_VU_20130208_A0218_10B_R_FR09.mz5",
      "/dors/bioinfo/zhanglab/tcga_colon_proteomics/TCGA-AA-3534-01A-22_W_VU_20130209_A0218_10B_R_FR10.mz5",
      "/dors/bioinfo/zhanglab/tcga_colon_proteomics/TCGA-AA-3534-01A-22_W_VU_20130209_A0218_10B_R_FR11.mz5",
      "/dors/bioinfo/zhanglab/tcga_colon_proteomics/TCGA-AA-3534-01A-22_W_VU_20130209_A0218_10B_R_FR12.mz5",
      "/dors/bioinfo/zhanglab/tcga_colon_proteomics/TCGA-AA-3534-01A-22_W_VU_20130209_A0218_10B_R_FR13.mz5",
      "/dors/bioinfo/zhanglab/tcga_colon_proteomics/TCGA-AA-3534-01A-22_W_VU_20130209_A0218_10B_R_FR14.mz5",
      "/dors/bioinfo/zhanglab/tcga_colon_proteomics/TCGA-AA-3534-01A-22_W_VU_20130209_A0218_10B_R_FR15.mz5"
    ],
    "TCGA-AA-3552-01A-22" => [
      "/dors/bioinfo/zhanglab/tcga_colon_proteomics/TCGA-AA-3552-01A-22_W_VU_20130125_A0218_9F_R_FR01.mz5",
      "/dors/bioinfo/zhanglab/tcga_colon_proteomics/TCGA-AA-3552-01A-22_W_VU_20130125_A0218_9F_R_FR02.mz5",
      "/dors/bioinfo/zhanglab/tcga_colon_proteomics/TCGA-AA-3552-01A-22_W_VU_20130125_A0218_9F_R_FR03.mz5",
      "/dors/bioinfo/zhanglab/tcga_colon_proteomics/TCGA-AA-3552-01A-22_W_VU_20130125_A0218_9F_R_FR04.mz5",
      "/dors/bioinfo/zhanglab/tcga_colon_proteomics/TCGA-AA-3552-01A-22_W_VU_20130125_A0218_9F_R_FR05.mz5",
      "/dors/bioinfo/zhanglab/tcga_colon_proteomics/TCGA-AA-3552-01A-22_W_VU_20130125_A0218_9F_R_FR06.mz5",
      "/dors/bioinfo/zhanglab/tcga_colon_proteomics/TCGA-AA-3552-01A-22_W_VU_20130126_A0218_9F_R_FR07.mz5",
      "/dors/bioinfo/zhanglab/tcga_colon_proteomics/TCGA-AA-3552-01A-22_W_VU_20130126_A0218_9F_R_FR08.mz5",
      "/dors/bioinfo/zhanglab/tcga_colon_proteomics/TCGA-AA-3552-01A-22_W_VU_20130126_A0218_9F_R_FR09.mz5",
      "/dors/bioinfo/zhanglab/tcga_colon_proteomics/TCGA-AA-3552-01A-22_W_VU_20130126_A0218_9F_R_FR10.mz5",
      "/dors/bioinfo/zhanglab/tcga_colon_proteomics/TCGA-AA-3552-01A-22_W_VU_20130126_A0218_9F_R_FR11.mz5",
      "/dors/bioinfo/zhanglab/tcga_colon_proteomics/TCGA-AA-3552-01A-22_W_VU_20130126_A0218_9F_R_FR12.mz5",
      "/dors/bioinfo/zhanglab/tcga_colon_proteomics/TCGA-AA-3552-01A-22_W_VU_20130126_A0218_9F_R_FR13.mz5",
      "/dors/bioinfo/zhanglab/tcga_colon_proteomics/TCGA-AA-3552-01A-22_W_VU_20130126_A0218_9F_R_FR14.mz5",
      "/dors/bioinfo/zhanglab/tcga_colon_proteomics/TCGA-AA-3552-01A-22_W_VU_20130126_A0218_9F_R_FR15.mz5"
    ],
  },
  msgf_target_buildsummary => {
    class           => "Proteomics::Summary::BuildSummary",
    perform         => 1,
    target_dir      => "${target_dir}/msgf_target_buildsummary",
    option          => "",
    source_ref      => "files",
    parameter_file  => $buildsummary_msgf_target_file,
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

#performConfig($config);
performTask( $config, "msgf_target_buildsummary" );

1;
