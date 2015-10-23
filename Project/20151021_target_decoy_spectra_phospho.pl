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
my $msgf_mod_file         = "/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/config/msgf_mods_itraq_phospho.txt";
my $msgf_option           = "-t 20ppm -ti 0,1 -tda 0 -m 3 -inst 1 -e 1 -protocol 3 -ntt 2 -n 2 -minLength 7 -maxLength 50 -addFeatures 1";
my $target_database       = "/scratch/cqs/shengq1/proteomics/20151001_target_decoy_spectra/database/humanRefSeq_Version54_with_tryp.fasta";
my $target_decoy_database = "/scratch/cqs/shengq1/proteomics/20151001_target_decoy_spectra/database/humanRefSeq_Version54_with_tryp_DECOY.fasta";

my $buildsummary_msgf_target_file                = "/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/config/buildsummary_msgf_target.param";
my $buildsummary_msgf_target_decoy_file          = "/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/config/buildsummary_msgf_target_decoy.param";
my $buildsummary_msgf_target_unique2_file        = "/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/config/buildsummary_msgf_target_unique2.param";
my $buildsummary_msgf_target_decoy_unique2_file  = "/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/config/buildsummary_msgf_target_decoy_unique2.param";
my $buildsummary_msgf_target_ratio2_file         = "/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/config/buildsummary_msgf_target_ratio2.param";
my $buildsummary_msgf_target_unique2_ratio2_file = "/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/config/buildsummary_msgf_target_unique2_ratio2.param";
my $buildsummary_bins                            = [ 1, 2, 3, 4, 5 ];

my $config = {
  general => { task_name => "phos" },
  files   => {
    "TCGA_13-1489_42-2590_36-2529_117C_P_PNNL_B2S5" => [
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_13-1489_42-2590_36-2529_117C_P_PNNL_B2S5_f01.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_13-1489_42-2590_36-2529_117C_P_PNNL_B2S5_f02.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_13-1489_42-2590_36-2529_117C_P_PNNL_B2S5_f03.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_13-1489_42-2590_36-2529_117C_P_PNNL_B2S5_f04.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_13-1489_42-2590_36-2529_117C_P_PNNL_B2S5_f05.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_13-1489_42-2590_36-2529_117C_P_PNNL_B2S5_f06.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_13-1489_42-2590_36-2529_117C_P_PNNL_B2S5_f07.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_13-1489_42-2590_36-2529_117C_P_PNNL_B2S5_f08.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_13-1489_42-2590_36-2529_117C_P_PNNL_B2S5_f09.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_13-1489_42-2590_36-2529_117C_P_PNNL_B2S5_f10.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_13-1489_42-2590_36-2529_117C_P_PNNL_B2S5_f11.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_13-1489_42-2590_36-2529_117C_P_PNNL_B2S5_f12.mgf"
    ],
    "TCGA_24-1428_24-1923_24-1563_117C_P_PNNL_B4S2" => [
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_24-1428_24-1923_24-1563_117C_P_PNNL_B4S2_f01.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_24-1428_24-1923_24-1563_117C_P_PNNL_B4S2_f02.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_24-1428_24-1923_24-1563_117C_P_PNNL_B4S2_f03.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_24-1428_24-1923_24-1563_117C_P_PNNL_B4S2_f04.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_24-1428_24-1923_24-1563_117C_P_PNNL_B4S2_f05.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_24-1428_24-1923_24-1563_117C_P_PNNL_B4S2_f06.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_24-1428_24-1923_24-1563_117C_P_PNNL_B4S2_f07.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_24-1428_24-1923_24-1563_117C_P_PNNL_B4S2_f08.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_24-1428_24-1923_24-1563_117C_P_PNNL_B4S2_f09.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_24-1428_24-1923_24-1563_117C_P_PNNL_B4S2_f10.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_24-1428_24-1923_24-1563_117C_P_PNNL_B4S2_f11.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_24-1428_24-1923_24-1563_117C_P_PNNL_B4S2_f12.mgf"
    ],
    "TCGA_24-1430_24-1466_24-1552_117C_P_PNNL_B4S3" => [
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_24-1430_24-1466_24-1552_117C_P_PNNL_B4S3_f01.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_24-1430_24-1466_24-1552_117C_P_PNNL_B4S3_f02.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_24-1430_24-1466_24-1552_117C_P_PNNL_B4S3_f03.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_24-1430_24-1466_24-1552_117C_P_PNNL_B4S3_f04.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_24-1430_24-1466_24-1552_117C_P_PNNL_B4S3_f05.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_24-1430_24-1466_24-1552_117C_P_PNNL_B4S3_f06.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_24-1430_24-1466_24-1552_117C_P_PNNL_B4S3_f07.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_24-1430_24-1466_24-1552_117C_P_PNNL_B4S3_f08.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_24-1430_24-1466_24-1552_117C_P_PNNL_B4S3_f09.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_24-1430_24-1466_24-1552_117C_P_PNNL_B4S3_f10.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_24-1430_24-1466_24-1552_117C_P_PNNL_B4S3_f11.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_24-1430_24-1466_24-1552_117C_P_PNNL_B4S3_f12.mgf"
    ],
    "TCGA_24-1467_29-2432_25-1321_117C_P_PNNL_B1S2" => [
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_24-1467_29-2432_25-1321_117C_P_PNNL_B1S2_f01.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_24-1467_29-2432_25-1321_117C_P_PNNL_B1S2_f02.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_24-1467_29-2432_25-1321_117C_P_PNNL_B1S2_f03.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_24-1467_29-2432_25-1321_117C_P_PNNL_B1S2_f04.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_24-1467_29-2432_25-1321_117C_P_PNNL_B1S2_f05.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_24-1467_29-2432_25-1321_117C_P_PNNL_B1S2_f06.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_24-1467_29-2432_25-1321_117C_P_PNNL_B1S2_f07.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_24-1467_29-2432_25-1321_117C_P_PNNL_B1S2_f08.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_24-1467_29-2432_25-1321_117C_P_PNNL_B1S2_f09.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_24-1467_29-2432_25-1321_117C_P_PNNL_B1S2_f10.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_24-1467_29-2432_25-1321_117C_P_PNNL_B1S2_f11.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_24-1467_29-2432_25-1321_117C_P_PNNL_B1S2_f12.mgf"
    ],
    "TCGA_24-1545_24-2288_29-1693_117C_P_PNNL_B1S5" => [
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_24-1545_24-2288_29-1693_117C_P_PNNL_B1S5_f01.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_24-1545_24-2288_29-1693_117C_P_PNNL_B1S5_f02.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_24-1545_24-2288_29-1693_117C_P_PNNL_B1S5_f03.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_24-1545_24-2288_29-1693_117C_P_PNNL_B1S5_f04.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_24-1545_24-2288_29-1693_117C_P_PNNL_B1S5_f05.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_24-1545_24-2288_29-1693_117C_P_PNNL_B1S5_f06.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_24-1545_24-2288_29-1693_117C_P_PNNL_B1S5_f07.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_24-1545_24-2288_29-1693_117C_P_PNNL_B1S5_f08.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_24-1545_24-2288_29-1693_117C_P_PNNL_B1S5_f09.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_24-1545_24-2288_29-1693_117C_P_PNNL_B1S5_f10.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_24-1545_24-2288_29-1693_117C_P_PNNL_B1S5_f11.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_24-1545_24-2288_29-1693_117C_P_PNNL_B1S5_f12.mgf"
    ],
    "TCGA_24-1556_24-2267_24-2260_117C_P_PNNL_B3S5" => [
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_24-1556_24-2267_24-2260_117C_P_PNNL_B3S5_f01.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_24-1556_24-2267_24-2260_117C_P_PNNL_B3S5_f02.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_24-1556_24-2267_24-2260_117C_P_PNNL_B3S5_f03.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_24-1556_24-2267_24-2260_117C_P_PNNL_B3S5_f04.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_24-1556_24-2267_24-2260_117C_P_PNNL_B3S5_f05.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_24-1556_24-2267_24-2260_117C_P_PNNL_B3S5_f06.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_24-1556_24-2267_24-2260_117C_P_PNNL_B3S5_f07.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_24-1556_24-2267_24-2260_117C_P_PNNL_B3S5_f08.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_24-1556_24-2267_24-2260_117C_P_PNNL_B3S5_f09.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_24-1556_24-2267_24-2260_117C_P_PNNL_B3S5_f10.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_24-1556_24-2267_24-2260_117C_P_PNNL_B3S5_f11.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_24-1556_24-2267_24-2260_117C_P_PNNL_B3S5_f12.mgf"
    ],
    "TCGA_24-1604_24-2020_25-1320_117C_P_PNNL_B2S3" => [
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_24-1604_24-2020_25-1320_117C_P_PNNL_B2S3_f01.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_24-1604_24-2020_25-1320_117C_P_PNNL_B2S3_f02.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_24-1604_24-2020_25-1320_117C_P_PNNL_B2S3_f03.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_24-1604_24-2020_25-1320_117C_P_PNNL_B2S3_f04.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_24-1604_24-2020_25-1320_117C_P_PNNL_B2S3_f05.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_24-1604_24-2020_25-1320_117C_P_PNNL_B2S3_f06.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_24-1604_24-2020_25-1320_117C_P_PNNL_B2S3_f07.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_24-1604_24-2020_25-1320_117C_P_PNNL_B2S3_f08.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_24-1604_24-2020_25-1320_117C_P_PNNL_B2S3_f09.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_24-1604_24-2020_25-1320_117C_P_PNNL_B2S3_f10.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_24-1604_24-2020_25-1320_117C_P_PNNL_B2S3_f11.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_24-1604_24-2020_25-1320_117C_P_PNNL_B2S3_f12.mgf"
    ],
    "TCGA_24-2023_24-1603_61-2612_117C_P_PNNL_B1S3" => [
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_24-2023_24-1603_61-2612_117C_P_PNNL_B1S3_f01.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_24-2023_24-1603_61-2612_117C_P_PNNL_B1S3_f02.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_24-2023_24-1603_61-2612_117C_P_PNNL_B1S3_f03.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_24-2023_24-1603_61-2612_117C_P_PNNL_B1S3_f04.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_24-2023_24-1603_61-2612_117C_P_PNNL_B1S3_f05.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_24-2023_24-1603_61-2612_117C_P_PNNL_B1S3_f06.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_24-2023_24-1603_61-2612_117C_P_PNNL_B1S3_f07.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_24-2023_24-1603_61-2612_117C_P_PNNL_B1S3_f08.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_24-2023_24-1603_61-2612_117C_P_PNNL_B1S3_f09.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_24-2023_24-1603_61-2612_117C_P_PNNL_B1S3_f10.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_24-2023_24-1603_61-2612_117C_P_PNNL_B1S3_f11.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_24-2023_24-1603_61-2612_117C_P_PNNL_B1S3_f12.mgf"
    ],
    "TCGA_25-1316_24-1550_24-1555_117C_P_PNNL_B1S1" => [
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_25-1316_24-1550_24-1555_117C_P_PNNL_B1S1_f01.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_25-1316_24-1550_24-1555_117C_P_PNNL_B1S1_f02.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_25-1316_24-1550_24-1555_117C_P_PNNL_B1S1_f03.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_25-1316_24-1550_24-1555_117C_P_PNNL_B1S1_f04.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_25-1316_24-1550_24-1555_117C_P_PNNL_B1S1_f05.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_25-1316_24-1550_24-1555_117C_P_PNNL_B1S1_f06.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_25-1316_24-1550_24-1555_117C_P_PNNL_B1S1_f07.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_25-1316_24-1550_24-1555_117C_P_PNNL_B1S1_f08.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_25-1316_24-1550_24-1555_117C_P_PNNL_B1S1_f09.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_25-1316_24-1550_24-1555_117C_P_PNNL_B1S1_f10.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_25-1316_24-1550_24-1555_117C_P_PNNL_B1S1_f11.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_25-1316_24-1550_24-1555_117C_P_PNNL_B1S1_f12.mgf"
    ],
    "TCGA_25-1319_24-2289_42-2588_117C_P_PNNL_B2S6" => [
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_25-1319_24-2289_42-2588_117C_P_PNNL_B2S6_f01.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_25-1319_24-2289_42-2588_117C_P_PNNL_B2S6_f02.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_25-1319_24-2289_42-2588_117C_P_PNNL_B2S6_f03.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_25-1319_24-2289_42-2588_117C_P_PNNL_B2S6_f04.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_25-1319_24-2289_42-2588_117C_P_PNNL_B2S6_f05.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_25-1319_24-2289_42-2588_117C_P_PNNL_B2S6_f06.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_25-1319_24-2289_42-2588_117C_P_PNNL_B2S6_f07.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_25-1319_24-2289_42-2588_117C_P_PNNL_B2S6_f08.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_25-1319_24-2289_42-2588_117C_P_PNNL_B2S6_f09.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_25-1319_24-2289_42-2588_117C_P_PNNL_B2S6_f10.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_25-1319_24-2289_42-2588_117C_P_PNNL_B2S6_f11.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_25-1319_24-2289_42-2588_117C_P_PNNL_B2S6_f12.mgf"
    ],
    "TCGA_25-1628_13-1494_24-1104_117C_P_PNNL_B1S4" => [
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_25-1628_13-1494_24-1104_117C_P_PNNL_B1S4_f01.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_25-1628_13-1494_24-1104_117C_P_PNNL_B1S4_f02.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_25-1628_13-1494_24-1104_117C_P_PNNL_B1S4_f03.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_25-1628_13-1494_24-1104_117C_P_PNNL_B1S4_f04.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_25-1628_13-1494_24-1104_117C_P_PNNL_B1S4_f05.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_25-1628_13-1494_24-1104_117C_P_PNNL_B1S4_f06.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_25-1628_13-1494_24-1104_117C_P_PNNL_B1S4_f07.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_25-1628_13-1494_24-1104_117C_P_PNNL_B1S4_f08.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_25-1628_13-1494_24-1104_117C_P_PNNL_B1S4_f09.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_25-1628_13-1494_24-1104_117C_P_PNNL_B1S4_f10.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_25-1628_13-1494_24-1104_117C_P_PNNL_B1S4_f11.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_25-1628_13-1494_24-1104_117C_P_PNNL_B1S4_f12.mgf"
    ],
    "TCGA_25-1630_13-1499_29-1775_117C_P_PNNL_B4S5" => [
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_25-1630_13-1499_29-1775_117C_P_PNNL_B4S5_f01.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_25-1630_13-1499_29-1775_117C_P_PNNL_B4S5_f02.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_25-1630_13-1499_29-1775_117C_P_PNNL_B4S5_f03.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_25-1630_13-1499_29-1775_117C_P_PNNL_B4S5_f04.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_25-1630_13-1499_29-1775_117C_P_PNNL_B4S5_f05.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_25-1630_13-1499_29-1775_117C_P_PNNL_B4S5_f06.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_25-1630_13-1499_29-1775_117C_P_PNNL_B4S5_f07.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_25-1630_13-1499_29-1775_117C_P_PNNL_B4S5_f08.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_25-1630_13-1499_29-1775_117C_P_PNNL_B4S5_f09.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_25-1630_13-1499_29-1775_117C_P_PNNL_B4S5_f10.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_25-1630_13-1499_29-1775_117C_P_PNNL_B4S5_f11.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_25-1630_13-1499_29-1775_117C_P_PNNL_B4S5_f12.mgf"
    ],
    "TCGA_25-2400_30-1866_25-2409_117C_P_PNNL_B3S6" => [
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_25-2400_30-1866_25-2409_117C_P_PNNL_B3S6_f01.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_25-2400_30-1866_25-2409_117C_P_PNNL_B3S6_f02.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_25-2400_30-1866_25-2409_117C_P_PNNL_B3S6_f03.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_25-2400_30-1866_25-2409_117C_P_PNNL_B3S6_f04.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_25-2400_30-1866_25-2409_117C_P_PNNL_B3S6_f05.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_25-2400_30-1866_25-2409_117C_P_PNNL_B3S6_f06.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_25-2400_30-1866_25-2409_117C_P_PNNL_B3S6_f07.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_25-2400_30-1866_25-2409_117C_P_PNNL_B3S6_f08.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_25-2400_30-1866_25-2409_117C_P_PNNL_B3S6_f09.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_25-2400_30-1866_25-2409_117C_P_PNNL_B3S6_f10.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_25-2400_30-1866_25-2409_117C_P_PNNL_B3S6_f11.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_25-2400_30-1866_25-2409_117C_P_PNNL_B3S6_f12.mgf"
    ],
    "TCGA_29-1688_24-1435_24-1562_117C_P_PNNL_B3S1" => [
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_29-1688_24-1435_24-1562_117C_P_PNNL_B3S1_f01.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_29-1688_24-1435_24-1562_117C_P_PNNL_B3S1_f02.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_29-1688_24-1435_24-1562_117C_P_PNNL_B3S1_f03.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_29-1688_24-1435_24-1562_117C_P_PNNL_B3S1_f04.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_29-1688_24-1435_24-1562_117C_P_PNNL_B3S1_f05.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_29-1688_24-1435_24-1562_117C_P_PNNL_B3S1_f06.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_29-1688_24-1435_24-1562_117C_P_PNNL_B3S1_f07.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_29-1688_24-1435_24-1562_117C_P_PNNL_B3S1_f08.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_29-1688_24-1435_24-1562_117C_P_PNNL_B3S1_f09.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_29-1688_24-1435_24-1562_117C_P_PNNL_B3S1_f10.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_29-1688_24-1435_24-1562_117C_P_PNNL_B3S1_f11.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_29-1688_24-1435_24-1562_117C_P_PNNL_B3S1_f12.mgf"
    ],
    "TCGA_29-1762_13-2071_36-2543_117C_P_PNNL_B1S6" => [
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_29-1762_13-2071_36-2543_117C_P_PNNL_B1S6_f01.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_29-1762_13-2071_36-2543_117C_P_PNNL_B1S6_f02.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_29-1762_13-2071_36-2543_117C_P_PNNL_B1S6_f03.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_29-1762_13-2071_36-2543_117C_P_PNNL_B1S6_f04.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_29-1762_13-2071_36-2543_117C_P_PNNL_B1S6_f05.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_29-1762_13-2071_36-2543_117C_P_PNNL_B1S6_f06.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_29-1762_13-2071_36-2543_117C_P_PNNL_B1S6_f07.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_29-1762_13-2071_36-2543_117C_P_PNNL_B1S6_f08.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_29-1762_13-2071_36-2543_117C_P_PNNL_B1S6_f09.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_29-1762_13-2071_36-2543_117C_P_PNNL_B1S6_f10.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_29-1762_13-2071_36-2543_117C_P_PNNL_B1S6_f11.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_29-1762_13-2071_36-2543_117C_P_PNNL_B1S6_f12.mgf"
    ],
    "TCGA_29-1763_61-1741_24-2024_117C_P_PNNL_B4S1" => [
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_29-1763_61-1741_24-2024_117C_P_PNNL_B4S1_f01.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_29-1763_61-1741_24-2024_117C_P_PNNL_B4S1_f02.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_29-1763_61-1741_24-2024_117C_P_PNNL_B4S1_f03.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_29-1763_61-1741_24-2024_117C_P_PNNL_B4S1_f04.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_29-1763_61-1741_24-2024_117C_P_PNNL_B4S1_f05.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_29-1763_61-1741_24-2024_117C_P_PNNL_B4S1_f06.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_29-1763_61-1741_24-2024_117C_P_PNNL_B4S1_f07.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_29-1763_61-1741_24-2024_117C_P_PNNL_B4S1_f08.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_29-1763_61-1741_24-2024_117C_P_PNNL_B4S1_f09.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_29-1763_61-1741_24-2024_117C_P_PNNL_B4S1_f10.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_29-1763_61-1741_24-2024_117C_P_PNNL_B4S1_f11.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_29-1763_61-1741_24-2024_117C_P_PNNL_B4S1_f12.mgf"
    ],
    "TCGA_29-1774_36-1576_24-2030_117C_P_PNNL_B2S4" => [
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_29-1774_36-1576_24-2030_117C_P_PNNL_B2S4_f01.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_29-1774_36-1576_24-2030_117C_P_PNNL_B2S4_f02.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_29-1774_36-1576_24-2030_117C_P_PNNL_B2S4_f03.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_29-1774_36-1576_24-2030_117C_P_PNNL_B2S4_f04.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_29-1774_36-1576_24-2030_117C_P_PNNL_B2S4_f05.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_29-1774_36-1576_24-2030_117C_P_PNNL_B2S4_f06.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_29-1774_36-1576_24-2030_117C_P_PNNL_B2S4_f07.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_29-1774_36-1576_24-2030_117C_P_PNNL_B2S4_f08.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_29-1774_36-1576_24-2030_117C_P_PNNL_B2S4_f09.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_29-1774_36-1576_24-2030_117C_P_PNNL_B2S4_f10.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_29-1774_36-1576_24-2030_117C_P_PNNL_B2S4_f11.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_29-1774_36-1576_24-2030_117C_P_PNNL_B2S4_f12.mgf"
    ],
    "TCGA_29-1777_61-2008_24-1551_117C_P_PNNL_B2S1" => [
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_29-1777_61-2008_24-1551_117C_P_PNNL_B2S1_f01.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_29-1777_61-2008_24-1551_117C_P_PNNL_B2S1_f02.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_29-1777_61-2008_24-1551_117C_P_PNNL_B2S1_f03.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_29-1777_61-2008_24-1551_117C_P_PNNL_B2S1_f04.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_29-1777_61-2008_24-1551_117C_P_PNNL_B2S1_f05.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_29-1777_61-2008_24-1551_117C_P_PNNL_B2S1_f06.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_29-1777_61-2008_24-1551_117C_P_PNNL_B2S1_f07.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_29-1777_61-2008_24-1551_117C_P_PNNL_B2S1_f08.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_29-1777_61-2008_24-1551_117C_P_PNNL_B2S1_f09.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_29-1777_61-2008_24-1551_117C_P_PNNL_B2S1_f10.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_29-1777_61-2008_24-1551_117C_P_PNNL_B2S1_f11.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_29-1777_61-2008_24-1551_117C_P_PNNL_B2S1_f12.mgf"
    ],
    "TCGA_29-2414_61-1995_13-1484_117C_P_PNNL_B3S2" => [
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_29-2414_61-1995_13-1484_117C_P_PNNL_B3S2_f01.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_29-2414_61-1995_13-1484_117C_P_PNNL_B3S2_f02.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_29-2414_61-1995_13-1484_117C_P_PNNL_B3S2_f03.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_29-2414_61-1995_13-1484_117C_P_PNNL_B3S2_f04.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_29-2414_61-1995_13-1484_117C_P_PNNL_B3S2_f05.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_29-2414_61-1995_13-1484_117C_P_PNNL_B3S2_f06.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_29-2414_61-1995_13-1484_117C_P_PNNL_B3S2_f07.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_29-2414_61-1995_13-1484_117C_P_PNNL_B3S2_f08.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_29-2414_61-1995_13-1484_117C_P_PNNL_B3S2_f09.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_29-2414_61-1995_13-1484_117C_P_PNNL_B3S2_f10.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_29-2414_61-1995_13-1484_117C_P_PNNL_B3S2_f11.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_29-2414_61-1995_13-1484_117C_P_PNNL_B3S2_f12.mgf"
    ],
    "TCGA_36-1581_23-1123_24-1553_117C_P_PNNL_B2S2" => [
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_36-1581_23-1123_24-1553_117C_P_PNNL_B2S2_f01.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_36-1581_23-1123_24-1553_117C_P_PNNL_B2S2_f02.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_36-1581_23-1123_24-1553_117C_P_PNNL_B2S2_f03.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_36-1581_23-1123_24-1553_117C_P_PNNL_B2S2_f04.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_36-1581_23-1123_24-1553_117C_P_PNNL_B2S2_f05.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_36-1581_23-1123_24-1553_117C_P_PNNL_B2S2_f06.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_36-1581_23-1123_24-1553_117C_P_PNNL_B2S2_f07.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_36-1581_23-1123_24-1553_117C_P_PNNL_B2S2_f08.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_36-1581_23-1123_24-1553_117C_P_PNNL_B2S2_f09.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_36-1581_23-1123_24-1553_117C_P_PNNL_B2S2_f10.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_36-1581_23-1123_24-1553_117C_P_PNNL_B2S2_f11.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_36-1581_23-1123_24-1553_117C_P_PNNL_B2S2_f12.mgf"
    ],
    "TCGA_36-2530_24-1103_29-1776_117C_P_PNNL_B3S3" => [
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_36-2530_24-1103_29-1776_117C_P_PNNL_B3S3_f01.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_36-2530_24-1103_29-1776_117C_P_PNNL_B3S3_f02.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_36-2530_24-1103_29-1776_117C_P_PNNL_B3S3_f03.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_36-2530_24-1103_29-1776_117C_P_PNNL_B3S3_f04.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_36-2530_24-1103_29-1776_117C_P_PNNL_B3S3_f05.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_36-2530_24-1103_29-1776_117C_P_PNNL_B3S3_f06.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_36-2530_24-1103_29-1776_117C_P_PNNL_B3S3_f07.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_36-2530_24-1103_29-1776_117C_P_PNNL_B3S3_f08.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_36-2530_24-1103_29-1776_117C_P_PNNL_B3S3_f09.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_36-2530_24-1103_29-1776_117C_P_PNNL_B3S3_f10.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_36-2530_24-1103_29-1776_117C_P_PNNL_B3S3_f11.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_36-2530_24-1103_29-1776_117C_P_PNNL_B3S3_f12.mgf"
    ],
    "TCGA_61-1907_23-1124_25-1635_117C_P_PNNL_B3S4" => [
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_61-1907_23-1124_25-1635_117C_P_PNNL_B3S4_f01.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_61-1907_23-1124_25-1635_117C_P_PNNL_B3S4_f02.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_61-1907_23-1124_25-1635_117C_P_PNNL_B3S4_f03.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_61-1907_23-1124_25-1635_117C_P_PNNL_B3S4_f04.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_61-1907_23-1124_25-1635_117C_P_PNNL_B3S4_f05.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_61-1907_23-1124_25-1635_117C_P_PNNL_B3S4_f06.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_61-1907_23-1124_25-1635_117C_P_PNNL_B3S4_f07.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_61-1907_23-1124_25-1635_117C_P_PNNL_B3S4_f08.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_61-1907_23-1124_25-1635_117C_P_PNNL_B3S4_f09.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_61-1907_23-1124_25-1635_117C_P_PNNL_B3S4_f10.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_61-1907_23-1124_25-1635_117C_P_PNNL_B3S4_f11.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_61-1907_23-1124_25-1635_117C_P_PNNL_B3S4_f12.mgf"
    ],
    "TCGA_61-2095_61-1919_25-2404_117C_P_PNNL_B4S4" => [
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_61-2095_61-1919_25-2404_117C_P_PNNL_B4S4_f01.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_61-2095_61-1919_25-2404_117C_P_PNNL_B4S4_f02.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_61-2095_61-1919_25-2404_117C_P_PNNL_B4S4_f03.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_61-2095_61-1919_25-2404_117C_P_PNNL_B4S4_f04.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_61-2095_61-1919_25-2404_117C_P_PNNL_B4S4_f05.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_61-2095_61-1919_25-2404_117C_P_PNNL_B4S4_f06.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_61-2095_61-1919_25-2404_117C_P_PNNL_B4S4_f07.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_61-2095_61-1919_25-2404_117C_P_PNNL_B4S4_f08.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_61-2095_61-1919_25-2404_117C_P_PNNL_B4S4_f09.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_61-2095_61-1919_25-2404_117C_P_PNNL_B4S4_f10.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_61-2095_61-1919_25-2404_117C_P_PNNL_B4S4_f11.mgf",
      "/gpfs21/scratch/cqs/shengq1/proteomics/20151021_target_decoy_spectra_phospho/phosphorylation_data/TCGA_61-2095_61-1919_25-2404_117C_P_PNNL_B4S4_f12.mgf"
    ],
  },
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
    source     => { T1_individual => [ "shift_precursor", "msgf_target", "msgf_target_psm", "msgf_target_decoy", "msgf_target_decoy_psm" ], },
    sh_direct  => 0,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  }
};

performConfig($config);

1;
