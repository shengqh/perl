#!/usr/bin/perl
use strict;
use warnings;

use CQS::QC;
use CQS::DNASeq;
use CQS::RNASeq;
use CQS::CQSTools;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;
use CQS::Cutadapt;

my $root;
my $cqs_tools;
my $hsa_gffs;
my $rno_gffs;
my $mmu_gffs;

if ( is_linux() ) {
  $root      = "/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2";
  $cqs_tools = "/home/shengq1/cqstools/CQS.Tools.exe";
  $hsa_gffs  = "/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/smrnapipeline/hsa_tableL.bed";
  $rno_gffs  = "/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/smrnapipeline/rno_tableL.bed";
  $mmu_gffs  = "/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/smrnapipeline/mmu_tableL.bed";
}
else {
  $root      = "d:/temp";
  $cqs_tools = "E:/sqh/programs/csharp/OmicsLabCSharp/CQS.Tools/bin/Release/CQS.Tools.exe";
  $hsa_gffs  = "H:/shengquanhu/projects/vangard/VANGARD00055_guoyan_mirna/hsa.gff3";
  $rno_gffs  = "H:/shengquanhu/projects/vangard/VANGARD00055_guoyan_mirna/rno.gff3";
  $mmu_gffs  = "H:/shengquanhu/projects/vangard/VANGARD00055_guoyan_mirna/mmu.gff3";
}
my $target_dir = create_directory_or_die($root);

my $target_rat_dir   = create_directory_or_die( $target_dir . "/rat" );
my $target_human_dir = create_directory_or_die( $target_dir . "/human" );
my $target_mouse_dir = create_directory_or_die( $target_dir . "/mouse" );

my $email     = "quanhu.sheng\@vanderbilt.edu";
my $task_name = "VANGARD00055";

my $bowtie1_option = "-a -m 20 --best --strata -v 3 -l 12 -p 8";

my $bowtie1_rat_index   = "/data/cqs/shengq1/reference/rn4/bowtie1_index/rn4";
my $bowtie1_human_index = "/data/cqs/guoy1/reference/hg19/bowtie_index/hg19";
my $bowtie1_mouse_index = "/data/cqs/guoy1/reference/mm9/bowtie_index/mm9";

#shrimp2 gmapper set mirna mode
#static int
#set_mode_from_string(char const * s) {
#  if (!strcmp(s, "mirna")) {
#    mode_mirna = true;
#
#    //load_default_mirna_seeds();
#
#    Hflag = true;
#    gapless_sw = true;
#    anchor_width = 0;
#    a_gap_open_score = -255;
#    b_gap_open_score = -255;
#    hash_filter_calls = false;
#    match_mode = 1;
#    window_len = 100.0;
#    Gflag = false;
#    compute_mapping_qualities = false;
#
#    return 1;
#  } else {
#    return 0;
#  }
#}
my $shrimp2_option              = "-Q -N 8 -o 1 --qv-offset 33";
my $shrimp2_rat_miRBase_index   = "/data/cqs/shengq1/reference/miRBase19/shrimp_index_2.2.3/rno.mature.dna-ls";
my $shrimp2_human_miRBase_index = "/data/cqs/shengq1/reference/miRBase19/shrimp_index_2.2.3/hsa.mature.dna-ls";
my $shrimp2_mouse_miRBase_index = "/data/cqs/shengq1/reference/miRBase19/shrimp_index_2.2.3/mmu.mature.dna-ls";

my $rat = {
  source => {
    "2516-01" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/2516-KCV-1_1.fastq"],
    "2516-02" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/2516-KCV-2_1.fastq"],
    "2516-03" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/2516-KCV-3_1.fastq"],
    "2516-04" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/2516-KCV-4_1.fastq"],
    "2516-05" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/2516-KCV-5_1.fastq"],
    "2516-06" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/2516-KCV-6_1.fastq"],
    "2516-07" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/2516-KCV-7_1.fastq"],
    "2516-08" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/2516-KCV-8_1.fastq"],
    "2516-09" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/2516-KCV-9_1.fastq"],
    "2570-01" => ["/autofs/blue_sequencer/Runs/projects/2570-KCV/2013-06-19/2570-KCV-1_1.fastq.gz"],
    "2570-02" => ["/autofs/blue_sequencer/Runs/projects/2570-KCV/2013-06-19/2570-KCV-2_1.fastq.gz"],
    "2570-03" => ["/autofs/blue_sequencer/Runs/projects/2570-KCV/2013-06-19/2570-KCV-3_1.fastq.gz"],
    "2570-04" => ["/autofs/blue_sequencer/Runs/projects/2570-KCV/2013-06-19/2570-KCV-4_1.fastq.gz"],
    "2570-05" => ["/autofs/blue_sequencer/Runs/projects/2570-KCV/2013-06-19/2570-KCV-5_1.fastq.gz"],
    "2570-06" => ["/autofs/blue_sequencer/Runs/projects/2570-KCV/2013-06-19/2570-KCV-6_1.fastq.gz"],
    "2570-07" => ["/autofs/blue_sequencer/Runs/projects/2570-KCV/2013-06-19/2570-KCV-7_1.fastq.gz"],
    "2570-08" => ["/autofs/blue_sequencer/Runs/projects/2570-KCV/2013-06-19/2570-KCV-8_1.fastq.gz"],
    "2570-09" => ["/autofs/blue_sequencer/Runs/projects/2570-KCV/2013-06-19/2570-KCV-9_1.fastq.gz"],
    "2570-10" => ["/autofs/blue_sequencer/Runs/projects/2570-KCV/2013-06-19/2570-KCV-10_1.fastq.gz"],
    "2570-11" => ["/autofs/blue_sequencer/Runs/projects/2570-KCV/2013-06-19/2570-KCV-11_1.fastq.gz"],
    "2570-12" => ["/autofs/blue_sequencer/Runs/projects/2570-KCV/2013-06-19/2570-KCV-12_1.fastq.gz"],
    "2570-13" => ["/autofs/blue_sequencer/Runs/projects/2570-KCV/2013-06-19/2570-KCV-13_1.fastq.gz"],
    "2570-14" => ["/autofs/blue_sequencer/Runs/projects/2570-KCV/2013-06-19/2570-KCV-14_1.fastq.gz"],
    "2570-15" => ["/autofs/blue_sequencer/Runs/projects/2570-KCV/2013-06-19/2570-KCV-15_1.fastq.gz"],
    "2570-16" => ["/autofs/blue_sequencer/Runs/projects/2570-KCV/2013-06-19/2570-KCV-16_1.fastq.gz"],
    "2570-17" => ["/autofs/blue_sequencer/Runs/projects/2570-KCV/2013-06-19/2570-KCV-17_1.fastq.gz"],
    "2570-18" => ["/autofs/blue_sequencer/Runs/projects/2570-KCV/2013-06-19/2570-KCV-18_1.fastq.gz"],
  },
  coordinate    => $rno_gffs,
  bowtie1_index => $bowtie1_rat_index,
  shrimp2_index => $shrimp2_rat_miRBase_index,
  target_dir    => $target_rat_dir,
  task_name     => $task_name . "_rat",
};
my $human = {
  source => {
    "2516-10"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/2516-KCV-10_1.fastq"],
    "2516-11"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/2516-KCV-11_1.fastq"],
    "2516-12"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/2516-KCV-12_1.fastq"],
    "2516-13"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/2516-KCV-13_1.fastq"],
    "2516-14"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/2516-KCV-14_1.fastq"],
    "2516-15"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/2516-KCV-15_1.fastq"],
    "2516-16"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/2516-KCV-16_1.fastq"],
    "2516-17"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/2516-KCV-17_1.fastq"],
    "2516-18"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/2516-KCV-18_1.fastq"],
    "2516-19"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/2516-KCV-19_1.fastq"],
    "2516-20"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/2516-KCV-20_1.fastq"],
    "2516-21"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/2516-KCV-21_1.fastq"],
    "2516-22"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/2516-KCV-22_1.fastq"],
    "2516-23"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/2516-KCV-23_1.fastq"],
    "2516-24"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/2516-KCV-24_1.fastq"],
    "2516-25"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/2516-KCV-25_1.fastq"],
    "2516-26"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/2516-KCV-26_1.fastq"],
    "2516-27"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/2516-KCV-27_1.fastq"],
    "2516-28"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/2516-KCV-28_1.fastq"],
    "2516-29"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/2516-KCV-29_1.fastq"],
    "2516-30"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/2516-KCV-30_1.fastq"],
    "2516-31"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/2516-KCV-31_1.fastq"],
    "2516-32"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/2516-KCV-32_1.fastq"],
    "2516-33"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/2516-KCV-33_1.fastq"],
    "2516-34"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/2516-KCV-34_1.fastq"],
    "2516-35"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/2516-KCV-35_1.fastq"],
    "2516-36"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/2516-KCV-36_1.fastq"],
    "2516-37"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/2516-KCV-37_1.fastq"],
    "2516-38"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/2516-KCV-38_1.fastq"],
    "2516-39"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/2516-KCV-39_1.fastq"],
    "2516-40"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/2516-KCV-40_1.fastq"],
    "2516-41"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/2516-KCV-41_1.fastq"],
    "2516-42"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/2516-KCV-42_1.fastq"],
    "2516-43"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/2516-KCV-43_1.fastq"],
    "2516-44"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/2516-KCV-44_1.fastq"],
    "2516-45"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/2516-KCV-45_1.fastq"],
    "2516-46"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/2516-KCV-46_1.fastq"],
    "2516-47"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/2516-KCV-47_1.fastq"],
    "2516-48"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/2516-KCV-48_1.fastq"],
    "KCV2_1N2" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/KCV2_1N2_GCCAAT_L003_R1_001.fastq"],
    "KCV2_1N3" => [
      "/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/KCV2_1N3_CAGATC_L003_R1_001.fastq",
      "/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/KCV2_1N3_CAGATC_L003_R1_002.fastq"
    ],
    "KCV2_1N4" => [
      "/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/KCV2_1N4_ACTTGA_L003_R1_001.fastq",
      "/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/KCV2_1N4_ACTTGA_L003_R1_002.fastq"
    ],
    "KCV2_1N5" => [
      "/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/KCV2_1N5_GATCAG_L003_R1_001.fastq",
      "/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/KCV2_1N5_GATCAG_L003_R1_002.fastq"
    ],
    "KCV2_1N6" => [
      "/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/KCV2_1N6_TAGCTT_L003_R1_001.fastq",
      "/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/KCV2_1N6_TAGCTT_L003_R1_002.fastq"
    ],
    "KCV2_2N1" => [
      "/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/KCV2_2N1_GGCTAC_L003_R1_001.fastq",
      "/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/KCV2_2N1_GGCTAC_L003_R1_002.fastq",
      "/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/KCV2_2N1_GGCTAC_L003_R1_003.fastq",
      "/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/KCV2_2N1_GGCTAC_L003_R1_004.fastq",
      "/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/KCV2_2N1_GGCTAC_L003_R1_005.fastq"
    ],
    "KCV2_2N2" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/KCV2_2N2_GCCGCG_L003_R1_001.fastq"],
    "KCV2_2N3" => [
      "/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/KCV2_2N3_CTTGTA_L003_R1_001.fastq",
      "/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/KCV2_2N3_CTTGTA_L003_R1_002.fastq",
      "/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/KCV2_2N3_CTTGTA_L003_R1_003.fastq",
      "/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/KCV2_2N3_CTTGTA_L003_R1_004.fastq",
      "/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/KCV2_2N3_CTTGTA_L003_R1_005.fastq",
      "/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/KCV2_2N3_CTTGTA_L003_R1_006.fastq",
      "/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/KCV2_2N3_CTTGTA_L003_R1_007.fastq",
      "/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/KCV2_2N3_CTTGTA_L003_R1_008.fastq",
      "/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/KCV2_2N3_CTTGTA_L003_R1_009.fastq",
      "/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/KCV2_2N3_CTTGTA_L003_R1_010.fastq"
    ],
    "KCV2_2N4"       => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/KCV2_2N4_GCCTTA_L003_R1_001.fastq"],
    "KCV2_2N5"       => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/KCV2_2N5_GCTCCA_L003_R1_001.fastq"],
    "KCV3_1C2"       => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/KCV3_1C2_GGCACA_L004_R1_001.fastq"],
    "KCV3_1C3"       => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/KCV3_1C3_GGCCTG_L004_R1_001.fastq"],
    "KCV3_1C4"       => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/KCV3_1C4_TCTACC_L004_R1_001.fastq"],
    "KCV3_1C5"       => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/KCV3_1C5_TGAAGT_L004_R1_001.fastq"],
    "KCV3_1C6"       => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/KCV3_1C6_TGCCAT_L004_R1_001.fastq"],
    "KCV3_2C1"       => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/KCV3_2C1_TGCTGG_L004_R1_001.fastq"],
    "KCV3_2C2"       => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/KCV3_2C2_TGGCGC_L004_R1_001.fastq"],
    "KCV3_2C3"       => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/KCV3_2C3_TTCGAA_L004_R1_001.fastq"],
    "Sample1"        => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/Sample1_12.fastq"],
    "Sample2"        => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/Sample2_12.fastq"],
    "Sample3"        => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/Sample3_12.fastq"],
    "Sample4"        => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/Sample4_12.fastq"],
    "Sample5"        => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/Sample5_12.fastq"],
    "2570-KCV-01-19" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/mouseLiverTransplant1_GTGAAA_L003_R1_001.fastq.gz"],
    "2570-KCV-01-20" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/mouseLiverTransplant2_GTGGCC_L003_R1_001.fastq.gz"],
    "2570-KCV-01-21" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/mouseLiverTransplant3_GTTTCG_L003_R1_001.fastq.gz"],
    "2570-KCV-01-22" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/mouseLiverTransplant4_CGTACG_L003_R1_001.fastq.gz"],
    "2570-KCV-01-23" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/mouseLiverControl1_GAGTGG_L003_R1_001.fastq.gz"],
    "2570-KCV-01-24" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/mouseLiverControl2_GGTAGC_L003_R1_001.fastq.gz"],
    "2570-KCV-01-25" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/mouseLiverControl3_ACTGAT_L003_R1_001.fastq.gz"],
    "2570-KCV-01-26" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/mouseLiverControl4_ATGAGC_L003_R1_001.fastq.gz"],
    "2570-KCV-01-27" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/mouseLiverControl5_ATTCCT_L003_R1_001.fastq.gz"],
    "2571-KCV-1-31"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/HD01_20_CACGAT_L005_R1_001.fastq.gz"],
    "2571-KCV-1-32"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/HD01_21_CACTCA_L005_R1_001.fastq.gz"],
    "2571-KCV-1-33"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/HD01_22_CAGGCG_L005_R1_001.fastq.gz"],
    "2571-KCV-1-34"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/HD01_23_CATGGC_L005_R1_001.fastq.gz"],
    "2571-KCV-1-35"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/HD01_24_CATTTT_L005_R1_001.fastq.gz"],
    "2571-KCV-1-36"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/HD01_25_CCAACA_L005_R1_001.fastq.gz"],
    "2571-KCV-1-37"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/HD01_26_CGGAAT_L005_R1_001.fastq.gz"],
    "2571-KCV-1-38"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/HD01_27_CTAGCT_L005_R1_001.fastq.gz"],
    "2571-KCV-1-39"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/HD01_28_CTATAC_L005_R1_001.fastq.gz"],
    "2571-KCV-1-40"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/NC02_20_CTCAGA_L005_R1_001.fastq.gz"],
    "2571-KCV-1-41"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/NC02_21_GACGAC_L005_R1_001.fastq.gz"],
    "2571-KCV-1-42"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/NC02_22_TAATCG_L005_R1_001.fastq.gz"],
    "2571-KCV-1-43"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/NC02_23_TACAGC_L005_R1_001.fastq.gz"],
    "2571-KCV-1-44"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/NC02_24_TATAAT_L005_R1_001.fastq.gz"],
    "2571-KCV-1-45"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/NC02_25_TCATTC_L005_R1_001.fastq.gz"],
    "2571-KCV-1-46"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/NC02_26_TCCCGA_L005_R1_001.fastq.gz"],
    "2571-KCV-1-47"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/NC02_27_TCGAAG_L005_R1_001.fastq.gz"],
    "2571-KCV-1-48"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/NC02_28_TCGGCA_L005_R1_001.fastq.gz"],
    "2572-KCV-1-19"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/CH003_GTGAAA_L002_R1_001.fastq.gz"],
    "2572-KCV-1-20"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/CH002SS_GTGGCC_L002_R1_001.fastq.gz"],
    "2572-KCV-1-21"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/CH001_GTTTCG_L002_R1_001.fastq.gz"],
    "2572-KCV-1-22"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/PO14s_CGTACG_L002_R1_001.fastq.gz"],
    "2572-KCV-1-23"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/PO16D_GGTAGC_L002_R1_001.fastq.gz"],
    "2572-KCV-1-24"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/PO23F_GAGTGG_L002_R1_001.fastq.gz"],
    "2572-KCV-1-25"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/VK_ACTGAT_L002_R1_001.fastq.gz"],
    "2572-KCV-1-26"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/CC_ATGAGC_L002_R1_001.fastq.gz"],
    "2572-KCV-1-27"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/AE_ATTCCT_L002_R1_001.fastq.gz"],
    "2572-KCV-1-28"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/CHL001_CAAAAG_L002_R1_001.fastq.gz"],
    "2572-KCV-1-29"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/CHL002_CAACTA_L002_R1_001.fastq.gz"],
    "2572-KCV-1-30"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/CHL003_CACCGG_L002_R1_001.fastq.gz"],
  },
  coordinate    => $hsa_gffs,
  bowtie1_index => $bowtie1_human_index,
  shrimp2_index => $shrimp2_human_miRBase_index,
  target_dir    => $target_human_dir,
  task_name     => $task_name . "_human",
};

my $mouse = {
  source => {
    "2570-KCV-01-19" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/mouseLiverTransplant1_GTGAAA_L003_R1_001.fastq.gz"],
    "2570-KCV-01-20" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/mouseLiverTransplant2_GTGGCC_L003_R1_001.fastq.gz"],
    "2570-KCV-01-21" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/mouseLiverTransplant3_GTTTCG_L003_R1_001.fastq.gz"],
    "2570-KCV-01-22" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/mouseLiverTransplant4_CGTACG_L003_R1_001.fastq.gz"],
    "2570-KCV-01-23" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/mouseLiverControl1_GAGTGG_L003_R1_001.fastq.gz"],
    "2570-KCV-01-24" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/mouseLiverControl2_GGTAGC_L003_R1_001.fastq.gz"],
    "2570-KCV-01-25" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/mouseLiverControl3_ACTGAT_L003_R1_001.fastq.gz"],
    "2570-KCV-01-26" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/mouseLiverControl4_ATGAGC_L003_R1_001.fastq.gz"],
    "2570-KCV-01-27" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/mouseLiverControl5_ATTCCT_L003_R1_001.fastq.gz"],
  },
  coordinate    => $mmu_gffs,
  bowtie1_index => $bowtie1_mouse_index,
  shrimp2_index => $shrimp2_mouse_miRBase_index,
  target_dir    => $target_mouse_dir,
  task_name     => $task_name . "_mouse",
};

my @defs = ( $rat, $human, $mouse );

#my @defs = ($mouse);
foreach my $def (@defs) {
  my $cur_target_dir = $def->{target_dir};
  my $config         = {
    general  => { "task_name" => $def->{task_name}, },
    cutadapt => {
      class      => "Cutadapt",
      perform    => 0,
      target_dir => "${target_dir}/cutadapt",
      option     => "-O 10",
      source     => $def->{source},
      adaptor    => "TGGAATTCTCGGGTGCCAAGG",
      extension  => "_clipped.fastq",
      sh_direct  => 1,
      pbs        => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "24",
        "mem"      => "20gb"
      },
    },
    fastqlen => {
      class      => "FastqLen",
      perform    => 0,
      target_dir => "${target_dir}/fastqlen",
      option     => "",
      source_ref => "cutadapt",
      cqstools   => $cqs_tools,
      sh_direct  => 1,
      pbs        => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "24",
        "mem"      => "20gb"
      },
    },
    cutadapt_len => {
      class      => "Cutadapt",
      perform    => 0,
      target_dir => "${target_dir}/cutadapt_len",
      option     => "-O 10 -m 12 -M 49",
      source     => $def->{source},
      adaptor    => "TGGAATTCTCGGGTGCCAAGG",
      extension  => "_clipped.fastq",
      sh_direct  => 1,
      pbs        => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "24",
        "mem"      => "20gb"
      },
    },
    identical => {
      class      => "FastqIdentical",
      perform    => 0,
      target_dir => "${target_dir}/identical",
      option     => "",
      source_ref => "cutadapt_len",
      cqstools   => $cqs_tools,
      extension  => "_clipped_identical.fastq",
      sh_direct  => 1,
      pbs        => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "24",
        "mem"      => "20gb"
      },
    },
    bowtie1_genome_cutadapt_topN => {
      class         => "Bowtie1",
      perform       => 1,
      target_dir    => "${cur_target_dir}/topN_bowtie1_genome_cutadapt",
      option        => $bowtie1_option,
      source_ref    => [ "identical", ".fastq\$" ],
      bowtie1_index => $def->{bowtie1_index},
      samonly       => 0,
      sh_direct     => 0,
      pbs           => {
        "email"    => $email,
        "nodes"    => "1:ppn=8",
        "walltime" => "72",
        "mem"      => "40gb"
      },
    },
    mirna_count_bowtie1_genome_cutadapt_topN => {
      class           => "MirnaCount",
      perform         => 1,
      target_dir      => "${cur_target_dir}/topN_bowtie1_genome_cutadapt_count",
      option          => "-d -t 4",
      source_ref      => "bowtie1_genome_cutadapt_topN",
      fastq_files_ref => "identical",
      seqcount_ref    => [ "identical", ".dupcount\$" ],
      cqs_tools       => $cqs_tools,
      gff_file        => $def->{coordinate},
      sh_direct       => 0,
      pbs             => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "72",
        "mem"      => "40gb"
      },
    },
    shrimp2_miRBase_bowtie1_genome_cutadapt_topN => {
      class         => "Shrimp2",
      perform       => 0,
      target_dir    => "${cur_target_dir}/topN_bowtie1_genome_cutadapt_shrimp2_miRBase",
      option        => $shrimp2_option,
      source_ref    => [ "mirna_count_bowtie1_genome_cutadapt_topN", ".fastq\$" ],
      shrimp2_index => $def->{shrimp2_index},
      is_mirna      => 1,
      output_bam    => 1,
      sh_direct     => 0,
      pbs           => {
        "email"    => $email,
        "nodes"    => "1:ppn=8",
        "walltime" => "720",
        "mem"      => "40gb"
      },
    },
  };

  performConfig($config);
}

1;

