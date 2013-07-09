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
  $root      = "/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna";
  $cqs_tools = "/home/shengq1/cqstools/CQS.Tools.exe";
  $hsa_gffs  = "/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/smrnapipeline/hsa_tableL.bed";
  $rno_gffs  = "/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/smrnapipeline/rno_tableL.bed";
  $mmu_gffs  = "/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/smrnapipeline/rno_tableL.bed";
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

my $bowtie2_option_topN = "-D 20 -R 3 -N 1 -L 12 -i S,1,0.50 --gbar 50 --rdg 1000,1000 --rfg 1000,1000 -k 20 -p 8";

my $bowtie2_rat_index   = "/data/cqs/shengq1/reference/rn4/bowtie2_index/rn4";
my $bowtie2_human_index = "/data/cqs/guoy1/reference/hg19/bowtie2_index/hg19";
my $bowtie2_mouse_index = "/data/cqs/guoy1/reference/mm10/bowtie2_index/mm10";

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
    "2516-01" => ["${root}/data/2516-KCV-1_1.fastq"],
    "2516-02" => ["${root}/data/2516-KCV-2_1.fastq"],
    "2516-03" => ["${root}/data/2516-KCV-3_1.fastq"],
    "2516-04" => ["${root}/data/2516-KCV-4_1.fastq"],
    "2516-05" => ["${root}/data/2516-KCV-5_1.fastq"],
    "2516-06" => ["${root}/data/2516-KCV-6_1.fastq"],
    "2516-07" => ["${root}/data/2516-KCV-7_1.fastq"],
    "2516-08" => ["${root}/data/2516-KCV-8_1.fastq"],
    "2516-09" => ["${root}/data/2516-KCV-9_1.fastq"],
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
  bowtie2_index => $bowtie2_rat_index,
  shrimp2_index => $shrimp2_rat_miRBase_index,
  target_dir    => $target_rat_dir,
  task_name     => $task_name . "_rat",
};
my $human = {
  source => {
    "2516-10"  => ["${root}/data/2516-KCV-10_1.fastq"],
    "2516-11"  => ["${root}/data/2516-KCV-11_1.fastq"],
    "2516-12"  => ["${root}/data/2516-KCV-12_1.fastq"],
    "2516-13"  => ["${root}/data/2516-KCV-13_1.fastq"],
    "2516-14"  => ["${root}/data/2516-KCV-14_1.fastq"],
    "2516-15"  => ["${root}/data/2516-KCV-15_1.fastq"],
    "2516-16"  => ["${root}/data/2516-KCV-16_1.fastq"],
    "2516-17"  => ["${root}/data/2516-KCV-17_1.fastq"],
    "2516-18"  => ["${root}/data/2516-KCV-18_1.fastq"],
    "2516-19"  => ["${root}/data/2516-KCV-19_1.fastq"],
    "2516-20"  => ["${root}/data/2516-KCV-20_1.fastq"],
    "2516-21"  => ["${root}/data/2516-KCV-21_1.fastq"],
    "2516-22"  => ["${root}/data/2516-KCV-22_1.fastq"],
    "2516-23"  => ["${root}/data/2516-KCV-23_1.fastq"],
    "2516-24"  => ["${root}/data/2516-KCV-24_1.fastq"],
    "2516-25"  => ["${root}/data/2516-KCV-25_1.fastq"],
    "2516-26"  => ["${root}/data/2516-KCV-26_1.fastq"],
    "2516-27"  => ["${root}/data/2516-KCV-27_1.fastq"],
    "2516-28"  => ["${root}/data/2516-KCV-28_1.fastq"],
    "2516-29"  => ["${root}/data/2516-KCV-29_1.fastq"],
    "2516-30"  => ["${root}/data/2516-KCV-30_1.fastq"],
    "2516-31"  => ["${root}/data/2516-KCV-31_1.fastq"],
    "2516-32"  => ["${root}/data/2516-KCV-32_1.fastq"],
    "2516-33"  => ["${root}/data/2516-KCV-33_1.fastq"],
    "2516-34"  => ["${root}/data/2516-KCV-34_1.fastq"],
    "2516-35"  => ["${root}/data/2516-KCV-35_1.fastq"],
    "2516-36"  => ["${root}/data/2516-KCV-36_1.fastq"],
    "2516-37"  => ["${root}/data/2516-KCV-37_1.fastq"],
    "2516-38"  => ["${root}/data/2516-KCV-38_1.fastq"],
    "2516-39"  => ["${root}/data/2516-KCV-39_1.fastq"],
    "2516-40"  => ["${root}/data/2516-KCV-40_1.fastq"],
    "2516-41"  => ["${root}/data/2516-KCV-41_1.fastq"],
    "2516-42"  => ["${root}/data/2516-KCV-42_1.fastq"],
    "2516-43"  => ["${root}/data/2516-KCV-43_1.fastq"],
    "2516-44"  => ["${root}/data/2516-KCV-44_1.fastq"],
    "2516-45"  => ["${root}/data/2516-KCV-45_1.fastq"],
    "2516-46"  => ["${root}/data/2516-KCV-46_1.fastq"],
    "2516-47"  => ["${root}/data/2516-KCV-47_1.fastq"],
    "2516-48"  => ["${root}/data/2516-KCV-48_1.fastq"],
    "KCV2_1N2" => ["${root}/data/KCV2_1N2_GCCAAT_L003_R1_001.fastq"],
    "KCV2_1N3" => [ "${root}/data/KCV2_1N3_CAGATC_L003_R1_001.fastq", "${root}/data/KCV2_1N3_CAGATC_L003_R1_002.fastq" ],
    "KCV2_1N4" => [ "${root}/data/KCV2_1N4_ACTTGA_L003_R1_001.fastq", "${root}/data/KCV2_1N4_ACTTGA_L003_R1_002.fastq" ],
    "KCV2_1N5" => [ "${root}/data/KCV2_1N5_GATCAG_L003_R1_001.fastq", "${root}/data/KCV2_1N5_GATCAG_L003_R1_002.fastq" ],
    "KCV2_1N6" => [ "${root}/data/KCV2_1N6_TAGCTT_L003_R1_001.fastq", "${root}/data/KCV2_1N6_TAGCTT_L003_R1_002.fastq" ],
    "KCV2_2N1" => [
      "${root}/data/KCV2_2N1_GGCTAC_L003_R1_001.fastq", "${root}/data/KCV2_2N1_GGCTAC_L003_R1_002.fastq",
      "${root}/data/KCV2_2N1_GGCTAC_L003_R1_003.fastq", "${root}/data/KCV2_2N1_GGCTAC_L003_R1_004.fastq",
      "${root}/data/KCV2_2N1_GGCTAC_L003_R1_005.fastq"
    ],
    "KCV2_2N2" => ["${root}/data/KCV2_2N2_GCCGCG_L003_R1_001.fastq"],
    "KCV2_2N3" => [
      "${root}/data/KCV2_2N3_CTTGTA_L003_R1_001.fastq", "${root}/data/KCV2_2N3_CTTGTA_L003_R1_002.fastq",
      "${root}/data/KCV2_2N3_CTTGTA_L003_R1_003.fastq", "${root}/data/KCV2_2N3_CTTGTA_L003_R1_004.fastq",
      "${root}/data/KCV2_2N3_CTTGTA_L003_R1_005.fastq", "${root}/data/KCV2_2N3_CTTGTA_L003_R1_006.fastq",
      "${root}/data/KCV2_2N3_CTTGTA_L003_R1_007.fastq", "${root}/data/KCV2_2N3_CTTGTA_L003_R1_008.fastq",
      "${root}/data/KCV2_2N3_CTTGTA_L003_R1_009.fastq", "${root}/data/KCV2_2N3_CTTGTA_L003_R1_010.fastq"
    ],
    "KCV2_2N4" => ["${root}/data/KCV2_2N4_GCCTTA_L003_R1_001.fastq"],
    "KCV2_2N5" => ["${root}/data/KCV2_2N5_GCTCCA_L003_R1_001.fastq"],
    "KCV3_1C2" => ["${root}/data/KCV3_1C2_GGCACA_L004_R1_001.fastq"],
    "KCV3_1C3" => ["${root}/data/KCV3_1C3_GGCCTG_L004_R1_001.fastq"],
    "KCV3_1C4" => ["${root}/data/KCV3_1C4_TCTACC_L004_R1_001.fastq"],
    "KCV3_1C5" => ["${root}/data/KCV3_1C5_TGAAGT_L004_R1_001.fastq"],
    "KCV3_1C6" => ["${root}/data/KCV3_1C6_TGCCAT_L004_R1_001.fastq"],
    "KCV3_2C1" => ["${root}/data/KCV3_2C1_TGCTGG_L004_R1_001.fastq"],
    "KCV3_2C2" => ["${root}/data/KCV3_2C2_TGGCGC_L004_R1_001.fastq"],
    "KCV3_2C3" => ["${root}/data/KCV3_2C3_TTCGAA_L004_R1_001.fastq"],
    "Sample1"  => ["${root}/data/Sample1_12.fastq"],
    "Sample2"  => ["${root}/data/Sample2_12.fastq"],
    "Sample3"  => ["${root}/data/Sample3_12.fastq"],
    "Sample4"  => ["${root}/data/Sample4_12.fastq"],
    "Sample5"  => ["${root}/data/Sample5_12.fastq"],
  },
  coordinate    => $hsa_gffs,
  bowtie2_index => $bowtie2_human_index,
  shrimp2_index => $shrimp2_human_miRBase_index,
  target_dir    => $target_human_dir,
  task_name     => $task_name . "_human",
};

my $mouse = {
  source => {
    "Sample_AE"                 => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/AE_ATTCCT_L002_R1_001.fastq.gz"],
    "Sample_CC"                 => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/CC_ATGAGC_L002_R1_001.fastq.gz"],
    "Sample_CH001"              => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/CH001_GTTTCG_L002_R1_001.fastq.gz"],
    "Sample_CH002SS"            => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/CH002SS_GTGGCC_L002_R1_001.fastq.gz"],
    "Sample_CH003"              => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/CH003_GTGAAA_L002_R1_001.fastq.gz"],
    "Sample_CHL001"             => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/CHL001_CAAAAG_L002_R1_001.fastq.gz"],
    "Sample_CHL002"             => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/CHL002_CAACTA_L002_R1_001.fastq.gz"],
    "Sample_CHL003"             => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/CHL003_CACCGG_L002_R1_001.fastq.gz"],
    "Sample_HD01_20"            => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/HD01_20_CACGAT_L005_R1_001.fastq.gz"],
    "Sample_HD01_21"            => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/HD01_21_CACTCA_L005_R1_001.fastq.gz"],
    "Sample_HD01_22"            => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/HD01_22_CAGGCG_L005_R1_001.fastq.gz"],
    "Sample_HD01_23"            => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/HD01_23_CATGGC_L005_R1_001.fastq.gz"],
    "Sample_HD01_24"            => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/HD01_24_CATTTT_L005_R1_001.fastq.gz"],
    "Sample_HD01_25"            => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/HD01_25_CCAACA_L005_R1_001.fastq.gz"],
    "Sample_HD01_26"            => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/HD01_26_CGGAAT_L005_R1_001.fastq.gz"],
    "Sample_HD01_27"            => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/HD01_27_CTAGCT_L005_R1_001.fastq.gz"],
    "Sample_HD01_28"            => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/HD01_28_CTATAC_L005_R1_001.fastq.gz"],
    "Sample_mouseLiverControl1" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/mouseLiverControl1_GAGTGG_L003_R1_001.fastq.gz"],
    "Sample_mouseLiverControl2" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/mouseLiverControl2_GGTAGC_L003_R1_001.fastq.gz"],
    "Sample_mouseLiverControl3" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/mouseLiverControl3_ACTGAT_L003_R1_001.fastq.gz"],
    "Sample_mouseLiverControl4" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/mouseLiverControl4_ATGAGC_L003_R1_001.fastq.gz"],
    "Sample_mouseLiverControl5" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/mouseLiverControl5_ATTCCT_L003_R1_001.fastq.gz"],
  },
  coordinate    => $mmu_gffs,
  bowtie2_index => $bowtie2_mouse_index,
  shrimp2_index => $shrimp2_mouse_miRBase_index,
  target_dir    => $target_mouse_dir,
  task_name     => $task_name . "_mouse",
};

#my @defs = ( $rat, $human, $mouse );
my @defs = ( $mouse );
foreach my $def (@defs) {
  my $config = {
    general  => { "task_name" => $def->{task_name}, },
    cutadapt => {
      class      => "Cutadapt",
      perform    => 1,
      target_dir => "${target_dir}/cutadapt",
      option     => "",
      source     => $def->{source},
      adaptor    => "TGGAATTCTCGGGTGCCAAGG",
      extension  => "_clipped.fastq",
      sh_direct  => 0,
      pbs        => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "24",
        "mem"      => "20gb"
      },
    },
    fastqlen => {
      class      => "FastqLen",
      perform    => 1,
      target_dir => "${target_dir}/fastqlen",
      option     => "",
      source_ref => "cutadapt",
      cqstools   => $cqs_tools,
      sh_direct  => 0,
      pbs        => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "24",
        "mem"      => "20gb"
      },
    },
    cutadapt_len => {
      class      => "Cutadapt",
      perform    => 1,
      target_dir => "${target_dir}/cutadapt_len",
      option     => "-m 12 -M 49",
      source     => $def->{source},
      adaptor    => "TGGAATTCTCGGGTGCCAAGG",
      extension  => "_clipped.fastq",
      sh_direct  => 0,
      pbs        => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "24",
        "mem"      => "20gb"
      },
    },
    bowtie2_genome_cutadapt_topN => {
      class         => "Bowtie2",
      perform       => 1,
      target_dir    => "${target_rat_dir}/topN_bowtie2_genome_cutadapt",
      option        => $bowtie2_option_topN,
      source_ref    => "cutadapt_len",
      bowtie2_index => $def->{bowtie2_index},
      samonly       => 0,
      sh_direct     => 0,
      pbs           => {
        "email"    => $email,
        "nodes"    => "1:ppn=8",
        "walltime" => "72",
        "mem"      => "40gb"
      },
    },
    mirna_count_bowtie2_genome_cutadapt_topN => {
      class        => "MirnaCount",
      perform      => 1,
      target_dir   => "${target_rat_dir}/topN_bowtie2_genome_cutadapt_count",
      option       => "",
      source_ref   => "bowtie2_genome_cutadapt_topN",
      cqs_tools    => $cqs_tools,
      gff_file     => $def->{coordinate},
      fasta_format => 0,
      sh_direct    => 0,
      pbs          => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "24",
        "mem"      => "40gb"
      },
    },
    shrimp2_miRBase_bowtie2_genome_cutadapt_topN => {
      class         => "Shrimp2",
      perform       => 1,
      target_dir    => "${target_rat_dir}/topN_bowtie2_genome_cutadapt_shrimp2_miRBase",
      option        => $shrimp2_option,
      source_ref    => [ "mirna_count_bowtie2_genome_cutadapt_topN", ".fastq\$" ],
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
