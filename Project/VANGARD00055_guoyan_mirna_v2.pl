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
my $cqstools;
my $hsa_gffs;
my $hsa_trna_gffs;
my $rno_gffs;
my $rno_trna_gffs;
my $mmu_gffs;
my $mmu_trna_gffs;

if ( is_linux() ) {
  $root          = "/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2";
  $cqstools      = "/home/shengq1/cqstools/CQS.Tools.exe";
  $hsa_gffs      = "/data/cqs/shengq1/reference/miRBase20/hsa.gff3";
  $hsa_trna_gffs = "/data/cqs/shengq1/reference/trna/hg19_tRNA.bed";
  $rno_gffs      = "/data/cqs/shengq1/reference/miRBase20/rno.gff3";
  $rno_trna_gffs = "/data/cqs/shengq1/reference/trna/rn4_tRNA.bed";
  $mmu_gffs      = "/data/cqs/shengq1/reference/miRBase20/mmu.gff3";
  $mmu_trna_gffs = "/data/cqs/shengq1/reference/trna/mm10_tRNA.bed";
}
else {
  $root     = "d:/temp";
  $cqstools = "E:/sqh/programs/csharp/OmicsLabCSharp/CQS.Tools/bin/Release/CQS.Tools.exe";
  $hsa_gffs = "H:/shengquanhu/projects/vangard/VANGARD00055_guoyan_mirna/hsa.gff3";
  $rno_gffs = "H:/shengquanhu/projects/vangard/VANGARD00055_guoyan_mirna/rno.gff3";
  $mmu_gffs = "H:/shengquanhu/projects/vangard/VANGARD00055_guoyan_mirna/mmu.gff3";
}
my $target_dir = create_directory_or_die($root);

my $target_rat_dir   = create_directory_or_die( $target_dir . "/rat" );
my $target_human_dir = create_directory_or_die( $target_dir . "/human" );
my $target_mouse_dir = create_directory_or_die( $target_dir . "/mouse" );

my $email     = "quanhu.sheng\@vanderbilt.edu";
my $task_name = "VANGARD00055";

my $samtools = "/home/shengq1/local/bin/samtools/samtools";

my $bowtie1_option_pm       = "-a -m 100 --best --strata -v 0 -l 12 -p 8";
my $bowtie1_option_1mm      = "-a -m 100 --best --strata -v 1 -l 12 -p 8";
my $bowtie1_option_1mm_trim = "-a -m 100 --best --strata -v 1 -l 12 -p 8 --trim5 2 --trim3 3";
my $bowtie1_option_3mm      = "-a -m 100 --best --strata -v 3 -l 12 -p 8";

my $bowtie1_rat_index   = "/data/cqs/shengq1/reference/rn4/bowtie1_index/rn4";
my $bowtie1_human_index = "/data/cqs/guoy1/reference/hg19/bowtie_index/hg19";

#my $bowtie1_mouse_index = "/data/cqs/guoy1/reference/mm9/bowtie_index/mm9";
my $bowtie1_mouse_index = "/data/cqs/shengq1/reference/mm10/bowtie_index/mm10";

my $bowtie2_option       = "-D 100 -R 3 -N 1 -L 8 -i S,1,0.50 --gbar 50 --rdg 1000,1000 --rfg 1000,1000 -k 20 -p 8";
my $bowtie2_local_option = "$bowtie2_option --local";

my $bowtie2_rat_index   = "/data/cqs/shengq1/reference/rn4/bowtie2_index/rn4";
my $bowtie2_human_index = "/data/cqs/guoy1/reference/hg19/bowtie2_index/hg19";

#my $bowtie2_mouse_index = "/data/cqs/guoy1/reference/mm9/bowtie2_index/mm9";
my $bowtie2_mouse_index = "/data/cqs/guoy1/reference/mm10/bowtie2_index/mm10";

my $mirnacount_option          = "-s";                                                #ignore score
my $trnacount_option           = "-s";                                                #ignore score
my $mirna_overlap_count_option = "-s --gtf_key miRNA";
my $fasta_file                 = "/data/cqs/shengq1/reference/miRBase20/mature.fa";

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

my $bwa_option       = "-o 0 -l 8 -n 3 -t 8";
my $bwa_hsammu_fasta = "/data/cqs/shengq1/reference/hg19mm9/bwa_0.7.4_index/hg19mm9.fa";
my $hsammu_gffs      = "/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/smrnapipeline/hsa_mmu_tableL.bed";
my $bwa_clip_option  = "-o 2 -e 3 -l 8 -n 3 -t 8";

my $rat = {
  notidentical_search => {
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
  source => {
    "NIH_Rat_HDL_01"            => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/20130925_NIH_Rat_HDL/Vickers_Rat_HDL_1_ATCACG_L001_R1.fastq.gz"],
    "NIH_Rat_HDL_02"            => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/20130925_NIH_Rat_HDL/Vickers_Rat_HDL_2_CGATGT_L001_R1.fastq.gz"],
    "NIH_Rat_HDL_03"            => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/20130925_NIH_Rat_HDL/Vickers_Rat_HDL_3_TTAGGC_L001_R1.fastq.gz"],
    "NIH_Rat_HDL_04"            => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/20130925_NIH_Rat_HDL/Vickers_Rat_HDL_4_TGACCA_L001_R1.fastq.gz"],
    "NIH_Rat_HDL_05"            => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/20130925_NIH_Rat_HDL/Vickers_Rat_HDL_5_ACAGTG_L001_R1.fastq.gz"],
    "NIH_Rat_HDL_06"            => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/20130925_NIH_Rat_HDL/Vickers_Rat_HDL_6_GCCAAT_L001_R1.fastq.gz"],
    "NIH_Rat_HDL_07"            => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/20130925_NIH_Rat_HDL/Vickers_Rat_HDL_7_CAGATC_L001_R1.fastq.gz"],
    "NIH_Rat_HDL_08"            => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/20130925_NIH_Rat_HDL/Vickers_Rat_HDL_8_ACTTGA_L001_R1.fastq.gz"],
    "NIH_Rat_HDL_09"            => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/20130925_NIH_Rat_HDL/Vickers_Rat_HDL_9_GATCAG_L001_R1.fastq.gz"],
    "NIH_Rat_HDL_10"            => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/20130925_NIH_Rat_HDL/Vickers_Rat_HDL_10_TAGCTT_L001_R1.fastq.gz"],
    "NIH_Rat_HDL_11"            => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/20130925_NIH_Rat_HDL/Vickers_Rat_HDL_11_GGCTAC_L001_R1.fastq.gz"],
    "2516-01"                   => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/2516-KCV-1_1.fastq"],
    "2516-02"                   => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/2516-KCV-2_1.fastq"],
    "2516-03"                   => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/2516-KCV-3_1.fastq"],
    "2516-04"                   => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/2516-KCV-4_1.fastq"],
    "2516-05"                   => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/2516-KCV-5_1.fastq"],
    "2516-06"                   => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/2516-KCV-6_1.fastq"],
    "2516-07"                   => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/2516-KCV-7_1.fastq"],
    "2516-08"                   => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/2516-KCV-8_1.fastq"],
    "2516-09"                   => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/2516-KCV-9_1.fastq"],
    "2570-01"                   => ["/autofs/blue_sequencer/Runs/projects/2570-KCV/2013-06-19/2570-KCV-1_1.fastq.gz"],
    "2570-02"                   => ["/autofs/blue_sequencer/Runs/projects/2570-KCV/2013-06-19/2570-KCV-2_1.fastq.gz"],
    "2570-03"                   => ["/autofs/blue_sequencer/Runs/projects/2570-KCV/2013-06-19/2570-KCV-3_1.fastq.gz"],
    "2570-04"                   => ["/autofs/blue_sequencer/Runs/projects/2570-KCV/2013-06-19/2570-KCV-4_1.fastq.gz"],
    "2570-05"                   => ["/autofs/blue_sequencer/Runs/projects/2570-KCV/2013-06-19/2570-KCV-5_1.fastq.gz"],
    "2570-06"                   => ["/autofs/blue_sequencer/Runs/projects/2570-KCV/2013-06-19/2570-KCV-6_1.fastq.gz"],
    "2570-07"                   => ["/autofs/blue_sequencer/Runs/projects/2570-KCV/2013-06-19/2570-KCV-7_1.fastq.gz"],
    "2570-08"                   => ["/autofs/blue_sequencer/Runs/projects/2570-KCV/2013-06-19/2570-KCV-8_1.fastq.gz"],
    "2570-09"                   => ["/autofs/blue_sequencer/Runs/projects/2570-KCV/2013-06-19/2570-KCV-9_1.fastq.gz"],
    "2570-10"                   => ["/autofs/blue_sequencer/Runs/projects/2570-KCV/2013-06-19/2570-KCV-10_1.fastq.gz"],
    "2570-11"                   => ["/autofs/blue_sequencer/Runs/projects/2570-KCV/2013-06-19/2570-KCV-11_1.fastq.gz"],
    "2570-12"                   => ["/autofs/blue_sequencer/Runs/projects/2570-KCV/2013-06-19/2570-KCV-12_1.fastq.gz"],
    "2570-13"                   => ["/autofs/blue_sequencer/Runs/projects/2570-KCV/2013-06-19/2570-KCV-13_1.fastq.gz"],
    "2570-14"                   => ["/autofs/blue_sequencer/Runs/projects/2570-KCV/2013-06-19/2570-KCV-14_1.fastq.gz"],
    "2570-15"                   => ["/autofs/blue_sequencer/Runs/projects/2570-KCV/2013-06-19/2570-KCV-15_1.fastq.gz"],
    "2570-16"                   => ["/autofs/blue_sequencer/Runs/projects/2570-KCV/2013-06-19/2570-KCV-16_1.fastq.gz"],
    "2570-17"                   => ["/autofs/blue_sequencer/Runs/projects/2570-KCV/2013-06-19/2570-KCV-17_1.fastq.gz"],
    "2570-18"                   => ["/autofs/blue_sequencer/Runs/projects/2570-KCV/2013-06-19/2570-KCV-18_1.fastq.gz"],
    "2261-ASK-102_1"            => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/2261-ASK-102_1_sequence.txt.gz"],
    "2261-ASK-103_1"            => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/2261-ASK-103_1_sequence.txt.gz"],
    "2261-ASK-57_1"             => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/2261-ASK-57_1_sequence.txt.gz"],
    "2261-ASK-59_1"             => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/2261-ASK-59_1_sequence.txt.gz"],
    "2261-ASK-61_1"             => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/2261-ASK-61_1_sequence.txt.gz"],
    "2261-ASK-62_1"             => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/2261-ASK-62_1_sequence.txt.gz"],
    "2261-ASK-63_1"             => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/2261-ASK-63_1_sequence.txt.gz"],
    "2261-ASK-66_1"             => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/2261-ASK-66_1_sequence.txt.gz"],
    "2261-ASK-68_1"             => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/2261-ASK-68_1_sequence.txt.gz"],
    "2261-ASK-71_1"             => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/2261-ASK-71_1_sequence.txt.gz"],
    "2261-ASK-82_1"             => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/2261-ASK-82_1_sequence.txt.gz"],
    "2261-ASK-93_1"             => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/2261-ASK-93_1_sequence.txt.gz"],
    "2261-ASK-94_1"             => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/2261-ASK-94_1_sequence.txt.gz"],
    "2261-ASK-95_1"             => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/2261-ASK-95_1_sequence.txt.gz"],
    "2261-ASK-96_1"             => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/2261-ASK-96_1_sequence.txt.gz"],
    "2288-RDB-54_1"             => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/2288-RDB-54_1_sequence.txt.gz"],
    "2288-RDB-55_1"             => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/2288-RDB-55_1_sequence.txt.gz"],
    "2288-RDB-56_1"             => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/2288-RDB-56_1_sequence.txt.gz"],
    "2288-RDB-57_1"             => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/2288-RDB-57_1_sequence.txt.gz"],
    "2288-RDB-58_1"             => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/2288-RDB-58_1_sequence.txt.gz"],
    "2288-RDB-59_1"             => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/2288-RDB-59_1_sequence.txt.gz"],
    "2288-RDB-60_1"             => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/2288-RDB-60_1_sequence.txt.gz"],
    "2288-RDB-61_1"             => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/2288-RDB-61_1_sequence.txt.gz"],
    "2634-SH-1_1"               => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/2634-SH-1_1_sequence.txt.gz"],
    "2634-SH-2_1"               => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/2634-SH-2_1_sequence.txt.gz"],
    "2634-SH-3_1"               => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/2634-SH-3_1_sequence.txt.gz"],
    "2634-SH-4_1"               => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/2634-SH-4_1_sequence.txt.gz"],
    "2634-SH-5_1"               => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/2634-SH-5_1_sequence.txt.gz"],
    "2634-SH-6_1"               => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/2634-SH-6_1_sequence.txt.gz"],
    "2646-ASK-1_1"              => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/2646-ASK-1_1_sequence.txt.gz"],
    "2646-ASK-2_1"              => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/2646-ASK-2_1_sequence.txt.gz"],
    "2646-ASK-3_1"              => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/2646-ASK-3_1_sequence.txt.gz"],
    "2646-ASK-4_1"              => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/2646-ASK-4_1_sequence.txt.gz"],
    "2646-ASK-5_1"              => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/2646-ASK-5_1_sequence.txt.gz"],
    "2661-AR-1_1"               => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/2661-AR-1_1_sequence.txt.gz"],
    "2661-AR-2_1"               => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/2661-AR-2_1_sequence.txt.gz"],
    "2661-AR-3_1"               => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/2661-AR-3_1_sequence.txt.gz"],
    "2661-AR-4_1"               => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/2661-AR-4_1_sequence.txt.gz"],
    "Colesevelam-3251.1_CGATGT" => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/Colesevelam-3251.1_CGATGT_L003_R1_001.fastq.gz"],
    "Colesevelam-3251.2_ATCACG" => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/Colesevelam-3251.2_ATCACG_L003_R1_001.fastq.gz"],
    "Colesevelam-3258.1_ACAGTG" => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/Colesevelam-3258.1_ACAGTG_L003_R1_001.fastq.gz"],
    "Colesevelam-3258.2_TGACCA" => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/Colesevelam-3258.2_TGACCA_L003_R1_001.fastq.gz"],
    "Colesevelam-3269.1_GCCAAT" => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/Colesevelam-3269.1_GCCAAT_L003_R1_001.fastq.gz"],
    "Colesevelam-3269.2_TTAGGC" => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/Colesevelam-3269.2_TTAGGC_L003_R1_001.fastq.gz"],
    "Lean-3183.1_GATCAG"        => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/Lean-3183.1_GATCAG_L003_R1_001.fastq.gz"],
    "Lean-3185.1_CAGATC"        => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/Lean-3185.1_CAGATC_L003_R1_001.fastq.gz"],
    "Lean-3185.2_ACTTGA"        => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/Lean-3185.2_ACTTGA_L003_R1_001.fastq.gz"],
    "Lean-KR-4073.1_TAGCTT"     => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/Lean-KR-4073.1_TAGCTT_L003_R1_001.fastq.gz"],
    "Lean-KR-4073.2_CTTGTA"     => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/Lean-KR-4073.2_CTTGTA_L003_R1_001.fastq.gz"],
    "Lean-KR-4074.2_AGTCAA"     => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/Lean-KR-4074.2_AGTCAA_L003_R1_001.fastq.gz"],
    "Lean-KR-4074.3_GGCTAC"     => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/Lean-KR-4074.3_GGCTAC_L003_R1_001.fastq.gz"],
    "Vehicle-3259.1_CCGTCC"     => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/Vehicle-3259.1_CCGTCC_L003_R1_001.fastq.gz"],
    "Vehicle-3259.2_ATGTCA"     => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/Vehicle-3259.2_ATGTCA_L003_R1_001.fastq.gz"],
    "Vehicle-3263.1_GTGAAA"     => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/Vehicle-3263.1_GTGAAA_L003_R1_001.fastq.gz"],
    "Vehicle-3263.2_AGTTCC"     => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/Vehicle-3263.2_AGTTCC_L003_R1_001.fastq.gz"],
    "Vehicle-3271.1_GTCCGC"     => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/Vehicle-3271.1_GTCCGC_L003_R1_001.fastq.gz"],
    "Vehicle-3271.2_GTAGAG"     => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/Vehicle-3271.2_GTAGAG_L003_R1_001.fastq.gz"],
  },
  coordinate      => $rno_gffs,
  trna_coordinate => $rno_trna_gffs,
  bowtie1_index   => $bowtie1_rat_index,
  bowtie2_index   => $bowtie2_rat_index,
  shrimp2_index   => $shrimp2_rat_miRBase_index,
  target_dir      => $target_rat_dir,
  task_name       => $task_name . "_rat",
};
my $human = {
  notidentical_search => {
    "2570-KCV-01-19"     => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/mouseLiverTransplant1_GTGAAA_L003_R1_001.fastq.gz"],
    "2570-KCV-01-20"     => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/mouseLiverTransplant2_GTGGCC_L003_R1_001.fastq.gz"],
    "2570-KCV-01-21"     => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/mouseLiverTransplant3_GTTTCG_L003_R1_001.fastq.gz"],
    "2570-KCV-01-22"     => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/mouseLiverTransplant4_CGTACG_L003_R1_001.fastq.gz"],
    "2570-KCV-01-23"     => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/mouseLiverControl1_GAGTGG_L003_R1_001.fastq.gz"],
    "2570-KCV-01-24"     => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/mouseLiverControl2_GGTAGC_L003_R1_001.fastq.gz"],
    "2570-KCV-01-25"     => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/mouseLiverControl3_ACTGAT_L003_R1_001.fastq.gz"],
    "2570-KCV-01-26"     => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/mouseLiverControl4_ATGAGC_L003_R1_001.fastq.gz"],
    "2570-KCV-01-27"     => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/mouseLiverControl5_ATTCCT_L003_R1_001.fastq.gz"],
    "01-018-Post_CTTGTA" => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/01-018-Post_CTTGTA_L005_R1_001.fastq.gz"],
    "01-018-Pre_GGCTAC"  => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/01-018-Pre_GGCTAC_L005_R1_001.fastq.gz"],
    "01-031-Post_GGTAGC" => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/01-031-Post_GGTAGC_L005_R1_001.fastq.gz"],
    "01-031-Pre_GAGTGG"  => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/01-031-Pre_GAGTGG_L005_R1_001.fastq.gz"],
    "01-061-Post_ATGAGC" => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/01-061-Post_ATGAGC_L004_R1_001.fastq.gz"],
    "01-061-Pre_ACTGAT"  => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/01-061-Pre_ACTGAT_L004_R1_001.fastq.gz"],
    "01-28-Post_CGATGT"  => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/01-28-Post_CGATGT_L004_R1_001.fastq.gz"],
    "01-28-Pre_ATCACG"   => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/01-28-Pre_ATCACG_L004_R1_001.fastq.gz"],
    "01-29-Pre_TGACCA"   => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/01-29-Pre_TGACCA_L005_R1_001.fastq.gz"],
    "01-29-Pre_TTAGGC"   => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/01-29-Pre_TTAGGC_L005_R1_001.fastq.gz"],
    "01-36-Post_GCCAAT"  => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/01-36-Post_GCCAAT_L004_R1_001.fastq.gz"],
    "01-36-Pre_ACAGTG"   => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/01-36-Pre_ACAGTG_L004_R1_001.fastq.gz"],
    "03-007-Post_CAAAAG" => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/03-007-Post_CAAAAG_L005_R1_001.fastq.gz"],
    "03-007-Pre_ATTCCT"  => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/03-007-Pre_ATTCCT_L005_R1_001.fastq.gz"],
    "03-011-Post_CACCGG" => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/03-011-Post_CACCGG_L004_R1_001.fastq.gz"],
    "03-011-Pre_CAACTA"  => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/03-011-Pre_CAACTA_L004_R1_001.fastq.gz"],
    "03-015-Post_CACTCA" => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/03-015-Post_CACTCA_L005_R1_001.fastq.gz"],
    "03-015-Pre_CACGAT"  => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/03-015-Pre_CACGAT_L005_R1_001.fastq.gz"],
    "03-018-Post_CGTACG" => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/03-018-Post_CGTACG_L004_R1_001.fastq.gz"],
    "03-018-Pre_GTTTCG"  => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/03-018-Pre_GTTTCG_L004_R1_001.fastq.gz"],
    "03-026-Post_AGTTCC" => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/03-026-Post_AGTTCC_L004_R1_001.fastq.gz"],
    "03-026-Pre_AGTCAA"  => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/03-026-Pre_AGTCAA_L004_R1_001.fastq.gz"],
    "03-031-Post_CATGGC" => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/03-031-Post_CATGGC_L004_R1_001.fastq.gz"],
    "03-031-Pre_CAGGCG"  => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/03-031-Pre_CAGGCG_L004_R1_001.fastq.gz"],
    "03-033-Post_CCAACA" => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/03-033-Post_CCAACA_L005_R1_001.fastq.gz"],
    "03-033-Pre_CATTTT"  => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/03-033-Pre_CATTTT_L005_R1_001.fastq.gz"],
    "03-036-Post_CCGTCC" => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/03-036-Post_CCGTCC_L005_R1_001.fastq.gz"],
    "03-036-Pre_ATGTCA"  => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/03-036-Pre_ATGTCA_L005_R1_001.fastq.gz"],
    "03-047-Post_CTAGCT" => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/03-047-Post_CTAGCT_L004_R1_001.fastq.gz"],
    "03-047-Pre_CGGAAT"  => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/03-047-Pre_CGGAAT_L004_R1_001.fastq.gz"],
    "03-049-Post_CTCAGA" => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/03-049-Post_CTCAGA_L005_R1_001.fastq.gz"],
    "03-049-Pre_CTATAC"  => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/03-049-Pre_CTATAC_L005_R1_001.fastq.gz"],
    "03-063-Post_GTGGCC" => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/03-063-Post_GTGGCC_L005_R1_001.fastq.gz"],
    "03-063-Pre_GTGAAA"  => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/03-063-Pre_GTGAAA_L005_R1_001.fastq.gz"],
    "03-065-Post_GTCCGC" => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/03-065-Post_GTCCGC_L004_R1_001.fastq.gz"],
    "03-065-Pre_GTAGAG"  => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/03-065-Pre_GTAGAG_L004_R1_001.fastq.gz"],
    "03-16-Post_ACTTGA"  => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/03-16-Post_ACTTGA_L005_R1_001.fastq.gz"],
    "03-16-Pre_CAGATC"   => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/03-16-Pre_CAGATC_L005_R1_001.fastq.gz"],
    "03-17-Post_TAGCTT"  => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/03-17-Post_TAGCTT_L004_R1_001.fastq.gz"],
    "03-17-Pre_GATCAG"   => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/03-17-Pre_GATCAG_L004_R1_001.fastq.gz"],
  },
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
    "KCV2_2N4"           => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/KCV2_2N4_GCCTTA_L003_R1_001.fastq"],
    "KCV2_2N5"           => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/KCV2_2N5_GCTCCA_L003_R1_001.fastq"],
    "KCV3_1C2"           => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/KCV3_1C2_GGCACA_L004_R1_001.fastq"],
    "KCV3_1C3"           => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/KCV3_1C3_GGCCTG_L004_R1_001.fastq"],
    "KCV3_1C4"           => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/KCV3_1C4_TCTACC_L004_R1_001.fastq"],
    "KCV3_1C5"           => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/KCV3_1C5_TGAAGT_L004_R1_001.fastq"],
    "KCV3_1C6"           => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/KCV3_1C6_TGCCAT_L004_R1_001.fastq"],
    "KCV3_2C1"           => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/KCV3_2C1_TGCTGG_L004_R1_001.fastq"],
    "KCV3_2C2"           => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/KCV3_2C2_TGGCGC_L004_R1_001.fastq"],
    "KCV3_2C3"           => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/KCV3_2C3_TTCGAA_L004_R1_001.fastq"],
    "Sample1"            => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/Sample1_12.fastq"],
    "Sample2"            => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/Sample2_12.fastq"],
    "Sample3"            => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/Sample3_12.fastq"],
    "Sample4"            => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/Sample4_12.fastq"],
    "Sample5"            => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/Sample5_12.fastq"],
    "2570-KCV-01-19"     => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/mouseLiverTransplant1_GTGAAA_L003_R1_001.fastq.gz"],
    "2570-KCV-01-20"     => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/mouseLiverTransplant2_GTGGCC_L003_R1_001.fastq.gz"],
    "2570-KCV-01-21"     => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/mouseLiverTransplant3_GTTTCG_L003_R1_001.fastq.gz"],
    "2570-KCV-01-22"     => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/mouseLiverTransplant4_CGTACG_L003_R1_001.fastq.gz"],
    "2570-KCV-01-23"     => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/mouseLiverControl1_GAGTGG_L003_R1_001.fastq.gz"],
    "2570-KCV-01-24"     => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/mouseLiverControl2_GGTAGC_L003_R1_001.fastq.gz"],
    "2570-KCV-01-25"     => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/mouseLiverControl3_ACTGAT_L003_R1_001.fastq.gz"],
    "2570-KCV-01-26"     => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/mouseLiverControl4_ATGAGC_L003_R1_001.fastq.gz"],
    "2570-KCV-01-27"     => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/mouseLiverControl5_ATTCCT_L003_R1_001.fastq.gz"],
    "2571-KCV-1-31"      => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/HD01_20_CACGAT_L005_R1_001.fastq.gz"],
    "2571-KCV-1-32"      => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/HD01_21_CACTCA_L005_R1_001.fastq.gz"],
    "2571-KCV-1-33"      => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/HD01_22_CAGGCG_L005_R1_001.fastq.gz"],
    "2571-KCV-1-34"      => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/HD01_23_CATGGC_L005_R1_001.fastq.gz"],
    "2571-KCV-1-35"      => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/HD01_24_CATTTT_L005_R1_001.fastq.gz"],
    "2571-KCV-1-36"      => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/HD01_25_CCAACA_L005_R1_001.fastq.gz"],
    "2571-KCV-1-37"      => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/HD01_26_CGGAAT_L005_R1_001.fastq.gz"],
    "2571-KCV-1-38"      => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/HD01_27_CTAGCT_L005_R1_001.fastq.gz"],
    "2571-KCV-1-39"      => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/HD01_28_CTATAC_L005_R1_001.fastq.gz"],
    "2571-KCV-1-40"      => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/NC02_20_CTCAGA_L005_R1_001.fastq.gz"],
    "2571-KCV-1-41"      => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/NC02_21_GACGAC_L005_R1_001.fastq.gz"],
    "2571-KCV-1-42"      => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/NC02_22_TAATCG_L005_R1_001.fastq.gz"],
    "2571-KCV-1-43"      => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/NC02_23_TACAGC_L005_R1_001.fastq.gz"],
    "2571-KCV-1-44"      => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/NC02_24_TATAAT_L005_R1_001.fastq.gz"],
    "2571-KCV-1-45"      => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/NC02_25_TCATTC_L005_R1_001.fastq.gz"],
    "2571-KCV-1-46"      => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/NC02_26_TCCCGA_L005_R1_001.fastq.gz"],
    "2571-KCV-1-47"      => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/NC02_27_TCGAAG_L005_R1_001.fastq.gz"],
    "2571-KCV-1-48"      => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/NC02_28_TCGGCA_L005_R1_001.fastq.gz"],
    "2572-KCV-1-19"      => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/CH003_GTGAAA_L002_R1_001.fastq.gz"],
    "2572-KCV-1-20"      => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/CH002SS_GTGGCC_L002_R1_001.fastq.gz"],
    "2572-KCV-1-21"      => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/CH001_GTTTCG_L002_R1_001.fastq.gz"],
    "2572-KCV-1-22"      => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/PO14s_CGTACG_L002_R1_001.fastq.gz"],
    "2572-KCV-1-23"      => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/PO16D_GGTAGC_L002_R1_001.fastq.gz"],
    "2572-KCV-1-24"      => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/PO23F_GAGTGG_L002_R1_001.fastq.gz"],
    "2572-KCV-1-25"      => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/VK_ACTGAT_L002_R1_001.fastq.gz"],
    "2572-KCV-1-26"      => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/CC_ATGAGC_L002_R1_001.fastq.gz"],
    "2572-KCV-1-27"      => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/AE_ATTCCT_L002_R1_001.fastq.gz"],
    "2572-KCV-1-28"      => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/CHL001_CAAAAG_L002_R1_001.fastq.gz"],
    "2572-KCV-1-29"      => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/CHL002_CAACTA_L002_R1_001.fastq.gz"],
    "2572-KCV-1-30"      => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/data/VickersTemp/CHL003_CACCGG_L002_R1_001.fastq.gz"],
    "01-018-Post_CTTGTA" => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/01-018-Post_CTTGTA_L005_R1_001.fastq.gz"],
    "01-018-Pre_GGCTAC"  => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/01-018-Pre_GGCTAC_L005_R1_001.fastq.gz"],
    "01-031-Post_GGTAGC" => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/01-031-Post_GGTAGC_L005_R1_001.fastq.gz"],
    "01-031-Pre_GAGTGG"  => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/01-031-Pre_GAGTGG_L005_R1_001.fastq.gz"],
    "01-061-Post_ATGAGC" => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/01-061-Post_ATGAGC_L004_R1_001.fastq.gz"],
    "01-061-Pre_ACTGAT"  => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/01-061-Pre_ACTGAT_L004_R1_001.fastq.gz"],
    "01-28-Post_CGATGT"  => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/01-28-Post_CGATGT_L004_R1_001.fastq.gz"],
    "01-28-Pre_ATCACG"   => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/01-28-Pre_ATCACG_L004_R1_001.fastq.gz"],
    "01-29-Pre_TGACCA"   => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/01-29-Pre_TGACCA_L005_R1_001.fastq.gz"],
    "01-29-Pre_TTAGGC"   => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/01-29-Pre_TTAGGC_L005_R1_001.fastq.gz"],
    "01-36-Post_GCCAAT"  => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/01-36-Post_GCCAAT_L004_R1_001.fastq.gz"],
    "01-36-Pre_ACAGTG"   => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/01-36-Pre_ACAGTG_L004_R1_001.fastq.gz"],
    "03-007-Post_CAAAAG" => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/03-007-Post_CAAAAG_L005_R1_001.fastq.gz"],
    "03-007-Pre_ATTCCT"  => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/03-007-Pre_ATTCCT_L005_R1_001.fastq.gz"],
    "03-011-Post_CACCGG" => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/03-011-Post_CACCGG_L004_R1_001.fastq.gz"],
    "03-011-Pre_CAACTA"  => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/03-011-Pre_CAACTA_L004_R1_001.fastq.gz"],
    "03-015-Post_CACTCA" => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/03-015-Post_CACTCA_L005_R1_001.fastq.gz"],
    "03-015-Pre_CACGAT"  => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/03-015-Pre_CACGAT_L005_R1_001.fastq.gz"],
    "03-018-Post_CGTACG" => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/03-018-Post_CGTACG_L004_R1_001.fastq.gz"],
    "03-018-Pre_GTTTCG"  => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/03-018-Pre_GTTTCG_L004_R1_001.fastq.gz"],
    "03-026-Post_AGTTCC" => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/03-026-Post_AGTTCC_L004_R1_001.fastq.gz"],
    "03-026-Pre_AGTCAA"  => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/03-026-Pre_AGTCAA_L004_R1_001.fastq.gz"],
    "03-031-Post_CATGGC" => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/03-031-Post_CATGGC_L004_R1_001.fastq.gz"],
    "03-031-Pre_CAGGCG"  => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/03-031-Pre_CAGGCG_L004_R1_001.fastq.gz"],
    "03-033-Post_CCAACA" => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/03-033-Post_CCAACA_L005_R1_001.fastq.gz"],
    "03-033-Pre_CATTTT"  => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/03-033-Pre_CATTTT_L005_R1_001.fastq.gz"],
    "03-036-Post_CCGTCC" => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/03-036-Post_CCGTCC_L005_R1_001.fastq.gz"],
    "03-036-Pre_ATGTCA"  => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/03-036-Pre_ATGTCA_L005_R1_001.fastq.gz"],
    "03-047-Post_CTAGCT" => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/03-047-Post_CTAGCT_L004_R1_001.fastq.gz"],
    "03-047-Pre_CGGAAT"  => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/03-047-Pre_CGGAAT_L004_R1_001.fastq.gz"],
    "03-049-Post_CTCAGA" => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/03-049-Post_CTCAGA_L005_R1_001.fastq.gz"],
    "03-049-Pre_CTATAC"  => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/03-049-Pre_CTATAC_L005_R1_001.fastq.gz"],
    "03-063-Post_GTGGCC" => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/03-063-Post_GTGGCC_L005_R1_001.fastq.gz"],
    "03-063-Pre_GTGAAA"  => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/03-063-Pre_GTGAAA_L005_R1_001.fastq.gz"],
    "03-065-Post_GTCCGC" => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/03-065-Post_GTCCGC_L004_R1_001.fastq.gz"],
    "03-065-Pre_GTAGAG"  => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/03-065-Pre_GTAGAG_L004_R1_001.fastq.gz"],
    "03-16-Post_ACTTGA"  => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/03-16-Post_ACTTGA_L005_R1_001.fastq.gz"],
    "03-16-Pre_CAGATC"   => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/03-16-Pre_CAGATC_L005_R1_001.fastq.gz"],
    "03-17-Post_TAGCTT"  => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/03-17-Post_TAGCTT_L004_R1_001.fastq.gz"],
    "03-17-Pre_GATCAG"   => ["/autofs/blue_sequencer/Runs/130823_SN508_0279_AD2BAFACXX/publish/03-17-Pre_GATCAG_L004_R1_001.fastq.gz"],
  },
  coordinate      => $hsa_gffs,
  trna_coordinate => $hsa_trna_gffs,
  bowtie1_index   => $bowtie1_human_index,
  bowtie2_index   => $bowtie2_human_index,
  shrimp2_index   => $shrimp2_human_miRBase_index,
  target_dir      => $target_human_dir,
  task_name       => $task_name . "_human",
};

my $mouse = {
  notidentical_search => {
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
  coordinate      => $mmu_gffs,
  trna_coordinate => $mmu_trna_gffs,
  bowtie1_index   => $bowtie1_mouse_index,
  bowtie2_index   => $bowtie2_mouse_index,
  shrimp2_index   => $shrimp2_mouse_miRBase_index,
  target_dir      => $target_mouse_dir,
  task_name       => $task_name . "_mouse",
};

#my @defs = ( $rat, $human, $mouse );

my @defs = ($rat);
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
      cqstools   => $cqstools,
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
      cqstools   => $cqstools,
      extension  => "_clipped_identical.fastq",
      sh_direct  => 1,
      pbs        => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "24",
        "mem"      => "20gb"
      },
    },

    #perfect match
    bowtie1_genome_cutadapt_topN_pm => {
      class         => "Bowtie1",
      perform       => 0,
      target_dir    => "${cur_target_dir}/topN_bowtie1_genome_cutadapt_pm",
      option        => $bowtie1_option_pm,
      source_ref    => [ "identical", ".fastq\$" ],
      bowtie1_index => $def->{bowtie1_index},
      samonly       => 0,
      sh_direct     => 1,
      pbs           => {
        "email"    => $email,
        "nodes"    => "1:ppn=8",
        "walltime" => "72",
        "mem"      => "40gb"
      },
    },
    miRNA_overlap_count_bowtie1_genome_cutadapt_topN_pm => {
      class           => "CQSMappedCount",
      perform         => 0,
      target_dir      => "${cur_target_dir}/topN_bowtie1_genome_cutadapt_pm_count_miRNA_overlap",
      option          => $mirna_overlap_count_option,
      source_ref      => "bowtie1_genome_cutadapt_topN_pm",
      fastq_files_ref => "identical",
      seqcount_ref    => [ "identical", ".dupcount\$" ],
      cqs_tools       => $cqstools,
      gff_file        => $def->{coordinate},
      samtools        => $samtools,
      sh_direct       => 1,
      pbs             => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "72",
        "mem"      => "40gb"
      },
    },
    miRNA_overlap_position_pm => {
      class      => "CQSMappedPosition",
      perform    => 0,
      target_dir => "${cur_target_dir}/topN_bowtie1_genome_cutadapt_pm_count_miRNA_overlap_position",
      option     => "-o " . $def->{task_name} . "_miRNA.position",
      source_ref => "miRNA_overlap_count_bowtie1_genome_cutadapt_topN_pm",
      cqs_tools  => $cqstools,
      sh_direct  => 1,
      pbs        => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "10",
        "mem"      => "10gb"
      },
    },
    tRNA_overlap_count_bowtie1_genome_cutadapt_topN_pm => {
      class           => "CQSMappedCount",
      perform         => 0,
      target_dir      => "${cur_target_dir}/topN_bowtie1_genome_cutadapt_pm_count_tRNA",
      option          => $trnacount_option,
      source_ref      => "bowtie1_genome_cutadapt_topN_pm",
      fastq_files_ref => "identical",
      seqcount_ref    => [ "identical", ".dupcount\$" ],
      cqs_tools       => $cqstools,
      gff_file        => $def->{trna_coordinate},
      samtools        => $samtools,
      sh_direct       => 1,
      pbs             => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "72",
        "mem"      => "40gb"
      },
    },
    tRNA_position_pm => {
      class      => "CQSMappedPosition",
      perform    => 0,
      target_dir => "${cur_target_dir}/topN_bowtie1_genome_cutadapt_pm_count_tRNA_position",
      option     => "-o " . $def->{task_name} . "_tRNA.position",
      source_ref => "tRNA_overlap_count_bowtie1_genome_cutadapt_topN_pm",
      cqs_tools  => $cqstools,
      sh_direct  => 1,
      pbs        => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "10",
        "mem"      => "10gb"
      },
    },

    #1 mismatch notidentical search
    cutadapt_len_fake => {
      class      => "Cutadapt",
      perform    => 0,
      target_dir => "${target_dir}/cutadapt_len",
      option     => "-O 10 -m 12 -M 49",
      source     => $def->{notidentical_search},
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
    bowtie1_genome_cutadapt_topN_1mm_notidentical => {
      class         => "Bowtie1",
      perform       => 0,
      target_dir    => "${cur_target_dir}/topN_bowtie1_genome_cutadapt_1mm_notidentical",
      option        => $bowtie1_option_1mm,
      source_ref    => "cutadapt_len_fake",
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

    #1 mismatch notidentical search
    bowtie1_genome_cutadapt_topN_1mm => {
      class         => "Bowtie1",
      perform       => 0,
      target_dir    => "${cur_target_dir}/topN_bowtie1_genome_cutadapt_1mm",
      option        => $bowtie1_option_1mm,
      source_ref    => [ "identical", ".fastq\$" ],
      bowtie1_index => $def->{bowtie1_index},
      samonly       => 0,
      sh_direct     => 1,
      pbs           => {
        "email"    => $email,
        "nodes"    => "1:ppn=8",
        "walltime" => "72",
        "mem"      => "40gb"
      },
    },
    mirna_1mm_count => {
      class           => "MirnaCount",
      perform         => 0,
      target_dir      => "${cur_target_dir}/topN_bowtie1_genome_cutadapt_1mm_count_miRNA",
      option          => $mirnacount_option,
      source_ref      => "bowtie1_genome_cutadapt_topN_1mm",
      fastq_files_ref => "identical",
      seqcount_ref    => [ "identical", ".dupcount\$" ],
      cqs_tools       => $cqstools,
      gff_file        => $def->{coordinate},
      fasta_file      => $fasta_file,
      samtools        => $samtools,
      sh_direct       => 1,
      pbs             => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "72",
        "mem"      => "40gb"
      },
    },
    miRNA_1mm_table => {
      class      => "CQSMirnaTable",
      perform    => 0,
      target_dir => "${target_dir}/summary_miRNA",
      option     => "-o " . $def->{task_name} . "_miRNA.count",
      source_ref => "mirna_1mm_count",
      cqs_tools  => $cqstools,
      sh_direct  => 1,
      pbs        => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "10",
        "mem"      => "10gb"
      },
    },
    miRNA_1mm_count_overlap => {
      class           => "CQSMappedCount",
      perform         => 1,
      target_dir      => "${cur_target_dir}/topN_bowtie1_genome_cutadapt_1mm_count_miRNA_overlap",
      option          => $mirna_overlap_count_option,
      source_ref      => "bowtie1_genome_cutadapt_topN_1mm",
      fastq_files_ref => "identical",
      seqcount_ref    => [ "identical", ".dupcount\$" ],
      cqs_tools       => $cqstools,
      gff_file        => $def->{coordinate},
      samtools        => $samtools,
      sh_direct       => 1,
      pbs             => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "72",
        "mem"      => "40gb"
      },
    },
    miRNA_1mm_overlap_table => {
      class      => "CQSMappedTable",
      perform    => 0,
      target_dir => "${target_dir}/summary_miRNA_overlap",
      option     => "-i 1 -v 2 -o " . $def->{task_name} . "_miRNA_overlap.count",
      source_ref => "miRNA_1mm_count_overlap",
      cqs_tools  => $cqstools,
      sh_direct  => 1,
      pbs        => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "10",
        "mem"      => "10gb"
      },
    },
    miRNA_1mm_overlap_position => {
      class      => "CQSMappedPosition",
      perform    => 0,
      target_dir => "${cur_target_dir}/topN_bowtie1_genome_cutadapt_1mm_count_miRNA_overlap_position",
      option     => "-o " . $def->{task_name} . "_miRNA.position",
      source_ref => "miRNA_1mm_count_overlap",
      cqs_tools  => $cqstools,
      sh_direct  => 1,
      pbs        => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "10",
        "mem"      => "10gb"
      },
    },
    tRNA_1mm_count => {
      class           => "CQSMappedCount",
      perform         => 0,
      target_dir      => "${cur_target_dir}/topN_bowtie1_genome_cutadapt_1mm_count_tRNA",
      option          => $trnacount_option,
      source_ref      => "bowtie1_genome_cutadapt_topN_1mm",
      fastq_files_ref => "identical",
      seqcount_ref    => [ "identical", ".dupcount\$" ],
      cqs_tools       => $cqstools,
      gff_file        => $def->{trna_coordinate},
      samtools        => $samtools,
      sh_direct       => 1,
      pbs             => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "72",
        "mem"      => "40gb"
      },
    },
    tRNA_1mm_table => {
      class      => "CQSMappedTable",
      perform    => 0,
      target_dir => "${target_dir}/summary_tRNA",
      option     => "-i 1 -v 2 -o " . $def->{task_name} . "_tRNA.count",
      source_ref => "tRNA_1mm_count",
      cqs_tools  => $cqstools,
      sh_direct  => 1,
      pbs        => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "10",
        "mem"      => "10gb"
      },
    },
    tRNA_1mm_position => {
      class      => "CQSMappedPosition",
      perform    => 0,
      target_dir => "${cur_target_dir}/topN_bowtie1_genome_cutadapt_1mm_count_tRNA_position",
      option     => "-o " . $def->{task_name} . "_tRNA.position",
      source_ref => "tRNA_1mm_count",
      cqs_tools  => $cqstools,
      sh_direct  => 1,
      pbs        => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "10",
        "mem"      => "10gb"
      },
    },

    #    cqs_pileup_bowtie1_genome_cutadapt_topN_1mm_miRNA => {
    #      class        => "CQSPileup",
    #      perform      => 0,
    #      target_dir   => "${cur_target_dir}/topN_bowtie1_genome_cutadapt_1mm_cqspileup_miRNA",
    #      option       => "--export_igv",
    #      source_ref   => "bowtie1_genome_cutadapt_topN_1mm",
    #      seqcount_ref => [ "identical", ".dupcount\$" ],
    #      cqs_tools    => $cqstools,
    #      gff_file     => $def->{coordinate},
    #      samtools     => $samtools,
    #      sh_direct    => 1,
    #      pbs          => {
    #        "email"    => $email,
    #        "nodes"    => "1:ppn=1",
    #        "walltime" => "72",
    #        "mem"      => "40gb"
    #      },
    #    },
    #    cqs_pileup_bowtie1_genome_cutadapt_topN_1mm_tRNA => {
    #      class        => "CQSPileup",
    #      perform      => 0,
    #      target_dir   => "${cur_target_dir}/topN_bowtie1_genome_cutadapt_1mm_cqspileup_tRNA",
    #      option       => "--export_igv",
    #      source_ref   => "bowtie1_genome_cutadapt_topN_1mm",
    #      seqcount_ref => [ "identical", ".dupcount\$" ],
    #      cqs_tools    => $cqstools,
    #      gff_file     => $def->{trna_coordinate},
    #      samtools     => $samtools,
    #      sh_direct    => 1,
    #      pbs          => {
    #        "email"    => $email,
    #        "nodes"    => "1:ppn=1",
    #        "walltime" => "72",
    #        "mem"      => "40gb"
    #      },
    #    },
    #    bowtie2_mirna_count_bowtie1_genome_cutadapt_topN => {
    #      class         => "Bowtie2",
    #      perform       => 0,
    #      target_dir    => "${cur_target_dir}/topN_bowtie1_genome_cutadapt_count_bowtie2",
    #      option        => $bowtie2_local_option,
    #      source_ref    => [ "mirna_count_bowtie1_genome_cutadapt_topN", ".fastq\$" ],
    #      bowtie2_index => $def->{bowtie2_index},
    #      sh_direct     => 1,
    #      pbs           => {
    #        "email"    => $email,
    #        "nodes"    => "1:ppn=8",
    #        "walltime" => "24",
    #        "mem"      => "20gb"
    #      },
    #    },
    #    shrimp2_miRBase_bowtie1_genome_cutadapt_topN => {
    #      class         => "Shrimp2",
    #      perform       => 0,
    #      target_dir    => "${cur_target_dir}/topN_bowtie1_genome_cutadapt_shrimp2_miRBase",
    #      option        => $shrimp2_option,
    #      source_ref    => [ "mirna_count_bowtie1_genome_cutadapt_topN", ".fastq\$" ],
    #      shrimp2_index => $def->{shrimp2_index},
    #      is_mirna      => 1,
    #      output_bam    => 1,
    #      sh_direct     => 0,
    #      pbs           => {
    #        "email"    => $email,
    #        "nodes"    => "1:ppn=8",
    #        "walltime" => "720",
    #        "mem"      => "40gb"
    #      },
    #    },
    #    bowtie1_genome_cutadapt_topN_pm => {
    #      class         => "Bowtie1",
    #      perform       => 0,
    #      target_dir    => "${cur_target_dir}/topN_bowtie1_genome_cutadapt_pm",
    #      option        => $bowtie1_option_pm,
    #      source_ref    => [ "identical", ".fastq\$" ],
    #      bowtie1_index => $def->{bowtie1_index},
    #      samonly       => 0,
    #      sh_direct     => 1,
    #      pbs           => {
    #        "email"    => $email,
    #        "nodes"    => "1:ppn=8",
    #        "walltime" => "72",
    #        "mem"      => "40gb"
    #      },
    #    },
    #    mirna_count_bowtie1_genome_cutadapt_topN_pm => {
    #      class           => "MirnaCount",
    #      perform         => 0,
    #      target_dir      => "${cur_target_dir}/topN_bowtie1_genome_cutadapt_count_pm",
    #      option          => $mirnacount_option,
    #      source_ref      => "bowtie1_genome_cutadapt_topN_pm",
    #      fastq_files_ref => "identical",
    #      seqcount_ref    => [ "identical", ".dupcount\$" ],
    #      cqs_tools       => $cqstools,
    #      gff_file        => $def->{coordinate},
    #      samtools        => $samtools,
    #      sh_direct       => 1,
    #      pbs             => {
    #        "email"    => $email,
    #        "nodes"    => "1:ppn=1",
    #        "walltime" => "72",
    #        "mem"      => "40gb"
    #      },
    #    },
    #    bowtie1_genome_cutadapt_topN_pm_unmatched => {
    #      class         => "Bowtie1",
    #      perform       => 0,
    #      target_dir    => "${cur_target_dir}/topN_bowtie1_genome_cutadapt_pm_unmatched",
    #      option        => $bowtie1_option_pm,
    #      source_ref    => [ "mirna_count_bowtie1_genome_cutadapt_topN_pm", ".fastq\$" ],
    #      bowtie1_index => $human->{bowtie1_index},
    #      samonly       => 0,
    #      sh_direct     => 1,
    #      pbs           => {
    #        "email"    => $email,
    #        "nodes"    => "1:ppn=8",
    #        "walltime" => "72",
    #        "mem"      => "40gb"
    #      },
    #    },
    #    mirna_count_bowtie1_genome_cutadapt_topN_pm_unmatched => {
    #      class           => "MirnaCount",
    #      perform         => 0,
    #      target_dir      => "${cur_target_dir}/topN_bowtie1_genome_cutadapt_count_pm_unmatched",
    #      option          => $mirnacount_option,
    #      source_ref      => "bowtie1_genome_cutadapt_topN_pm_unmatched",
    #      fastq_files_ref => "identical",
    #      seqcount_ref    => [ "identical", ".dupcount\$" ],
    #      cqs_tools       => $cqstools,
    #      gff_file        => $human->{coordinate},
    #      samtools        => $samtools,
    #      sh_direct       => 1,
    #      pbs             => {
    #        "email"    => $email,
    #        "nodes"    => "1:ppn=1",
    #        "walltime" => "72",
    #        "mem"      => "40gb"
    #      },
    #    },
    #    bwa_genome_cutadapt_topN => {
    #      class      => "BWA",
    #      perform    => 0,
    #      target_dir => "${cur_target_dir}/topN_bwa_genome_cutadapt",
    #      option     => $bwa_option,
    #      source_ref => [ "identical", ".fastq\$" ],
    #      fasta_file => $def->{bwa_fasta},
    #      samonly    => 0,
    #      sh_direct  => 1,
    #      pbs        => {
    #        "email"    => $email,
    #        "nodes"    => "1:ppn=8",
    #        "walltime" => "72",
    #        "mem"      => "40gb"
    #      },
    #    },
    #    bwa_genome_cutadapt_softclip_topN => {
    #      class      => "BWA",
    #      perform    => 0,
    #      target_dir => "${cur_target_dir}/topN_bwa_genome_softclip_cutadapt",
    #      option     => $bwa_clip_option,
    #      source_ref => [ "identical", ".fastq\$" ],
    #      fasta_file => $def->{bwa_fasta},
    #      samonly    => 0,
    #      sh_direct  => 1,
    #      pbs        => {
    #        "email"    => $email,
    #        "nodes"    => "1:ppn=8",
    #        "walltime" => "72",
    #        "mem"      => "40gb"
    #      },
    #    },
    #    mirna_count_bwa_genome_cutadapt_topN => {
    #      class           => "MirnaCount",
    #      perform         => 0,
    #      target_dir      => "${cur_target_dir}/topN_bwa_genome_cutadapt_count",
    #      option          => $mirnacount_option . " -e 3",
    #      source_ref      => "bwa_genome_cutadapt_topN",
    #      fastq_files_ref => "identical",
    #      seqcount_ref    => [ "identical", ".dupcount\$" ],
    #      cqs_tools       => $cqstools,
    #      gff_file        => $hsammu_gffs,
    #      samtools        => $samtools,
    #      sh_direct       => 1,
    #      pbs             => {
    #        "email"    => $email,
    #        "nodes"    => "1:ppn=1",
    #        "walltime" => "72",
    #        "mem"      => "40gb"
    #      },
    #    },
    #    bowtie1_genome_cutadapt_topN_3mm => {
    #      class         => "Bowtie1",
    #      perform       => 0,
    #      target_dir    => "${cur_target_dir}/topN_bowtie1_genome_cutadapt_3mm",
    #      option        => $bowtie1_option_3mm,
    #      source_ref    => [ "identical", ".fastq\$" ],
    #      bowtie1_index => $def->{bowtie1_index},
    #      samonly       => 0,
    #      sh_direct     => 1,
    #      pbs           => {
    #        "email"    => $email,
    #        "nodes"    => "1:ppn=8",
    #        "walltime" => "72",
    #        "mem"      => "40gb"
    #      },
    #    },

  };

  performConfig($config);
#  performTrace($config);
}

1;

