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

my $root     = create_directory_or_die("/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2");
my $cqstools = "/home/shengq1/cqstools/CQS.Tools.exe";

my $hg19_mrna_gff     = "/data/cqs/shengq1/reference/miRBase20/hsa.gff3";
my $hg19_trna_bed     = "/data/cqs/guoy1/reference/smallrna/hg19_tRNA_ucsc_ensembl.bed";
my $hg19_trna_fasta   = "/data/cqs/guoy1/reference/smallrna/hg19_tRNA_ucsc_ensembl.bed.fa";
my $hg19_smallrna_bed = "/data/cqs/guoy1/reference/smallrna/hg19_smallRNA_ucsc_ensembl.bed";

my $mm10_mrna_gff     = "/data/cqs/shengq1/reference/miRBase20/mmu.gff3";
my $mm10_trna_bed     = "/data/cqs/guoy1/reference/smallrna/mm10_tRNA_ucsc_ensembl.bed";
my $mm10_trna_fasta   = "/data/cqs/guoy1/reference/smallrna/mm10_tRNA_ucsc_ensembl.bed.fa";
my $mm10_smallrna_bed = "/data/cqs/guoy1/reference/smallrna/mm10_smallRNA_ucsc_ensembl.bed";

my $rn4_mrna_gff     = "/data/cqs/shengq1/reference/miRBase20/rno.gff3";
my $rn4_trna_bed     = "/data/cqs/guoy1/reference/smallrna/rn4_tRNA_ucsc_ensembl.bed";
my $rn4_trna_fasta   = "/data/cqs/guoy1/reference/smallrna/rn4_tRNA_ucsc_ensembl.bed.fa";
my $rn4_smallrna_bed = "/data/cqs/guoy1/reference/smallrna/rn4_smallRNA_ucsc_ensembl.bed";

my $target_dir = $root;

my $target_rat_dir   = create_directory_or_die( $target_dir . "/rat" );
my $target_human_dir = create_directory_or_die( $target_dir . "/human" );
my $target_mouse_dir = create_directory_or_die( $target_dir . "/mouse" );

my $email     = "quanhu.sheng\@vanderbilt.edu";
my $task_name = "VANGARD00055";

my $samtools = "/home/shengq1/local/bin/samtools/samtools";

my $bowtie1_option_pm       = "-a -m 100 --best --strata -v 0 -l 12 -p 8";
my $bowtie1_option_1mm      = "-a -m 100 --best --strata -v 1 -l 12 -p 8";
my $bowtie1_option_1mm_trim = "-a -m 100 --best --strata -v 1 -l 12 -p 8 --trim3 3";
my $bowtie1_option_3mm      = "-a -m 100 --best --strata -v 3 -l 12 -p 8";

my $bowtie1_rat_index   = "/data/cqs/shengq1/reference/rn4/bowtie1_index/rn4";
my $bowtie1_human_index = "/data/cqs/guoy1/reference/hg19/bowtie_index_hg19_rCRS_1.0.0/hg19_rCRS";
my $bowtie1_mouse_index = "/data/cqs/shengq1/reference/mm10/bowtie_index/mm10";

my $mirnacount_option          = "-s";                                                    #ignore score
my $trnacount_option           = "--length --sequence";
my $mirna_overlap_count_option = "-s --gtf_key miRNA";
my $mirna_fasta                = "/data/cqs/shengq1/reference/miRBase20/mature.dna.fa";

my $rat = {
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
    "2516-01"                   => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/2516-KCV-1_1.fastq.gz"],
    "2516-02"                   => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/2516-KCV-2_1.fastq.gz"],
    "2516-03"                   => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/2516-KCV-3_1.fastq.gz"],
    "2516-04"                   => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/2516-KCV-4_1.fastq.gz"],
    "2516-05"                   => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/2516-KCV-5_1.fastq.gz"],
    "2516-06"                   => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/2516-KCV-6_1.fastq.gz"],
    "2516-07"                   => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/2516-KCV-7_1.fastq.gz"],
    "2516-08"                   => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/2516-KCV-8_1.fastq.gz"],
    "2516-09"                   => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/2516-KCV-9_1.fastq.gz"],
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
  coordinate          => $rn4_mrna_gff,
  trna_coordinate     => $rn4_trna_bed,
  trna_fasta          => $rn4_trna_fasta,
  smallrna_coordinate => $rn4_smallrna_bed,
  bowtie1_index       => $bowtie1_rat_index,
  target_dir          => $target_rat_dir,
  task_name           => $task_name . "_rat",
  groups              => {
    "Islets_diabetes" => [ "2516-01", "2516-02", "2516-03", "2516-04", "2516-05", "2516-06", "2516-07", "2516-08", "2516-09" ],
    "Liver_diabetes"  => [
      "2570-01", "2570-02", "2570-03", "2570-04", "2570-05", "2570-06", "2570-07", "2570-08", "2570-09", "2570-10",
      "2570-11", "2570-12", "2570-13", "2570-14", "2570-15", "2570-16", "2570-17", "2570-18"
    ],
    "HDL_diabetes" => [
      "Colesevelam-3251.1_CGATGT", "Colesevelam-3251.2_ATCACG", "Colesevelam-3258.1_ACAGTG", "Colesevelam-3258.2_TGACCA", "Colesevelam-3269.1_GCCAAT", "Colesevelam-3269.2_TTAGGC",
      "Lean-3183.1_GATCAG",        "Lean-3185.1_CAGATC",        "Lean-3185.2_ACTTGA",        "Lean-KR-4073.1_TAGCTT",     "Lean-KR-4073.2_CTTGTA",     "Lean-KR-4074.2_AGTCAA",
      "Lean-KR-4074.3_GGCTAC",     "Vehicle-3259.1_CCGTCC",     "Vehicle-3259.2_ATGTCA",     "Vehicle-3263.1_GTGAAA",     "Vehicle-3263.2_AGTTCC",     "Vehicle-3271.1_GTCCGC",
      "Vehicle-3271.2_GTAGAG"
    ],
    "NIH_HDL" => [
      "NIH_Rat_HDL_01", "NIH_Rat_HDL_02", "NIH_Rat_HDL_03", "NIH_Rat_HDL_04", "NIH_Rat_HDL_05", "NIH_Rat_HDL_06",
      "NIH_Rat_HDL_07", "NIH_Rat_HDL_08", "NIH_Rat_HDL_09", "NIH_Rat_HDL_10", "NIH_Rat_HDL_11"
    ]
  }
};

my $transplant = [ "2570-KCV-01-19", "2570-KCV-01-20", "2570-KCV-01-21", "2570-KCV-01-22", "2570-KCV-01-23", "2570-KCV-01-24", "2570-KCV-01-25", "2570-KCV-01-26", "2570-KCV-01-27" ];

my $human = {
  source => {
    "2516-10"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/2516-KCV-10_1.fastq.gz"],
    "2516-11"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/2516-KCV-11_1.fastq.gz"],
    "2516-12"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/2516-KCV-12_1.fastq.gz"],
    "2516-13"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/2516-KCV-13_1.fastq.gz"],
    "2516-14"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/2516-KCV-14_1.fastq.gz"],
    "2516-15"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/2516-KCV-15_1.fastq.gz"],
    "2516-16"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/2516-KCV-16_1.fastq.gz"],
    "2516-17"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/2516-KCV-17_1.fastq.gz"],
    "2516-18"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/2516-KCV-18_1.fastq.gz"],
    "2516-19"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/2516-KCV-19_1.fastq.gz"],
    "2516-20"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/2516-KCV-20_1.fastq.gz"],
    "2516-21"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/2516-KCV-21_1.fastq.gz"],
    "2516-22"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/2516-KCV-22_1.fastq.gz"],
    "2516-23"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/2516-KCV-23_1.fastq.gz"],
    "2516-24"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/2516-KCV-24_1.fastq.gz"],
    "2516-25"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/2516-KCV-25_1.fastq.gz"],
    "2516-26"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/2516-KCV-26_1.fastq.gz"],
    "2516-27"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/2516-KCV-27_1.fastq.gz"],
    "2516-28"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/2516-KCV-28_1.fastq.gz"],
    "2516-29"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/2516-KCV-29_1.fastq.gz"],
    "2516-30"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/2516-KCV-30_1.fastq.gz"],
    "2516-31"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/2516-KCV-31_1.fastq.gz"],
    "2516-32"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/2516-KCV-32_1.fastq.gz"],
    "2516-33"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/2516-KCV-33_1.fastq.gz"],
    "2516-34"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/2516-KCV-34_1.fastq.gz"],
    "2516-35"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/2516-KCV-35_1.fastq.gz"],
    "2516-36"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/2516-KCV-36_1.fastq.gz"],
    "2516-37"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/2516-KCV-37_1.fastq.gz"],
    "2516-38"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/2516-KCV-38_1.fastq.gz"],
    "2516-39"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/2516-KCV-39_1.fastq.gz"],
    "2516-40"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/2516-KCV-40_1.fastq.gz"],
    "2516-41"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/2516-KCV-41_1.fastq.gz"],
    "2516-42"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/2516-KCV-42_1.fastq.gz"],
    "2516-43"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/2516-KCV-43_1.fastq.gz"],
    "2516-44"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/2516-KCV-44_1.fastq.gz"],
    "2516-45"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/2516-KCV-45_1.fastq.gz"],
    "2516-46"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/2516-KCV-46_1.fastq.gz"],
    "2516-47"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/2516-KCV-47_1.fastq.gz"],
    "2516-48"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/2516-KCV-48_1.fastq.gz"],
    "KCV2_1N2" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/KCV2_1N2_GCCAAT_L003_R1_001.fastq.gz"],
    "KCV2_1N3" => [
      "/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/KCV2_1N3_CAGATC_L003_R1_001.fastq",
      "/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/KCV2_1N3_CAGATC_L003_R1_002.fastq"
    ],
    "KCV2_1N4" => [
      "/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/KCV2_1N4_ACTTGA_L003_R1_001.fastq",
      "/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/KCV2_1N4_ACTTGA_L003_R1_002.fastq"
    ],
    "KCV2_1N5" => [
      "/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/KCV2_1N5_GATCAG_L003_R1_001.fastq",
      "/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/KCV2_1N5_GATCAG_L003_R1_002.fastq"
    ],
    "KCV2_1N6" => [
      "/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/KCV2_1N6_TAGCTT_L003_R1_001.fastq",
      "/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/KCV2_1N6_TAGCTT_L003_R1_002.fastq"
    ],
    "KCV2_2N1" => [
      "/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/KCV2_2N1_GGCTAC_L003_R1_001.fastq",
      "/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/KCV2_2N1_GGCTAC_L003_R1_002.fastq",
      "/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/KCV2_2N1_GGCTAC_L003_R1_003.fastq",
      "/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/KCV2_2N1_GGCTAC_L003_R1_004.fastq",
      "/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/KCV2_2N1_GGCTAC_L003_R1_005.fastq"
    ],
    "KCV2_2N2" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/KCV2_2N2_GCCGCG_L003_R1_001.fastq.gz"],
    "KCV2_2N3" => [
      "/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/KCV2_2N3_CTTGTA_L003_R1_001.fastq",
      "/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/KCV2_2N3_CTTGTA_L003_R1_002.fastq",
      "/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/KCV2_2N3_CTTGTA_L003_R1_003.fastq",
      "/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/KCV2_2N3_CTTGTA_L003_R1_004.fastq",
      "/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/KCV2_2N3_CTTGTA_L003_R1_005.fastq",
      "/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/KCV2_2N3_CTTGTA_L003_R1_006.fastq",
      "/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/KCV2_2N3_CTTGTA_L003_R1_007.fastq",
      "/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/KCV2_2N3_CTTGTA_L003_R1_008.fastq",
      "/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/KCV2_2N3_CTTGTA_L003_R1_009.fastq",
      "/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/KCV2_2N3_CTTGTA_L003_R1_010.fastq"
    ],
    "KCV2_2N4"           => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/KCV2_2N4_GCCTTA_L003_R1_001.fastq.gz"],
    "KCV2_2N5"           => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/KCV2_2N5_GCTCCA_L003_R1_001.fastq.gz"],
    "KCV3_1C2"           => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/KCV3_1C2_GGCACA_L004_R1_001.fastq.gz"],
    "KCV3_1C3"           => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/KCV3_1C3_GGCCTG_L004_R1_001.fastq.gz"],
    "KCV3_1C4"           => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/KCV3_1C4_TCTACC_L004_R1_001.fastq.gz"],
    "KCV3_1C5"           => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/KCV3_1C5_TGAAGT_L004_R1_001.fastq.gz"],
    "KCV3_1C6"           => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/KCV3_1C6_TGCCAT_L004_R1_001.fastq.gz"],
    "KCV3_2C1"           => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/KCV3_2C1_TGCTGG_L004_R1_001.fastq.gz"],
    "KCV3_2C2"           => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/KCV3_2C2_TGGCGC_L004_R1_001.fastq.gz"],
    "KCV3_2C3"           => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/KCV3_2C3_TTCGAA_L004_R1_001.fastq.gz"],
    "Sample1"            => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/Sample1_12.fastq.gz"],
    "Sample2"            => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/Sample2_12.fastq.gz"],
    "Sample3"            => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/Sample3_12.fastq.gz"],
    "Sample4"            => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/Sample4_12.fastq.gz"],
    "Sample5"            => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/Sample5_12.fastq.gz"],
    "2570-KCV-01-19"     => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/VickersTemp/mouseLiverTransplant1_GTGAAA_L003_R1_001.fastq.gz"],
    "2570-KCV-01-20"     => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/VickersTemp/mouseLiverTransplant2_GTGGCC_L003_R1_001.fastq.gz"],
    "2570-KCV-01-21"     => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/VickersTemp/mouseLiverTransplant3_GTTTCG_L003_R1_001.fastq.gz"],
    "2570-KCV-01-22"     => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/VickersTemp/mouseLiverTransplant4_CGTACG_L003_R1_001.fastq.gz"],
    "2570-KCV-01-23"     => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/VickersTemp/mouseLiverControl1_GAGTGG_L003_R1_001.fastq.gz"],
    "2570-KCV-01-24"     => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/VickersTemp/mouseLiverControl2_GGTAGC_L003_R1_001.fastq.gz"],
    "2570-KCV-01-25"     => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/VickersTemp/mouseLiverControl3_ACTGAT_L003_R1_001.fastq.gz"],
    "2570-KCV-01-26"     => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/VickersTemp/mouseLiverControl4_ATGAGC_L003_R1_001.fastq.gz"],
    "2570-KCV-01-27"     => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/VickersTemp/mouseLiverControl5_ATTCCT_L003_R1_001.fastq.gz"],
    "2571-KCV-1-31"      => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/VickersTemp/HD01_20_CACGAT_L005_R1_001.fastq.gz"],
    "2571-KCV-1-32"      => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/VickersTemp/HD01_21_CACTCA_L005_R1_001.fastq.gz"],
    "2571-KCV-1-33"      => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/VickersTemp/HD01_22_CAGGCG_L005_R1_001.fastq.gz"],
    "2571-KCV-1-34"      => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/VickersTemp/HD01_23_CATGGC_L005_R1_001.fastq.gz"],
    "2571-KCV-1-35"      => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/VickersTemp/HD01_24_CATTTT_L005_R1_001.fastq.gz"],
    "2571-KCV-1-36"      => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/VickersTemp/HD01_25_CCAACA_L005_R1_001.fastq.gz"],
    "2571-KCV-1-37"      => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/VickersTemp/HD01_26_CGGAAT_L005_R1_001.fastq.gz"],
    "2571-KCV-1-38"      => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/VickersTemp/HD01_27_CTAGCT_L005_R1_001.fastq.gz"],
    "2571-KCV-1-39"      => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/VickersTemp/HD01_28_CTATAC_L005_R1_001.fastq.gz"],
    "2571-KCV-1-40"      => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/VickersTemp/NC02_20_CTCAGA_L005_R1_001.fastq.gz"],
    "2571-KCV-1-41"      => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/VickersTemp/NC02_21_GACGAC_L005_R1_001.fastq.gz"],
    "2571-KCV-1-42"      => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/VickersTemp/NC02_22_TAATCG_L005_R1_001.fastq.gz"],
    "2571-KCV-1-43"      => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/VickersTemp/NC02_23_TACAGC_L005_R1_001.fastq.gz"],
    "2571-KCV-1-44"      => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/VickersTemp/NC02_24_TATAAT_L005_R1_001.fastq.gz"],
    "2571-KCV-1-45"      => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/VickersTemp/NC02_25_TCATTC_L005_R1_001.fastq.gz"],
    "2571-KCV-1-46"      => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/VickersTemp/NC02_26_TCCCGA_L005_R1_001.fastq.gz"],
    "2571-KCV-1-47"      => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/VickersTemp/NC02_27_TCGAAG_L005_R1_001.fastq.gz"],
    "2571-KCV-1-48"      => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/VickersTemp/NC02_28_TCGGCA_L005_R1_001.fastq.gz"],
    "2572-KCV-1-19"      => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/VickersTemp/CH003_GTGAAA_L002_R1_001.fastq.gz"],
    "2572-KCV-1-20"      => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/VickersTemp/CH002SS_GTGGCC_L002_R1_001.fastq.gz"],
    "2572-KCV-1-21"      => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/VickersTemp/CH001_GTTTCG_L002_R1_001.fastq.gz"],
    "2572-KCV-1-22"      => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/VickersTemp/PO14s_CGTACG_L002_R1_001.fastq.gz"],
    "2572-KCV-1-23"      => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/VickersTemp/PO16D_GGTAGC_L002_R1_001.fastq.gz"],
    "2572-KCV-1-24"      => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/VickersTemp/PO23F_GAGTGG_L002_R1_001.fastq.gz"],
    "2572-KCV-1-25"      => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/VickersTemp/VK_ACTGAT_L002_R1_001.fastq.gz"],
    "2572-KCV-1-26"      => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/VickersTemp/CC_ATGAGC_L002_R1_001.fastq.gz"],
    "2572-KCV-1-27"      => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/VickersTemp/AE_ATTCCT_L002_R1_001.fastq.gz"],
    "2572-KCV-1-28"      => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/VickersTemp/CHL001_CAAAAG_L002_R1_001.fastq.gz"],
    "2572-KCV-1-29"      => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/VickersTemp/CHL002_CAACTA_L002_R1_001.fastq.gz"],
    "2572-KCV-1-30"      => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/VickersTemp/CHL003_CACCGG_L002_R1_001.fastq.gz"],
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
  coordinate          => $hg19_mrna_gff,
  trna_coordinate     => $hg19_trna_bed,
  trna_fasta          => $hg19_trna_fasta,
  smallrna_coordinate => $hg19_smallrna_bed,
  bowtie1_index       => $bowtie1_human_index,
  target_dir          => $target_human_dir,
  task_name           => $task_name . "_human",
  groups              => {
    "Human_HDL_N3" =>
      [ "2516-10", "2516-11", "2516-12", "2516-13", "2516-14", "2516-15", "2516-16", "2516-17", "2516-18", "2516-19", "2516-20", "2516-21", "2516-22", "2516-23", "2516-24", "2516-47", "2516-48" ],
    "CKD_HDL_1" => [
      "2516-25", "2516-26", "2516-27", "2516-28", "2516-29", "2516-30", "2516-31", "2516-32", "2516-33", "2516-34", "2516-35", "2516-36",
      "2516-37", "2516-38", "2516-39", "2516-40", "2516-41", "2516-42", "2516-43", "2516-44", "2516-45", "2516-46"
    ],
    "TransplantLiver" => $transplant,
    "CKD_HDL_2"       => [
      "2571-KCV-1-31", "2571-KCV-1-32", "2571-KCV-1-33", "2571-KCV-1-34", "2571-KCV-1-35", "2571-KCV-1-36", "2571-KCV-1-37", "2571-KCV-1-38", "2571-KCV-1-39", "2571-KCV-1-40",
      "2571-KCV-1-41", "2571-KCV-1-42", "2571-KCV-1-43", "2571-KCV-1-44", "2571-KCV-1-45", "2571-KCV-1-46", "2571-KCV-1-47", "2571-KCV-1-48"
    ],
    "RA"   => [ "2572-KCV-1-19", "2572-KCV-1-20", "2572-KCV-1-21", "2572-KCV-1-22", "2572-KCV-1-23", "2572-KCV-1-24" ],
    "SLE"  => [ "2572-KCV-1-25", "2572-KCV-1-26", "2572-KCV-1-27", "2572-KCV-1-28", "2572-KCV-1-29", "2572-KCV-1-30" ],
    "Nash" => [
      "KCV2_1N2", "KCV2_1N3", "KCV2_1N4", "KCV2_1N5", "KCV2_1N6", "KCV2_2N1", "KCV2_2N2", "KCV2_2N3", "KCV2_2N4", "KCV2_2N5",
      "KCV3_1C2", "KCV3_1C3", "KCV3_1C4", "KCV3_1C5", "KCV3_1C6", "KCV3_2C1", "KCV3_2C2", "KCV3_2C3"
    ],
    "Plasma_fractions" => [ "Sample1", "Sample2", "Sample3", "Sample4", "Sample5", ],
    "Human_HDL_N10"    => [
      "01-018-Post_CTTGTA", "01-018-Pre_GGCTAC", "01-031-Post_GGTAGC", "01-031-Pre_GAGTGG", "01-061-Post_ATGAGC", "01-061-Pre_ACTGAT", "01-28-Post_CGATGT",  "01-28-Pre_ATCACG",
      "01-29-Pre_TGACCA",   "01-29-Pre_TTAGGC",  "01-36-Post_GCCAAT",  "01-36-Pre_ACAGTG",  "03-007-Post_CAAAAG", "03-007-Pre_ATTCCT", "03-011-Post_CACCGG", "03-011-Pre_CAACTA",
      "03-015-Post_CACTCA", "03-015-Pre_CACGAT", "03-018-Post_CGTACG", "03-018-Pre_GTTTCG", "03-026-Post_AGTTCC", "03-026-Pre_AGTCAA", "03-031-Post_CATGGC", "03-031-Pre_CAGGCG",
      "03-033-Post_CCAACA", "03-033-Pre_CATTTT", "03-036-Post_CCGTCC", "03-036-Pre_ATGTCA", "03-047-Post_CTAGCT", "03-047-Pre_CGGAAT", "03-049-Post_CTCAGA", "03-049-Pre_CTATAC",
      "03-063-Post_GTGGCC", "03-063-Pre_GTGAAA", "03-065-Post_GTCCGC", "03-065-Pre_GTAGAG", "03-16-Post_ACTTGA",  "03-16-Pre_CAGATC",  "03-17-Post_TAGCTT",  "03-17-Pre_GATCAG"
    ],
  }
};

my $mouse = {
  source => {
    "2570-KCV-01-19" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/VickersTemp/mouseLiverTransplant1_GTGAAA_L003_R1_001.fastq.gz"],
    "2570-KCV-01-20" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/VickersTemp/mouseLiverTransplant2_GTGGCC_L003_R1_001.fastq.gz"],
    "2570-KCV-01-21" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/VickersTemp/mouseLiverTransplant3_GTTTCG_L003_R1_001.fastq.gz"],
    "2570-KCV-01-22" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/VickersTemp/mouseLiverTransplant4_CGTACG_L003_R1_001.fastq.gz"],
    "2570-KCV-01-23" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/VickersTemp/mouseLiverControl1_GAGTGG_L003_R1_001.fastq.gz"],
    "2570-KCV-01-24" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/VickersTemp/mouseLiverControl2_GGTAGC_L003_R1_001.fastq.gz"],
    "2570-KCV-01-25" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/VickersTemp/mouseLiverControl3_ACTGAT_L003_R1_001.fastq.gz"],
    "2570-KCV-01-26" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/VickersTemp/mouseLiverControl4_ATGAGC_L003_R1_001.fastq.gz"],
    "2570-KCV-01-27" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/VickersTemp/mouseLiverControl5_ATTCCT_L003_R1_001.fastq.gz"],
  },
  coordinate          => $mm10_mrna_gff,
  trna_coordinate     => $mm10_trna_bed,
  trna_fasta          => $mm10_trna_fasta,
  smallrna_coordinate => $mm10_smallrna_bed,
  bowtie1_index       => $bowtie1_mouse_index,
  target_dir          => $target_mouse_dir,
  task_name           => $task_name . "_mouse",
  groups              => { "TransplantLiver" => $transplant }
};

my @defs = ( $rat, $human, $mouse );

#my @defs = ($human);
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
    mirna_pm_count => {
      class           => "MirnaCount",
      perform         => 0,
      target_dir      => "${cur_target_dir}/topN_bowtie1_genome_cutadapt_pm_count_miRNA",
      option          => $mirnacount_option,
      source_ref      => "bowtie1_genome_cutadapt_topN_pm",
      fastq_files_ref => "identical",
      seqcount_ref    => [ "identical", ".dupcount\$" ],
      cqs_tools       => $cqstools,
      gff_file        => $def->{coordinate},
      fasta_file      => $mirna_fasta,
      samtools        => $samtools,
      sh_direct       => 1,
      pbs             => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "72",
        "mem"      => "40gb"
      },
    },
    miRNA_pm_table => {
      class      => "CQSMirnaTable",
      perform    => 1,
      target_dir => "${target_dir}/summary_pm",
      option     => "",
      source_ref => "mirna_pm_count",
      groups     => $def->{groups},
      cqs_tools  => $cqstools,
      sh_direct  => 1,
      prefix     => "miRNA_pm_",
      pbs        => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "10",
        "mem"      => "10gb"
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
      fasta_file      => $mirna_fasta,
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
      fasta_file      => $def->{trna_fasta},
      samtools        => $samtools,
      sh_direct       => 1,
      pbs             => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "72",
        "mem"      => "40gb"
      },
    },

    #for distinct count only
    fake_human_tRNA_overlap_count_bowtie1_genome_cutadapt_topN_pm => {
      class           => "CQSMappedCount",
      perform         => 0,
      target_dir      => $human->{target_dir} . "/topN_bowtie1_genome_cutadapt_pm_count_tRNA",
      option          => $trnacount_option,
      source_ref      => "bowtie1_genome_cutadapt_topN_pm",
      fastq_files_ref => "identical",
      seqcount_ref    => [ "identical", ".dupcount\$" ],
      cqs_tools       => $cqstools,
      gff_file        => $human->{trna_coordinate},
      fasta_file      => $human->{trna_fasta},
      samtools        => $samtools,
      sh_direct       => 1,
      pbs             => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "72",
        "mem"      => "40gb"
      },
    },
    tRNA_overlap_count_bowtie1_genome_cutadapt_topN_pm_distinct => {
      class         => "CQS::MappedDistinct",
      perform       => 0,
      target_dir    => "${target_dir}/pm_tRNA_distinct",
      option        => "",
      source_ref    => [ "tRNA_overlap_count_bowtie1_genome_cutadapt_topN_pm", ".mapped.xml" ],
      second_ref    => [ "fake_human_tRNA_overlap_count_bowtie1_genome_cutadapt_topN_pm", ".mapped.xml" ],
      cqs_tools     => $cqstools,
      first_suffix  => "mouse_",
      second_suffix => "human_",
      sh_direct     => 1,
      pbs           => {
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
    bowtie1_genome_cutadapt_topN_1mm_notidentical => {
      class         => "Bowtie1",
      perform       => 0,
      target_dir    => "${cur_target_dir}/topN_bowtie1_genome_cutadapt_1mm_notidentical",
      option        => $bowtie1_option_1mm,
      source_ref    => "cutadapt_len",
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

    #1 mismatch search
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
    bowtie1_genome_cutadapt_topN_1mm_trim => {
      class         => "Bowtie1",
      perform       => 0,
      target_dir    => "${cur_target_dir}/topN_bowtie1_genome_cutadapt_1mm_trim",
      option        => $bowtie1_option_1mm . " -5 2 -3 3 ",
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
      fasta_file      => $mirna_fasta,
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
      perform    => 1,
      target_dir => "${target_dir}/summary_1mm",
      option     => "",
      source_ref => "mirna_1mm_count",
      cqs_tools  => $cqstools,
      groups     => $def->{groups},
      prefix     => "miRNA_1mm_",
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
      perform         => 0,
      target_dir      => "${cur_target_dir}/topN_bowtie1_genome_cutadapt_1mm_count_miRNA_overlap",
      option          => $mirna_overlap_count_option,
      source_ref      => "bowtie1_genome_cutadapt_topN_1mm",
      fastq_files_ref => "identical",
      seqcount_ref    => [ "identical", ".dupcount\$" ],
      cqs_tools       => $cqstools,
      gff_file        => $def->{coordinate},
      fasta_file      => $mirna_fasta,
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
      perform    => 1,
      target_dir => "${target_dir}/summary_1mm_overlap",
      option     => "",
      source_ref => "miRNA_1mm_count_overlap",
      groups     => $def->{groups},
      cqs_tools  => $cqstools,
      prefix     => "miRNA_overlap_1mm_",
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
      fasta_file      => $def->{trna_fasta},
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
      perform    => 1,
      target_dir => "${target_dir}/summary",
      option     => "",
      source_ref => [ "tRNA_1mm_count", ".xml" ],
      groups     => $def->{groups},
      cqs_tools  => $cqstools,
      prefix     => "tRNA_1mm_",
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
    smallRNA_1mm_count => {
      class           => "CQSMappedCount",
      perform         => 0,
      target_dir      => "${cur_target_dir}/topN_bowtie1_genome_cutadapt_1mm_count_smallRNA",
      option          => $trnacount_option,
      source_ref      => "bowtie1_genome_cutadapt_topN_1mm",
      fastq_files_ref => "identical",
      seqcount_ref    => [ "identical", ".dupcount\$" ],
      cqs_tools       => $cqstools,
      gff_file        => $def->{smallrna_coordinate},
      samtools        => $samtools,
      sh_direct       => 1,
      pbs             => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "72",
        "mem"      => "40gb"
      },
    },
    smallRNA_1mm_category => {
      class           => "CQSSmallRNACategory",
      perform         => 0,
      target_dir      => "${cur_target_dir}/topN_bowtie1_genome_cutadapt_1mm_count_smallRNA_category",
      option          => "",
      source_ref      => [ "smallRNA_1mm_count", ".mapped.xml\$" ],
      mirna_count_ref => [ "mirna_1mm_count", ".mapped.xml\$" ],
      groups          => $def->{groups},
      cqs_tools       => $cqstools,
      sh_direct       => 1,
      pbs             => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "72",
        "mem"      => "40gb"
      },
    },
  };

  #  performConfig($config, "" , 1);
  performConfig($config);

  if ( $def == $mouse ) {
    #performTask( $config, "tRNA_overlap_count_bowtie1_genome_cutadapt_topN_pm_distinct" );
  }
}

my $w83_bowtie1_index = "/scratch/cqs/shengq1/references/gingivalis_W83/bowtie_1.0.0_index/Gingivalis_W83";
my $w83_gtf_index     = "/scratch/cqs/shengq1/references/gingivalis_W83/Gingivalis_W83.gtf";

my $target_w83_dir = create_directory_or_die( $target_dir . "/w83" );

my $w83config = {
  general  => { "task_name" => "w83", },
  bamfiles => {
    "2572-KCV-1-19" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm/result/2572-KCV-1-19/2572-KCV-1-19.bam"],
    "2572-KCV-1-20" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm/result/2572-KCV-1-20/2572-KCV-1-20.bam"],
    "2572-KCV-1-21" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm/result/2572-KCV-1-21/2572-KCV-1-21.bam"],
    "2572-KCV-1-22" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm/result/2572-KCV-1-22/2572-KCV-1-22.bam"],
    "2572-KCV-1-23" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm/result/2572-KCV-1-23/2572-KCV-1-23.bam"],
    "2572-KCV-1-24" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm/result/2572-KCV-1-24/2572-KCV-1-24.bam"],
    "2572-KCV-1-25" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm/result/2572-KCV-1-25/2572-KCV-1-25.bam"],
    "2572-KCV-1-26" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm/result/2572-KCV-1-26/2572-KCV-1-26.bam"],
    "2572-KCV-1-27" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm/result/2572-KCV-1-27/2572-KCV-1-27.bam"],
    "2572-KCV-1-28" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm/result/2572-KCV-1-28/2572-KCV-1-28.bam"],
    "2572-KCV-1-29" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm/result/2572-KCV-1-29/2572-KCV-1-29.bam"],
    "2572-KCV-1-30" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm/result/2572-KCV-1-30/2572-KCV-1-30.bam"],
  },
  countfiles => {
    "2572-KCV-1-19" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/identical/result/2572-KCV-1-19_clipped_identical.dupcount"],
    "2572-KCV-1-20" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/identical/result/2572-KCV-1-20_clipped_identical.dupcount"],
    "2572-KCV-1-21" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/identical/result/2572-KCV-1-21_clipped_identical.dupcount"],
    "2572-KCV-1-22" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/identical/result/2572-KCV-1-22_clipped_identical.dupcount"],
    "2572-KCV-1-23" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/identical/result/2572-KCV-1-23_clipped_identical.dupcount"],
    "2572-KCV-1-24" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/identical/result/2572-KCV-1-24_clipped_identical.dupcount"],
    "2572-KCV-1-25" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/identical/result/2572-KCV-1-25_clipped_identical.dupcount"],
    "2572-KCV-1-26" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/identical/result/2572-KCV-1-26_clipped_identical.dupcount"],
    "2572-KCV-1-27" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/identical/result/2572-KCV-1-27_clipped_identical.dupcount"],
    "2572-KCV-1-28" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/identical/result/2572-KCV-1-28_clipped_identical.dupcount"],
    "2572-KCV-1-29" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/identical/result/2572-KCV-1-29_clipped_identical.dupcount"],
    "2572-KCV-1-30" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/identical/result/2572-KCV-1-30_clipped_identical.dupcount"],
  },
  bam2fastq => {
    class               => "Bam2Fastq",
    perform             => 1,
    target_dir          => "${target_w83_dir}/bam2fastq",
    option              => "",
    source_ref          => "bamfiles",
    cqstools            => $cqstools,
    ispaired            => 0,
    unmapped_only       => 1,
    sort_before_convert => 0,
    sort_thread         => 12,
    unzipped            => 1,
    sh_direct           => 1,
    pbs                 => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "20gb"
    },
  },
  bowtie1_genome_cutadapt_topN_1mm => {
    class         => "Bowtie1",
    perform       => 1,
    target_dir    => "${target_w83_dir}/topN_bowtie1_genome_cutadapt_1mm",
    option        => $bowtie1_option_1mm,
    source_ref    => "bam2fastq",
    bowtie1_index => $w83_bowtie1_index,
    samonly       => 0,
    sh_direct     => 1,
    pbs           => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  count_1mm => {
    class           => "CQSMappedCount",
    perform         => 1,
    target_dir      => "${target_w83_dir}/bowtie1_genome_cutadapt_1mm_count",
    option          => "",
    source_ref      => "bowtie1_genome_cutadapt_topN_1mm",
    fastq_files_ref => "bam2fastq",
    seqcount_ref    => "countfiles",
    cqs_tools       => $cqstools,
    gff_file        => $w83_gtf_index,
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
    perform    => 1,
    target_dir => "${target_w83_dir}/summary",
    option     => "",
    source_ref => [ "count_1mm", ".xml" ],
    cqs_tools  => $cqstools,
    prefix     => "gingivalis_W83_1mm_",
    sh_direct  => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "10",
      "mem"      => "10gb"
    },
  },
};

#performConfig($w83config);

#
#my $parclip_files = {
#  "Parclip_01" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/PARCLIP/Vickers_Parclip_1_ATCACG_L002_R1.fastq.gz"],
#  "Parclip_02" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/PARCLIP/Vickers_Parclip_2_CGATGT_L002_R1.fastq.gz"],
#  "Parclip_03" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/PARCLIP/Vickers_Parclip_3_TTAGGC_L002_R1.fastq.gz"],
#  "Parclip_04" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/PARCLIP/Vickers_Parclip_4_TGACCA_L002_R1.fastq.gz"],
#  "Parclip_05" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/PARCLIP/Vickers_Parclip_5_ACAGTG_L002_R1.fastq.gz"],
#  "Parclip_06" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/PARCLIP/Vickers_Parclip_6_GCCAAT_L002_R1.fastq.gz"],
#  "Parclip_07" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/PARCLIP/Vickers_Parclip_7_CAGATC_L002_R1.fastq.gz"],
#  "Parclip_08" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/data/PARCLIP/Vickers_Parclip_8_ACTTGA_L002_R1.fastq.gz"]
#};
#
#my $target_parclip_dir = create_directory_or_die( $target_dir . "/parclip" );
#
#my $parclip_config = {
#  general      => { "task_name" => "parclip", },
#  fastqfiles   => $parclip_files,
#  cutadapt_len => {
#    class      => "Cutadapt",
#    perform    => 0,
#    target_dir => "${target_parclip_dir}/cutadapt_len",
#    option     => "-O 10 -m 12",
#    source_ref => "fastqfiles",
#    adaptor    => "TGGAATTCTCGGGTGCCAAGG",
#    extension  => "_clipped.fastq",
#    sh_direct  => 0,
#    gzipped    => 1,
#    pbs        => {
#      "email"    => $email,
#      "nodes"    => "1:ppn=1",
#      "walltime" => "24",
#      "mem"      => "20gb"
#    },
#  },
#  fastqlen => {
#    class      => "FastqLen",
#    perform    => 0,
#    target_dir => "${target_parclip_dir}/fastqlen",
#    option     => "",
#    source_ref => "cutadapt_len",
#    cqstools   => $cqstools,
#    sh_direct  => 1,
#    pbs        => {
#      "email"    => $email,
#      "nodes"    => "1:ppn=1",
#      "walltime" => "24",
#      "mem"      => "20gb"
#    },
#  },
#  bowtie1out => {
#    class         => "Bowtie1",
#    perform       => 0,
#    target_dir    => "${target_parclip_dir}/bowtie1out",
#    option        => "-v 2 -m 10 --best --strata",
#    source_ref    => [ "cutadapt_len", "fastq.gz\$" ],
#    bowtie1_index => "/data/cqs/guoy1/reference/hg19chr/bowtie1_index/hg19",
#    samformat     => 0,
#    samonly       => 0,
#    sh_direct     => 0,
#    pbs           => {
#      "email"    => $email,
#      "nodes"    => "1:ppn=1",
#      "walltime" => "72",
#      "mem"      => "20gb"
#    },
#  },
#  bowtie1bam => {
#    class         => "Bowtie1",
#    perform       => 0,
#    target_dir    => "${target_parclip_dir}/bowtie1bam",
#    option        => "-v 2 -m 10 --best --strata",
#    source_ref    => [ "cutadapt_len", "fastq.gz\$" ],
#    bowtie1_index => "/data/cqs/guoy1/reference/hg19chr/bowtie1_index/hg19",
#    samformat     => 1,
#    samonly       => 0,
#    sh_direct     => 0,
#    pbs           => {
#      "email"    => $email,
#      "nodes"    => "1:ppn=1",
#      "walltime" => "72",
#      "mem"      => "20gb"
#    },
#  },
#  PARalyzer => {
#    class      => "ParClip::PARalyzer",
#    perform    => 0,
#    target_dir => "${target_parclip_dir}/paralyzer",
#    option     => "",
#    source_ref => "bowtie1out",
#    genome2bit => "/data/cqs/guoy1/reference/hg19chr/hg19.2bit",
#    mirna_db   => "/data/cqs/shengq1/reference/miRBase20/hsa.mature.dna.db",
#    sh_direct  => 1,
#    pbs        => {
#      "email"    => $email,
#      "nodes"    => "1:ppn=1",
#      "walltime" => "72",
#      "mem"      => "20gb"
#    },
#  },
#  annotation => {
#    class            => "CQS::ParalyzerClusterAnnotator",
#    perform          => 1,
#    target_dir       => "${target_parclip_dir}/paralyzer",
#    option           => "-f miRNA",
#    source_ref       => [ "PARalyzer", ".cluster.csv" ],
#    cqstools         => $cqstools,
#    coordinate_files => [ $hg19_mrna_gff, $hg19_trna_bed ],
#    sh_direct        => 1,
#    pbs              => {
#      "email"    => $email,
#      "nodes"    => "1:ppn=1",
#      "walltime" => "72",
#      "mem"      => "20gb"
#    },
#  },
#};
#
#performConfig($parclip_config);

1;
