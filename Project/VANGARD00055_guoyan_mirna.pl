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

if ( is_linux() ) {
  $root      = "/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna";
  $cqs_tools = "/home/shengq1/cqstools/CQS.Tools.exe";
  $hsa_gffs  = "/data/cqs/shengq1/reference/gff3/hsa.gff3";
  $rno_gffs  = "/data/cqs/shengq1/reference/gff3/rno.gff3";
}
else {
  $root      = "d:/temp";
  $cqs_tools = "E:/sqh/programs/csharp/OmicsLabCSharp/CQS.Tools/bin/Release/CQS.Tools.exe";
  $hsa_gffs  = "H:/shengquanhu/projects/vangard/VANGARD00055_guoyan_mirna/hsa.gff3";
  $rno_gffs  = "H:/shengquanhu/projects/vangard/VANGARD00055_guoyan_mirna/rno.gff3";
}
my $target_dir = create_directory_or_die($root);

my $target_rat_dir   = create_directory_or_die( $target_dir . "/rat" );
my $target_human_dir = create_directory_or_die( $target_dir . "/human" );

my $email     = "quanhu.sheng\@vanderbilt.edu";
my $task_name = "VANGARD00055";

my $bwa_option             = "-l 8 -n 1 -o 0";
my $bwa_option_wholegenome = $bwa_option . " -t 8";
my $option_samse_mirna     = "";

my $bowtie2_option_top1 = "-D 20 -R 3 -N 1 -L 12 -i S,1,0.50 --gbar 50 --rdg 1000,1000 --rfg 1000,1000 -p 8";
my $bowtie2_option_topN = "-D 20 -R 3 -N 1 -L 12 -i S,1,0.50 --gbar 50 --rdg 1000,1000 --rfg 1000,1000 -k 20 -p 8";

my $bowtie2_rat_index   = "/data/cqs/shengq1/reference/rn4/bowtie2_index/rn4";
my $bowtie2_human_index = "/data/cqs/guoy1/reference/hg19/bowtie2_index/hg19";

my $bowtie1_option_topN = "-v 1 -n 1 -l 12 -k 100 --best --strata -p 8";

my $novoalign_option = "-l 15 -t 30 -r Random -m";

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
my $shrimp2_option      = "-Q -N 8 -o 1 --qv-offset 33";
my $shrimp2_rat_index   = "/data/cqs/shengq1/reference/rn4/shrimp2_index_ls_mirna/rn4-ls";
my $shrimp2_human_index = "/data/cqs/shengq1/reference/hg19/shrimp2_index_ls_mirna/hg19_chr-ls";
my $shrimp2_rat_miRBase_index   = "/data/cqs/shengq1/reference/miRBase19/shrimp_index_2.2.3/rno.mature.dna-ls";
my $shrimp2_human_miRBase_index = "/data/cqs/shengq1/reference/miRBase19/shrimp_index_2.2.3/hsa.mature.dna-ls";

my $rat = {
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
};
my $human = {
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
};

my $config_rat = {
  general  => { "task_name" => $task_name . "_rat", },
  cutadapt => {
    class      => "Cutadapt",
    perform    => 0,
    target_dir => "${target_dir}/cutadapt",
    option     => "",
    source     => $rat,
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
    option     => "-m 12 -M 49",
    source     => $rat,
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
  bowtie2_genome_cutadapt_top1 => {
    class         => "Bowtie2",
    perform       => 0,
    target_dir    => "${target_rat_dir}/top1_bowtie2_genome_cutadapt",
    option        => $bowtie2_option_top1,
    source_ref    => "cutadapt_len",
    bowtie2_index => $bowtie2_rat_index,
    samonly       => 0,
    sh_direct     => 0,
    pbs           => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "24",
      "mem"      => "20gb"
    },
  },
  mirna_count_bowtie2_genome_cutadapt_top1 => {
    class        => "MirnaCount",
    perform      => 0,
    target_dir   => "${target_rat_dir}/top1_bowtie2_genome_cutadapt_count",
    option       => "",
    source_ref   => "bowtie2_genome_cutadapt_top1",
    cqs_tools    => $cqs_tools,
    gff_file     => $rno_gffs,
    fasta_format => 0,
    sh_direct    => 1,
    pbs          => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "24",
      "mem"      => "40gb"
    },
  },
  shrimp2_bowtie2_genome_cutadapt_top1 => {
    class         => "Shrimp2",
    perform       => 0,
    target_dir    => "${target_rat_dir}/top1_bowtie2_genome_cutadapt_shrimp2",
    option        => $shrimp2_option,
    source_ref    => [ "mirna_count_bowtie2_genome_cutadapt_top1", ".fastq\$" ],
    shrimp2_index => $shrimp2_rat_index,
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
  bowtie2_genome_cutadapt_topN => {
    class         => "Bowtie2",
    perform       => 0,
    target_dir    => "${target_rat_dir}/topN_bowtie2_genome_cutadapt",
    option        => $bowtie2_option_topN,
    source_ref    => "cutadapt_len",
    bowtie2_index => $bowtie2_rat_index,
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
    perform      => 0,
    target_dir   => "${target_rat_dir}/topN_bowtie2_genome_cutadapt_count",
    option       => "",
    source_ref   => "bowtie2_genome_cutadapt_topN",
    cqs_tools    => $cqs_tools,
    gff_file     => $rno_gffs,
    fasta_format => 0,
    sh_direct    => 0,
    pbs          => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "24",
      "mem"      => "40gb"
    },
  },
  mirna_count_bowtie2_genome_cutadapt_topN_tRNA => {
    class        => "MirnaCount",
    perform      => 1,
    target_dir   => "${target_rat_dir}/topN_bowtie2_genome_cutadapt_count_tRNA",
    option       => "",
    source_ref   => "bowtie2_genome_cutadapt_topN",
    cqs_tools    => $cqs_tools,
    gff_file     => "/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/smrnapipeline/rn4_tRNA.bed",
    fasta_format => 0,
    sh_direct    => 1,
    pbs          => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "24",
      "mem"      => "40gb"
    },
  },
  shrimp2_bowtie2_genome_cutadapt_topN => {
    class         => "Shrimp2",
    perform       => 0,
    target_dir    => "${target_rat_dir}/topN_bowtie2_genome_cutadapt_shrimp2_miRBase",
    option        => $shrimp2_option,
    source_ref    => [ "mirna_count_bowtie2_genome_cutadapt_topN", ".fastq\$" ],
    shrimp2_index => $shrimp2_rat_miRBase_index,
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
  bowtie2_genome_identical_topN => {
    class         => "Bowtie2",
    perform       => 0,
    target_dir    => "${target_rat_dir}/topN_bowtie2_genome_identical",
    option        => $bowtie2_option_topN,
    source_ref    => [ "identical", ".fastq\$" ],
    bowtie2_index => $bowtie2_rat_index,
    samonly       => 0,
    sh_direct     => 0,
    pbs           => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "24",
      "mem"      => "20gb"
    },
  },
  mirna_count_bowtie2_genome_identical_topN => {
    class        => "MirnaCount",
    perform      => 0,
    target_dir   => "${target_rat_dir}/topN_bowtie2_genome_identical_count",
    option       => "",
    source_ref   => "bowtie2_genome_identical_topN",
    cqs_tools    => $cqs_tools,
    gff_file     => $rno_gffs,
    seqcount_ref => [ "identical", ".dupcount\$" ],
    fasta_format => 0,
    sh_direct    => 1,
    pbs          => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "24",
      "mem"      => "40gb"
    },
  },
  shrimp2_bowtie2_genome_identical_topN => {
    class         => "Shrimp2",
    perform       => 0,
    target_dir    => "${target_rat_dir}/topN_bowtie2_genome_identical_shrimp2",
    option        => $shrimp2_option,
    source_ref    => [ "mirna_count_bowtie2_genome_identical_topN", ".fastq\$" ],
    shrimp2_index => $shrimp2_rat_index,
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

  #  bowtie1_genome_cutadapt => {
  #    class         => "Bowtie1",
  #    perform       => 0,
  #    target_dir    => "${target_rat_dir}/bowtie1_genome_cutadapt",
  #    option        => $bowtie1_option_topN,
  #    source_ref    => "cutadapt",
  #    bowtie1_index => "/data/cqs/shengq1/reference/rn4/bowtie1_index/rn4",
  #    fasta_file    => "/data/cqs/shengq1/reference/rn4/bowtie1_index/rn4.fa",
  #    samonly       => 0,
  #    pbs           => {
  #      "email"    => $email,
  #      "nodes"    => "1:ppn=8",
  #      "walltime" => "24",
  #      "mem"      => "20gb"
  #    },
  #  },
  #  bowtie1_genome_cutadapt_mirna_count => {
  #    class        => "MirnaCount",
  #    perform      => 0,
  #    target_dir   => "${target_rat_dir}/bowtie1_genome_cutadapt",
  #    option       => "",
  #    source_ref   => "bowtie1_genome_cutadapt",
  #    cqs_tools    => $cqs_tools,
  #    gff_file     => $rno_gffs,
  #    fasta_format => 0,
  #    pbs          => {
  #      "email"    => $email,
  #      "nodes"    => "1:ppn=1",
  #      "walltime" => "24",
  #      "mem"      => "40gb"
  #    },
  #  },
  #  bowtie1_genome_identical => {
  #    class         => "Bowtie1",
  #    perform       => 0,
  #    target_dir    => "${target_rat_dir}/bowtie1_genome_identical",
  #    option        => $bowtie1_option_topN,
  #    source_ref    => [ "identical", ".fastq\$" ],
  #    bowtie1_index => "/data/cqs/shengq1/reference/rn4/bowtie1_index/rn4",
  #    fasta_file    => "/data/cqs/shengq1/reference/rn4/bowtie1_index/rn4.fa",
  #    samonly       => 0,
  #    pbs           => {
  #      "email"    => $email,
  #      "nodes"    => "1:ppn=8",
  #      "walltime" => "24",
  #      "mem"      => "20gb"
  #    },
  #  },
  #  bowtie1_genome_identical_mirna_count => {
  #    class        => "MirnaCount",
  #    perform      => 0,
  #    target_dir   => "${target_rat_dir}/bowtie1_genome_identical",
  #    option       => "",
  #    source_ref   => "bowtie1_genome_identical",
  #    cqs_tools    => $cqs_tools,
  #    gff_file     => $rno_gffs,
  #    fasta_format => 0,
  #    pbs          => {
  #      "email"    => $email,
  #      "nodes"    => "1:ppn=1",
  #      "walltime" => "24",
  #      "mem"      => "40gb"
  #    },
  #  },

  #  bwa_mature => {
  #    target_dir   => "${target_rat_dir}/bwa_miRBase_species",
  #    option       => $bwa_option,
  #    option_samse => $option_samse_mirna,
  #    source_ref   => "fastqfiles",
  #    fasta_file   => "/data/cqs/shengq1/reference/miRBase19/rno.mature.dna.fa",
  #    pbs          => {
  #      "email"    => $email,
  #      "nodes"    => "1:ppn=1",
  #      "walltime" => "24",
  #      "mem"      => "20gb"
  #    },
  #  },
  #  bwa => {
  #    target_dir      => "${target_rat_dir}/bwa_genome",
  #    option          => $bwa_option_wholegenome,
  #    option_samse    => "",
  #    source_ref      => "fastqfiles",
  #    fasta_file      => "/data/cqs/shengq1/reference/rn4/bwa_index/rn4.fa",
  #    estimate_insert => 0,
  #    source_ref      => "fastqfiles",
  #    pbs             => {
  #      "email"    => $email,
  #      "nodes"    => "1:ppn=8",
  #      "walltime" => "24",
  #      "mem"      => "20gb"
  #    },
  #  },
  #  mirna_count_bwa => {
  #    target_dir   => "${target_rat_dir}/bwa_genome",
  #    option       => "",
  #    source_ref   => "bwa",
  #    cqs_tools    => $cqs_tools,
  #    gff_file     => $rno_gffs,
  #    fasta_format => 0,
  #    pbs          => {
  #      "email"    => $email,
  #      "nodes"    => "1:ppn=1",
  #      "walltime" => "24",
  #      "mem"      => "40gb"
  #    },
  #  },
  #  mirna_count_bowtie1 => {
  #    target_dir   => "${target_rat_dir}/bowtie1_genome",
  #    option       => "",
  #    source_ref   => "bowtie1",
  #    cqs_tools    => $cqs_tools,
  #    gff_file     => $rno_gffs,
  #    fasta_format => 0,
  #    pbs          => {
  #      "email"    => $email,
  #      "nodes"    => "1:ppn=1",
  #      "walltime" => "24",
  #      "mem"      => "40gb"
  #    },
  #  },
  #  bowtie2_mature => {
  #    target_dir    => "${target_rat_dir}/bowtie2_mature",
  #    option        => $bowtie2_option_topN,
  #    source_ref    => "unmappedfiles",
  #    bowtie2_index => "/data/cqs/shengq1/reference/miRBase19/rno.mature.dna",
  #    fasta_file    => "/data/cqs/shengq1/reference/miRBase19/rno.mature.dna.fa",
  #    samonly       => 0,
  #    pbs           => {
  #      "email"    => $email,
  #      "nodes"    => "1:ppn=8",
  #      "walltime" => "24",
  #      "mem"      => "20gb"
  #    },
  #  },
};

my $config_human = {
  general  => { task_name => $task_name . "_human" },
  cutadapt => {
    class      => "Cutadapt",
    perform    => 0,
    target_dir => "${target_dir}/cutadapt",
    option     => "",
    source     => $human,
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
    option     => "-m 12 -M 49",
    source     => $human,
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
  bowtie2_genome_cutadapt_top1 => {
    class         => "Bowtie2",
    perform       => 0,
    target_dir    => "${target_human_dir}/top1_bowtie2_genome_cutadapt",
    option        => $bowtie2_option_top1,
    source_ref    => "cutadapt_len",
    bowtie2_index => $bowtie2_human_index,
    samonly       => 0,
    sh_direct     => 1,
    pbs           => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "24",
      "mem"      => "20gb"
    },
  },
  mirna_count_bowtie2_genome_cutadapt_top1 => {
    class        => "MirnaCount",
    perform      => 0,
    target_dir   => "${target_human_dir}/top1_bowtie2_genome_cutadapt_count",
    option       => "",
    source_ref   => "bowtie2_genome_cutadapt_top1",
    cqs_tools    => $cqs_tools,
    gff_file     => $hsa_gffs,
    fasta_format => 0,
    sh_direct    => 1,
    pbs          => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "24",
      "mem"      => "40gb"
    },
  },
  shrimp2_bowtie2_genome_cutadapt_top1 => {
    class         => "Shrimp2",
    perform       => 0,
    target_dir    => "${target_human_dir}/top1_bowtie2_genome_cutadapt_shrimp2",
    option        => $shrimp2_option,
    source_ref    => [ "mirna_count_bowtie2_genome_cutadapt_top1", ".fastq\$" ],
    shrimp2_index => $shrimp2_human_index,
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
  bowtie2_genome_cutadapt_topN => {
    class         => "Bowtie2",
    perform       => 0,
    target_dir    => "${target_human_dir}/topN_bowtie2_genome_cutadapt",
    option        => $bowtie2_option_topN,
    source_ref    => "cutadapt_len",
    bowtie2_index => $bowtie2_human_index,
    samonly       => 0,
    sh_direct     => 1,
    pbs           => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "24",
      "mem"      => "20gb"
    },
  },
  mirna_count_bowtie2_genome_cutadapt_topN => {
    class        => "MirnaCount",
    perform      => 0,
    target_dir   => "${target_human_dir}/topN_bowtie2_genome_cutadapt_count",
    option       => "",
    source_ref   => "bowtie2_genome_cutadapt_topN",
    cqs_tools    => $cqs_tools,
    gff_file     => $hsa_gffs,
    fasta_format => 0,
    sh_direct    => 1,
    pbs          => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "24",
      "mem"      => "40gb"
    },
  },
  mirna_count_bowtie2_genome_cutadapt_topN_tRNA => {
    class        => "MirnaCount",
    perform      => 1,
    target_dir   => "${target_human_dir}/topN_bowtie2_genome_cutadapt_count_tRNA",
    option       => "",
    source_ref   => "bowtie2_genome_cutadapt_topN",
    cqs_tools    => $cqs_tools,
    gff_file     => "/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna/smrnapipeline/hg19_tRNA.bed",
    fasta_format => 0,
    sh_direct    => 1,
    pbs          => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "24",
      "mem"      => "40gb"
    },
  },
  shrimp2_bowtie2_genome_cutadapt_topN => {
    class         => "Shrimp2",
    perform       => 0,
    target_dir    => "${target_human_dir}/topN_bowtie2_genome_cutadapt_shrimp2_miRBase",
    option        => $shrimp2_option,
    source_ref    => [ "mirna_count_bowtie2_genome_cutadapt_topN", ".fastq\$" ],
    shrimp2_index => $shrimp2_human_miRBase_index,
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
  bowtie2_genome_identical_topN => {
    class         => "Bowtie2",
    perform       => 0,
    target_dir    => "${target_human_dir}/topN_bowtie2_genome_identical",
    option        => $bowtie2_option_topN,
    source_ref    => [ "identical", ".fastq\$" ],
    bowtie2_index => $bowtie2_human_index,
    samonly       => 0,
    sh_direct     => 1,
    pbs           => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "24",
      "mem"      => "20gb"
    },
  },
  mirna_count_bowtie2_genome_identical_topN => {
    class        => "MirnaCount",
    perform      => 0,
    target_dir   => "${target_human_dir}/topN_bowtie2_genome_identical_count",
    option       => "",
    source_ref   => "bowtie2_genome_identical_topN",
    seqcount_ref => [ "identical", ".dupcount\$" ],
    cqs_tools    => $cqs_tools,
    gff_file     => $hsa_gffs,
    fasta_format => 0,
    sh_direct    => 1,
    pbs          => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "24",
      "mem"      => "40gb"
    },
  },
  shrimp2_bowtie2_genome_identical_topN => {
    class         => "Shrimp2",
    perform       => 0,
    target_dir    => "${target_human_dir}/topN_bowtie2_genome_identical_shrimp2",
    option        => $shrimp2_option,
    source_ref    => [ "mirna_count_bowtie2_genome_identical_topN", ".fastq\$" ],
    shrimp2_index => $shrimp2_human_index,
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

  #  bowtie1_cutadapt => {
  #    class         => "Bowtie1",
  #    target_dir    => "${target_human_dir}/bowtie1_genome_cutadapt",
  #    option        => $bowtie1_option_topN,
  #    source_ref    => "cutadapt",
  #    bowtie1_index => "/data/cqs/guoy1/reference/hg19/bowtie_index/hg19",
  #    fasta_file    => "/data/cqs/guoy1/reference/hg19/bowtie_index/hg19.fa",
  #    sh_direct     => 0,
  #    samonly       => 0,
  #    pbs           => {
  #      "email"    => $email,
  #      "nodes"    => "1:ppn=8",
  #      "walltime" => "24",
  #      "mem"      => "20gb"
  #    },
  #  },
  #  bowtie1_identical => {
  #    class         => "Bowtie1",
  #    target_dir    => "${target_human_dir}/bowtie1_genome_identical",
  #    option        => $bowtie1_option_topN,
  #    source_ref    => [ "identical", ".fastq\$" ],
  #    bowtie1_index => "/data/cqs/guoy1/reference/hg19/bowtie_index/hg19",
  #    fasta_file    => "/data/cqs/guoy1/reference/hg19/bowtie_index/hg19.fa",
  #    sh_direct     => 0,
  #    samonly       => 0,
  #    pbs           => {
  #      "email"    => $email,
  #      "nodes"    => "1:ppn=8",
  #      "walltime" => "24",
  #      "mem"      => "20gb"
  #    },
  #  },
  #  mirna_count_bowtie1_identical => {
  #    class        => "MirnaCount",
  #    target_dir   => "${target_human_dir}/bowtie1_genome_identical",
  #    option       => "",
  #    source_ref   => "bowtie1_identical",
  #    seqcount_ref => [ "identical", ".count\$" ],
  #    cqs_tools    => $cqs_tools,
  #    gff_file     => $hsa_gffs,
  #    sh_direct    => 1,
  #    fasta_format => 0,
  #    pbs          => {
  #      "email"    => $email,
  #      "nodes"    => "1:ppn=1",
  #      "walltime" => "24",
  #      "mem"      => "20gb"
  #    },
  #  },

  #  bwa_mature => {
  #    target_dir   => "${target_human_dir}/bwa_miRBase_species",
  #    option       => $bwa_option,
  #    option_samse => $option_samse_mirna,
  #    source_ref   => "fastqfiles",
  #    fasta_file   => "/data/cqs/shengq1/reference/miRBase19/hsa.mature.dna.fa",
  #    pbs          => {
  #      "email"    => $email,
  #      "nodes"    => "1:ppn=1",
  #      "walltime" => "24",
  #      "mem"      => "20gb"
  #    },
  #  },
  #  bwa => {
  #    target_dir      => "${target_human_dir}/bwa_genome",
  #    option          => $bwa_option_wholegenome,
  #    option_samse    => "",
  #    source_ref      => "fastqfiles",
  #    fasta_file      => "/data/cqs/shengq1/reference/hg19/hg19_chr.fa",
  #    estimate_insert => 0,
  #    source_ref      => "fastqfiles",
  #    sh_direct       => 1,
  #    pbs             => {
  #      "email"    => $email,
  #      "nodes"    => "1:ppn=8",
  #      "walltime" => "24",
  #      "mem"      => "20gb"
  #    },
  #  },
  #  mirna_count_bwa => {
  #    target_dir   => "${target_human_dir}/bwa_genome",
  #    option       => "",
  #    source_ref   => "bwa",
  #    cqs_tools    => $cqs_tools,
  #    gff_file     => $hsa_gffs,
  #    sh_direct    => 1,
  #    fasta_format => 0,
  #    pbs          => {
  #      "email"    => $email,
  #      "nodes"    => "1:ppn=1",
  #      "walltime" => "24",
  #      "mem"      => "20gb"
  #    },
  #  },
  #  bowtie2_mature => {
  #    class =>"Bowtiw2",
  #    perform => 0,
  #    target_dir    => "${target_human_dir}/bowtie2_mature",
  #    option        => $bowtie2_option_topN,
  #    source_ref    => "fastqfiles",
  #    bowtie2_index => "/data/cqs/shengq1/reference/miRBase19/hsa.mature.dna",
  #    samonly       => 0,
  #    sh_direct     => 1,
  #    pbs           => {
  #      "email"    => $email,
  #      "nodes"    => "1:ppn=8",
  #      "walltime" => "24",
  #      "mem"      => "20gb"
  #    },
  #  },
};

#my $cutadapt = new CQS::Cutadapt();
#my $all_cutadapt = { %{ $cutadapt->result( $config_rat, "cutadapt" ) }, %{ $cutadapt->result( $config_human, "cutadapt" ) } };
#foreach my $k ( sort keys %{$all_cutadapt} ) {
#  print "$k => @{$all_cutadapt->{$k}}\n";
#}
#my $config_mirna = {
#  general => {
#    path_file => "/home/shengq1/local/bin/path.txt",
#    task_name => $task_name . "_mirna"
#  },
#  fastqfiles => $all_cutadapt,
#  bwa_mature => {
#    target_dir   => "${target_dir}/bwa_miRBase_mature",
#    option       => $bwa_option,
#    option_samse => $option_samse_mirna,
#    source_ref   => "fastqfiles",
#    fasta_file   => "/data/cqs/shengq1/reference/miRBase19/mature.dna.fa",
#    pbs          => {
#      "email"    => $email,
#      "nodes"    => "1:ppn=1",
#      "walltime" => "24",
#      "mem"      => "20gb"
#    },
#  },
#  bwa_hairpin => {
#    target_dir   => "${target_dir}/bwa_miRBase_hairpin",
#    option       => $bwa_option,
#    option_samse => $option_samse_mirna,
#    source_ref   => "fastqfiles",
#    fasta_file   => "/data/cqs/shengq1/reference/miRBase19/hairpin.dna.fa",
#    pbs          => {
#      "email"    => $email,
#      "nodes"    => "1:ppn=1",
#      "walltime" => "24",
#      "mem"      => "20gb"
#    },
#  },
#  bwa_illumina => {
#    target_dir   => "${target_dir}/bwa_illumina_miRNA",
#    option       => $bwa_option,
#    option_samse => $option_samse_mirna,
#    source_ref   => "fastqfiles",
#    fasta_file   => "/data/cqs/shengq1/reference/miRNA_illumina/mir.fa",
#    pbs          => {
#      "email"    => $email,
#      "nodes"    => "1:ppn=1",
#      "walltime" => "24",
#      "mem"      => "20gb"
#    },
#  },
#  novoalign_mature => {
#    target_dir   => "${target_dir}/novoalign_miRBase_mature",
#    option       => $novoalign_option,
#    option_samse => "",
#    source_ref   => "fastqfiles",
#    novoindex    => "/data/cqs/shengq1/reference/miRBase19/mature.dna.nix",
#    pbs          => {
#      "email"    => $email,
#      "nodes"    => "1:ppn=1",
#      "walltime" => "24",
#      "mem"      => "20gb"
#    },
#  },
#  novoalign_hairpin => {
#    target_dir   => "${target_dir}/novoalign_miRBase_hairpin",
#    option       => $novoalign_option,
#    option_samse => "",
#    source_ref   => "fastqfiles",
#    novoindex    => "/data/cqs/shengq1/reference/miRBase19/hairpin.dna.nix",
#    pbs          => {
#      "email"    => $email,
#      "nodes"    => "1:ppn=1",
#      "walltime" => "24",
#      "mem"      => "20gb"
#    },
#  },
#  novoalign_illumina => {
#    target_dir   => "${target_dir}/novoalign_illumina",
#    option       => $novoalign_option,
#    option_samse => "",
#    source_ref   => "fastqfiles",
#    novoindex    => "/data/cqs/shengq1/reference/miRNA_illumina/mir.nix",
#    pbs          => {
#      "email"    => $email,
#      "nodes"    => "1:ppn=1",
#      "walltime" => "24",
#      "mem"      => "20gb"
#    },
#  },
#
#};

#performConfig($config_rat);
#performConfig($config_human);

performTask($config_rat,"mirna_count_bowtie2_genome_cutadapt_topN");

#performConfig($config_rat, "^shrimp2", 1);
#performConfig($config_human, "^shrimp2", 1);

#bwa_by_pbs_single( $config_mirna, "bwa_mature" );
#bwa_by_pbs_single( $config_mirna, "bwa_hairpin" );
#bwa_by_pbs_single( $config_mirna, "bwa_illumina" );

#novoalign( $config_mirna, "novoalign_mature" );
#novoalign( $config_mirna, "novoalign_hairpin" );
#novoalign( $config_mirna, "novoalign_illumina" );

#bwa_by_pbs_single( $config_rat,   "bwa_mature" );
#bwa_by_pbs_single( $config_human, "bwa_mature" );

#bwa_by_pbs_single( $config_rat, "bwa" );
#bwa_by_pbs_single( $config_human, "bwa" );

1;
