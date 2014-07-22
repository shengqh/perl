#!/usr/bin/perl
use strict;
use warnings;

use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;

my $target_dir = create_directory_or_die("/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201405_smallRNA_2829_2570/");
my $cqstools   = "/home/shengq1/cqstools/CQS.Tools.exe";

my $email     = "quanhu.sheng\@vanderbilt.edu";
my $task_name = "Vicky_2829_2570";

my $config = {
  general         => { "task_name" => $task_name },
  miRNA_1mm_table => {
    class      => "CQSMirnaTable",
    perform    => 0,
    target_dir => "${target_dir}/topN_bowtie1_genome_cutadapt_1mm_count_miRNA_table",
    option     => "",
    source     => {
      "01-018-Post_CTTGTA" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_miRNA/result/01-018-Post_CTTGTA/01-018-Post_CTTGTA.bam.count"],
      "01-018-Pre_GGCTAC" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_miRNA/result/01-018-Pre_GGCTAC/01-018-Pre_GGCTAC.bam.count"],
      "01-031-Post_GGTAGC" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_miRNA/result/01-031-Post_GGTAGC/01-031-Post_GGTAGC.bam.count"],
      "01-031-Pre_GAGTGG" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_miRNA/result/01-031-Pre_GAGTGG/01-031-Pre_GAGTGG.bam.count"],
      "01-061-Post_ATGAGC" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_miRNA/result/01-061-Post_ATGAGC/01-061-Post_ATGAGC.bam.count"],
      "01-061-Pre_ACTGAT" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_miRNA/result/01-061-Pre_ACTGAT/01-061-Pre_ACTGAT.bam.count"],
      "01-28-Post_CGATGT" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_miRNA/result/01-28-Post_CGATGT/01-28-Post_CGATGT.bam.count"],
      "01-28-Pre_ATCACG"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_miRNA/result/01-28-Pre_ATCACG/01-28-Pre_ATCACG.bam.count"],
      "01-29-Pre_TGACCA"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_miRNA/result/01-29-Pre_TGACCA/01-29-Pre_TGACCA.bam.count"],
      "01-29-Pre_TTAGGC"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_miRNA/result/01-29-Pre_TTAGGC/01-29-Pre_TTAGGC.bam.count"],
      "01-36-Post_GCCAAT" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_miRNA/result/01-36-Post_GCCAAT/01-36-Post_GCCAAT.bam.count"],
      "01-36-Pre_ACAGTG" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_miRNA/result/01-36-Pre_ACAGTG/01-36-Pre_ACAGTG.bam.count"],
      "03-007-Post_CAAAAG" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_miRNA/result/03-007-Post_CAAAAG/03-007-Post_CAAAAG.bam.count"],
      "03-007-Pre_ATTCCT" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_miRNA/result/03-007-Pre_ATTCCT/03-007-Pre_ATTCCT.bam.count"],
      "03-011-Post_CACCGG" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_miRNA/result/03-011-Post_CACCGG/03-011-Post_CACCGG.bam.count"],
      "03-011-Pre_CAACTA" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_miRNA/result/03-011-Pre_CAACTA/03-011-Pre_CAACTA.bam.count"],
      "03-015-Post_CACTCA" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_miRNA/result/03-015-Post_CACTCA/03-015-Post_CACTCA.bam.count"],
      "03-015-Pre_CACGAT" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_miRNA/result/03-015-Pre_CACGAT/03-015-Pre_CACGAT.bam.count"],
      "03-018-Post_CGTACG" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_miRNA/result/03-018-Post_CGTACG/03-018-Post_CGTACG.bam.count"],
      "03-018-Pre_GTTTCG" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_miRNA/result/03-018-Pre_GTTTCG/03-018-Pre_GTTTCG.bam.count"],
      "03-026-Post_AGTTCC" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_miRNA/result/03-026-Post_AGTTCC/03-026-Post_AGTTCC.bam.count"],
      "03-026-Pre_AGTCAA" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_miRNA/result/03-026-Pre_AGTCAA/03-026-Pre_AGTCAA.bam.count"],
      "03-031-Post_CATGGC" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_miRNA/result/03-031-Post_CATGGC/03-031-Post_CATGGC.bam.count"],
      "03-031-Pre_CAGGCG" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_miRNA/result/03-031-Pre_CAGGCG/03-031-Pre_CAGGCG.bam.count"],
      "03-033-Post_CCAACA" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_miRNA/result/03-033-Post_CCAACA/03-033-Post_CCAACA.bam.count"],
      "03-033-Pre_CATTTT" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_miRNA/result/03-033-Pre_CATTTT/03-033-Pre_CATTTT.bam.count"],
      "03-036-Post_CCGTCC" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_miRNA/result/03-036-Post_CCGTCC/03-036-Post_CCGTCC.bam.count"],
      "03-036-Pre_ATGTCA" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_miRNA/result/03-036-Pre_ATGTCA/03-036-Pre_ATGTCA.bam.count"],
      "03-047-Post_CTAGCT" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_miRNA/result/03-047-Post_CTAGCT/03-047-Post_CTAGCT.bam.count"],
      "03-047-Pre_CGGAAT" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_miRNA/result/03-047-Pre_CGGAAT/03-047-Pre_CGGAAT.bam.count"],
      "03-049-Post_CTCAGA" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_miRNA/result/03-049-Post_CTCAGA/03-049-Post_CTCAGA.bam.count"],
      "03-049-Pre_CTATAC" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_miRNA/result/03-049-Pre_CTATAC/03-049-Pre_CTATAC.bam.count"],
      "03-063-Post_GTGGCC" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_miRNA/result/03-063-Post_GTGGCC/03-063-Post_GTGGCC.bam.count"],
      "03-063-Pre_GTGAAA" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_miRNA/result/03-063-Pre_GTGAAA/03-063-Pre_GTGAAA.bam.count"],
      "03-065-Post_GTCCGC" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_miRNA/result/03-065-Post_GTCCGC/03-065-Post_GTCCGC.bam.count"],
      "03-065-Pre_GTAGAG" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_miRNA/result/03-065-Pre_GTAGAG/03-065-Pre_GTAGAG.bam.count"],
      "03-16-Post_ACTTGA" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_miRNA/result/03-16-Post_ACTTGA/03-16-Post_ACTTGA.bam.count"],
      "03-16-Pre_CAGATC"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_miRNA/result/03-16-Pre_CAGATC/03-16-Pre_CAGATC.bam.count"],
      "03-17-Post_TAGCTT" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_miRNA/result/03-17-Post_TAGCTT/03-17-Post_TAGCTT.bam.count"],
      "03-17-Pre_GATCAG" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_miRNA/result/03-17-Pre_GATCAG/03-17-Pre_GATCAG.bam.count"],
      "2829-KCV-1A"      => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201405_smallRNA_2829/topN_bowtie1_genome_cutadapt_1mm_count_miRNA/result/2829-KCV-1A/2829-KCV-1A.bam.count"],
      "2829-KCV-1B"      => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201405_smallRNA_2829/topN_bowtie1_genome_cutadapt_1mm_count_miRNA/result/2829-KCV-1B/2829-KCV-1B.bam.count"],
      "2829-KCV-1C"      => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201405_smallRNA_2829/topN_bowtie1_genome_cutadapt_1mm_count_miRNA/result/2829-KCV-1C/2829-KCV-1C.bam.count"],
      "2829-KCV-1D"      => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201405_smallRNA_2829/topN_bowtie1_genome_cutadapt_1mm_count_miRNA/result/2829-KCV-1D/2829-KCV-1D.bam.count"],
      "2829-KCV-1E"      => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201405_smallRNA_2829/topN_bowtie1_genome_cutadapt_1mm_count_miRNA/result/2829-KCV-1E/2829-KCV-1E.bam.count"],
      "2829-KCV-1F"      => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201405_smallRNA_2829/topN_bowtie1_genome_cutadapt_1mm_count_miRNA/result/2829-KCV-1F/2829-KCV-1F.bam.count"],
      "2829-KCV-1G"      => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201405_smallRNA_2829/topN_bowtie1_genome_cutadapt_1mm_count_miRNA/result/2829-KCV-1G/2829-KCV-1G.bam.count"],
      "2829-KCV-1H"      => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201405_smallRNA_2829/topN_bowtie1_genome_cutadapt_1mm_count_miRNA/result/2829-KCV-1H/2829-KCV-1H.bam.count"],
      "2829-KCV-1I"      => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201405_smallRNA_2829/topN_bowtie1_genome_cutadapt_1mm_count_miRNA/result/2829-KCV-1I/2829-KCV-1I.bam.count"],
      "2829-KCV-1J"      => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201405_smallRNA_2829/topN_bowtie1_genome_cutadapt_1mm_count_miRNA/result/2829-KCV-1J/2829-KCV-1J.bam.count"],
    },
    cqs_tools => $cqstools,
    prefix    => "miRNA_1mm_",
    sh_direct => 1,
    pbs       => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "10",
      "mem"      => "10gb"
    },
  },
  tRNA_1mm_table => {
    class      => "CQSMappedTable",
    perform    => 0,
    target_dir => "${target_dir}/topN_bowtie1_genome_cutadapt_1mm_count_tRNA_table",
    option     => "",
    source     => {
      "01-018-Post_CTTGTA" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_tRNA/result/01-018-Post_CTTGTA/01-018-Post_CTTGTA.bam.count.mapped.xml"],
      "01-018-Pre_GGCTAC" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_tRNA/result/01-018-Pre_GGCTAC/01-018-Pre_GGCTAC.bam.count.mapped.xml"],
      "01-031-Post_GGTAGC" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_tRNA/result/01-031-Post_GGTAGC/01-031-Post_GGTAGC.bam.count.mapped.xml"],
      "01-031-Pre_GAGTGG" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_tRNA/result/01-031-Pre_GAGTGG/01-031-Pre_GAGTGG.bam.count.mapped.xml"],
      "01-061-Post_ATGAGC" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_tRNA/result/01-061-Post_ATGAGC/01-061-Post_ATGAGC.bam.count.mapped.xml"],
      "01-061-Pre_ACTGAT" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_tRNA/result/01-061-Pre_ACTGAT/01-061-Pre_ACTGAT.bam.count.mapped.xml"],
      "01-28-Post_CGATGT" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_tRNA/result/01-28-Post_CGATGT/01-28-Post_CGATGT.bam.count.mapped.xml"],
      "01-28-Pre_ATCACG" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_tRNA/result/01-28-Pre_ATCACG/01-28-Pre_ATCACG.bam.count.mapped.xml"],
      "01-29-Pre_TGACCA" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_tRNA/result/01-29-Pre_TGACCA/01-29-Pre_TGACCA.bam.count.mapped.xml"],
      "01-29-Pre_TTAGGC" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_tRNA/result/01-29-Pre_TTAGGC/01-29-Pre_TTAGGC.bam.count.mapped.xml"],
      "01-36-Post_GCCAAT" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_tRNA/result/01-36-Post_GCCAAT/01-36-Post_GCCAAT.bam.count.mapped.xml"],
      "01-36-Pre_ACAGTG" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_tRNA/result/01-36-Pre_ACAGTG/01-36-Pre_ACAGTG.bam.count.mapped.xml"],
      "03-007-Post_CAAAAG" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_tRNA/result/03-007-Post_CAAAAG/03-007-Post_CAAAAG.bam.count.mapped.xml"],
      "03-007-Pre_ATTCCT" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_tRNA/result/03-007-Pre_ATTCCT/03-007-Pre_ATTCCT.bam.count.mapped.xml"],
      "03-011-Post_CACCGG" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_tRNA/result/03-011-Post_CACCGG/03-011-Post_CACCGG.bam.count.mapped.xml"],
      "03-011-Pre_CAACTA" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_tRNA/result/03-011-Pre_CAACTA/03-011-Pre_CAACTA.bam.count.mapped.xml"],
      "03-015-Post_CACTCA" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_tRNA/result/03-015-Post_CACTCA/03-015-Post_CACTCA.bam.count.mapped.xml"],
      "03-015-Pre_CACGAT" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_tRNA/result/03-015-Pre_CACGAT/03-015-Pre_CACGAT.bam.count.mapped.xml"],
      "03-018-Post_CGTACG" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_tRNA/result/03-018-Post_CGTACG/03-018-Post_CGTACG.bam.count.mapped.xml"],
      "03-018-Pre_GTTTCG" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_tRNA/result/03-018-Pre_GTTTCG/03-018-Pre_GTTTCG.bam.count.mapped.xml"],
      "03-026-Post_AGTTCC" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_tRNA/result/03-026-Post_AGTTCC/03-026-Post_AGTTCC.bam.count.mapped.xml"],
      "03-026-Pre_AGTCAA" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_tRNA/result/03-026-Pre_AGTCAA/03-026-Pre_AGTCAA.bam.count.mapped.xml"],
      "03-031-Post_CATGGC" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_tRNA/result/03-031-Post_CATGGC/03-031-Post_CATGGC.bam.count.mapped.xml"],
      "03-031-Pre_CAGGCG" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_tRNA/result/03-031-Pre_CAGGCG/03-031-Pre_CAGGCG.bam.count.mapped.xml"],
      "03-033-Post_CCAACA" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_tRNA/result/03-033-Post_CCAACA/03-033-Post_CCAACA.bam.count.mapped.xml"],
      "03-033-Pre_CATTTT" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_tRNA/result/03-033-Pre_CATTTT/03-033-Pre_CATTTT.bam.count.mapped.xml"],
      "03-036-Post_CCGTCC" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_tRNA/result/03-036-Post_CCGTCC/03-036-Post_CCGTCC.bam.count.mapped.xml"],
      "03-036-Pre_ATGTCA" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_tRNA/result/03-036-Pre_ATGTCA/03-036-Pre_ATGTCA.bam.count.mapped.xml"],
      "03-047-Post_CTAGCT" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_tRNA/result/03-047-Post_CTAGCT/03-047-Post_CTAGCT.bam.count.mapped.xml"],
      "03-047-Pre_CGGAAT" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_tRNA/result/03-047-Pre_CGGAAT/03-047-Pre_CGGAAT.bam.count.mapped.xml"],
      "03-049-Post_CTCAGA" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_tRNA/result/03-049-Post_CTCAGA/03-049-Post_CTCAGA.bam.count.mapped.xml"],
      "03-049-Pre_CTATAC" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_tRNA/result/03-049-Pre_CTATAC/03-049-Pre_CTATAC.bam.count.mapped.xml"],
      "03-063-Post_GTGGCC" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_tRNA/result/03-063-Post_GTGGCC/03-063-Post_GTGGCC.bam.count.mapped.xml"],
      "03-063-Pre_GTGAAA" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_tRNA/result/03-063-Pre_GTGAAA/03-063-Pre_GTGAAA.bam.count.mapped.xml"],
      "03-065-Post_GTCCGC" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_tRNA/result/03-065-Post_GTCCGC/03-065-Post_GTCCGC.bam.count.mapped.xml"],
      "03-065-Pre_GTAGAG" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_tRNA/result/03-065-Pre_GTAGAG/03-065-Pre_GTAGAG.bam.count.mapped.xml"],
      "03-16-Post_ACTTGA" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_tRNA/result/03-16-Post_ACTTGA/03-16-Post_ACTTGA.bam.count.mapped.xml"],
      "03-16-Pre_CAGATC" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_tRNA/result/03-16-Pre_CAGATC/03-16-Pre_CAGATC.bam.count.mapped.xml"],
      "03-17-Post_TAGCTT" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_tRNA/result/03-17-Post_TAGCTT/03-17-Post_TAGCTT.bam.count.mapped.xml"],
      "03-17-Pre_GATCAG" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_tRNA/result/03-17-Pre_GATCAG/03-17-Pre_GATCAG.bam.count.mapped.xml"],
      "2829-KCV-1A" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201405_smallRNA_2829/topN_bowtie1_genome_cutadapt_1mm_count_tRNA/result/2829-KCV-1A/2829-KCV-1A.bam.count.mapped.xml"],
      "2829-KCV-1B" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201405_smallRNA_2829/topN_bowtie1_genome_cutadapt_1mm_count_tRNA/result/2829-KCV-1B/2829-KCV-1B.bam.count.mapped.xml"],
      "2829-KCV-1C" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201405_smallRNA_2829/topN_bowtie1_genome_cutadapt_1mm_count_tRNA/result/2829-KCV-1C/2829-KCV-1C.bam.count.mapped.xml"],
      "2829-KCV-1D" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201405_smallRNA_2829/topN_bowtie1_genome_cutadapt_1mm_count_tRNA/result/2829-KCV-1D/2829-KCV-1D.bam.count.mapped.xml"],
      "2829-KCV-1E" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201405_smallRNA_2829/topN_bowtie1_genome_cutadapt_1mm_count_tRNA/result/2829-KCV-1E/2829-KCV-1E.bam.count.mapped.xml"],
      "2829-KCV-1F" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201405_smallRNA_2829/topN_bowtie1_genome_cutadapt_1mm_count_tRNA/result/2829-KCV-1F/2829-KCV-1F.bam.count.mapped.xml"],
      "2829-KCV-1G" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201405_smallRNA_2829/topN_bowtie1_genome_cutadapt_1mm_count_tRNA/result/2829-KCV-1G/2829-KCV-1G.bam.count.mapped.xml"],
      "2829-KCV-1H" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201405_smallRNA_2829/topN_bowtie1_genome_cutadapt_1mm_count_tRNA/result/2829-KCV-1H/2829-KCV-1H.bam.count.mapped.xml"],
      "2829-KCV-1I" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201405_smallRNA_2829/topN_bowtie1_genome_cutadapt_1mm_count_tRNA/result/2829-KCV-1I/2829-KCV-1I.bam.count.mapped.xml"],
      "2829-KCV-1J" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201405_smallRNA_2829/topN_bowtie1_genome_cutadapt_1mm_count_tRNA/result/2829-KCV-1J/2829-KCV-1J.bam.count.mapped.xml"],
    },
    cqs_tools => $cqstools,
    prefix    => "tRNA_1mm_",
    sh_direct => 1,
    pbs       => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "10",
      "mem"      => "10gb"
    },
  },
  smallRNA_1mm_table => {
    class      => "CQSMappedTable",
    perform    => 1,
    target_dir => "${target_dir}/topN_bowtie1_genome_cutadapt_1mm_count_smallRNA_table",
    option     => "",
    source     => {
      "01-018-Post_CTTGTA" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_smallRNA/result/01-018-Post_CTTGTA/01-018-Post_CTTGTA.bam.count.mapped.xml"],
      "01-018-Pre_GGCTAC" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_smallRNA/result/01-018-Pre_GGCTAC/01-018-Pre_GGCTAC.bam.count.mapped.xml"],
      "01-031-Post_GGTAGC" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_smallRNA/result/01-031-Post_GGTAGC/01-031-Post_GGTAGC.bam.count.mapped.xml"],
      "01-031-Pre_GAGTGG" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_smallRNA/result/01-031-Pre_GAGTGG/01-031-Pre_GAGTGG.bam.count.mapped.xml"],
      "01-061-Post_ATGAGC" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_smallRNA/result/01-061-Post_ATGAGC/01-061-Post_ATGAGC.bam.count.mapped.xml"],
      "01-061-Pre_ACTGAT" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_smallRNA/result/01-061-Pre_ACTGAT/01-061-Pre_ACTGAT.bam.count.mapped.xml"],
      "01-28-Post_CGATGT" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_smallRNA/result/01-28-Post_CGATGT/01-28-Post_CGATGT.bam.count.mapped.xml"],
      "01-28-Pre_ATCACG" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_smallRNA/result/01-28-Pre_ATCACG/01-28-Pre_ATCACG.bam.count.mapped.xml"],
      "01-29-Pre_TGACCA" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_smallRNA/result/01-29-Pre_TGACCA/01-29-Pre_TGACCA.bam.count.mapped.xml"],
      "01-29-Pre_TTAGGC" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_smallRNA/result/01-29-Pre_TTAGGC/01-29-Pre_TTAGGC.bam.count.mapped.xml"],
      "01-36-Post_GCCAAT" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_smallRNA/result/01-36-Post_GCCAAT/01-36-Post_GCCAAT.bam.count.mapped.xml"],
      "01-36-Pre_ACAGTG" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_smallRNA/result/01-36-Pre_ACAGTG/01-36-Pre_ACAGTG.bam.count.mapped.xml"],
      "03-007-Post_CAAAAG" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_smallRNA/result/03-007-Post_CAAAAG/03-007-Post_CAAAAG.bam.count.mapped.xml"],
      "03-007-Pre_ATTCCT" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_smallRNA/result/03-007-Pre_ATTCCT/03-007-Pre_ATTCCT.bam.count.mapped.xml"],
      "03-011-Post_CACCGG" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_smallRNA/result/03-011-Post_CACCGG/03-011-Post_CACCGG.bam.count.mapped.xml"],
      "03-011-Pre_CAACTA" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_smallRNA/result/03-011-Pre_CAACTA/03-011-Pre_CAACTA.bam.count.mapped.xml"],
      "03-015-Post_CACTCA" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_smallRNA/result/03-015-Post_CACTCA/03-015-Post_CACTCA.bam.count.mapped.xml"],
      "03-015-Pre_CACGAT" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_smallRNA/result/03-015-Pre_CACGAT/03-015-Pre_CACGAT.bam.count.mapped.xml"],
      "03-018-Post_CGTACG" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_smallRNA/result/03-018-Post_CGTACG/03-018-Post_CGTACG.bam.count.mapped.xml"],
      "03-018-Pre_GTTTCG" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_smallRNA/result/03-018-Pre_GTTTCG/03-018-Pre_GTTTCG.bam.count.mapped.xml"],
      "03-026-Post_AGTTCC" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_smallRNA/result/03-026-Post_AGTTCC/03-026-Post_AGTTCC.bam.count.mapped.xml"],
      "03-026-Pre_AGTCAA" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_smallRNA/result/03-026-Pre_AGTCAA/03-026-Pre_AGTCAA.bam.count.mapped.xml"],
      "03-031-Post_CATGGC" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_smallRNA/result/03-031-Post_CATGGC/03-031-Post_CATGGC.bam.count.mapped.xml"],
      "03-031-Pre_CAGGCG" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_smallRNA/result/03-031-Pre_CAGGCG/03-031-Pre_CAGGCG.bam.count.mapped.xml"],
      "03-033-Post_CCAACA" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_smallRNA/result/03-033-Post_CCAACA/03-033-Post_CCAACA.bam.count.mapped.xml"],
      "03-033-Pre_CATTTT" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_smallRNA/result/03-033-Pre_CATTTT/03-033-Pre_CATTTT.bam.count.mapped.xml"],
      "03-036-Post_CCGTCC" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_smallRNA/result/03-036-Post_CCGTCC/03-036-Post_CCGTCC.bam.count.mapped.xml"],
      "03-036-Pre_ATGTCA" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_smallRNA/result/03-036-Pre_ATGTCA/03-036-Pre_ATGTCA.bam.count.mapped.xml"],
      "03-047-Post_CTAGCT" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_smallRNA/result/03-047-Post_CTAGCT/03-047-Post_CTAGCT.bam.count.mapped.xml"],
      "03-047-Pre_CGGAAT" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_smallRNA/result/03-047-Pre_CGGAAT/03-047-Pre_CGGAAT.bam.count.mapped.xml"],
      "03-049-Post_CTCAGA" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_smallRNA/result/03-049-Post_CTCAGA/03-049-Post_CTCAGA.bam.count.mapped.xml"],
      "03-049-Pre_CTATAC" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_smallRNA/result/03-049-Pre_CTATAC/03-049-Pre_CTATAC.bam.count.mapped.xml"],
      "03-063-Post_GTGGCC" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_smallRNA/result/03-063-Post_GTGGCC/03-063-Post_GTGGCC.bam.count.mapped.xml"],
      "03-063-Pre_GTGAAA" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_smallRNA/result/03-063-Pre_GTGAAA/03-063-Pre_GTGAAA.bam.count.mapped.xml"],
      "03-065-Post_GTCCGC" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_smallRNA/result/03-065-Post_GTCCGC/03-065-Post_GTCCGC.bam.count.mapped.xml"],
      "03-065-Pre_GTAGAG" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_smallRNA/result/03-065-Pre_GTAGAG/03-065-Pre_GTAGAG.bam.count.mapped.xml"],
      "03-16-Post_ACTTGA" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_smallRNA/result/03-16-Post_ACTTGA/03-16-Post_ACTTGA.bam.count.mapped.xml"],
      "03-16-Pre_CAGATC" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_smallRNA/result/03-16-Pre_CAGATC/03-16-Pre_CAGATC.bam.count.mapped.xml"],
      "03-17-Post_TAGCTT" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_smallRNA/result/03-17-Post_TAGCTT/03-17-Post_TAGCTT.bam.count.mapped.xml"],
      "03-17-Pre_GATCAG" =>
        ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm_count_smallRNA/result/03-17-Pre_GATCAG/03-17-Pre_GATCAG.bam.count.mapped.xml"],
      "2829-KCV-1A" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201405_smallRNA_2829/topN_bowtie1_genome_cutadapt_1mm_count_smallRNA/result/2829-KCV-1A/2829-KCV-1A.bam.count.mapped.xml"],
      "2829-KCV-1B" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201405_smallRNA_2829/topN_bowtie1_genome_cutadapt_1mm_count_smallRNA/result/2829-KCV-1B/2829-KCV-1B.bam.count.mapped.xml"],
      "2829-KCV-1C" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201405_smallRNA_2829/topN_bowtie1_genome_cutadapt_1mm_count_smallRNA/result/2829-KCV-1C/2829-KCV-1C.bam.count.mapped.xml"],
      "2829-KCV-1D" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201405_smallRNA_2829/topN_bowtie1_genome_cutadapt_1mm_count_smallRNA/result/2829-KCV-1D/2829-KCV-1D.bam.count.mapped.xml"],
      "2829-KCV-1E" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201405_smallRNA_2829/topN_bowtie1_genome_cutadapt_1mm_count_smallRNA/result/2829-KCV-1E/2829-KCV-1E.bam.count.mapped.xml"],
      "2829-KCV-1F" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201405_smallRNA_2829/topN_bowtie1_genome_cutadapt_1mm_count_smallRNA/result/2829-KCV-1F/2829-KCV-1F.bam.count.mapped.xml"],
      "2829-KCV-1G" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201405_smallRNA_2829/topN_bowtie1_genome_cutadapt_1mm_count_smallRNA/result/2829-KCV-1G/2829-KCV-1G.bam.count.mapped.xml"],
      "2829-KCV-1H" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201405_smallRNA_2829/topN_bowtie1_genome_cutadapt_1mm_count_smallRNA/result/2829-KCV-1H/2829-KCV-1H.bam.count.mapped.xml"],
      "2829-KCV-1I" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201405_smallRNA_2829/topN_bowtie1_genome_cutadapt_1mm_count_smallRNA/result/2829-KCV-1I/2829-KCV-1I.bam.count.mapped.xml"],
      "2829-KCV-1J" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201405_smallRNA_2829/topN_bowtie1_genome_cutadapt_1mm_count_smallRNA/result/2829-KCV-1J/2829-KCV-1J.bam.count.mapped.xml"],
    },
    cqs_tools => $cqstools,
    prefix    => "smallRNA_1mm_",
    sh_direct => 1,
    pbs       => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "10",
      "mem"      => "10gb"
    },
  },
  groups => {
    "2829" => [ "2829-KCV-1A", "2829-KCV-1B", "2829-KCV-1C", "2829-KCV-1D", "2829-KCV-1E", "2829-KCV-1F", "2829-KCV-1G", "2829-KCV-1H", "2829-KCV-1I", "2829-KCV-1J" ],

    #cqstools file_def -i . -r -f .+Pre.+xml$ -n \(.+-Pre\) -g \(Pre\)
    "2570PRE" => [
      "01-018-Pre", "01-031-Pre", "01-061-Pre", "01-28-Pre",  "01-29-Pre",  "01-36-Pre",  "03-007-Pre", "03-011-Pre", "03-015-Pre", "03-018-Pre",
      "03-026-Pre", "03-031-Pre", "03-033-Pre", "03-036-Pre", "03-047-Pre", "03-049-Pre", "03-063-Pre", "03-065-Pre", "03-16-Pre",  "03-17-Pre"
    ],
  },
  pairs => {
    "2819_VS_2570PRE" => { groups => [ "2829", "2570PRE" ], },
  },
  deseq2 => {
    class         => "DESeq2",
    perform       => 1,
    target_dir    => "${target_dir}/deseq2",
    option        => "",
    source_ref    => "pairs",
    groups_ref    => "groups",
    countfile_ref => "smallRNA_1mm_table",
    sh_direct     => 1,
    pbs           => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "10",
      "mem"      => "10gb"
    },
  },
};

performConfig($config);

1;
