#!/usr/bin/perl
use strict;
use warnings;

use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;

my $target_dir = create_directory_or_die("/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201405_smallRNA_2829_2570/");
my $cqstools   = "/home/shengq1/cqstools/CQS.Tools.exe";
my $samtools   = "/home/shengq1/local/bin/samtools/samtools";

my $email     = "quanhu.sheng\@vanderbilt.edu";
my $task_name = "Vicky_2829_2570";

my $bowtie1_option_pm = "-a -m 100 --best --strata -v 0 -l 12 -p 8";

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
    perform    => 0,
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

    #cqstools file_def -i . -r -f .+Pre.+xml$ -n \(.+\).bam -g \(Pre\)
    "2570PRE" => [
      "01-018-Pre_GGCTAC", "01-031-Pre_GAGTGG", "01-061-Pre_ACTGAT", "01-28-Pre_ATCACG",  "01-29-Pre_TGACCA",  "01-29-Pre_TTAGGC",  "01-36-Pre_ACAGTG",  "03-007-Pre_ATTCCT",
      "03-011-Pre_CAACTA", "03-015-Pre_CACGAT", "03-018-Pre_GTTTCG", "03-026-Pre_AGTCAA", "03-031-Pre_CAGGCG", "03-033-Pre_CATTTT", "03-036-Pre_ATGTCA", "03-047-Pre_CGGAAT",
      "03-049-Pre_CTATAC", "03-063-Pre_GTGAAA", "03-065-Pre_GTAGAG", "03-16-Pre_CAGATC",  "03-17-Pre_GATCAG"
    ],
    "2570POST_Placebo" => [
      "01-28-Post_CGATGT",  "01-36-Post_GCCAAT",  "03-16-Post_ACTTGA", "03-17-Post_TAGCTT", "03-018-Post_CGTACG", "03-026-Post_AGTTCC",
      "03-036-Post_CCGTCC", "03-063-Post_GTGGCC", "03-065-Post_GTCCGC"
    ],
    "2570POST_Colesevelam" => [
      "01-018-Post_CTTGTA", "01-031-Post_GGTAGC", "01-061-Post_ATGAGC", "03-007-Post_CAAAAG", "03-011-Post_CACCGG", "03-015-Post_CACTCA",
      "03-031-Post_CATGGC", "03-033-Post_CCAACA", "03-047-Post_CTAGCT", "03-049-Post_CTCAGA"
    ],
  },
  pairs => {
    "2819_VS_2570PRE" => { groups => [ "2829", "2570PRE" ], },

    #"2819_VS_2570"    => { groups => [ "2829", "2570PRE", "2570POST_Placebo", "2570POST_Colesevelam" ], },
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

  #2 perfect match search
  chromosome_count => {
    class      => "CQS::CQSChromosomeCount",
    perform    => 1,
    target_dir => "${target_dir}/chromoeome_count",
    option     => "",
    source     => {
      "01-18-Post"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/topN_bowtie1_genome_cutadapt_miRbase_pm/result/01-18-Post/01-18-Post.bam"],
      "01-18-Pre"   => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/topN_bowtie1_genome_cutadapt_miRbase_pm/result/01-18-Pre/01-18-Pre.bam"],
      "01-28-Post"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/topN_bowtie1_genome_cutadapt_miRbase_pm/result/01-28-Post/01-28-Post.bam"],
      "01-28-Pre"   => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/topN_bowtie1_genome_cutadapt_miRbase_pm/result/01-28-Pre/01-28-Pre.bam"],
      "01-29-Post"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/topN_bowtie1_genome_cutadapt_miRbase_pm/result/01-29-Post/01-29-Post.bam"],
      "01-29-Pre"   => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/topN_bowtie1_genome_cutadapt_miRbase_pm/result/01-29-Pre/01-29-Pre.bam"],
      "01-31-Post"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/topN_bowtie1_genome_cutadapt_miRbase_pm/result/01-31-Post/01-31-Post.bam"],
      "01-31-Pre"   => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/topN_bowtie1_genome_cutadapt_miRbase_pm/result/01-31-Pre/01-31-Pre.bam"],
      "01-36-Post"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/topN_bowtie1_genome_cutadapt_miRbase_pm/result/01-36-Post/01-36-Post.bam"],
      "01-36-Pre"   => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/topN_bowtie1_genome_cutadapt_miRbase_pm/result/01-36-Pre/01-36-Pre.bam"],
      "01-61-Post"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/topN_bowtie1_genome_cutadapt_miRbase_pm/result/01-61-Post/01-61-Post.bam"],
      "01-61-Pre"   => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/topN_bowtie1_genome_cutadapt_miRbase_pm/result/01-61-Pre/01-61-Pre.bam"],
      "03-07-Post"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/topN_bowtie1_genome_cutadapt_miRbase_pm/result/03-07-Post/03-07-Post.bam"],
      "03-07-Pre"   => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/topN_bowtie1_genome_cutadapt_miRbase_pm/result/03-07-Pre/03-07-Pre.bam"],
      "03-11-Post"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/topN_bowtie1_genome_cutadapt_miRbase_pm/result/03-11-Post/03-11-Post.bam"],
      "03-11-Pre"   => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/topN_bowtie1_genome_cutadapt_miRbase_pm/result/03-11-Pre/03-11-Pre.bam"],
      "03-15-Post"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/topN_bowtie1_genome_cutadapt_miRbase_pm/result/03-15-Post/03-15-Post.bam"],
      "03-15-Pre"   => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/topN_bowtie1_genome_cutadapt_miRbase_pm/result/03-15-Pre/03-15-Pre.bam"],
      "03-16-Post"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/topN_bowtie1_genome_cutadapt_miRbase_pm/result/03-16-Post/03-16-Post.bam"],
      "03-16-Pre"   => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/topN_bowtie1_genome_cutadapt_miRbase_pm/result/03-16-Pre/03-16-Pre.bam"],
      "03-17-Post"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/topN_bowtie1_genome_cutadapt_miRbase_pm/result/03-17-Post/03-17-Post.bam"],
      "03-17-Pre"   => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/topN_bowtie1_genome_cutadapt_miRbase_pm/result/03-17-Pre/03-17-Pre.bam"],
      "03-18-Post"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/topN_bowtie1_genome_cutadapt_miRbase_pm/result/03-18-Post/03-18-Post.bam"],
      "03-18-Pre"   => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/topN_bowtie1_genome_cutadapt_miRbase_pm/result/03-18-Pre/03-18-Pre.bam"],
      "03-26-Post"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/topN_bowtie1_genome_cutadapt_miRbase_pm/result/03-26-Post/03-26-Post.bam"],
      "03-26-Pre"   => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/topN_bowtie1_genome_cutadapt_miRbase_pm/result/03-26-Pre/03-26-Pre.bam"],
      "03-31-Post"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/topN_bowtie1_genome_cutadapt_miRbase_pm/result/03-31-Post/03-31-Post.bam"],
      "03-31-Pre"   => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/topN_bowtie1_genome_cutadapt_miRbase_pm/result/03-31-Pre/03-31-Pre.bam"],
      "03-33-Post"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/topN_bowtie1_genome_cutadapt_miRbase_pm/result/03-33-Post/03-33-Post.bam"],
      "03-33-Pre"   => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/topN_bowtie1_genome_cutadapt_miRbase_pm/result/03-33-Pre/03-33-Pre.bam"],
      "03-36-Post"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/topN_bowtie1_genome_cutadapt_miRbase_pm/result/03-36-Post/03-36-Post.bam"],
      "03-36-Pre"   => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/topN_bowtie1_genome_cutadapt_miRbase_pm/result/03-36-Pre/03-36-Pre.bam"],
      "03-47-Post"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/topN_bowtie1_genome_cutadapt_miRbase_pm/result/03-47-Post/03-47-Post.bam"],
      "03-47-Pre"   => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/topN_bowtie1_genome_cutadapt_miRbase_pm/result/03-47-Pre/03-47-Pre.bam"],
      "03-49-Post"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/topN_bowtie1_genome_cutadapt_miRbase_pm/result/03-49-Post/03-49-Post.bam"],
      "03-49-Pre"   => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/topN_bowtie1_genome_cutadapt_miRbase_pm/result/03-49-Pre/03-49-Pre.bam"],
      "03-63-Post"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/topN_bowtie1_genome_cutadapt_miRbase_pm/result/03-63-Post/03-63-Post.bam"],
      "03-63-Pre"   => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/topN_bowtie1_genome_cutadapt_miRbase_pm/result/03-63-Pre/03-63-Pre.bam"],
      "03-65-Post"  => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/topN_bowtie1_genome_cutadapt_miRbase_pm/result/03-65-Post/03-65-Post.bam"],
      "03-65-Pre"   => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/topN_bowtie1_genome_cutadapt_miRbase_pm/result/03-65-Pre/03-65-Pre.bam"],
      "2829-KCV-1A" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201405_smallRNA_2829/topN_bowtie1_genome_cutadapt_miRbase_pm/result/2829-KCV-1A/2829-KCV-1A.bam"],
      "2829-KCV-1B" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201405_smallRNA_2829/topN_bowtie1_genome_cutadapt_miRbase_pm/result/2829-KCV-1B/2829-KCV-1B.bam"],
      "2829-KCV-1C" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201405_smallRNA_2829/topN_bowtie1_genome_cutadapt_miRbase_pm/result/2829-KCV-1C/2829-KCV-1C.bam"],
      "2829-KCV-1D" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201405_smallRNA_2829/topN_bowtie1_genome_cutadapt_miRbase_pm/result/2829-KCV-1D/2829-KCV-1D.bam"],
      "2829-KCV-1E" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201405_smallRNA_2829/topN_bowtie1_genome_cutadapt_miRbase_pm/result/2829-KCV-1E/2829-KCV-1E.bam"],
      "2829-KCV-1F" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201405_smallRNA_2829/topN_bowtie1_genome_cutadapt_miRbase_pm/result/2829-KCV-1F/2829-KCV-1F.bam"],
      "2829-KCV-1G" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201405_smallRNA_2829/topN_bowtie1_genome_cutadapt_miRbase_pm/result/2829-KCV-1G/2829-KCV-1G.bam"],
      "2829-KCV-1H" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201405_smallRNA_2829/topN_bowtie1_genome_cutadapt_miRbase_pm/result/2829-KCV-1H/2829-KCV-1H.bam"],
      "2829-KCV-1I" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201405_smallRNA_2829/topN_bowtie1_genome_cutadapt_miRbase_pm/result/2829-KCV-1I/2829-KCV-1I.bam"],
      "2829-KCV-1J" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201405_smallRNA_2829/topN_bowtie1_genome_cutadapt_miRbase_pm/result/2829-KCV-1J/2829-KCV-1J.bam"],
    },
    seqcount => {
      "01-18-Post"                    => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/identical/result/01-18-Post_clipped_identical.fastq.dupcount"],
      "01-18-Pre"                     => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/identical/result/01-18-Pre_clipped_identical.fastq.dupcount"],
      "01-28-Post"                    => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/identical/result/01-28-Post_clipped_identical.fastq.dupcount"],
      "01-28-Pre"                     => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/identical/result/01-28-Pre_clipped_identical.fastq.dupcount"],
      "01-29-Post"                    => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/identical/result/01-29-Post_clipped_identical.fastq.dupcount"],
      "01-29-Pre"                     => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/identical/result/01-29-Pre_clipped_identical.fastq.dupcount"],
      "01-31-Post"                    => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/identical/result/01-31-Post_clipped_identical.fastq.dupcount"],
      "01-31-Pre"                     => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/identical/result/01-31-Pre_clipped_identical.fastq.dupcount"],
      "01-36-Post"                    => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/identical/result/01-36-Post_clipped_identical.fastq.dupcount"],
      "01-36-Pre"                     => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/identical/result/01-36-Pre_clipped_identical.fastq.dupcount"],
      "01-61-Post"                    => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/identical/result/01-61-Post_clipped_identical.fastq.dupcount"],
      "01-61-Pre"                     => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/identical/result/01-61-Pre_clipped_identical.fastq.dupcount"],
      "03-07-Post"                    => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/identical/result/03-07-Post_clipped_identical.fastq.dupcount"],
      "03-07-Pre"                     => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/identical/result/03-07-Pre_clipped_identical.fastq.dupcount"],
      "03-11-Post"                    => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/identical/result/03-11-Post_clipped_identical.fastq.dupcount"],
      "03-11-Pre"                     => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/identical/result/03-11-Pre_clipped_identical.fastq.dupcount"],
      "03-15-Post"                    => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/identical/result/03-15-Post_clipped_identical.fastq.dupcount"],
      "03-15-Pre"                     => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/identical/result/03-15-Pre_clipped_identical.fastq.dupcount"],
      "03-16-Post"                    => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/identical/result/03-16-Post_clipped_identical.fastq.dupcount"],
      "03-16-Pre"                     => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/identical/result/03-16-Pre_clipped_identical.fastq.dupcount"],
      "03-17-Post"                    => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/identical/result/03-17-Post_clipped_identical.fastq.dupcount"],
      "03-17-Pre"                     => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/identical/result/03-17-Pre_clipped_identical.fastq.dupcount"],
      "03-18-Post"                    => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/identical/result/03-18-Post_clipped_identical.fastq.dupcount"],
      "03-18-Pre"                     => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/identical/result/03-18-Pre_clipped_identical.fastq.dupcount"],
      "03-26-Post"                    => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/identical/result/03-26-Post_clipped_identical.fastq.dupcount"],
      "03-26-Pre"                     => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/identical/result/03-26-Pre_clipped_identical.fastq.dupcount"],
      "03-31-Post"                    => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/identical/result/03-31-Post_clipped_identical.fastq.dupcount"],
      "03-31-Pre"                     => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/identical/result/03-31-Pre_clipped_identical.fastq.dupcount"],
      "03-33-Post"                    => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/identical/result/03-33-Post_clipped_identical.fastq.dupcount"],
      "03-33-Pre"                     => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/identical/result/03-33-Pre_clipped_identical.fastq.dupcount"],
      "03-36-Post"                    => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/identical/result/03-36-Post_clipped_identical.fastq.dupcount"],
      "03-36-Pre"                     => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/identical/result/03-36-Pre_clipped_identical.fastq.dupcount"],
      "03-47-Post"                    => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/identical/result/03-47-Post_clipped_identical.fastq.dupcount"],
      "03-47-Pre"                     => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/identical/result/03-47-Pre_clipped_identical.fastq.dupcount"],
      "03-49-Post"                    => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/identical/result/03-49-Post_clipped_identical.fastq.dupcount"],
      "03-49-Pre"                     => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/identical/result/03-49-Pre_clipped_identical.fastq.dupcount"],
      "03-63-Post"                    => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/identical/result/03-63-Post_clipped_identical.fastq.dupcount"],
      "03-63-Pre"                     => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/identical/result/03-63-Pre_clipped_identical.fastq.dupcount"],
      "03-65-Post"                    => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/identical/result/03-65-Post_clipped_identical.fastq.dupcount"],
      "03-65-Pre"                     => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201408_smallRNA_Human_HDL_N10/identical/result/03-65-Pre_clipped_identical.fastq.dupcount"],
      "2829-KCV-1A_clipped_identical" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201405_smallRNA_2829/identical/result/2829-KCV-1A_clipped_identical.dupcount"],
      "2829-KCV-1B_clipped_identical" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201405_smallRNA_2829/identical/result/2829-KCV-1B_clipped_identical.dupcount"],
      "2829-KCV-1C_clipped_identical" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201405_smallRNA_2829/identical/result/2829-KCV-1C_clipped_identical.dupcount"],
      "2829-KCV-1D_clipped_identical" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201405_smallRNA_2829/identical/result/2829-KCV-1D_clipped_identical.dupcount"],
      "2829-KCV-1E_clipped_identical" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201405_smallRNA_2829/identical/result/2829-KCV-1E_clipped_identical.dupcount"],
      "2829-KCV-1F_clipped_identical" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201405_smallRNA_2829/identical/result/2829-KCV-1F_clipped_identical.dupcount"],
      "2829-KCV-1G_clipped_identical" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201405_smallRNA_2829/identical/result/2829-KCV-1G_clipped_identical.dupcount"],
      "2829-KCV-1H_clipped_identical" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201405_smallRNA_2829/identical/result/2829-KCV-1H_clipped_identical.dupcount"],
      "2829-KCV-1I_clipped_identical" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201405_smallRNA_2829/identical/result/2829-KCV-1I_clipped_identical.dupcount"],
      "2829-KCV-1J_clipped_identical" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201405_smallRNA_2829/identical/result/2829-KCV-1J_clipped_identical.dupcount"],
    },
    cqs_tools => $cqstools,
    samtools  => $samtools,
    pbs       => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },

};

#performConfig($config);
performTask( $config, "chromosome_count" );

1;
