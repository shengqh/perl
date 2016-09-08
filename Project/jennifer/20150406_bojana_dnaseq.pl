#!/usr/bin/perl
use strict;
use warnings;

use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;

my $target_dir = "/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes";
my $cqstools   = "/home/shengq1/cqstools/CQS.Tools.exe";
my $email      = "quanhu.sheng\@vanderbilt.edu";

my $bwa_fasta      = "/scratch/cqs/shengq1/references/hg19_16569_MT/bwa_index_0.7.12/hg19_16569_MT.fa";
my $transcript_gtf = "/scratch/cqs/shengq1/references/ensembl_gtf/v75/Homo_sapiens.GRCh37.75.MT.gtf";
my $name_map_file  = "/scratch/cqs/shengq1/references/ensembl_gtf/v75/Homo_sapiens.GRCh37.75.MT.map";

my $dbsnp  = "/scratch/cqs/shengq1/references/dbsnp/human_GRCh37_v142_16569_MT.vcf";
my $hapmap = "/scratch/cqs/shengq1/references/gatk/b37/hapmap_3.3.b37.vcf";
my $omni   = "/scratch/cqs/shengq1/references/gatk/b37/1000G_omni2.5.b37.vcf";
my $g1000  = "/scratch/cqs/shengq1/references/gatk/b37/1000G_phase1.snps.high_confidence.b37.vcf";
my $mills  = "/scratch/cqs/shengq1/references/gatk/b37/Mills_and_1000G_gold_standard.indels.b37.vcf";

my $cosmic    = "/scratch/cqs/shengq1/references/cosmic/cosmic_v71_hg19_16569_M.vcf";
my $mutect    = "/home/shengq1/local/bin/mutect-1.1.7.jar";
my $affy_file = "/data/cqs/shengq1/reference/affy/HG-U133_Plus_2.na33.annot.csv";

my $annovar_param = "-protocol refGene,snp138,cosmic70 -operation g,f,f --remove";
my $annovar_db    = "/scratch/cqs/shengq1/references/annovar/humandb/";
my $gatk_jar      = "/home/shengq1/local/bin/GATK/GenomeAnalysisTK.jar";
my $picard_jar    = "/scratch/cqs/shengq1/local/bin/picard/picard.jar";
my $conifer       = "/home/shengq1/pylibs/bin/conifer.py";
my $covered_bed   = "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/0699701_Covered.bed";

my $qc3_perl = "/scratch/cqs/shengq1/local/bin/qc3/qc3.pl";

my $cluster = "slurm";

my $config = {
  general => { task_name => "tnbc_dnaseq" },

  files => {
    "SA001" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_1/1_CATCAAGT_L007_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_1/1_CATCAAGT_L007_R2_001.fastq.gz"
    ],
    "SA003BP" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_3BP/3BP_TCTTCACA_L008_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_3BP/3BP_TCTTCACA_L008_R2_001.fastq.gz"
    ],
    "SA004BP" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_4BP/4BP_TGAAGAGA_L008_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_4BP/4BP_TGAAGAGA_L008_R2_001.fastq.gz"
    ],
    "SA004P" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_4-P/4-P_TGGCTTCA_L008_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_4-P/4-P_TGGCTTCA_L008_R2_001.fastq.gz"
    ],
    "SA005" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_5/5_CAGATCTG_L007_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_5/5_CAGATCTG_L007_R2_001.fastq.gz"
    ],
    "SA006" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_6/6_ACATTGGC_L007_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_6/6_ACATTGGC_L007_R2_001.fastq.gz"
    ],
    "SA007" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_7/7_CTGTAGCC_L007_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_7/7_CTGTAGCC_L007_R2_001.fastq.gz"
    ],
    "SA007P" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_7-P/7-P_CCTAATCC_L007_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_7-P/7-P_CCTAATCC_L007_R2_001.fastq.gz"
    ],
    "SA009" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_9/9_CGACACAC_L008_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_9/9_CGACACAC_L008_R2_001.fastq.gz"
    ],
    "SA010" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_10/10_CCTCTATC_L008_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_10/10_CCTCTATC_L008_R2_001.fastq.gz"
    ],
    "SA011" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_11/11_ACCACTGT_L007_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_11/11_ACCACTGT_L007_R2_001.fastq.gz"
    ],
    "SA012" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_12/12_AGTGGTCA_L007_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_12/12_AGTGGTCA_L007_R2_001.fastq.gz"
    ],
    "SA014" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_14/14_AACGTGAT_L008_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_14/14_AACGTGAT_L008_R2_001.fastq.gz"
    ],
    "SA016" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_16/16_AAACATCG_L008_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_16/16_AAACATCG_L008_R2_001.fastq.gz"
    ],
    "SA017" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_17/17_ATGCCTAA_L007_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_17/17_ATGCCTAA_L007_R2_001.fastq.gz"
    ],
    "SA018" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_18/18_CGGATTGC_L008_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_18/18_CGGATTGC_L008_R2_001.fastq.gz"
    ],
    "SA019" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_19/19_TGGTGGTA_L008_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_19/19_TGGTGGTA_L008_R2_001.fastq.gz"
    ],
    "SA022" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_22/22_AGCCATGC_L008_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_22/22_AGCCATGC_L008_R2_001.fastq.gz"
    ],
    "SA022P" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_22-P/22-P_TTCACGCA_L008_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_22-P/22-P_TTCACGCA_L008_R2_001.fastq.gz"
    ],
    "SA023" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_23/23_AACTCACC_L008_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_23/23_AACTCACC_L008_R2_001.fastq.gz"
    ],
    "SA026" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_26/26_AAACATCG_L007_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_26/26_AAACATCG_L007_R2_001.fastq.gz"
    ],
    "SA026P" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_26-P/26-P_AAGAGATC_L008_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_26-P/26-P_AAGAGATC_L008_R2_001.fastq.gz"
    ],
    "SA029" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_29/29_AACGTGAT_L007_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_29/29_AACGTGAT_L007_R2_001.fastq.gz"
    ],
    "SA033" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_33/33_ACAAGCTA_L007_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_33/33_ACAAGCTA_L007_R2_001.fastq.gz"
    ],
    "SA036" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_36/36_CGCTGATC_L007_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_36/36_CGCTGATC_L007_R2_001.fastq.gz"
    ],
    "SA036P" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_36-P/36-P_AGGCTAAC_L008_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_36-P/36-P_AGGCTAAC_L008_R2_001.fastq.gz"
    ],
    "SA040" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_40/40_CAAGACTA_L007_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_40/40_CAAGACTA_L007_R2_001.fastq.gz"
    ],
    "SA040P" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_40-P/40-P_AAGGACAC_L008_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_40-P/40-P_AAGGACAC_L008_R2_001.fastq.gz"
    ],
    "SA041" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_41/41_AGTACAAG_L007_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_41/41_AGTACAAG_L007_R2_001.fastq.gz"
    ],
    "SA042" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_42/42_AACAACCA_L007_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_42/42_AACAACCA_L007_R2_001.fastq.gz"
    ],
    "SA043" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_43/43_ATAGCGAC_L007_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_43/43_ATAGCGAC_L007_R2_001.fastq.gz"
    ],
    "SA044" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_44/44_AACCGAGA_L007_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_44/44_AACCGAGA_L007_R2_001.fastq.gz"
    ],
    "SA045" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_45/45_CAATGGAA_L007_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_45/45_CAATGGAA_L007_R2_001.fastq.gz"
    ],
    "SA047" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_47/47_AACGCTTA_L007_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_47/47_AACGCTTA_L007_R2_001.fastq.gz"
    ],
    "SA048" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_48/48_GACAGTGC_L008_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_48/48_GACAGTGC_L008_R2_001.fastq.gz"
    ],
    "SA050" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_50/50_CACTTCGA_L007_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_50/50_CACTTCGA_L007_R2_001.fastq.gz"
    ],
    "SA051" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_51/51_CAGCGTTA_L007_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_51/51_CAGCGTTA_L007_R2_001.fastq.gz"
    ],
    "SA052" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_52/52_CATACCAA_L007_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_52/52_CATACCAA_L007_R2_001.fastq.gz"
    ],
    "SA054" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_54/54_AATCCGTC_L008_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_54/54_AATCCGTC_L008_R2_001.fastq.gz"
    ],
    "SA055" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_55/55_AAGACGGA_L007_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_55/55_AAGACGGA_L007_R2_001.fastq.gz"
    ],
    "SA056" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_56/56_AAGGTACA_L007_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_56/56_AAGGTACA_L007_R2_001.fastq.gz"
    ],
    "SA057" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_57/57_ACACAGAA_L007_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_57/57_ACACAGAA_L007_R2_001.fastq.gz"
    ],
    "SA060" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_60/60_CCAGTTCA_L007_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_60/60_CCAGTTCA_L007_R2_001.fastq.gz"
    ],
    "SA061" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_61/61_CTAAGGTC_L008_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_61/61_CTAAGGTC_L008_R2_001.fastq.gz"
    ],
    "SA062" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_62/62_ATCATTCC_L007_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_62/62_ATCATTCC_L007_R2_001.fastq.gz"
    ],
    "SA062P" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_62-P/62-P_AATGTTGC_L008_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_62-P/62-P_AATGTTGC_L008_R2_001.fastq.gz"
    ],
    "SA063" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_63/63_ACAGCAGA_L007_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_63/63_ACAGCAGA_L007_R2_001.fastq.gz"
    ],
    "SA064" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_64/64_CCGAAGTA_L007_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_64/64_CCGAAGTA_L007_R2_001.fastq.gz"
    ],
    "SA065" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_65/65_CAAGGAGC_L007_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_65/65_CAAGGAGC_L007_R2_001.fastq.gz"
    ],
    "SA066" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_66/66_ACCTCCAA_L007_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_66/66_ACCTCCAA_L007_R2_001.fastq.gz"
    ],
    "SA068" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_68/68_ACGCTCGA_L007_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_68/68_ACGCTCGA_L007_R2_001.fastq.gz"
    ],
    "SA068P" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_68-P/68-P_ACACGACC_L008_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_68-P/68-P_ACACGACC_L008_R2_001.fastq.gz"
    ],
    "SA069" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_69/69_ACGTATCA_L007_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_69/69_ACGTATCA_L007_R2_001.fastq.gz"
    ],
    "SA070" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_70/70_GAACAGGC_L008_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_70/70_GAACAGGC_L008_R2_001.fastq.gz"
    ],
    "SA071" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_71/71_ACTATGCA_L007_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_71/71_ACTATGCA_L007_R2_001.fastq.gz"
    ],
    "SA073" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_73/73_AGAGTCAA_L007_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_73/73_AGAGTCAA_L007_R2_001.fastq.gz"
    ],
    "SA074" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_74/74_AGATCGCA_L007_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_74/74_AGATCGCA_L007_R2_001.fastq.gz"
    ],
    "SA075" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_75/75_GAGTTAGC_L008_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_75/75_GAGTTAGC_L008_R2_001.fastq.gz"
    ],
    "SA077" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_77/77_AGCAGGAA_L007_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_77/77_AGCAGGAA_L007_R2_001.fastq.gz"
    ],
    "SA081" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_81/81_ATTGGCTC_L007_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_81/81_ATTGGCTC_L007_R2_001.fastq.gz"
    ],
    "SA082" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_82/82_AGTCACTA_L007_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_82/82_AGTCACTA_L007_R2_001.fastq.gz"
    ],
    "SA083" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_83/83_ATCCTGTA_L007_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_83/83_ATCCTGTA_L007_R2_001.fastq.gz"
    ],
    "SA084" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_84/84_ATTGAGGA_L007_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_84/84_ATTGAGGA_L007_R2_001.fastq.gz"
    ],
    "SA087BP" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_87BP/87BP_TGGAACAA_L008_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_87BP/87BP_TGGAACAA_L008_R2_001.fastq.gz"
    ],
    "SA089" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_89/89_CCGTGAGA_L007_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_89/89_CCGTGAGA_L007_R2_001.fastq.gz"
    ],
    "SA091" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_91/91_CACCTTAC_L007_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_91/91_CACCTTAC_L007_R2_001.fastq.gz"
    ],
    "SA092" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_92/92_CCTCCTGA_L007_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_92/92_CCTCCTGA_L007_R2_001.fastq.gz"
    ],
    "SA093" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_93/93_CGAACTTA_L007_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_93/93_CGAACTTA_L007_R2_001.fastq.gz"
    ],
    "SA094" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_94/94_CGACTGGA_L007_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_94/94_CGACTGGA_L007_R2_001.fastq.gz"
    ],
    "SA095" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_95/95_CAACCACA_L007_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_95/95_CAACCACA_L007_R2_001.fastq.gz"
    ],
    "SA097" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_97/97_CGCATACA_L008_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_97/97_CGCATACA_L008_R2_001.fastq.gz"
    ],
    "SA098" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_98/98_CTCAATGA_L008_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_98/98_CTCAATGA_L008_R2_001.fastq.gz"
    ],
    "SA100" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_100/100_GATGAATC_L008_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_100/100_GATGAATC_L008_R2_001.fastq.gz"
    ],
    "SA102" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_102/102_CCATCCTC_L007_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_102/102_CCATCCTC_L007_R2_001.fastq.gz"
    ],
    "SA102P" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_102-P/102-P_ACAGATTC_L008_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_102-P/102-P_ACAGATTC_L008_R2_001.fastq.gz"
    ],
    "SA103" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_103/103_CTGAGCCA_L008_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_103/103_CTGAGCCA_L008_R2_001.fastq.gz"
    ],
    "SA104" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_104/104_CTGGCATA_L008_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_104/104_CTGGCATA_L008_R2_001.fastq.gz"
    ],
    "SA105" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_105/105_GAATCTGA_L008_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_105/105_GAATCTGA_L008_R2_001.fastq.gz"
    ],
    "SA106" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_106/106_GACTAGTA_L008_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_106/106_GACTAGTA_L008_R2_001.fastq.gz"
    ],
    "SA106P" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_106-P/106-P_AGATGTAC_L008_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_106-P/106-P_AGATGTAC_L008_R2_001.fastq.gz"
    ],
    "SA107" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_107/107_CCGACAAC_L007_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_107/107_CCGACAAC_L007_R2_001.fastq.gz"
    ],
    "SA111" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_111/111_GAGCTGAA_L008_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_111/111_GAGCTGAA_L008_R2_001.fastq.gz"
    ],
    "SA112" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_112/112_GATAGACA_L008_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_112/112_GATAGACA_L008_R2_001.fastq.gz"
    ],
    "SA114" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_114/114_GCCACATA_L008_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_114/114_GCCACATA_L008_R2_001.fastq.gz"
    ],
    "SA118" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_118/118_GCGAGTAA_L008_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_118/118_GCGAGTAA_L008_R2_001.fastq.gz"
    ],
    "SA118P" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_118-P/118-P_AGCACCTC_L008_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_118-P/118-P_AGCACCTC_L008_R2_001.fastq.gz"
    ],
    "SA124" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_124/124_GCTAACGA_L008_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_124/124_GCTAACGA_L008_R2_001.fastq.gz"
    ],
    "SA125" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_125/125_GCTCGGTA_L008_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_125/125_GCTCGGTA_L008_R2_001.fastq.gz"
    ],
    "SA126" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_126/126_GGAGAACA_L008_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_126/126_GGAGAACA_L008_R2_001.fastq.gz"
    ],
    "SA127" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_127/127_GGTGCGAA_L008_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_127/127_GGTGCGAA_L008_R2_001.fastq.gz"
    ],
    "SA128" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_128/128_GTACGCAA_L008_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_128/128_GTACGCAA_L008_R2_001.fastq.gz"
    ],
    "SA130" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_130/130_GTCGTAGA_L008_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_130/130_GTCGTAGA_L008_R2_001.fastq.gz"
    ],
    "SA133" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_133/133_GTCTGTCA_L008_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_133/133_GTCTGTCA_L008_R2_001.fastq.gz"
    ],
    "SA136" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_136/136_GTGTTCTA_L008_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_136/136_GTGTTCTA_L008_R2_001.fastq.gz"
    ],
    "SA139" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_139/139_GCCAAGAC_L008_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_139/139_GCCAAGAC_L008_R2_001.fastq.gz"
    ],
    "SA142" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_142/142_TAGGATGA_L008_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_142/142_TAGGATGA_L008_R2_001.fastq.gz"
    ],
    "SA144" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_144/144_TATCAGCA_L008_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_144/144_TATCAGCA_L008_R2_001.fastq.gz"
    ],
    "SA145" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_145/145_TCCGTCTA_L008_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_145/145_TCCGTCTA_L008_R2_001.fastq.gz"
    ],
  },
  groups => {
    "SA007" => [ "SA007P", "SA007" ],
    "SA022" => [ "SA022P", "SA022" ],
    "SA026" => [ "SA026P", "SA026" ],
    "SA036" => [ "SA036P", "SA036" ],
    "SA040" => [ "SA040P", "SA040" ],
    "SA062" => [ "SA062P", "SA062" ],
    "SA068" => [ "SA068P", "SA068" ],
    "SA102" => [ "SA102P", "SA102" ],
    "SA106" => [ "SA106P", "SA106" ],
    "SA118" => [ "SA118P", "SA118" ],
  },

  fastqc => {
    class      => "QC::FastQC",
    perform    => 0,
    target_dir => "${target_dir}/fastqc",
    option     => "",
    source_ref => "files",
    sh_direct  => 1,
    cluster    => $cluster,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=2",
      "walltime" => "2",
      "mem"      => "40gb"
    },
  },
  fastqc_summary => {
    class      => "QC::FastQCSummary",
    perform    => 0,
    target_dir => "${target_dir}/fastqc",
    option     => "",
    cluster    => $cluster,
    cqstools   => $cqstools,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "2",
      "mem"      => "10gb"
    },
  },
  bwa => {
    class      => "Alignment::BWA",
    perform    => 1,
    target_dir => "${target_dir}/bwa",
    option     => "",
    bwa_index  => $bwa_fasta,
    source_ref => "files",
    picard_jar => $picard_jar,
    sh_direct  => 0,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  bwa_refine => {
    class       => "GATK::Refine",
    perform     => 1,
    target_dir  => "${target_dir}/bwa_refine",
    option      => "-Xmx40g",
    gatk_option => "--fix_misencoded_quality_scores",
    fasta_file  => $bwa_fasta,
    source_ref  => "bwa",
    vcf_files   => [$dbsnp],
    gatk_jar    => $gatk_jar,
    picard_jar  => $picard_jar,
    sh_direct   => 0,
    sorted      => 1,
    pbs         => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "240",
      "mem"      => "40gb"
    },
  },
  bwa_refine_hc_gvcf => {
    class       => "GATK::HaplotypeCallerGVCF",
    perform     => 1,
    target_dir  => "${target_dir}/bwa_refine_hc_gvcf",
    option      => "",
    source_ref  => "bwa_refine",
    java_option => "",
    fasta_file  => $bwa_fasta,
    dbsnp_vcf   => $dbsnp,
    gatk_jar    => $gatk_jar,
    extension   => ".gvcf",
    sh_direct   => 0,
    pbs         => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  bwa_refine_hc_gvcf_vqsr => {
    class       => "GATK::VariantFilterVQSR",
    perform     => 1,
    target_dir  => "${target_dir}/bwa_refine_hc_gvcf_vqsr",
    option      => "",
    source_ref  => "bwa_refine_hc_gvcf",
    java_option => "",
    fasta_file  => $bwa_fasta,
    dbsnp_vcf   => $dbsnp,
    hapmap_vcf  => $hapmap,
    omni_vcf    => $omni,
    g1000_vcf   => $g1000,
    mills_vcf   => $mills,
    gatk_jar    => $gatk_jar,
    cqstools    => $cqstools,
    sh_direct   => 1,
    pbs         => {
      "email"    => $email,
      "nodes"    => "1:ppn=24",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  qc3vcf => {
    class      => "QC::QC3vcf",
    perform    => 1,
    target_dir => "${target_dir}/bwa_refine_hc_gvcf_vqsr_qc3",
    option     => "",
    qc3_perl   => $qc3_perl,
    source_ref => [ "bwa_refine_hc_gvcf_vqsr", "snp" ],
    annovar_db => $annovar_db,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  bwa_refine_hc_gvcf_vqsr_annovar => {
    class      => "Annotation::Annovar",
    perform    => 1,
    target_dir => "${target_dir}/bwa_refine_hc_gvcf_vqsr_annovar",
    source_ref => [ "bwa_refine_hc_gvcf_vqsr", "snp" ],
    option     => $annovar_param,
    annovar_db => $annovar_db,
    buildver   => "hg19",
    sh_direct  => 1,
    isvcf      => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "10gb"
    },
  },
  conifer => {
    class       => "CNV::Conifer",
    perform     => 1,
    target_dir  => "${target_dir}/conifer",
    option      => "",
    source_ref  => "bwa",
    conifer     => $conifer,
    bedfile     => $covered_bed,
    isbamsorted => 1,
    sh_direct   => 1,
    pbs         => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "720",
      "mem"      => "10gb"
    },
  },
  cnmops => {
    class       => "CNV::cnMops",
    perform     => 1,
    target_dir  => "${target_dir}/cnmops",
    option      => "",
    source_ref  => "bwa",
    bedfile     => $covered_bed,
    pairmode    => "paired",
    isbamsorted => 1,
    sh_direct   => 1,
    pbs         => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "720",
      "mem"      => "40gb"
    },
  },
  muTect => {
    class        => "GATK::MuTect",
    perform      => 1,
    target_dir   => "${target_dir}/muTect",
    option       => "--min_qscore 20 --filter_reads_with_N_cigar",
    java_option  => "-Xmx40g",
    source_ref   => "bwa_refine",
    groups_ref   => "groups",
    fasta_file   => $bwa_fasta,
    cosmic_file  => $cosmic,
    dbsnp_file   => $dbsnp,
    bychromosome => 0,
    sh_direct    => 0,
    muTect_jar   => $mutect,
    pbs          => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "240",
      "mem"      => "40gb"
    },
  },
  annovar_muTect => {
    class      => "Annotation::Annovar",
    perform    => 1,
    target_dir => "${target_dir}/muTect",
    option     => $annovar_param,
    source_ref => [ "muTect", ".pass.vcf\$" ],
    annovar_db => $annovar_db,
    buildver   => "hg19",
    cqstools   => $cqstools,
    affy_file  => $affy_file,
    sh_direct  => 1,
    isvcf      => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "10gb"
    },
  },

  sequencetask => {
    class      => "CQS::SequenceTask",
    perform    => 1,
    target_dir => "${target_dir}/sequencetask",
    option     => "",
    source     => {
      step1 => ["fastqc"],
      step2 => ["fastqc_summary"],
      step3 => [ "bwa", "bwa_refine", "bwa_refine_hc_gvcf" ],
      step4 => ["bwa_refine_hc_gvcf_vqsr"],
      step5 => ["bwa_refine_hc_gvcf_vqsr_annovar"],
    },
    sh_direct => 0,
    pbs       => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
};

#performConfig($config);

#performTask( $config, "muTect" );

performTask( $config, "bwa_refine_hc_gvcf_vqsr_annovar" );

1;

