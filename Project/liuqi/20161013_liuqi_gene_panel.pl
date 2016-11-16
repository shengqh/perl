#!/usr/bin/perl
use strict;
use warnings;

use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;

my $target_dir = create_directory_or_die("/scratch/cqs/shengq1/dnaseq/20161013_liuqi_gene_panel");
my $cqstools   = "/home/shengq1/cqstools/cqstools.exe";
my $email      = "quanhu.sheng\@vanderbilt.edu";

my $bwa_fasta      = "/scratch/cqs/shengq1/references/gatk/b37/bwa_index_0.7.12/human_g1k_v37.fasta";
my $transcript_gtf = "/scratch/cqs/shengq1/references/ensembl_gtf/v75/Homo_sapiens.GRCh37.75.MT.gtf";
my $name_map_file  = "/scratch/cqs/shengq1/references/ensembl_gtf/v75/Homo_sapiens.GRCh37.75.MT.map";

my $dbsnp  = "/scratch/cqs/shengq1/references/dbsnp/human_GRCh37_v142_16569_MT.vcf";
my $hapmap = "/scratch/cqs/shengq1/references/gatk/b37/hapmap_3.3.b37.vcf";
my $omni   = "/scratch/cqs/shengq1/references/gatk/b37/1000G_omni2.5.b37.vcf";
my $g1000  = "/scratch/cqs/shengq1/references/gatk/b37/1000G_phase1.snps.high_confidence.b37.vcf";
my $mills  = "/scratch/cqs/shengq1/references/gatk/b37/Mills_and_1000G_gold_standard.indels.b37.vcf";

my $cosmic = "/scratch/cqs/shengq1/references/cosmic/cosmic_v71_hg19_16569_MT.vcf";

my $annovar_param = "-protocol refGene,avsnp147,cosmic70,exac03 -operation g,f,f,f";
my $annovar_db    = "/scratch/cqs/shengq1/references/annovar/humandb/";
my $gatk_jar      = "/home/shengq1/local/bin/GATK/GenomeAnalysisTK.jar";
my $picard_jar    = "/scratch/cqs/shengq1/local/bin/picard/picard.jar";
my $qc3_perl      = "/scratch/cqs/shengq1/local/bin/qc3/qc3.pl";

my $covered_bed = "/gpfs21/scratch/cqs/shengq1/dnaseq/20160829_liuqi_gene_panel/documents/TCPS_Labcorp_genes_20160829.bed";

my $gene_bed = "/gpfs21/scratch/cqs/shengq1/dnaseq/20160829_liuqi_gene_panel/documents/interesting_gene.bed";

my $cluster = "slurm";

my $genes = "APC ROBO2 KRAS NRAS BRAF TP53";

my $config = {
  general => { task_name => "Adenoma" },
  files   => {
    "TP0097_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233961/BHCHJNBBXX_DS-233961_AGTTCC_L005_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233961/BHCHJNBBXX_DS-233961_AGTTCC_L005_R2_001.fastq.gz"
    ],
    "TP0109_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233949/AH7YCFBBXX_DS-233949_CGATGT_L005_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233949/AH7YCFBBXX_DS-233949_CGATGT_L005_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233949/BHFLGLBBXX_DS-233949_CGATGT_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233949/BHFLGLBBXX_DS-233949_CGATGT_L007_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233949/BHFLGLBBXX_DS-233949_CGATGT_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233949/BHFLGLBBXX_DS-233949_CGATGT_L008_R2_001.fastq.gz"
    ],
    "TP0192_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233953/AH7YCFBBXX_DS-233953_GCCAAT_L005_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233953/AH7YCFBBXX_DS-233953_GCCAAT_L005_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233953/BHFLGLBBXX_DS-233953_GCCAAT_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233953/BHFLGLBBXX_DS-233953_GCCAAT_L007_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233953/BHFLGLBBXX_DS-233953_GCCAAT_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233953/BHFLGLBBXX_DS-233953_GCCAAT_L008_R2_001.fastq.gz"
    ],
    "TP0242_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233454/BHCHJNBBXX_DS-233454_ACAGTG_L002_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233454/BHCHJNBBXX_DS-233454_ACAGTG_L002_R2_001.fastq.gz"
    ],
    "TP0399_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233914/BHCHJNBBXX_DS-233914_TGACCA_L004_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233914/BHCHJNBBXX_DS-233914_TGACCA_L004_R2_001.fastq.gz"
    ],
    "TP0449_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234554/AH7YCFBBXX_DS-234554_GTTTCG_L001_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234554/AH7YCFBBXX_DS-234554_GTTTCG_L001_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234554/BHFG5FBBXX_DS-234554_GTTTCG_L002_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234554/BHFG5FBBXX_DS-234554_GTTTCG_L002_R2_001.fastq.gz"
    ],
    "TP0637_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233402/AH7W5TBBXX_DS-233402_CCGTCC_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233402/AH7W5TBBXX_DS-233402_CCGTCC_L008_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233402/AH7YCFBBXX_DS-233402_CCGTCC_L004_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233402/AH7YCFBBXX_DS-233402_CCGTCC_L004_R2_001.fastq.gz"
    ],
    "TP0779_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233879/BHCHJNBBXX_DS-233879_GTCCGC_L002_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233879/BHCHJNBBXX_DS-233879_GTCCGC_L002_R2_001.fastq.gz"
    ],
    "TP0919_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232719/AH7W5TBBXX_DS-232719_GTTTCG_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232719/AH7W5TBBXX_DS-232719_GTTTCG_L007_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232719/BH7WHCBBXX_DS-232719_GTTTCG_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232719/BH7WHCBBXX_DS-232719_GTTTCG_L007_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232719/BH7WHCBBXX_DS-232719_GTTTCG_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232719/BH7WHCBBXX_DS-232719_GTTTCG_L008_R2_001.fastq.gz"
    ],
    "TP0919_Recur_Base2" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232680/AH7W5TBBXX_DS-232680_CGATGT_L006_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232680/AH7W5TBBXX_DS-232680_CGATGT_L006_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232680/BH7WHCBBXX_DS-232680_CGATGT_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232680/BH7WHCBBXX_DS-232680_CGATGT_L007_R2_001.fastq.gz"
    ],
    "TP1089_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234525/BHCHJNBBXX_DS-234525_ACAGTG_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234525/BHCHJNBBXX_DS-234525_ACAGTG_L008_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234525/BHFG5FBBXX_DS-234525_ACAGTG_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234525/BHFG5FBBXX_DS-234525_ACAGTG_L008_R2_001.fastq.gz"
    ],
    "TP1095_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233927/BHCHJNBBXX_DS-233927_GTCCGC_L004_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233927/BHCHJNBBXX_DS-233927_GTCCGC_L004_R2_001.fastq.gz"
    ],
    "TP1095_Recur_Recur" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233920/BHCHJNBBXX_DS-233920_TAGCTT_L004_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233920/BHCHJNBBXX_DS-233920_TAGCTT_L004_R2_001.fastq.gz"
    ],
    "TP1145_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232685/AH7W5TBBXX_DS-232685_GCCAAT_L006_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232685/AH7W5TBBXX_DS-232685_GCCAAT_L006_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232685/BH7WHCBBXX_DS-232685_GCCAAT_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232685/BH7WHCBBXX_DS-232685_GCCAAT_L007_R2_001.fastq.gz"
    ],
    "TP1145_Recur_Recur" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232801/AH7W5TBBXX_DS-232801_ATCACG_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232801/AH7W5TBBXX_DS-232801_ATCACG_L008_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232801/AH7YCFBBXX_DS-232801_ATCACG_L004_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232801/AH7YCFBBXX_DS-232801_ATCACG_L004_R2_001.fastq.gz"
    ],
    "TP1166_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234533/AH7YCFBBXX_DS-234533_AGTCAA_L002_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234533/AH7YCFBBXX_DS-234533_AGTCAA_L002_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234533/BHFLGLBBXX_DS-234533_AGTCAA_L006_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234533/BHFLGLBBXX_DS-234533_AGTCAA_L006_R2_001.fastq.gz"
    ],
    "TP1166_Norecur_Base2" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234569/AH7YCFBBXX_DS-234569_GGCTAC_L006_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234569/AH7YCFBBXX_DS-234569_GGCTAC_L006_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234569/BHFLGLBBXX_DS-234569_GGCTAC_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234569/BHFLGLBBXX_DS-234569_GGCTAC_L007_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234569/BHFLGLBBXX_DS-234569_GGCTAC_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234569/BHFLGLBBXX_DS-234569_GGCTAC_L008_R2_001.fastq.gz"
    ],
    "TP1215_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233959/BHCHJNBBXX_DS-233959_CTTGTA_L005_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233959/BHCHJNBBXX_DS-233959_CTTGTA_L005_R2_001.fastq.gz"
    ],
    "TP1215_Recur_Base2" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233893/BHCHJNBBXX_DS-233893_CAGATC_L003_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233893/BHCHJNBBXX_DS-233893_CAGATC_L003_R2_001.fastq.gz"
    ],
    "TP1215_Recur_Recur" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233955/AH7YCFBBXX_DS-233955_ACTTGA_L005_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233955/AH7YCFBBXX_DS-233955_ACTTGA_L005_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233955/BHFLGLBBXX_DS-233955_ACTTGA_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233955/BHFLGLBBXX_DS-233955_ACTTGA_L007_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233955/BHFLGLBBXX_DS-233955_ACTTGA_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233955/BHFLGLBBXX_DS-233955_ACTTGA_L008_R2_001.fastq.gz"
    ],
    "TP1341_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233438/BHCHJNBBXX_DS-233438_AGTTCC_L001_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233438/BHCHJNBBXX_DS-233438_AGTTCC_L001_R2_001.fastq.gz"
    ],
    "TP1351_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233908/BHCHJNBBXX_DS-233908_GAGTGG_L003_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233908/BHCHJNBBXX_DS-233908_GAGTGG_L003_R2_001.fastq.gz"
    ],
    "TP1367_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233896/BH7WHCBBXX_DS-233896_TAGCTT_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233896/BH7WHCBBXX_DS-233896_TAGCTT_L008_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233896/BHCHJNBBXX_DS-233896_TAGCTT_L003_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233896/BHCHJNBBXX_DS-233896_TAGCTT_L003_R2_001.fastq.gz"
    ],
    "TP1367_Recur_Recur" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233903/BHCHJNBBXX_DS-233903_GTCCGC_L003_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233903/BHCHJNBBXX_DS-233903_GTCCGC_L003_R2_001.fastq.gz"
    ],
    "TP1515_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233984/BHCHJNBBXX_DS-233984_ACTGAT_L006_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233984/BHCHJNBBXX_DS-233984_ACTGAT_L006_R2_001.fastq.gz"
    ],
    "TP1622_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233439/BHCHJNBBXX_DS-233439_ATCACG_L002_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233439/BHCHJNBBXX_DS-233439_ATCACG_L002_R2_001.fastq.gz"
    ],
    "TP1626_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233888/BHCHJNBBXX_DS-233888_CGATGT_L003_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233888/BHCHJNBBXX_DS-233888_CGATGT_L003_R2_001.fastq.gz"
    ],
    "TP1626_Recur_Base2" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233947/BHCHJNBBXX_DS-233947_ATTCCT_L004_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233947/BHCHJNBBXX_DS-233947_ATTCCT_L004_R2_001.fastq.gz"
    ],
    "TP1626_Recur_Recur" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233965/AH7YCFBBXX_DS-233965_TGACCA_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233965/AH7YCFBBXX_DS-233965_TGACCA_L007_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233965/BHFLGLBBXX_DS-233965_TGACCA_L006_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233965/BHFLGLBBXX_DS-233965_TGACCA_L006_R2_001.fastq.gz"
    ],
    "TP1638_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233951/AH7YCFBBXX_DS-233951_TGACCA_L005_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233951/AH7YCFBBXX_DS-233951_TGACCA_L005_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233951/BHFLGLBBXX_DS-233951_TGACCA_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233951/BHFLGLBBXX_DS-233951_TGACCA_L007_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233951/BHFLGLBBXX_DS-233951_TGACCA_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233951/BHFLGLBBXX_DS-233951_TGACCA_L008_R2_001.fastq.gz"
    ],
    "TP1638_Norecur_Base2" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233980/BHCHJNBBXX_DS-233980_GTGGCC_L006_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233980/BHCHJNBBXX_DS-233980_GTGGCC_L006_R2_001.fastq.gz"
    ],
    "TP2263_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233957/BHCHJNBBXX_DS-233957_TAGCTT_L005_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233957/BHCHJNBBXX_DS-233957_TAGCTT_L005_R2_001.fastq.gz"
    ],
    "TP2293_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233921/BHCHJNBBXX_DS-233921_GGCTAC_L004_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233921/BHCHJNBBXX_DS-233921_GGCTAC_L004_R2_001.fastq.gz"
    ],
    "TP2479_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234520/BHCHJNBBXX_DS-234520_ATTCCT_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234520/BHCHJNBBXX_DS-234520_ATTCCT_L008_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234520/BHFG5FBBXX_DS-234520_ATTCCT_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234520/BHFG5FBBXX_DS-234520_ATTCCT_L008_R2_001.fastq.gz"
    ],
    "TP2613_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233474/BH7WHCBBXX_DS-233474_AGTCAA_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233474/BH7WHCBBXX_DS-233474_AGTCAA_L007_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233474/BHCHJNBBXX_DS-233474_AGTCAA_L002_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233474/BHCHJNBBXX_DS-233474_AGTCAA_L002_R2_001.fastq.gz"
    ],
    "TP2613_Recur_Base2" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233982/BHCHJNBBXX_DS-233982_CGTACG_L006_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233982/BHCHJNBBXX_DS-233982_CGTACG_L006_R2_001.fastq.gz"
    ],
    "TP2613_Recur_Recur" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233468/BH7WHCBBXX_DS-233468_GGCTAC_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233468/BH7WHCBBXX_DS-233468_GGCTAC_L007_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233468/BHCHJNBBXX_DS-233468_GGCTAC_L002_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233468/BHCHJNBBXX_DS-233468_GGCTAC_L002_R2_001.fastq.gz"
    ],
    "TP2794_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233395/AH7W5TBBXX_DS-233395_CTTGTA_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233395/AH7W5TBBXX_DS-233395_CTTGTA_L008_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233395/AH7YCFBBXX_DS-233395_CTTGTA_L004_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233395/AH7YCFBBXX_DS-233395_CTTGTA_L004_R2_001.fastq.gz"
    ],
    "TP2963_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233399/AH7W5TBBXX_DS-233399_AGTTCC_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233399/AH7W5TBBXX_DS-233399_AGTTCC_L008_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233399/AH7YCFBBXX_DS-233399_AGTTCC_L004_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233399/AH7YCFBBXX_DS-233399_AGTTCC_L004_R2_001.fastq.gz"
    ],
    "TP2963_Recur_Recur" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233406/BHCHJNBBXX_DS-233406_GTCCGC_L001_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233406/BHCHJNBBXX_DS-233406_GTCCGC_L001_R2_001.fastq.gz"
    ],
    "TP3017_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233944/BHCHJNBBXX_DS-233944_CGTACG_L004_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233944/BHCHJNBBXX_DS-233944_CGTACG_L004_R2_001.fastq.gz"
    ],
    "TP3017_Norecur_Base2" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233952/AH7YCFBBXX_DS-233952_ACAGTG_L005_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233952/AH7YCFBBXX_DS-233952_ACAGTG_L005_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233952/BHFLGLBBXX_DS-233952_ACAGTG_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233952/BHFLGLBBXX_DS-233952_ACAGTG_L007_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233952/BHFLGLBBXX_DS-233952_ACAGTG_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233952/BHFLGLBBXX_DS-233952_ACAGTG_L008_R2_001.fastq.gz"
    ],
    "TP3348_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233890/BHCHJNBBXX_DS-233890_TGACCA_L003_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233890/BHCHJNBBXX_DS-233890_TGACCA_L003_R2_001.fastq.gz"
    ],
    "TP3348_Recur_Recur" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233922/BHCHJNBBXX_DS-233922_CTTGTA_L004_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233922/BHCHJNBBXX_DS-233922_CTTGTA_L004_R2_001.fastq.gz"
    ],
    "TP3361_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233987/BHCHJNBBXX_DS-233987_CGATGT_L006_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233987/BHCHJNBBXX_DS-233987_CGATGT_L006_R2_001.fastq.gz"
    ],
    "TP3964_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232704/AH7YCFBBXX_DS-232704_GGCTAC_L003_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232704/AH7YCFBBXX_DS-232704_GGCTAC_L003_R2_001.fastq.gz"
    ],
    "TP3964_Recur_Base2" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232711/AH7YCFBBXX_DS-232711_AGTTCC_L003_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232711/AH7YCFBBXX_DS-232711_AGTTCC_L003_R2_001.fastq.gz"
    ],
    "TP4048_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232712/AH7YCFBBXX_DS-232712_ATGTCA_L003_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232712/AH7YCFBBXX_DS-232712_ATGTCA_L003_R2_001.fastq.gz"
    ],
    "TP4054_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234479/BHCHJNBBXX_DS-234479_CAGATC_L006_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234479/BHCHJNBBXX_DS-234479_CAGATC_L006_R2_001.fastq.gz"
    ],
    "TP4054_Recur_Base2" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234483/BHCHJNBBXX_DS-234483_GGCTAC_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234483/BHCHJNBBXX_DS-234483_GGCTAC_L007_R2_001.fastq.gz"
    ],
    "TP4054_Recur_Recur" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234529/AH7YCFBBXX_DS-234529_GATCAG_L002_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234529/AH7YCFBBXX_DS-234529_GATCAG_L002_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234529/BHFLGLBBXX_DS-234529_GATCAG_L006_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234529/BHFLGLBBXX_DS-234529_GATCAG_L006_R2_001.fastq.gz"
    ],
    "TP4113_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232784/AH7W5TBBXX_DS-232784_GGCTAC_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232784/AH7W5TBBXX_DS-232784_GGCTAC_L007_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232784/AH7YCFBBXX_DS-232784_GGCTAC_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232784/AH7YCFBBXX_DS-232784_GGCTAC_L008_R2_001.fastq.gz"
    ],
    "TP4129_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232796/AH7W5TBBXX_DS-232796_CGTACG_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232796/AH7W5TBBXX_DS-232796_CGTACG_L008_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232796/AH7YCFBBXX_DS-232796_CGTACG_L004_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232796/AH7YCFBBXX_DS-232796_CGTACG_L004_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232796/AH7YCFBBXX_DS-232796_CGTACG_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232796/AH7YCFBBXX_DS-232796_CGTACG_L008_R2_001.fastq.gz"
    ],
    "TP4190_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234487/BHCHJNBBXX_DS-234487_ATGTCA_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234487/BHCHJNBBXX_DS-234487_ATGTCA_L007_R2_001.fastq.gz"
    ],
    "TP4424_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233962/AH7YCFBBXX_DS-233962_ATCACG_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233962/AH7YCFBBXX_DS-233962_ATCACG_L007_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233962/BHFLGLBBXX_DS-233962_ATCACG_L006_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233962/BHFLGLBBXX_DS-233962_ATCACG_L006_R2_001.fastq.gz"
    ],
    "TP4830_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232737/AH7W5TBBXX_DS-232737_ACAGTG_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232737/AH7W5TBBXX_DS-232737_ACAGTG_L007_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232737/AH7YCFBBXX_DS-232737_ACAGTG_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232737/AH7YCFBBXX_DS-232737_ACAGTG_L008_R2_001.fastq.gz"
    ],
    "TP4830_Recur_Base2" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232716/AH7W5TBBXX_DS-232716_GTCCGC_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232716/AH7W5TBBXX_DS-232716_GTCCGC_L007_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232716/BH7WHCBBXX_DS-232716_GTCCGC_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232716/BH7WHCBBXX_DS-232716_GTCCGC_L007_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232716/BH7WHCBBXX_DS-232716_GTCCGC_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232716/BH7WHCBBXX_DS-232716_GTCCGC_L008_R2_001.fastq.gz"
    ],
    "TP4830_Recur_Recur" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232743/AH7W5TBBXX_DS-232743_GCCAAT_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232743/AH7W5TBBXX_DS-232743_GCCAAT_L007_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232743/AH7YCFBBXX_DS-232743_GCCAAT_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232743/AH7YCFBBXX_DS-232743_GCCAAT_L008_R2_001.fastq.gz"
    ],
    "TP4924_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234553/AH7YCFBBXX_DS-234553_GTGGCC_L001_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234553/AH7YCFBBXX_DS-234553_GTGGCC_L001_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234553/BHFG5FBBXX_DS-234553_GTGGCC_L002_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234553/BHFG5FBBXX_DS-234553_GTGGCC_L002_R2_001.fastq.gz"
    ],
    "TP4991_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232787/AH7W5TBBXX_DS-232787_AGTCAA_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232787/AH7W5TBBXX_DS-232787_AGTCAA_L007_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232787/AH7YCFBBXX_DS-232787_AGTCAA_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232787/AH7YCFBBXX_DS-232787_AGTCAA_L008_R2_001.fastq.gz"
    ],
    "TP5078_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233886/BHCHJNBBXX_DS-233886_ATTCCT_L002_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233886/BHCHJNBBXX_DS-233886_ATTCCT_L002_R2_001.fastq.gz"
    ],
    "TP5755_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234576/AH7YCFBBXX_DS-234576_GTGAAA_L006_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234576/AH7YCFBBXX_DS-234576_GTGAAA_L006_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234576/BHFLGLBBXX_DS-234576_GTGAAA_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234576/BHFLGLBBXX_DS-234576_GTGAAA_L007_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234576/BHFLGLBBXX_DS-234576_GTGAAA_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234576/BHFLGLBBXX_DS-234576_GTGAAA_L008_R2_001.fastq.gz"
    ],
    "TP5789_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233925/BHCHJNBBXX_DS-233925_ATGTCA_L004_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233925/BHCHJNBBXX_DS-233925_ATGTCA_L004_R2_001.fastq.gz"
    ],
    "TP5789_Recur_Base2" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233916/BHCHJNBBXX_DS-233916_GCCAAT_L004_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233916/BHCHJNBBXX_DS-233916_GCCAAT_L004_R2_001.fastq.gz"
    ],
    "TP5789_Recur_Recur" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233895/BH7WHCBBXX_DS-233895_GATCAG_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233895/BH7WHCBBXX_DS-233895_GATCAG_L008_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233895/BHCHJNBBXX_DS-233895_GATCAG_L003_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233895/BHCHJNBBXX_DS-233895_GATCAG_L003_R2_001.fastq.gz"
    ],
    "TP5987_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232797/AH7W5TBBXX_DS-232797_GAGTGG_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232797/AH7W5TBBXX_DS-232797_GAGTGG_L008_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232797/AH7YCFBBXX_DS-232797_GAGTGG_L004_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232797/AH7YCFBBXX_DS-232797_GAGTGG_L004_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232797/AH7YCFBBXX_DS-232797_GAGTGG_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232797/AH7YCFBBXX_DS-232797_GAGTGG_L008_R2_001.fastq.gz"
    ],
    "TP6106_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234538/AH7YCFBBXX_DS-234538_TGACCA_L002_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234538/AH7YCFBBXX_DS-234538_TGACCA_L002_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234538/BHFG5FBBXX_DS-234538_TGACCA_L001_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234538/BHFG5FBBXX_DS-234538_TGACCA_L001_R2_001.fastq.gz"
    ],
    "TP6106_Recur_Recur" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234510/BHCHJNBBXX_DS-234510_AGTTCC_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234510/BHCHJNBBXX_DS-234510_AGTTCC_L008_R2_001.fastq.gz"
    ],
    "TP6114_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233968/AH7YCFBBXX_DS-233968_CAGATC_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233968/AH7YCFBBXX_DS-233968_CAGATC_L007_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233968/BHFLGLBBXX_DS-233968_CAGATC_L006_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233968/BHFLGLBBXX_DS-233968_CAGATC_L006_R2_001.fastq.gz"
    ],
    "TP6114_Recur_Recur" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233899/BH7WHCBBXX_DS-233899_AGTCAA_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233899/BH7WHCBBXX_DS-233899_AGTCAA_L008_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233899/BHCHJNBBXX_DS-233899_AGTCAA_L003_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233899/BHCHJNBBXX_DS-233899_AGTCAA_L003_R2_001.fastq.gz"
    ],
    "TP6118_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233465/BH7WHCBBXX_DS-233465_GATCAG_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233465/BH7WHCBBXX_DS-233465_GATCAG_L007_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233465/BHCHJNBBXX_DS-233465_GATCAG_L002_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233465/BHCHJNBBXX_DS-233465_GATCAG_L002_R2_001.fastq.gz"
    ],
    "TP6120_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232736/AH7W5TBBXX_DS-232736_TGACCA_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232736/AH7W5TBBXX_DS-232736_TGACCA_L007_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232736/AH7YCFBBXX_DS-232736_TGACCA_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232736/AH7YCFBBXX_DS-232736_TGACCA_L008_R2_001.fastq.gz"
    ],
    "TP6120_Recur_Base2" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232729/AH7W5TBBXX_DS-232729_ATTCCT_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232729/AH7W5TBBXX_DS-232729_ATTCCT_L007_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232729/BH7WHCBBXX_DS-232729_ATTCCT_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232729/BH7WHCBBXX_DS-232729_ATTCCT_L007_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232729/BH7WHCBBXX_DS-232729_ATTCCT_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232729/BH7WHCBBXX_DS-232729_ATTCCT_L008_R2_001.fastq.gz"
    ],
    "TP6120_Recur_Recur" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232681/AH7W5TBBXX_DS-232681_TTAGGC_L006_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232681/AH7W5TBBXX_DS-232681_TTAGGC_L006_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232681/BH7WHCBBXX_DS-232681_TTAGGC_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232681/BH7WHCBBXX_DS-232681_TTAGGC_L007_R2_001.fastq.gz"
    ],
    "TP6161_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234542/AH7YCFBBXX_DS-234542_ACTTGA_L002_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234542/AH7YCFBBXX_DS-234542_ACTTGA_L002_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234542/BHFG5FBBXX_DS-234542_ACTTGA_L001_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234542/BHFG5FBBXX_DS-234542_ACTTGA_L001_R2_001.fastq.gz"
    ],
    "TP6242_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234513/BHCHJNBBXX_DS-234513_GTCCGC_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234513/BHCHJNBBXX_DS-234513_GTCCGC_L008_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234513/BHFG5FBBXX_DS-234513_GTCCGC_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234513/BHFG5FBBXX_DS-234513_GTCCGC_L008_R2_001.fastq.gz"
    ],
    "TP6318_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234571/AH7YCFBBXX_DS-234571_AGTCAA_L006_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234571/AH7YCFBBXX_DS-234571_AGTCAA_L006_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234571/BHFLGLBBXX_DS-234571_AGTCAA_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234571/BHFLGLBBXX_DS-234571_AGTCAA_L007_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234571/BHFLGLBBXX_DS-234571_AGTCAA_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234571/BHFLGLBBXX_DS-234571_AGTCAA_L008_R2_001.fastq.gz"
    ],
    "TP6489_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234531/AH7YCFBBXX_DS-234531_GGCTAC_L002_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234531/AH7YCFBBXX_DS-234531_GGCTAC_L002_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234531/BHFLGLBBXX_DS-234531_GGCTAC_L006_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234531/BHFLGLBBXX_DS-234531_GGCTAC_L006_R2_001.fastq.gz"
    ],
    "TP6556_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234505/BHCHJNBBXX_DS-234505_GATCAG_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234505/BHCHJNBBXX_DS-234505_GATCAG_L008_R2_001.fastq.gz"
    ],
    "TP6671_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232799/AH7W5TBBXX_DS-232799_ATTCCT_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232799/AH7W5TBBXX_DS-232799_ATTCCT_L008_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232799/AH7YCFBBXX_DS-232799_ATTCCT_L004_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232799/AH7YCFBBXX_DS-232799_ATTCCT_L004_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232799/AH7YCFBBXX_DS-232799_ATTCCT_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232799/AH7YCFBBXX_DS-232799_ATTCCT_L008_R2_001.fastq.gz"
    ],
    "TP6681_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234548/AH7YCFBBXX_DS-234548_AGTTCC_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234548/AH7YCFBBXX_DS-234548_AGTTCC_L007_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234548/BHFG5FBBXX_DS-234548_AGTTCC_L001_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234548/BHFG5FBBXX_DS-234548_AGTTCC_L001_R2_001.fastq.gz"
    ],
    "TP6933_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233966/AH7YCFBBXX_DS-233966_ACAGTG_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233966/AH7YCFBBXX_DS-233966_ACAGTG_L007_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233966/BHFLGLBBXX_DS-233966_ACAGTG_L006_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233966/BHFLGLBBXX_DS-233966_ACAGTG_L006_R2_001.fastq.gz"
    ],
    "TP6933_Recur_Base2" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234566/AH7YCFBBXX_DS-234566_ACTTGA_L001_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234566/AH7YCFBBXX_DS-234566_ACTTGA_L001_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234566/BHFG5FBBXX_DS-234566_ACTTGA_L002_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234566/BHFG5FBBXX_DS-234566_ACTTGA_L002_R2_001.fastq.gz"
    ],
    "TP6933_Recur_Recur" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234504/BHCHJNBBXX_DS-234504_ACTTGA_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234504/BHCHJNBBXX_DS-234504_ACTTGA_L007_R2_001.fastq.gz"
    ],
    "TP7022_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233431/BHCHJNBBXX_DS-233431_GATCAG_L001_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233431/BHCHJNBBXX_DS-233431_GATCAG_L001_R2_001.fastq.gz"
    ],
    "TP7084_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233884/BHCHJNBBXX_DS-233884_GAGTGG_L002_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233884/BHCHJNBBXX_DS-233884_GAGTGG_L002_R2_001.fastq.gz"
    ],
    "TP7084_Recur_Base2" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233930/BHCHJNBBXX_DS-233930_GTTTCG_L004_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233930/BHCHJNBBXX_DS-233930_GTTTCG_L004_R2_001.fastq.gz"
    ],
    "TP7084_Recur_Recur" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233923/BHCHJNBBXX_DS-233923_AGTCAA_L004_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233923/BHCHJNBBXX_DS-233923_AGTCAA_L004_R2_001.fastq.gz"
    ],
    "TP7123_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233887/BHCHJNBBXX_DS-233887_ATCACG_L003_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233887/BHCHJNBBXX_DS-233887_ATCACG_L003_R2_001.fastq.gz"
    ],
    "TP7138_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233926/BHCHJNBBXX_DS-233926_CCGTCC_L004_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233926/BHCHJNBBXX_DS-233926_CCGTCC_L004_R2_001.fastq.gz"
    ],
    "TP7335_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234526/BHCHJNBBXX_DS-234526_GCCAAT_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234526/BHCHJNBBXX_DS-234526_GCCAAT_L008_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234526/BHFG5FBBXX_DS-234526_GCCAAT_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234526/BHFG5FBBXX_DS-234526_GCCAAT_L008_R2_001.fastq.gz"
    ],
    "TP7335_Recur_Base2" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234555/AH7YCFBBXX_DS-234555_CGTACG_L001_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234555/AH7YCFBBXX_DS-234555_CGTACG_L001_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234555/BHFG5FBBXX_DS-234555_CGTACG_L002_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234555/BHFG5FBBXX_DS-234555_CGTACG_L002_R2_001.fastq.gz"
    ],
    "TP7335_Recur_Recur" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234532/AH7YCFBBXX_DS-234532_CTTGTA_L002_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234532/AH7YCFBBXX_DS-234532_CTTGTA_L002_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234532/BHFLGLBBXX_DS-234532_CTTGTA_L006_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234532/BHFLGLBBXX_DS-234532_CTTGTA_L006_R2_001.fastq.gz"
    ],
    "TP7376_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233885/BHCHJNBBXX_DS-233885_ACTGAT_L002_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233885/BHCHJNBBXX_DS-233885_ACTGAT_L002_R2_001.fastq.gz"
    ],
    "TP7415_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233898/BH7WHCBBXX_DS-233898_CTTGTA_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233898/BH7WHCBBXX_DS-233898_CTTGTA_L008_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233898/BHCHJNBBXX_DS-233898_CTTGTA_L003_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233898/BHCHJNBBXX_DS-233898_CTTGTA_L003_R2_001.fastq.gz"
    ],
    "TP7415_Recur_Recur" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233910/BHCHJNBBXX_DS-233910_ATTCCT_L003_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233910/BHCHJNBBXX_DS-233910_ATTCCT_L003_R2_001.fastq.gz"
    ],
    "TP7520_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232788/AH7W5TBBXX_DS-232788_AGTTCC_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232788/AH7W5TBBXX_DS-232788_AGTTCC_L007_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232788/AH7YCFBBXX_DS-232788_AGTTCC_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232788/AH7YCFBBXX_DS-232788_AGTTCC_L008_R2_001.fastq.gz"
    ],
    "TP7520_Recur_Base2" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232789/AH7W5TBBXX_DS-232789_ATGTCA_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232789/AH7W5TBBXX_DS-232789_ATGTCA_L007_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232789/AH7YCFBBXX_DS-232789_ATGTCA_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232789/AH7YCFBBXX_DS-232789_ATGTCA_L008_R2_001.fastq.gz"
    ],
    "TP7520_Recur_Recur" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232709/AH7YCFBBXX_DS-232709_CTTGTA_L003_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232709/AH7YCFBBXX_DS-232709_CTTGTA_L003_R2_001.fastq.gz"
    ],
    "TP7644_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232735/AH7W5TBBXX_DS-232735_TTAGGC_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232735/AH7W5TBBXX_DS-232735_TTAGGC_L007_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232735/AH7YCFBBXX_DS-232735_TTAGGC_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232735/AH7YCFBBXX_DS-232735_TTAGGC_L008_R2_001.fastq.gz"
    ],
    "TP7669_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233946/BHCHJNBBXX_DS-233946_ACTGAT_L004_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233946/BHCHJNBBXX_DS-233946_ACTGAT_L004_R2_001.fastq.gz"
    ],
    "TP7745_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232720/AH7W5TBBXX_DS-232720_CGTACG_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232720/AH7W5TBBXX_DS-232720_CGTACG_L007_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232720/BH7WHCBBXX_DS-232720_CGTACG_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232720/BH7WHCBBXX_DS-232720_CGTACG_L007_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232720/BH7WHCBBXX_DS-232720_CGTACG_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232720/BH7WHCBBXX_DS-232720_CGTACG_L008_R2_001.fastq.gz"
    ],
    "TP7745_Recur_Recur" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232732/AH7W5TBBXX_DS-232732_ATCACG_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232732/AH7W5TBBXX_DS-232732_ATCACG_L007_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232732/AH7YCFBBXX_DS-232732_ATCACG_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232732/AH7YCFBBXX_DS-232732_ATCACG_L008_R2_001.fastq.gz"
    ],
    "TP7818_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233472/BH7WHCBBXX_DS-233472_CTTGTA_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233472/BH7WHCBBXX_DS-233472_CTTGTA_L007_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233472/BHCHJNBBXX_DS-233472_CTTGTA_L002_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233472/BHCHJNBBXX_DS-233472_CTTGTA_L002_R2_001.fastq.gz"
    ],
    "TP8093_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234543/AH7YCFBBXX_DS-234543_GATCAG_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234543/AH7YCFBBXX_DS-234543_GATCAG_L007_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234543/BHFG5FBBXX_DS-234543_GATCAG_L001_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234543/BHFG5FBBXX_DS-234543_GATCAG_L001_R2_001.fastq.gz"
    ],
    "TP8106_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232780/AH7W5TBBXX_DS-232780_GATCAG_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232780/AH7W5TBBXX_DS-232780_GATCAG_L007_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232780/AH7YCFBBXX_DS-232780_GATCAG_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232780/AH7YCFBBXX_DS-232780_GATCAG_L008_R2_001.fastq.gz"
    ],
    "TP8106_Recur_Base2" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232713/AH7YCFBBXX_DS-232713_CCGTCC_L003_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232713/AH7YCFBBXX_DS-232713_CCGTCC_L003_R2_001.fastq.gz"
    ],
    "TP8106_Recur_Recur" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232792/AH7W5TBBXX_DS-232792_GTGAAA_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232792/AH7W5TBBXX_DS-232792_GTGAAA_L008_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232792/AH7YCFBBXX_DS-232792_GTGAAA_L004_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232792/AH7YCFBBXX_DS-232792_GTGAAA_L004_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232792/AH7YCFBBXX_DS-232792_GTGAAA_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232792/AH7YCFBBXX_DS-232792_GTGAAA_L008_R2_001.fastq.gz"
    ],
    "TP8117_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233880/BHCHJNBBXX_DS-233880_GTGAAA_L002_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233880/BHCHJNBBXX_DS-233880_GTGAAA_L002_R2_001.fastq.gz"
    ],
    "TP8117_Recur_Base2" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233900/BH7WHCBBXX_DS-233900_AGTTCC_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233900/BH7WHCBBXX_DS-233900_AGTTCC_L008_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233900/BHCHJNBBXX_DS-233900_AGTTCC_L003_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233900/BHCHJNBBXX_DS-233900_AGTTCC_L003_R2_001.fastq.gz"
    ],
    "TP8192_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232722/AH7W5TBBXX_DS-232722_GAGTGG_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232722/AH7W5TBBXX_DS-232722_GAGTGG_L007_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232722/BH7WHCBBXX_DS-232722_GAGTGG_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232722/BH7WHCBBXX_DS-232722_GAGTGG_L007_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232722/BH7WHCBBXX_DS-232722_GAGTGG_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232722/BH7WHCBBXX_DS-232722_GAGTGG_L008_R2_001.fastq.gz"
    ],
    "TP8213_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234550/AH7YCFBBXX_DS-234550_CCGTCC_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234550/AH7YCFBBXX_DS-234550_CCGTCC_L007_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234550/BHFG5FBBXX_DS-234550_CCGTCC_L001_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234550/BHFG5FBBXX_DS-234550_CCGTCC_L001_R2_001.fastq.gz"
    ],
    "TP8476_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233416/BHCHJNBBXX_DS-233416_ACTGAT_L001_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233416/BHCHJNBBXX_DS-233416_ACTGAT_L001_R2_001.fastq.gz"
    ],
    "TP8476_Recur_Recur" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233377/AH7W5TBBXX_DS-233377_TGACCA_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233377/AH7W5TBBXX_DS-233377_TGACCA_L008_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233377/AH7YCFBBXX_DS-233377_TGACCA_L004_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233377/AH7YCFBBXX_DS-233377_TGACCA_L004_R2_001.fastq.gz"
    ],
    "TP8708_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234511/BHCHJNBBXX_DS-234511_ATGTCA_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234511/BHCHJNBBXX_DS-234511_ATGTCA_L008_R2_001.fastq.gz"
    ],
    "TP8838_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233464/BHCHJNBBXX_DS-233464_ACTTGA_L002_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233464/BHCHJNBBXX_DS-233464_ACTTGA_L002_R2_001.fastq.gz"
    ],
    "TP8925_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233442/BHCHJNBBXX_DS-233442_CGATGT_L002_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233442/BHCHJNBBXX_DS-233442_CGATGT_L002_R2_001.fastq.gz"
    ],
    "TP8925_Recur_Base2" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233381/AH7W5TBBXX_DS-233381_GCCAAT_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233381/AH7W5TBBXX_DS-233381_GCCAAT_L008_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233381/AH7YCFBBXX_DS-233381_GCCAAT_L004_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233381/AH7YCFBBXX_DS-233381_GCCAAT_L004_R2_001.fastq.gz"
    ],
    "TP8941_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234575/AH7YCFBBXX_DS-234575_GTCCGC_L006_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234575/AH7YCFBBXX_DS-234575_GTCCGC_L006_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234575/BHFLGLBBXX_DS-234575_GTCCGC_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234575/BHFLGLBBXX_DS-234575_GTCCGC_L007_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234575/BHFLGLBBXX_DS-234575_GTCCGC_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234575/BHFLGLBBXX_DS-234575_GTCCGC_L008_R2_001.fastq.gz"
    ],
    "TP8943_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233413/BHCHJNBBXX_DS-233413_GAGTGG_L001_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233413/BHCHJNBBXX_DS-233413_GAGTGG_L001_R2_001.fastq.gz"
    ],
    "TP8943_Recur_Recur" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233382/AH7W5TBBXX_DS-233382_CAGATC_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233382/AH7W5TBBXX_DS-233382_CAGATC_L008_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233382/AH7YCFBBXX_DS-233382_CAGATC_L004_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233382/AH7YCFBBXX_DS-233382_CAGATC_L004_R2_001.fastq.gz"
    ],
    "TP9033_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234503/BHCHJNBBXX_DS-234503_CAGATC_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234503/BHCHJNBBXX_DS-234503_CAGATC_L007_R2_001.fastq.gz"
    ],
    "TP9033_Recur_Recur" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234499/BHCHJNBBXX_DS-234499_TTAGGC_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234499/BHCHJNBBXX_DS-234499_TTAGGC_L007_R2_001.fastq.gz"
    ],
    "TP9049_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234496/BHCHJNBBXX_DS-234496_ATTCCT_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234496/BHCHJNBBXX_DS-234496_ATTCCT_L007_R2_001.fastq.gz"
    ],
    "TP9121_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233904/BHCHJNBBXX_DS-233904_GTGAAA_L003_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233904/BHCHJNBBXX_DS-233904_GTGAAA_L003_R2_001.fastq.gz"
    ],
    "TP9121_Recur_Base2" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233918/BHCHJNBBXX_DS-233918_ACTTGA_L004_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233918/BHCHJNBBXX_DS-233918_ACTTGA_L004_R2_001.fastq.gz"
    ],
    "TP9125_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234517/BHCHJNBBXX_DS-234517_CGTACG_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234517/BHCHJNBBXX_DS-234517_CGTACG_L008_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234517/BHFG5FBBXX_DS-234517_CGTACG_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234517/BHFG5FBBXX_DS-234517_CGTACG_L008_R2_001.fastq.gz"
    ],
    "TP9125_Recur_Base2" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234522/BHCHJNBBXX_DS-234522_CGATGT_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234522/BHCHJNBBXX_DS-234522_CGATGT_L008_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234522/BHFG5FBBXX_DS-234522_CGATGT_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234522/BHFG5FBBXX_DS-234522_CGATGT_L008_R2_001.fastq.gz"
    ],
    "TP9206_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233419/BH7WHCBBXX_DS-233419_CGATGT_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233419/BH7WHCBBXX_DS-233419_CGATGT_L008_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233419/BHCHJNBBXX_DS-233419_CGATGT_L001_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233419/BHCHJNBBXX_DS-233419_CGATGT_L001_R2_001.fastq.gz"
    ],
    "TP9304_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233383/AH7W5TBBXX_DS-233383_ACTTGA_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233383/AH7W5TBBXX_DS-233383_ACTTGA_L008_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233383/AH7YCFBBXX_DS-233383_ACTTGA_L004_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233383/AH7YCFBBXX_DS-233383_ACTTGA_L004_R2_001.fastq.gz"
    ],
    "TP9304_Recur_Base2" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233437/BHCHJNBBXX_DS-233437_AGTCAA_L001_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233437/BHCHJNBBXX_DS-233437_AGTCAA_L001_R2_001.fastq.gz"
    ],
    "TP9311_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233433/BHCHJNBBXX_DS-233433_TAGCTT_L001_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233433/BHCHJNBBXX_DS-233433_TAGCTT_L001_R2_001.fastq.gz"
    ],
    "TP9400_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233911/BHCHJNBBXX_DS-233911_ATCACG_L004_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233911/BHCHJNBBXX_DS-233911_ATCACG_L004_R2_001.fastq.gz"
    ],
    "TP9418_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233423/BH7WHCBBXX_DS-233423_TGACCA_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233423/BH7WHCBBXX_DS-233423_TGACCA_L008_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233423/BHCHJNBBXX_DS-233423_TGACCA_L001_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233423/BHCHJNBBXX_DS-233423_TGACCA_L001_R2_001.fastq.gz"
    ],
    "TP9453_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233976/BHCHJNBBXX_DS-233976_ATGTCA_L006_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233976/BHCHJNBBXX_DS-233976_ATGTCA_L006_R2_001.fastq.gz"
    ],
    "TP9453_Recur_Recur" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233883/BHCHJNBBXX_DS-233883_CGTACG_L002_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233883/BHCHJNBBXX_DS-233883_CGTACG_L002_R2_001.fastq.gz"
    ],
    "TP9647_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234527/BHCHJNBBXX_DS-234527_CAGATC_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234527/BHCHJNBBXX_DS-234527_CAGATC_L008_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234527/BHFG5FBBXX_DS-234527_CAGATC_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234527/BHFG5FBBXX_DS-234527_CAGATC_L008_R2_001.fastq.gz"
    ],
    "TQ10159_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233974/BHCHJNBBXX_DS-233974_AGTCAA_L006_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233974/BHCHJNBBXX_DS-233974_AGTCAA_L006_R2_001.fastq.gz"
    ],
    "TQ10293_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234515/BHCHJNBBXX_DS-234515_GTGGCC_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234515/BHCHJNBBXX_DS-234515_GTGGCC_L008_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234515/BHFG5FBBXX_DS-234515_GTGGCC_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234515/BHFG5FBBXX_DS-234515_GTGGCC_L008_R2_001.fastq.gz"
    ],
    "TQ10377_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234498/BHCHJNBBXX_DS-234498_CGATGT_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234498/BHCHJNBBXX_DS-234498_CGATGT_L007_R2_001.fastq.gz"
    ],
    "TQ10377_Recur_Base2" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234482/BHCHJNBBXX_DS-234482_TAGCTT_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234482/BHCHJNBBXX_DS-234482_TAGCTT_L007_R2_001.fastq.gz"
    ],
    "TQ10377_Recur_Recur" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234516/BHCHJNBBXX_DS-234516_GTTTCG_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234516/BHCHJNBBXX_DS-234516_GTTTCG_L008_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234516/BHFG5FBBXX_DS-234516_GTTTCG_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234516/BHFG5FBBXX_DS-234516_GTTTCG_L008_R2_001.fastq.gz"
    ],
    "TQ10388_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234573/AH7YCFBBXX_DS-234573_ATGTCA_L006_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234573/AH7YCFBBXX_DS-234573_ATGTCA_L006_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234573/BHFLGLBBXX_DS-234573_ATGTCA_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234573/BHFLGLBBXX_DS-234573_ATGTCA_L007_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234573/BHFLGLBBXX_DS-234573_ATGTCA_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234573/BHFLGLBBXX_DS-234573_ATGTCA_L008_R2_001.fastq.gz"
    ],
    "TQ10409_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234558/AH7YCFBBXX_DS-234558_ATTCCT_L001_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234558/AH7YCFBBXX_DS-234558_ATTCCT_L001_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234558/BHFG5FBBXX_DS-234558_ATTCCT_L002_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234558/BHFG5FBBXX_DS-234558_ATTCCT_L002_R2_001.fastq.gz"
    ],
    "TQ10409_Recur_Recur" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234492/BHCHJNBBXX_DS-234492_GTTTCG_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234492/BHCHJNBBXX_DS-234492_GTTTCG_L007_R2_001.fastq.gz"
    ],
    "TQ10411_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234572/AH7YCFBBXX_DS-234572_AGTTCC_L006_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234572/AH7YCFBBXX_DS-234572_AGTTCC_L006_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234572/BHFLGLBBXX_DS-234572_AGTTCC_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234572/BHFLGLBBXX_DS-234572_AGTTCC_L007_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234572/BHFLGLBBXX_DS-234572_AGTTCC_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234572/BHFLGLBBXX_DS-234572_AGTTCC_L008_R2_001.fastq.gz"
    ],
    "TQ10690_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232783/AH7W5TBBXX_DS-232783_TAGCTT_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232783/AH7W5TBBXX_DS-232783_TAGCTT_L007_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232783/AH7YCFBBXX_DS-232783_TAGCTT_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232783/AH7YCFBBXX_DS-232783_TAGCTT_L008_R2_001.fastq.gz"
    ],
    "TQ10690_Recur_Recur" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232794/AH7W5TBBXX_DS-232794_GTTTCG_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232794/AH7W5TBBXX_DS-232794_GTTTCG_L008_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232794/AH7YCFBBXX_DS-232794_GTTTCG_L004_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232794/AH7YCFBBXX_DS-232794_GTTTCG_L004_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232794/AH7YCFBBXX_DS-232794_GTTTCG_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232794/AH7YCFBBXX_DS-232794_GTTTCG_L008_R2_001.fastq.gz"
    ],
    "TQ10754_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233986/BHCHJNBBXX_DS-233986_ATCACG_L006_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233986/BHCHJNBBXX_DS-233986_ATCACG_L006_R2_001.fastq.gz"
    ],
    "TQ10754_Recur_Base2" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233906/BHCHJNBBXX_DS-233906_GTTTCG_L003_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233906/BHCHJNBBXX_DS-233906_GTTTCG_L003_R2_001.fastq.gz"
    ],
    "TQ10754_Recur_Recur" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233989/BHCHJNBBXX_DS-233989_TGACCA_L006_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233989/BHCHJNBBXX_DS-233989_TGACCA_L006_R2_001.fastq.gz"
    ],
    "TQ10809_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233443/BHCHJNBBXX_DS-233443_TTAGGC_L002_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233443/BHCHJNBBXX_DS-233443_TTAGGC_L002_R2_001.fastq.gz"
    ],
    "TQ10809_Recur_Recur" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233435/BHCHJNBBXX_DS-233435_GGCTAC_L001_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233435/BHCHJNBBXX_DS-233435_GGCTAC_L001_R2_001.fastq.gz"
    ],
    "TQ11113_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234534/AH7YCFBBXX_DS-234534_AGTTCC_L002_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234534/AH7YCFBBXX_DS-234534_AGTTCC_L002_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234534/BHFLGLBBXX_DS-234534_AGTTCC_L006_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234534/BHFLGLBBXX_DS-234534_AGTTCC_L006_R2_001.fastq.gz"
    ],
    "TQ11113_Recur_Base2" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234506/BHCHJNBBXX_DS-234506_TAGCTT_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234506/BHCHJNBBXX_DS-234506_TAGCTT_L008_R2_001.fastq.gz"
    ],
    "TQ11568_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233950/AH7YCFBBXX_DS-233950_TTAGGC_L005_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233950/AH7YCFBBXX_DS-233950_TTAGGC_L005_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233950/BHFLGLBBXX_DS-233950_TTAGGC_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233950/BHFLGLBBXX_DS-233950_TTAGGC_L007_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233950/BHFLGLBBXX_DS-233950_TTAGGC_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233950/BHFLGLBBXX_DS-233950_TTAGGC_L008_R2_001.fastq.gz"
    ],
    "TQ11761_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232688/AH7W5TBBXX_DS-232688_CAGATC_L006_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232688/AH7W5TBBXX_DS-232688_CAGATC_L006_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232688/BH7WHCBBXX_DS-232688_CAGATC_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232688/BH7WHCBBXX_DS-232688_CAGATC_L007_R2_001.fastq.gz"
    ],
    "TQ11771_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232693/AH7W5TBBXX_DS-232693_ACTTGA_L006_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232693/AH7W5TBBXX_DS-232693_ACTTGA_L006_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232693/BH7WHCBBXX_DS-232693_ACTTGA_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232693/BH7WHCBBXX_DS-232693_ACTTGA_L007_R2_001.fastq.gz"
    ],
    "TQ11931_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234551/AH7YCFBBXX_DS-234551_GTCCGC_L001_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234551/AH7YCFBBXX_DS-234551_GTCCGC_L001_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234551/BHFG5FBBXX_DS-234551_GTCCGC_L002_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234551/BHFG5FBBXX_DS-234551_GTCCGC_L002_R2_001.fastq.gz"
    ],
    "TQ11953_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232683/AH7W5TBBXX_DS-232683_ACAGTG_L006_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232683/AH7W5TBBXX_DS-232683_ACAGTG_L006_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232683/BH7WHCBBXX_DS-232683_ACAGTG_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232683/BH7WHCBBXX_DS-232683_ACAGTG_L007_R2_001.fastq.gz"
    ],
    "TQ12041_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233476/BH7WHCBBXX_DS-233476_AGTTCC_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233476/BH7WHCBBXX_DS-233476_AGTTCC_L007_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233476/BHCHJNBBXX_DS-233476_AGTTCC_L002_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233476/BHCHJNBBXX_DS-233476_AGTTCC_L002_R2_001.fastq.gz"
    ],
    "TQ12151_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234494/BHCHJNBBXX_DS-234494_GAGTGG_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234494/BHCHJNBBXX_DS-234494_GAGTGG_L007_R2_001.fastq.gz"
    ],
    "TQ12184_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234489/BHCHJNBBXX_DS-234489_GTCCGC_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234489/BHCHJNBBXX_DS-234489_GTCCGC_L007_R2_001.fastq.gz"
    ],
    "TQ12184_Recur_Recur" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234497/BHCHJNBBXX_DS-234497_ATCACG_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234497/BHCHJNBBXX_DS-234497_ATCACG_L007_R2_001.fastq.gz"
    ],
    "TQ12220_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233426/BH7WHCBBXX_DS-233426_CAGATC_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233426/BH7WHCBBXX_DS-233426_CAGATC_L008_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233426/BHCHJNBBXX_DS-233426_CAGATC_L001_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233426/BHCHJNBBXX_DS-233426_CAGATC_L001_R2_001.fastq.gz"
    ],
    "TQ12279_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233420/BH7WHCBBXX_DS-233420_TTAGGC_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233420/BH7WHCBBXX_DS-233420_TTAGGC_L008_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233420/BHCHJNBBXX_DS-233420_TTAGGC_L001_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233420/BHCHJNBBXX_DS-233420_TTAGGC_L001_R2_001.fastq.gz"
    ],
    "TQ12279_Recur_Base2" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233379/AH7W5TBBXX_DS-233379_ACAGTG_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233379/AH7W5TBBXX_DS-233379_ACAGTG_L008_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233379/AH7YCFBBXX_DS-233379_ACAGTG_L004_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233379/AH7YCFBBXX_DS-233379_ACAGTG_L004_R2_001.fastq.gz"
    ],
    "TQ12391_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233460/BHCHJNBBXX_DS-233460_CAGATC_L002_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233460/BHCHJNBBXX_DS-233460_CAGATC_L002_R2_001.fastq.gz"
    ],
    "TQ12391_Recur_Recur" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233424/BH7WHCBBXX_DS-233424_ACAGTG_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233424/BH7WHCBBXX_DS-233424_ACAGTG_L008_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233424/BHCHJNBBXX_DS-233424_ACAGTG_L001_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233424/BHCHJNBBXX_DS-233424_ACAGTG_L001_R2_001.fastq.gz"
    ],
    "TQ12939_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232785/AH7W5TBBXX_DS-232785_CTTGTA_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232785/AH7W5TBBXX_DS-232785_CTTGTA_L007_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232785/AH7YCFBBXX_DS-232785_CTTGTA_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232785/AH7YCFBBXX_DS-232785_CTTGTA_L008_R2_001.fastq.gz"
    ],
    "TQ12939_Recur_Base2" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232803/AH7W5TBBXX_DS-232803_CGATGT_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232803/AH7W5TBBXX_DS-232803_CGATGT_L008_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232803/AH7YCFBBXX_DS-232803_CGATGT_L004_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232803/AH7YCFBBXX_DS-232803_CGATGT_L004_R2_001.fastq.gz"
    ],
    "TQ13267_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234488/BHCHJNBBXX_DS-234488_CCGTCC_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234488/BHCHJNBBXX_DS-234488_CCGTCC_L007_R2_001.fastq.gz"
    ],
    "TQ13507_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233891/BHCHJNBBXX_DS-233891_ACAGTG_L003_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233891/BHCHJNBBXX_DS-233891_ACAGTG_L003_R2_001.fastq.gz"
    ],
    "TQ13507_Recur_Base2" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233983/BHCHJNBBXX_DS-233983_GAGTGG_L006_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233983/BHCHJNBBXX_DS-233983_GAGTGG_L006_R2_001.fastq.gz"
    ],
    "TQ13514_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233971/BHCHJNBBXX_DS-233971_TAGCTT_L006_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233971/BHCHJNBBXX_DS-233971_TAGCTT_L006_R2_001.fastq.gz"
    ],
    "TQ13514_Recur_Recur" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233954/AH7YCFBBXX_DS-233954_CAGATC_L005_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233954/AH7YCFBBXX_DS-233954_CAGATC_L005_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233954/BHFLGLBBXX_DS-233954_CAGATC_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233954/BHFLGLBBXX_DS-233954_CAGATC_L007_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233954/BHFLGLBBXX_DS-233954_CAGATC_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233954/BHFLGLBBXX_DS-233954_CAGATC_L008_R2_001.fastq.gz"
    ],
    "TQ13732_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234537/AH7YCFBBXX_DS-234537_TTAGGC_L002_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234537/AH7YCFBBXX_DS-234537_TTAGGC_L002_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234537/BHFG5FBBXX_DS-234537_TTAGGC_L001_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234537/BHFG5FBBXX_DS-234537_TTAGGC_L001_R2_001.fastq.gz"
    ],
    "TQ13864_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234521/BHCHJNBBXX_DS-234521_ATCACG_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234521/BHCHJNBBXX_DS-234521_ATCACG_L008_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234521/BHFG5FBBXX_DS-234521_ATCACG_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234521/BHFG5FBBXX_DS-234521_ATCACG_L008_R2_001.fastq.gz"
    ],
    "TQ13864_Recur_Recur" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234477/BHCHJNBBXX_DS-234477_ACAGTG_L006_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234477/BHCHJNBBXX_DS-234477_ACAGTG_L006_R2_001.fastq.gz"
    ],
    "TQ13920_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234485/BHCHJNBBXX_DS-234485_AGTCAA_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234485/BHCHJNBBXX_DS-234485_AGTCAA_L007_R2_001.fastq.gz"
    ],
    "TQ14132_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232793/AH7W5TBBXX_DS-232793_GTGGCC_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232793/AH7W5TBBXX_DS-232793_GTGGCC_L008_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232793/AH7YCFBBXX_DS-232793_GTGGCC_L004_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232793/AH7YCFBBXX_DS-232793_GTGGCC_L004_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232793/AH7YCFBBXX_DS-232793_GTGGCC_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232793/AH7YCFBBXX_DS-232793_GTGGCC_L008_R2_001.fastq.gz"
    ],
    "TQ14236_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234564/AH7YCFBBXX_DS-234564_GCCAAT_L001_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234564/AH7YCFBBXX_DS-234564_GCCAAT_L001_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234564/BHFG5FBBXX_DS-234564_GCCAAT_L002_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234564/BHFG5FBBXX_DS-234564_GCCAAT_L002_R2_001.fastq.gz"
    ],
    "TQ14544_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233929/BHCHJNBBXX_DS-233929_GTGGCC_L004_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233929/BHCHJNBBXX_DS-233929_GTGGCC_L004_R2_001.fastq.gz"
    ],
    "TQ14625_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233425/BH7WHCBBXX_DS-233425_GCCAAT_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233425/BH7WHCBBXX_DS-233425_GCCAAT_L008_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233425/BHCHJNBBXX_DS-233425_GCCAAT_L001_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233425/BHCHJNBBXX_DS-233425_GCCAAT_L001_R2_001.fastq.gz"
    ],
    "TQ14677_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233407/BHCHJNBBXX_DS-233407_GTGAAA_L001_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233407/BHCHJNBBXX_DS-233407_GTGAAA_L001_R2_001.fastq.gz"
    ],
    "TQ14696_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233409/BHCHJNBBXX_DS-233409_GTTTCG_L001_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233409/BHCHJNBBXX_DS-233409_GTTTCG_L001_R2_001.fastq.gz"
    ],
    "TQ14696_Recur_Base2" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233389/AH7W5TBBXX_DS-233389_GGCTAC_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233389/AH7W5TBBXX_DS-233389_GGCTAC_L008_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233389/AH7YCFBBXX_DS-233389_GGCTAC_L004_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233389/AH7YCFBBXX_DS-233389_GGCTAC_L004_R2_001.fastq.gz"
    ],
    "TQ14696_Recur_Recur" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233458/BHCHJNBBXX_DS-233458_GCCAAT_L002_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233458/BHCHJNBBXX_DS-233458_GCCAAT_L002_R2_001.fastq.gz"
    ],
    "TQ14773_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232678/AH7W5TBBXX_DS-232678_ATCACG_L006_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232678/AH7W5TBBXX_DS-232678_ATCACG_L006_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232678/BH7WHCBBXX_DS-232678_ATCACG_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232678/BH7WHCBBXX_DS-232678_ATCACG_L007_R2_001.fastq.gz"
    ],
    "TQ14852_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234560/AH7YCFBBXX_DS-234560_CGATGT_L001_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234560/AH7YCFBBXX_DS-234560_CGATGT_L001_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234560/BHFG5FBBXX_DS-234560_CGATGT_L002_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234560/BHFG5FBBXX_DS-234560_CGATGT_L002_R2_001.fastq.gz"
    ],
    "TQ15022_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232746/AH7W5TBBXX_DS-232746_ACTTGA_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232746/AH7W5TBBXX_DS-232746_ACTTGA_L007_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232746/AH7YCFBBXX_DS-232746_ACTTGA_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232746/AH7YCFBBXX_DS-232746_ACTTGA_L008_R2_001.fastq.gz"
    ],
    "TQ15022_Recur_Recur" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232804/AH7W5TBBXX_DS-232804_TTAGGC_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232804/AH7W5TBBXX_DS-232804_TTAGGC_L008_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232804/AH7YCFBBXX_DS-232804_TTAGGC_L004_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232804/AH7YCFBBXX_DS-232804_TTAGGC_L004_R2_001.fastq.gz"
    ],
    "TQ15291_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234565/AH7YCFBBXX_DS-234565_CAGATC_L001_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234565/AH7YCFBBXX_DS-234565_CAGATC_L001_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234565/BHFG5FBBXX_DS-234565_CAGATC_L002_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234565/BHFG5FBBXX_DS-234565_CAGATC_L002_R2_001.fastq.gz"
    ],
    "TQ15291_Recur_Recur" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234574/AH7YCFBBXX_DS-234574_CCGTCC_L006_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234574/AH7YCFBBXX_DS-234574_CCGTCC_L006_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234574/BHFLGLBBXX_DS-234574_CCGTCC_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234574/BHFLGLBBXX_DS-234574_CCGTCC_L007_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234574/BHFLGLBBXX_DS-234574_CCGTCC_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234574/BHFLGLBBXX_DS-234574_CCGTCC_L008_R2_001.fastq.gz"
    ],
    "TQ15316_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233905/BHCHJNBBXX_DS-233905_GTGGCC_L003_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233905/BHCHJNBBXX_DS-233905_GTGGCC_L003_R2_001.fastq.gz"
    ],
    "TQ15338_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233985/BHCHJNBBXX_DS-233985_ATTCCT_L006_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233985/BHCHJNBBXX_DS-233985_ATTCCT_L006_R2_001.fastq.gz"
    ],
    "TQ15700_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233882/BHCHJNBBXX_DS-233882_GTTTCG_L002_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233882/BHCHJNBBXX_DS-233882_GTTTCG_L002_R2_001.fastq.gz"
    ],
    "TQ15811_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233894/BHCHJNBBXX_DS-233894_ACTTGA_L003_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233894/BHCHJNBBXX_DS-233894_ACTTGA_L003_R2_001.fastq.gz"
    ],
    "TQ15811_Norecur_Base2" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233948/AH7YCFBBXX_DS-233948_ATCACG_L005_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233948/AH7YCFBBXX_DS-233948_ATCACG_L005_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233948/BHFLGLBBXX_DS-233948_ATCACG_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233948/BHFLGLBBXX_DS-233948_ATCACG_L007_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233948/BHFLGLBBXX_DS-233948_ATCACG_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233948/BHFLGLBBXX_DS-233948_ATCACG_L008_R2_001.fastq.gz"
    ],
    "TQ15819_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232798/AH7W5TBBXX_DS-232798_ACTGAT_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232798/AH7W5TBBXX_DS-232798_ACTGAT_L008_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232798/AH7YCFBBXX_DS-232798_ACTGAT_L004_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232798/AH7YCFBBXX_DS-232798_ACTGAT_L004_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232798/AH7YCFBBXX_DS-232798_ACTGAT_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232798/AH7YCFBBXX_DS-232798_ACTGAT_L008_R2_001.fastq.gz"
    ],
    "TQ15819_Norecur_Base2" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232733/AH7W5TBBXX_DS-232733_CGATGT_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232733/AH7W5TBBXX_DS-232733_CGATGT_L007_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232733/AH7YCFBBXX_DS-232733_CGATGT_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232733/AH7YCFBBXX_DS-232733_CGATGT_L008_R2_001.fastq.gz"
    ],
    "TQ15923_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234508/BHCHJNBBXX_DS-234508_CTTGTA_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234508/BHCHJNBBXX_DS-234508_CTTGTA_L008_R2_001.fastq.gz"
    ],
    "TQ15923_Recur_Base2" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234514/BHCHJNBBXX_DS-234514_GTGAAA_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234514/BHCHJNBBXX_DS-234514_GTGAAA_L008_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234514/BHFG5FBBXX_DS-234514_GTGAAA_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234514/BHFG5FBBXX_DS-234514_GTGAAA_L008_R2_001.fastq.gz"
    ],
    "TQ15923_Recur_Recur" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234493/BHCHJNBBXX_DS-234493_CGTACG_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234493/BHCHJNBBXX_DS-234493_CGTACG_L007_R2_001.fastq.gz"
    ],
    "TQ15960_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233967/AH7YCFBBXX_DS-233967_GCCAAT_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233967/AH7YCFBBXX_DS-233967_GCCAAT_L007_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233967/BHFLGLBBXX_DS-233967_GCCAAT_L006_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233967/BHFLGLBBXX_DS-233967_GCCAAT_L006_R2_001.fastq.gz"
    ],
    "TQ15960_Recur_Base2" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233958/BHCHJNBBXX_DS-233958_GGCTAC_L005_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233958/BHCHJNBBXX_DS-233958_GGCTAC_L005_R2_001.fastq.gz"
    ],
    "TQ15960_Recur_Recur" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233969/AH7YCFBBXX_DS-233969_ACTTGA_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233969/AH7YCFBBXX_DS-233969_ACTTGA_L007_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233969/BHFLGLBBXX_DS-233969_ACTTGA_L006_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233969/BHFLGLBBXX_DS-233969_ACTTGA_L006_R2_001.fastq.gz"
    ],
    "TQ15971_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234518/BHCHJNBBXX_DS-234518_GAGTGG_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234518/BHCHJNBBXX_DS-234518_GAGTGG_L008_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234518/BHFG5FBBXX_DS-234518_GAGTGG_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234518/BHFG5FBBXX_DS-234518_GAGTGG_L008_R2_001.fastq.gz"
    ],
    "TQ15983_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233913/BHCHJNBBXX_DS-233913_TTAGGC_L004_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233913/BHCHJNBBXX_DS-233913_TTAGGC_L004_R2_001.fastq.gz"
    ],
    "TQ15983_Recur_Recur" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233917/BHCHJNBBXX_DS-233917_CAGATC_L004_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233917/BHCHJNBBXX_DS-233917_CAGATC_L004_R2_001.fastq.gz"
    ],
    "TQ16049_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232682/AH7W5TBBXX_DS-232682_TGACCA_L006_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232682/AH7W5TBBXX_DS-232682_TGACCA_L006_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232682/BH7WHCBBXX_DS-232682_TGACCA_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232682/BH7WHCBBXX_DS-232682_TGACCA_L007_R2_001.fastq.gz"
    ],
    "TQ16049_Recur_Recur" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232710/AH7YCFBBXX_DS-232710_AGTCAA_L003_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232710/AH7YCFBBXX_DS-232710_AGTCAA_L003_R2_001.fastq.gz"
    ],
    "TQ16093_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234478/BHCHJNBBXX_DS-234478_GCCAAT_L006_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234478/BHCHJNBBXX_DS-234478_GCCAAT_L006_R2_001.fastq.gz"
    ],
    "TQ16172_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234480/BHCHJNBBXX_DS-234480_ACTTGA_L006_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234480/BHCHJNBBXX_DS-234480_ACTTGA_L006_R2_001.fastq.gz"
    ],
    "TQ16172_Recur_Base2" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234528/BHCHJNBBXX_DS-234528_ACTTGA_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234528/BHCHJNBBXX_DS-234528_ACTTGA_L008_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234528/BHFG5FBBXX_DS-234528_ACTTGA_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234528/BHFG5FBBXX_DS-234528_ACTTGA_L008_R2_001.fastq.gz"
    ],
    "TQ16205_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234547/AH7YCFBBXX_DS-234547_AGTCAA_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234547/AH7YCFBBXX_DS-234547_AGTCAA_L007_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234547/BHFG5FBBXX_DS-234547_AGTCAA_L001_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234547/BHFG5FBBXX_DS-234547_AGTCAA_L001_R2_001.fastq.gz"
    ],
    "TQ16534_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233897/BH7WHCBBXX_DS-233897_GGCTAC_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233897/BH7WHCBBXX_DS-233897_GGCTAC_L008_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233897/BHCHJNBBXX_DS-233897_GGCTAC_L003_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233897/BHCHJNBBXX_DS-233897_GGCTAC_L003_R2_001.fastq.gz"
    ],
    "TQ16534_Recur_Base2" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233919/BHCHJNBBXX_DS-233919_GATCAG_L004_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233919/BHCHJNBBXX_DS-233919_GATCAG_L004_R2_001.fastq.gz"
    ],
    "TQ16561_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234552/AH7YCFBBXX_DS-234552_GTGAAA_L001_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234552/AH7YCFBBXX_DS-234552_GTGAAA_L001_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234552/BHFG5FBBXX_DS-234552_GTGAAA_L002_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234552/BHFG5FBBXX_DS-234552_GTGAAA_L002_R2_001.fastq.gz"
    ],
    "TQ16561_Recur_Base2" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234539/AH7YCFBBXX_DS-234539_ACAGTG_L002_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234539/AH7YCFBBXX_DS-234539_ACAGTG_L002_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234539/BHFG5FBBXX_DS-234539_ACAGTG_L001_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234539/BHFG5FBBXX_DS-234539_ACAGTG_L001_R2_001.fastq.gz"
    ],
    "TQ16683_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234490/BHCHJNBBXX_DS-234490_GTGAAA_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234490/BHCHJNBBXX_DS-234490_GTGAAA_L007_R2_001.fastq.gz"
    ],
    "TQ16683_Recur_Recur" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234570/AH7YCFBBXX_DS-234570_CTTGTA_L006_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234570/AH7YCFBBXX_DS-234570_CTTGTA_L006_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234570/BHFLGLBBXX_DS-234570_CTTGTA_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234570/BHFLGLBBXX_DS-234570_CTTGTA_L007_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234570/BHFLGLBBXX_DS-234570_CTTGTA_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234570/BHFLGLBBXX_DS-234570_CTTGTA_L008_R2_001.fastq.gz"
    ],
    "TQ17072_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233960/BHCHJNBBXX_DS-233960_AGTCAA_L005_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233960/BHCHJNBBXX_DS-233960_AGTCAA_L005_R2_001.fastq.gz"
    ],
    "TQ17302_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234563/AH7YCFBBXX_DS-234563_ACAGTG_L001_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234563/AH7YCFBBXX_DS-234563_ACAGTG_L001_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234563/BHFG5FBBXX_DS-234563_ACAGTG_L002_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234563/BHFG5FBBXX_DS-234563_ACAGTG_L002_R2_001.fastq.gz"
    ],
    "TQ17652_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232717/AH7W5TBBXX_DS-232717_GTGAAA_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232717/AH7W5TBBXX_DS-232717_GTGAAA_L007_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232717/BH7WHCBBXX_DS-232717_GTGAAA_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232717/BH7WHCBBXX_DS-232717_GTGAAA_L007_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232717/BH7WHCBBXX_DS-232717_GTGAAA_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232717/BH7WHCBBXX_DS-232717_GTGAAA_L008_R2_001.fastq.gz"
    ],
    "TQ18012_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234556/AH7YCFBBXX_DS-234556_GAGTGG_L001_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234556/AH7YCFBBXX_DS-234556_GAGTGG_L001_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234556/BHFG5FBBXX_DS-234556_GAGTGG_L002_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234556/BHFG5FBBXX_DS-234556_GAGTGG_L002_R2_001.fastq.gz"
    ],
    "TQ18136_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233924/BHCHJNBBXX_DS-233924_AGTTCC_L004_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233924/BHCHJNBBXX_DS-233924_AGTTCC_L004_R2_001.fastq.gz"
    ],
    "TQ18290_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234524/BHCHJNBBXX_DS-234524_TGACCA_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234524/BHCHJNBBXX_DS-234524_TGACCA_L008_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234524/BHFG5FBBXX_DS-234524_TGACCA_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234524/BHFG5FBBXX_DS-234524_TGACCA_L008_R2_001.fastq.gz"
    ],
    "TQ18556_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234545/AH7YCFBBXX_DS-234545_GGCTAC_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234545/AH7YCFBBXX_DS-234545_GGCTAC_L007_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234545/BHFG5FBBXX_DS-234545_GGCTAC_L001_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234545/BHFG5FBBXX_DS-234545_GGCTAC_L001_R2_001.fastq.gz"
    ],
    "TQ18586_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234567/AH7YCFBBXX_DS-234567_GATCAG_L006_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234567/AH7YCFBBXX_DS-234567_GATCAG_L006_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234567/BHFLGLBBXX_DS-234567_GATCAG_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234567/BHFLGLBBXX_DS-234567_GATCAG_L007_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234567/BHFLGLBBXX_DS-234567_GATCAG_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234567/BHFLGLBBXX_DS-234567_GATCAG_L008_R2_001.fastq.gz"
    ],
    "TQ18586_Recur_Recur" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234544/AH7YCFBBXX_DS-234544_TAGCTT_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234544/AH7YCFBBXX_DS-234544_TAGCTT_L007_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234544/BHFG5FBBXX_DS-234544_TAGCTT_L001_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234544/BHFG5FBBXX_DS-234544_TAGCTT_L001_R2_001.fastq.gz"
    ],
    "TQ18826_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233907/BHCHJNBBXX_DS-233907_CGTACG_L003_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233907/BHCHJNBBXX_DS-233907_CGTACG_L003_R2_001.fastq.gz"
    ],
    "TQ18826_Norecur_Base2" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233963/AH7YCFBBXX_DS-233963_CGATGT_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233963/AH7YCFBBXX_DS-233963_CGATGT_L007_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233963/BHFLGLBBXX_DS-233963_CGATGT_L006_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233963/BHFLGLBBXX_DS-233963_CGATGT_L006_R2_001.fastq.gz"
    ],
    "TQ18926_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233467/BH7WHCBBXX_DS-233467_TAGCTT_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233467/BH7WHCBBXX_DS-233467_TAGCTT_L007_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233467/BHCHJNBBXX_DS-233467_TAGCTT_L002_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233467/BHCHJNBBXX_DS-233467_TAGCTT_L002_R2_001.fastq.gz"
    ],
    "TQ19912_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234540/AH7YCFBBXX_DS-234540_GCCAAT_L002_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234540/AH7YCFBBXX_DS-234540_GCCAAT_L002_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234540/BHFG5FBBXX_DS-234540_GCCAAT_L001_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234540/BHFG5FBBXX_DS-234540_GCCAAT_L001_R2_001.fastq.gz"
    ],
    "TQ19912_Recur_Recur" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234530/AH7YCFBBXX_DS-234530_TAGCTT_L002_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234530/AH7YCFBBXX_DS-234530_TAGCTT_L002_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234530/BHFLGLBBXX_DS-234530_TAGCTT_L006_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234530/BHFLGLBBXX_DS-234530_TAGCTT_L006_R2_001.fastq.gz"
    ],
    "TQ19966_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234500/BHCHJNBBXX_DS-234500_TGACCA_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234500/BHCHJNBBXX_DS-234500_TGACCA_L007_R2_001.fastq.gz"
    ],
    "TQ20052_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232791/AH7W5TBBXX_DS-232791_GTCCGC_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232791/AH7W5TBBXX_DS-232791_GTCCGC_L008_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232791/AH7YCFBBXX_DS-232791_GTCCGC_L004_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232791/AH7YCFBBXX_DS-232791_GTCCGC_L004_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232791/AH7YCFBBXX_DS-232791_GTCCGC_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232791/AH7YCFBBXX_DS-232791_GTCCGC_L008_R2_001.fastq.gz"
    ],
    "TQ20052_Norecur_Base2" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232718/AH7W5TBBXX_DS-232718_GTGGCC_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232718/AH7W5TBBXX_DS-232718_GTGGCC_L007_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232718/BH7WHCBBXX_DS-232718_GTGGCC_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232718/BH7WHCBBXX_DS-232718_GTGGCC_L007_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232718/BH7WHCBBXX_DS-232718_GTGGCC_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232718/BH7WHCBBXX_DS-232718_GTGGCC_L008_R2_001.fastq.gz"
    ],
    "TQ20166_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234481/BHCHJNBBXX_DS-234481_GATCAG_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234481/BHCHJNBBXX_DS-234481_GATCAG_L007_R2_001.fastq.gz"
    ],
    "TQ20232_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233988/BHCHJNBBXX_DS-233988_TTAGGC_L006_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233988/BHCHJNBBXX_DS-233988_TTAGGC_L006_R2_001.fastq.gz"
    ],
    "TQ20232_Recur_Base2" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233970/BHCHJNBBXX_DS-233970_GATCAG_L006_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233970/BHCHJNBBXX_DS-233970_GATCAG_L006_R2_001.fastq.gz"
    ],
    "TQ20256_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233418/BH7WHCBBXX_DS-233418_ATCACG_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233418/BH7WHCBBXX_DS-233418_ATCACG_L008_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233418/BHCHJNBBXX_DS-233418_ATCACG_L001_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233418/BHCHJNBBXX_DS-233418_ATCACG_L001_R2_001.fastq.gz"
    ],
    "TQ20256_Norecur_Base2" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233384/AH7W5TBBXX_DS-233384_GATCAG_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233384/AH7W5TBBXX_DS-233384_GATCAG_L008_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233384/AH7YCFBBXX_DS-233384_GATCAG_L004_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233384/AH7YCFBBXX_DS-233384_GATCAG_L004_R2_001.fastq.gz"
    ],
    "TQ20572_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233945/BHCHJNBBXX_DS-233945_GAGTGG_L004_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233945/BHCHJNBBXX_DS-233945_GAGTGG_L004_R2_001.fastq.gz"
    ],
    "TQ20944_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233417/BHCHJNBBXX_DS-233417_ATTCCT_L001_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233417/BHCHJNBBXX_DS-233417_ATTCCT_L001_R2_001.fastq.gz"
    ],
    "TQ20944_Recur_Recur" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233401/AH7W5TBBXX_DS-233401_ATGTCA_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233401/AH7W5TBBXX_DS-233401_ATGTCA_L008_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233401/AH7YCFBBXX_DS-233401_ATGTCA_L004_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233401/AH7YCFBBXX_DS-233401_ATGTCA_L004_R2_001.fastq.gz"
    ],
    "TQ21638_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234484/BHCHJNBBXX_DS-234484_CTTGTA_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234484/BHCHJNBBXX_DS-234484_CTTGTA_L007_R2_001.fastq.gz"
    ],
    "TQ21852_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232701/AH7YCFBBXX_DS-232701_TAGCTT_L003_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232701/AH7YCFBBXX_DS-232701_TAGCTT_L003_R2_001.fastq.gz"
    ],
    "TQ21852_Norecur_Base2" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232745/AH7W5TBBXX_DS-232745_CAGATC_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232745/AH7W5TBBXX_DS-232745_CAGATC_L007_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232745/AH7YCFBBXX_DS-232745_CAGATC_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232745/AH7YCFBBXX_DS-232745_CAGATC_L008_R2_001.fastq.gz"
    ],
    "TQ21872_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234495/BHCHJNBBXX_DS-234495_ACTGAT_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234495/BHCHJNBBXX_DS-234495_ACTGAT_L007_R2_001.fastq.gz"
    ],
    "TQ22116_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234536/AH7YCFBBXX_DS-234536_CGATGT_L002_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234536/AH7YCFBBXX_DS-234536_CGATGT_L002_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234536/BHFG5FBBXX_DS-234536_CGATGT_L001_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234536/BHFG5FBBXX_DS-234536_CGATGT_L001_R2_001.fastq.gz"
    ],
    "TQ22646_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234491/BHCHJNBBXX_DS-234491_GTGGCC_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234491/BHCHJNBBXX_DS-234491_GTGGCC_L007_R2_001.fastq.gz"
    ],
    "TQ22646_Norecur_Base2" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234549/AH7YCFBBXX_DS-234549_ATGTCA_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234549/AH7YCFBBXX_DS-234549_ATGTCA_L007_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234549/BHFG5FBBXX_DS-234549_ATGTCA_L001_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234549/BHFG5FBBXX_DS-234549_ATGTCA_L001_R2_001.fastq.gz"
    ],
    "TQ22694_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234561/AH7YCFBBXX_DS-234561_TTAGGC_L001_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234561/AH7YCFBBXX_DS-234561_TTAGGC_L001_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234561/BHFG5FBBXX_DS-234561_TTAGGC_L002_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234561/BHFG5FBBXX_DS-234561_TTAGGC_L002_R2_001.fastq.gz"
    ],
    "TQ22876_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233412/BHCHJNBBXX_DS-233412_CGTACG_L001_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233412/BHCHJNBBXX_DS-233412_CGTACG_L001_R2_001.fastq.gz"
    ],
    "TQ22876_Recur_Recur" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233408/BHCHJNBBXX_DS-233408_GTGGCC_L001_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233408/BHCHJNBBXX_DS-233408_GTGGCC_L001_R2_001.fastq.gz"
    ],
    "TQ22888_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233981/BHCHJNBBXX_DS-233981_GTTTCG_L006_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233981/BHCHJNBBXX_DS-233981_GTTTCG_L006_R2_001.fastq.gz"
    ],
    "TQ23418_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233972/BHCHJNBBXX_DS-233972_GGCTAC_L006_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233972/BHCHJNBBXX_DS-233972_GGCTAC_L006_R2_001.fastq.gz"
    ],
    "TQ23570_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234546/AH7YCFBBXX_DS-234546_CTTGTA_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234546/AH7YCFBBXX_DS-234546_CTTGTA_L007_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234546/BHFG5FBBXX_DS-234546_CTTGTA_L001_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234546/BHFG5FBBXX_DS-234546_CTTGTA_L001_R2_001.fastq.gz"
    ],
    "TQ23758_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233979/BHCHJNBBXX_DS-233979_GTGAAA_L006_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233979/BHCHJNBBXX_DS-233979_GTGAAA_L006_R2_001.fastq.gz"
    ],
    "TQ23758_Norecur_Base2" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233956/BHCHJNBBXX_DS-233956_GATCAG_L005_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233956/BHCHJNBBXX_DS-233956_GATCAG_L005_R2_001.fastq.gz"
    ],
    "TQ23830_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234557/AH7YCFBBXX_DS-234557_ACTGAT_L001_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234557/AH7YCFBBXX_DS-234557_ACTGAT_L001_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234557/BHFG5FBBXX_DS-234557_ACTGAT_L002_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234557/BHFG5FBBXX_DS-234557_ACTGAT_L002_R2_001.fastq.gz"
    ],
    "TQ23830_Recur_Recur" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234568/AH7YCFBBXX_DS-234568_TAGCTT_L006_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234568/AH7YCFBBXX_DS-234568_TAGCTT_L006_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234568/BHFLGLBBXX_DS-234568_TAGCTT_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234568/BHFLGLBBXX_DS-234568_TAGCTT_L007_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234568/BHFLGLBBXX_DS-234568_TAGCTT_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234568/BHFLGLBBXX_DS-234568_TAGCTT_L008_R2_001.fastq.gz"
    ],
    "TQ24060_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233964/AH7YCFBBXX_DS-233964_TTAGGC_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233964/AH7YCFBBXX_DS-233964_TTAGGC_L007_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233964/BHFLGLBBXX_DS-233964_TTAGGC_L006_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233964/BHFLGLBBXX_DS-233964_TTAGGC_L006_R2_001.fastq.gz"
    ],
    "TQ24060_Recur_Base2" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233878/BH7WHCBBXX_DS-233878_CCGTCC_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233878/BH7WHCBBXX_DS-233878_CCGTCC_L007_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233878/BHCHJNBBXX_DS-233878_CCGTCC_L002_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233878/BHCHJNBBXX_DS-233878_CCGTCC_L002_R2_001.fastq.gz"
    ],
    "TQ24060_Recur_Recur" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233977/BHCHJNBBXX_DS-233977_CCGTCC_L006_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233977/BHCHJNBBXX_DS-233977_CCGTCC_L006_R2_001.fastq.gz"
    ],
    "TQ24080_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234501/BHCHJNBBXX_DS-234501_ACAGTG_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234501/BHCHJNBBXX_DS-234501_ACAGTG_L007_R2_001.fastq.gz"
    ],
    "TQ24080_Recur_Recur" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234486/BHCHJNBBXX_DS-234486_AGTTCC_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234486/BHCHJNBBXX_DS-234486_AGTTCC_L007_R2_001.fastq.gz"
    ],
    "TQ24378_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233902/BH7WHCBBXX_DS-233902_CCGTCC_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233902/BH7WHCBBXX_DS-233902_CCGTCC_L008_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233902/BHCHJNBBXX_DS-233902_CCGTCC_L003_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233902/BHCHJNBBXX_DS-233902_CCGTCC_L003_R2_001.fastq.gz"
    ],
    "TQ24378_Norecur_Base2" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233912/BHCHJNBBXX_DS-233912_CGATGT_L004_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233912/BHCHJNBBXX_DS-233912_CGATGT_L004_R2_001.fastq.gz"
    ],
    "TQ24496_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234519/BHCHJNBBXX_DS-234519_ACTGAT_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234519/BHCHJNBBXX_DS-234519_ACTGAT_L008_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234519/BHFG5FBBXX_DS-234519_ACTGAT_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234519/BHFG5FBBXX_DS-234519_ACTGAT_L008_R2_001.fastq.gz"
    ],
    "TQ24496_Norecur_Base2" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234523/BHCHJNBBXX_DS-234523_TTAGGC_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234523/BHCHJNBBXX_DS-234523_TTAGGC_L008_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234523/BHFG5FBBXX_DS-234523_TTAGGC_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234523/BHFG5FBBXX_DS-234523_TTAGGC_L008_R2_001.fastq.gz"
    ],
    "TQ24562_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232700/AH7YCFBBXX_DS-232700_GATCAG_L003_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232700/AH7YCFBBXX_DS-232700_GATCAG_L003_R2_001.fastq.gz"
    ],
    "TQ24562_Norecur_Base2" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232727/AH7W5TBBXX_DS-232727_ACTGAT_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232727/AH7W5TBBXX_DS-232727_ACTGAT_L007_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232727/BH7WHCBBXX_DS-232727_ACTGAT_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232727/BH7WHCBBXX_DS-232727_ACTGAT_L007_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232727/BH7WHCBBXX_DS-232727_ACTGAT_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232727/BH7WHCBBXX_DS-232727_ACTGAT_L008_R2_001.fastq.gz"
    ],
    "TQ24808_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233397/AH7W5TBBXX_DS-233397_AGTCAA_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233397/AH7W5TBBXX_DS-233397_AGTCAA_L008_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233397/AH7YCFBBXX_DS-233397_AGTCAA_L004_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233397/AH7YCFBBXX_DS-233397_AGTCAA_L004_R2_001.fastq.gz"
    ],
    "TQ25102_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234559/AH7YCFBBXX_DS-234559_ATCACG_L001_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234559/AH7YCFBBXX_DS-234559_ATCACG_L001_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234559/BHFG5FBBXX_DS-234559_ATCACG_L002_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234559/BHFG5FBBXX_DS-234559_ATCACG_L002_R2_001.fastq.gz"
    ],
    "TQ25634_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234541/AH7YCFBBXX_DS-234541_CAGATC_L002_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234541/AH7YCFBBXX_DS-234541_CAGATC_L002_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234541/BHFG5FBBXX_DS-234541_CAGATC_L001_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234541/BHFG5FBBXX_DS-234541_CAGATC_L001_R2_001.fastq.gz"
    ],
    "TQ25780_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233909/BHCHJNBBXX_DS-233909_ACTGAT_L003_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233909/BHCHJNBBXX_DS-233909_ACTGAT_L003_R2_001.fastq.gz"
    ],
    "TQ26030_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233449/BHCHJNBBXX_DS-233449_TGACCA_L002_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233449/BHCHJNBBXX_DS-233449_TGACCA_L002_R2_001.fastq.gz"
    ],
    "TQ26108_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234562/AH7YCFBBXX_DS-234562_TGACCA_L001_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234562/AH7YCFBBXX_DS-234562_TGACCA_L001_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234562/BHFG5FBBXX_DS-234562_TGACCA_L002_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234562/BHFG5FBBXX_DS-234562_TGACCA_L002_R2_001.fastq.gz"
    ],
    "TQ26142_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234502/BHCHJNBBXX_DS-234502_GCCAAT_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234502/BHCHJNBBXX_DS-234502_GCCAAT_L007_R2_001.fastq.gz"
    ],
    "TQ27614_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233973/BHCHJNBBXX_DS-233973_CTTGTA_L006_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233973/BHCHJNBBXX_DS-233973_CTTGTA_L006_R2_001.fastq.gz"
    ],
    "TQ27818_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234535/AH7YCFBBXX_DS-234535_ATCACG_L002_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234535/AH7YCFBBXX_DS-234535_ATCACG_L002_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234535/BHFG5FBBXX_DS-234535_ATCACG_L001_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234535/BHFG5FBBXX_DS-234535_ATCACG_L001_R2_001.fastq.gz"
    ],
    "TQ28186_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233436/BHCHJNBBXX_DS-233436_CTTGTA_L001_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233436/BHCHJNBBXX_DS-233436_CTTGTA_L001_R2_001.fastq.gz"
    ],
    "TQ28186_Recur_Recur" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233385/AH7W5TBBXX_DS-233385_TAGCTT_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233385/AH7W5TBBXX_DS-233385_TAGCTT_L008_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233385/AH7YCFBBXX_DS-233385_TAGCTT_L004_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233385/AH7YCFBBXX_DS-233385_TAGCTT_L004_R2_001.fastq.gz"
    ],
    "TQ29934_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232790/AH7W5TBBXX_DS-232790_CCGTCC_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232790/AH7W5TBBXX_DS-232790_CCGTCC_L007_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232790/AH7YCFBBXX_DS-232790_CCGTCC_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-232790/AH7YCFBBXX_DS-232790_CCGTCC_L008_R2_001.fastq.gz"
    ],
    "TQ29992_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233889/BHCHJNBBXX_DS-233889_TTAGGC_L003_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233889/BHCHJNBBXX_DS-233889_TTAGGC_L003_R2_001.fastq.gz"
    ],
    "TQ30516_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233978/BHCHJNBBXX_DS-233978_GTCCGC_L006_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233978/BHCHJNBBXX_DS-233978_GTCCGC_L006_R2_001.fastq.gz"
    ],
    "TQ31150_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233892/BHCHJNBBXX_DS-233892_GCCAAT_L003_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233892/BHCHJNBBXX_DS-233892_GCCAAT_L003_R2_001.fastq.gz"
    ],
    "TQ31150_Norecur_Base2" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233975/BHCHJNBBXX_DS-233975_AGTTCC_L006_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233975/BHCHJNBBXX_DS-233975_AGTTCC_L006_R2_001.fastq.gz"
    ],
    "TQ31266_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233901/BH7WHCBBXX_DS-233901_ATGTCA_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233901/BH7WHCBBXX_DS-233901_ATGTCA_L008_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233901/BHCHJNBBXX_DS-233901_ATGTCA_L003_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233901/BHCHJNBBXX_DS-233901_ATGTCA_L003_R2_001.fastq.gz"
    ],
    "TQ33220_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234507/BHCHJNBBXX_DS-234507_GGCTAC_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234507/BHCHJNBBXX_DS-234507_GGCTAC_L008_R2_001.fastq.gz"
    ],
    "TQ33722_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233877/BH7WHCBBXX_DS-233877_ATGTCA_L007_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233877/BH7WHCBBXX_DS-233877_ATGTCA_L007_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233877/BHCHJNBBXX_DS-233877_ATGTCA_L002_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233877/BHCHJNBBXX_DS-233877_ATGTCA_L002_R2_001.fastq.gz"
    ],
    "TQ33796_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233915/BHCHJNBBXX_DS-233915_ACAGTG_L004_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233915/BHCHJNBBXX_DS-233915_ACAGTG_L004_R2_001.fastq.gz"
    ],
    "TQ34432_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233881/BHCHJNBBXX_DS-233881_GTGGCC_L002_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233881/BHCHJNBBXX_DS-233881_GTGGCC_L002_R2_001.fastq.gz"
    ],
    "TQ34432_Norecur_Base2" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233928/BHCHJNBBXX_DS-233928_GTGAAA_L004_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233928/BHCHJNBBXX_DS-233928_GTGAAA_L004_R2_001.fastq.gz"
    ],
    "TQ35642_Recur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234509/BHCHJNBBXX_DS-234509_AGTCAA_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234509/BHCHJNBBXX_DS-234509_AGTCAA_L008_R2_001.fastq.gz"
    ],
    "TQ35642_Recur_Base2" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234512/BHCHJNBBXX_DS-234512_CCGTCC_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-234512/BHCHJNBBXX_DS-234512_CCGTCC_L008_R2_001.fastq.gz"
    ],
    "TQ37680_Norecur_Base1" => [
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233429/BH7WHCBBXX_DS-233429_ACTTGA_L008_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233429/BH7WHCBBXX_DS-233429_ACTTGA_L008_R2_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233429/BHCHJNBBXX_DS-233429_ACTTGA_L001_R1_001.fastq.gz",
      "/data/cqs/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233429/BHCHJNBBXX_DS-233429_ACTTGA_L001_R2_001.fastq.gz"
    ],
  },
  merge_fastq => {
    class      => "Format::MergeFastq",
    perform    => 1,
    target_dir => "${target_dir}/merge_fastq",
    option     => "",
    source_ref => "files",
    sh_direct  => 1,
    is_paired  => 1,
    cluster    => $cluster,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=2",
      "walltime" => "2",
      "mem"      => "10gb"
    },
  },
  fastqc => {
    class      => "QC::FastQC",
    perform    => 1,
    target_dir => "${target_dir}/fastqc",
    option     => "",
    source_ref => "merge_fastq",
    sh_direct  => 1,
    cluster    => $cluster,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=2",
      "walltime" => "2",
      "mem"      => "10gb"
    },
  },
  fastqc_summary => {
    class      => "QC::FastQCSummary",
    perform    => 1,
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
    source_ref => "merge_fastq",
    picard_jar => $picard_jar,
    sh_direct  => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "24",
      "mem"      => "40gb"
    },
  },
  bwa_refine => {
    class      => "GATK::Refine",
    perform    => 1,
    target_dir => "${target_dir}/bwa_refine",
    option     => "-Xmx40g",

    #gatk_option => "--fix_misencoded_quality_scores",
    gatk_option              => "",
    fasta_file               => $bwa_fasta,
    source_ref               => "bwa",
    vcf_files                => [ $dbsnp, $mills ],
    gatk_jar                 => $gatk_jar,
    picard_jar               => $picard_jar,
    sh_direct                => 0,
    slim_print_reads         => 1,
    use_self_slim_method     => 1,
    samtools_baq_calibration => 0,
    sorted                   => 1,
    pbs                      => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "2",
      "mem"      => "40gb"
    },
  },
  bwa_refine_genes_bam => {
    class      => "Samtools::View",
    perform    => 1,
    target_dir => "${target_dir}/bwa_refine_genes_bam",
    option     => "-h -b -L $gene_bed",
    extension  => ".bam",
    source_ref => "bwa_refine",
    sh_direct  => 1,
    sorted     => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "24",
      "mem"      => "10gb"
    },
  },
  qc3 => {
    class             => "QC::QC3bam",
    perform           => 0,
    target_dir        => "${target_dir}/qc3bam",
    option            => "",
    target_region_bed => $covered_bed,
    transcript_gtf    => $transcript_gtf,
    qc3_perl          => $qc3_perl,
    source_ref        => "bwa",
    pbs               => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  qc3_refine => {
    class             => "QC::QC3bam",
    perform           => 0,
    target_dir        => "${target_dir}/qc3bam_refine",
    option            => "",
    target_region_bed => $covered_bed,
    transcript_gtf    => $transcript_gtf,
    qc3_perl          => $qc3_perl,
    source_ref        => "bwa_refine",
    pbs               => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  bwa_refine_hc_gvcf => {
    class         => "GATK::HaplotypeCaller",
    perform       => 1,
    target_dir    => "${target_dir}/bwa_refine_hc_gvcf",
    option        => "",
    source_ref    => "bwa_refine",
    java_option   => "",
    fasta_file    => $bwa_fasta,
    gatk_jar      => $gatk_jar,
    bed_file      => $covered_bed,
    extension     => ".g.vcf",
    by_chromosome => 0,                                    #since we have the bed file, we cannot use by_chromosome.
    gvcf          => 1,
    sh_direct     => 0,
    pbs           => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  bwa_refine_hc_gvcf_vqsr => {
    class       => "GATK::VariantFilter",
    perform     => 1,
    target_dir  => "${target_dir}/bwa_refine_hc_gvcf_vqsr",
    option      => "",
    vqsr_mode   => 1,
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
    bed_file    => $covered_bed,
    sh_direct   => 1,
    pbs         => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  bwa_refine_hc_gvcf_vqsr_annovar => {
    class      => "Annotation::Annovar",
    perform    => 1,
    target_dir => "${target_dir}/bwa_refine_hc_gvcf_vqsr_annovar",
    source_ref => "bwa_refine_hc_gvcf_vqsr",
    option     => $annovar_param,
    annovar_db => $annovar_db,
    buildver   => "hg19",
    sh_direct  => 1,
    isvcf      => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "2",
      "mem"      => "10gb"
    },
  },
  bwa_refine_hc_gvcf_vqsr_annovar_filter => {
    class               => "Annotation::FilterAnnovar",
    perform             => 1,
    target_dir          => "${target_dir}/bwa_refine_hc_gvcf_vqsr_annovar_filter",
    source_ref          => "bwa_refine_hc_gvcf_vqsr_annovar",
    option              => "",
    sh_direct           => 1,
    maximum_exac_values => "0.001,0.01,0.1,1.0",
    pbs                 => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "2",
      "mem"      => "10gb"
    },
  },
  bwa_refine_hc_gvcf_vqsr_annovar_filter_cbioportal => {
    class               => "Annotation::CBioPortal",
    perform             => 1,
    target_dir          => "${target_dir}/bwa_refine_hc_gvcf_vqsr_annovar_filter",
    source_ref          => ["bwa_refine_hc_gvcf_vqsr_annovar_filter", ".SNP.tsv"],
    option              => "",
    sample_name_pattern => "cur_",
    sample_name_suffix  => "",
    gene_names          => $genes,
    sh_direct           => 1,
    pbs                 => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "2",
      "mem"      => "10gb"
    },
  },
  bwa_refine_hc_gvcf_vqsr_annovar_filter_oncoprint => {
    class               => "Annotation::OncoPrint",
    perform             => 1,
    target_dir          => "${target_dir}/bwa_refine_hc_gvcf_vqsr_annovar_filter",
    source_ref          => ["bwa_refine_hc_gvcf_vqsr_annovar_filter", ".SNP.tsv"],
    picture_width       => 6000,
    picture_height      => 1000,
    sample_name_pattern => "cur_",
    sample_name_suffix  => "",
    gene_names          => $genes,
    sh_direct           => 1,
    pbs                 => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "2",
      "mem"      => "10gb"
    },
  },
  bwa_refine_hc_gvcf_vqsr_annovar_filter_base1 => {
    class               => "Annotation::FilterAnnovar",
    perform             => 1,
    target_dir          => "${target_dir}/bwa_refine_hc_gvcf_vqsr_annovar_filter_base1",
    source_ref          => "bwa_refine_hc_gvcf_vqsr_annovar",
    option              => "",
    sample_name_pattern => "_Base1",
    sample_name_suffix  => ".base1",
    sh_direct           => 1,
    maximum_exac_values => "0.001,0.01,0.1,1.0",
    pbs                 => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "2",
      "mem"      => "10gb"
    },
  },
  bwa_refine_hc_gvcf_vqsr_annovar_filter_base1_cbioportal => {
    class               => "Annotation::CBioPortal",
    perform             => 1,
    target_dir          => "${target_dir}/bwa_refine_hc_gvcf_vqsr_annovar_filter_base1",
    source_ref          => ["bwa_refine_hc_gvcf_vqsr_annovar_filter_base1", ".SNP.tsv"],
    option              => "",
    sample_name_pattern => "_Base1",
    sample_name_suffix  => "",
    gene_names          => $genes,
    sh_direct           => 1,
    pbs                 => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "2",
      "mem"      => "10gb"
    },
  },
  bwa_refine_hc_gvcf_vqsr_annovar_filter_base1_oncoprint => {
    class               => "Annotation::OncoPrint",
    perform             => 1,
    target_dir          => "${target_dir}/bwa_refine_hc_gvcf_vqsr_annovar_filter_base1",
    source_ref          => ["bwa_refine_hc_gvcf_vqsr_annovar_filter_base1", ".SNP.tsv"],
    picture_width       => 4000,
    picture_height      => 1000,
    sample_name_pattern => "_Base1",
    sample_name_suffix  => "",
    gene_names          => "APC ROBO2 KRAS NRAS BRAF TP53",
    sh_direct           => 1,
    pbs                 => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "2",
      "mem"      => "10gb"
    },
  },
  bwa_refine_hc_gvcf_hardfilter => {
    class       => "GATK::VariantFilter",
    perform     => 0,
    target_dir  => "${target_dir}/bwa_refine_hc_gvcf_hardfilter",
    option      => "",
    vqsr_mode   => 0,
    source_ref  => "bwa_refine_hc_gvcf",
    java_option => "",
    fasta_file  => $bwa_fasta,
    dbsnp_vcf   => $dbsnp,
    gatk_jar    => $gatk_jar,
    cqstools    => $cqstools,
    bed_file    => $covered_bed,
    is_rna      => 0,
    sh_direct   => 1,
    pbs         => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  bwa_refine_hc_gvcf_hardfilter_annovar => {
    class      => "Annotation::Annovar",
    perform    => 0,
    target_dir => "${target_dir}/bwa_refine_hc_gvcf_hardfilter_annovar",
    source_ref => "bwa_refine_hc_gvcf_hardfilter",
    option     => $annovar_param,
    annovar_db => $annovar_db,
    buildver   => "hg19",
    sh_direct  => 1,
    isvcf      => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "2",
      "mem"      => "10gb"
    },
  },
  sequencetask => {
    class      => "CQS::SequenceTask",
    perform    => 1,
    target_dir => "${target_dir}/sequencetask",
    option     => "",
    source     => {
      step1 => [ "merge_fastq", "fastqc", "bwa", "bwa_refine", "bwa_refine_hc_gvcf" ],
      step2 => [
        "fastqc_summary",

        # "qc3", "qc3_refine",
        "bwa_refine_hc_gvcf_vqsr", "bwa_refine_hc_gvcf_vqsr_annovar", 
        "bwa_refine_hc_gvcf_vqsr_annovar_filter", 
        "bwa_refine_hc_gvcf_vqsr_annovar_filter_oncoprint",
        "bwa_refine_hc_gvcf_vqsr_annovar_filter_cbioportal",
        "bwa_refine_hc_gvcf_vqsr_annovar_filter_base1", 
        "bwa_refine_hc_gvcf_vqsr_annovar_filter_base1_oncoprint",
        "bwa_refine_hc_gvcf_vqsr_annovar_filter_base1_cbioportal"

          #, "bwa_refine_hc_gvcf_hardfilter", "bwa_refine_hc_gvcf_hardfilter_annovar",
      ],
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

performConfig($config);

1;

