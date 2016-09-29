#!/usr/bin/perl
use strict;
use warnings;

use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;
use Hash::Merge qw( merge );

my $target_dir = "/scratch/cqs/shengq1/temp/20160921_3199_targetDnaSeqCNV";
my $cqstools   = "/home/shengq1/cqstools/CQS.Tools.exe";
my $email      = "shilin.zhao\@vanderbilt.edu";

#genome file
my $bwa_fasta_hg38 =
"/scratch/cqs/zhaos/reference/hg38/hg38bundle/bwa_index_0.7.12/Homo_sapiens_assembly38.fasta";

#GATK realign
my $dbsnp_hg38 =
"/scratch/cqs/zhaos/reference/hg38/hg38bundle/Homo_sapiens_assembly38.dbsnp138.vcf";
my $mills_g1000_indel_hg38 =
"/scratch/cqs/zhaos/reference/hg38/hg38bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf";
my $g1000_phase1_indel_hg38 =
"/scratch/cqs/zhaos/reference/hg38/hg38bundle/Homo_sapiens_assembly38.known_indels.vcf";

#target region
my $covered_bed_hg20 =
"/scratch/cqs/zhaos/Pierre/3199_bedfile/131212_HG19_cancer_TS_EZ_TiledOnly_95pToHg38.bed";
my $covered_bed_hg38_padding250_CNV =
"/scratch/cqs/zhaos/Pierre/3199_bedfile/131212_HG19_cancer_TS_EZ_TiledOnly_95pToHg38.bed.forCNV.pad250.bed";

my $gatk_jar = "/home/zhaos/bin/GenomeAnalysisTK36.jar";    #GATK version 36
my $gatk4_jar =
  "/scratch/cqs/zhaos/bin/gatk-protected/build/libs/gatk-protected.jar"
  ;    #GATK version 4 for CNV

my $picard_jar = "/scratch/cqs/shengq1/local/bin/picard/picard.jar";

my $gslLibraryFile =
  "/usr/local/gsl/latest/x86_64/gcc46/nonet/lib/libgslcblas.so";
my $hdfViewFolder = "/scratch/cqs/zhaos/bin/HDFView/lib/linux/";
my $PanelOfNormal =
"/scratch/cqs/zhaos/Pierre/20160920_3199_gatk4_CNV/test/Normal.bed.All.coverage.tsv.PON";

my $config = {
	general => { task_name => "target3199CNV" },

	normalBamFiles => {
		"Normal1" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160825_3502_NormalTargetDNASeq_hg38/bwa/result/Normal1/Normal1.bam"
		],
		"Normal10" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160825_3502_NormalTargetDNASeq_hg38/bwa/result/Normal10/Normal10.bam"
		],
		"Normal2" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160825_3502_NormalTargetDNASeq_hg38/bwa/result/Normal2/Normal2.bam"
		],
		"Normal3" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160825_3502_NormalTargetDNASeq_hg38/bwa/result/Normal3/Normal3.bam"
		],
		"Normal4" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160825_3502_NormalTargetDNASeq_hg38/bwa/result/Normal4/Normal4.bam"
		],
		"Normal5" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160825_3502_NormalTargetDNASeq_hg38/bwa/result/Normal5/Normal5.bam"
		],
		"Normal6" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160825_3502_NormalTargetDNASeq_hg38/bwa/result/Normal6/Normal6.bam"
		],
		"Normal7" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160825_3502_NormalTargetDNASeq_hg38/bwa/result/Normal7/Normal7.bam"
		],
		"Normal8" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160825_3502_NormalTargetDNASeq_hg38/bwa/result/Normal8/Normal8.bam"
		],
		"Normal9" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160825_3502_NormalTargetDNASeq_hg38/bwa/result/Normal9/Normal9.bam"
		],

	},
	sampleBamFiles => {
		"T102" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T102/T102.bam"
		],
		"T103" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T103/T103.bam"
		],
		"T104" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T104/T104.bam"
		],
		"T105" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T105/T105.bam"
		],
		"T106" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T106/T106.bam"
		],
		"T107" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T107/T107.bam"
		],
		"T108" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T108/T108.bam"
		],
		"T109" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T109/T109.bam"
		],
		"T11" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T11/T11.bam"
		],
		"T110" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T110/T110.bam"
		],
		"T111" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T111/T111.bam"
		],
		"T112" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T112/T112.bam"
		],
		"T113" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T113/T113.bam"
		],
		"T114" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T114/T114.bam"
		],
		"T115" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T115/T115.bam"
		],
		"T116" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T116/T116.bam"
		],
		"T117" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T117/T117.bam"
		],
		"T118" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T118/T118.bam"
		],
		"T119" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T119/T119.bam"
		],
		"T12_Rep1" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T12_Rep1/T12_Rep1.bam"
		],
		"T12_Rep2" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T12_Rep2/T12_Rep2.bam"
		],
		"T120" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T120/T120.bam"
		],
		"T121" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T121/T121.bam"
		],
		"T122" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T122/T122.bam"
		],
		"T123" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T123/T123.bam"
		],
		"T124" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T124/T124.bam"
		],
		"T125" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T125/T125.bam"
		],
		"T126" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T126/T126.bam"
		],
		"T127" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T127/T127.bam"
		],
		"T129" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T129/T129.bam"
		],
		"T130" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T130/T130.bam"
		],
		"T131" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T131/T131.bam"
		],
		"T14" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T14/T14.bam"
		],
		"T15_Rep1" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T15_Rep1/T15_Rep1.bam"
		],
		"T15_Rep2" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T15_Rep2/T15_Rep2.bam"
		],
		"T19" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T19/T19.bam"
		],
		"T2" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T2/T2.bam"
		],
		"T20" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T20/T20.bam"
		],
		"T21" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T21/T21.bam"
		],
		"T24_Rep1" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T24_Rep1/T24_Rep1.bam"
		],
		"T24_Rep2" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T24_Rep2/T24_Rep2.bam"
		],
		"T26" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T26/T26.bam"
		],
		"T27" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T27/T27.bam"
		],
		"T28" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T28/T28.bam"
		],
		"T29" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T29/T29.bam"
		],
		"T30" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T30/T30.bam"
		],
		"T31" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T31/T31.bam"
		],
		"T32" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T32/T32.bam"
		],
		"T33" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T33/T33.bam"
		],
		"T34" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T34/T34.bam"
		],
		"T35" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T35/T35.bam"
		],
		"T36" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T36/T36.bam"
		],
		"T37" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T37/T37.bam"
		],
		"T38" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T38/T38.bam"
		],
		"T39" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T39/T39.bam"
		],
		"T4" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T4/T4.bam"
		],
		"T40" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T40/T40.bam"
		],
		"T41" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T41/T41.bam"
		],
		"T43" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T43/T43.bam"
		],
		"T44" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T44/T44.bam"
		],
		"T45" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T45/T45.bam"
		],
		"T46" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T46/T46.bam"
		],
		"T47" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T47/T47.bam"
		],
		"T48" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T48/T48.bam"
		],
		"T49" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T49/T49.bam"
		],
		"T5_Rep1" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T5_Rep1/T5_Rep1.bam"
		],
		"T5_Rep2" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T5_Rep2/T5_Rep2.bam"
		],
		"T50" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T50/T50.bam"
		],
		"T51" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T51/T51.bam"
		],
		"T52" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T52/T52.bam"
		],
		"T53" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T53/T53.bam"
		],
		"T54" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T54/T54.bam"
		],
		"T55" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T55/T55.bam"
		],
		"T57" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T57/T57.bam"
		],
		"T58" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T58/T58.bam"
		],
		"T59" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T59/T59.bam"
		],
		"T6" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T6/T6.bam"
		],
		"T60" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T60/T60.bam"
		],
		"T61" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T61/T61.bam"
		],
		"T62" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T62/T62.bam"
		],
		"T63" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T63/T63.bam"
		],
		"T64" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T64/T64.bam"
		],
		"T65" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T65/T65.bam"
		],
		"T66" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T66/T66.bam"
		],
		"T67" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T67/T67.bam"
		],
		"T68" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T68/T68.bam"
		],
		"T69" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T69/T69.bam"
		],
		"T70" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T70/T70.bam"
		],
		"T71" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T71/T71.bam"
		],
		"T73" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T73/T73.bam"
		],
		"T74" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T74/T74.bam"
		],
		"T76" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T76/T76.bam"
		],
		"T79" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T79/T79.bam"
		],
		"T8" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T8/T8.bam"
		],
		"T80" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T80/T80.bam"
		],
		"T81" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T81/T81.bam"
		],
		"T82" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T82/T82.bam"
		],
		"T83" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T83/T83.bam"
		],
		"T84" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T84/T84.bam"
		],
		"T85" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T85/T85.bam"
		],
		"T86" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T86/T86.bam"
		],
		"T87" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T87/T87.bam"
		],
		"T88" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T88/T88.bam"
		],
		"T89" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T89/T89.bam"
		],
		"T9" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T9/T9.bam"
		],
		"T90" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T90/T90.bam"
		],
		"T91" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T91/T91.bam"
		],
		"T92" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T92/T92.bam"
		],
		"T93" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T93/T93.bam"
		],
		"T95" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T95/T95.bam"
		],
		"T96" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T96/T96.bam"
		],
		"T97" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T97/T97.bam"
		],
		"T99" => [
"/gpfs21/scratch/cqs/zhaos/Pierre/20160301_3199_targetDNASeq_hg38/bwa/result/T99/T99.bam"
		],

	},

	bwa_refine_keepDup => {
		class       => "GATK::Refine",
		perform     => 1,
		target_dir  => "${target_dir}/bwa_refine_keepDup",
		option      => "-Xmx40g",
		gatk_option => "",
		fasta_file  => $bwa_fasta_hg38,
		source_ref  => "normalBamFiles",
		indel_vcf_files =>
		  [ $mills_g1000_indel_hg38, $g1000_phase1_indel_hg38 ],
		known_vcf_files  => [$dbsnp_hg38],
		bed_file         => $covered_bed_hg20,
		interval_padding => 250,                 #CNV need to padding 250 bp
		gatk_jar         => $gatk_jar,
		picard_jar       => $picard_jar,
		sh_direct        => 0,
		sorted           => 1,
		slim_print_reads => 0,
		remove_duplicate => 0,
		mark_duplicate   => 1,
		pbs              => {
			"email"    => $email,
			"nodes"    => "1:ppn=8",
			"walltime" => "6",
			"mem"      => "40gb"
		},
	},
	bwa_refine_keepDup_samples => {
		class       => "GATK::Refine",
		perform     => 1,
		target_dir  => "${target_dir}/bwa_refine_keepDup_samples",
		option      => "-Xmx40g",
		gatk_option => "",
		fasta_file  => $bwa_fasta_hg38,
		source_ref  => "sampleBamFiles",
		indel_vcf_files =>
		  [ $mills_g1000_indel_hg38, $g1000_phase1_indel_hg38 ],
		known_vcf_files  => [$dbsnp_hg38],
		bed_file         => $covered_bed_hg20,
		interval_padding => 250,                 #CNV need to padding 250 bp
		gatk_jar         => $gatk_jar,
		picard_jar       => $picard_jar,
		sh_direct        => 0,
		sorted           => 1,
		slim_print_reads => 0,
		remove_duplicate => 0,
		mark_duplicate   => 1,
		pbs              => {
			"email"    => $email,
			"nodes"    => "1:ppn=8",
			"walltime" => "6",
			"mem"      => "40gb"
		},
	},

	gatk4_cnv_PON => {
		class          => "GATK4::CNV",
		perform        => 1,
		target_dir     => "${target_dir}/gatk4_cnv_PON",
		option         => "-Xmx40g",
		gatk_option    => "",
		fasta_file     => $bwa_fasta_hg38,
		source_ref     => "bwa_refine_keepDup",
		gslLibraryFile => $gslLibraryFile,
		hdfViewFolder  => $hdfViewFolder,
		bed_file       => $covered_bed_hg38_padding250_CNV,
		gatk_jar       => $gatk4_jar,
		makePON        => 1,
		sh_direct      => 0,
		pbs            => {
			"email"    => $email,
			"nodes"    => "1:ppn=8",
			"walltime" => "6",
			"mem"      => "40gb"
		},
	},

	gatk4_cnv_samples => {
		class          => "GATK4::CNV",
		perform        => 1,
		target_dir     => "${target_dir}/gatk4_cnv_samples",
		option         => "-Xmx40g",
		gatk_option    => "",
		fasta_file     => $bwa_fasta_hg38,
		source_ref     => "bwa_refine_keepDup_samples",
		gslLibraryFile => $gslLibraryFile,
		hdfViewFolder  => $hdfViewFolder,
		bed_file       => $covered_bed_hg38_padding250_CNV,
		gatk_jar       => $gatk4_jar,
		PanelOfNormal_ref  => ["gatk4_cnv_PON","PON"],
		sh_direct      => 0,
		pbs            => {
			"email"    => $email,
			"nodes"    => "1:ppn=8",
			"walltime" => "6",
			"mem"      => "40gb"
		},
	},
};

#performConfig($config);
performTask($config, "gatk4_cnv_samples");

1;
