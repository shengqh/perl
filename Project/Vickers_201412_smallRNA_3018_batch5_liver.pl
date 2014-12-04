#!/usr/bin/perl
use strict;
use warnings;

use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;

my $root = create_directory_or_die(
	"/scratch/cqs/shengq1/vickers/vickers_201412_smallRNA_3018_batch5_liver/");
my $cqstools = "/home/shengq1/cqstools/CQS.Tools.exe";

my $email = "quanhu.sheng\@vanderbilt.edu";

my $samtools           = "/home/shengq1/local/bin/samtools/samtools";
my $bowtie1_option_1mm = "-a -m 100 --best --strata -v 1 -p 8";

my $mirnacount_option          = "-s";                    #ignore score
my $trnacount_option           = "--length --sequence";
my $mirna_overlap_count_option = "-s --gtf_key miRNA";
my $mirna_fasta = "/data/cqs/shengq1/reference/miRBase20/mature.dna.fa";

my $bowtie1_option_pm = "-a -m 100 --best --strata -v 0 -p 8";

my $def = {
	task_name        => "3018_b2",
	mirna_coordinate => "/data/cqs/shengq1/reference/miRBase20/hsa.gff3",
	trna_coordinate =>
	  "/data/cqs/guoy1/reference/smallrna/hg19_tRNA_ucsc_ensembl.bed",
	trna_fasta =>
	  "/data/cqs/guoy1/reference/smallrna/hg19_tRNA_ucsc_ensembl.bed.fa",
	smallrna_coordinate =>
	  "/data/cqs/guoy1/reference/smallrna/hg19_smallRNA_ucsc_ensembl.bed",
	bowtie1_index =>
	  "/data/cqs/guoy1/reference/hg19/bowtie_index_hg19_rCRS_1.0.0/hg19_rCRS",
	target_dir => $root,
	files      => {
		"3018-KCV-5-1" => [
"/gpfs21/scratch/cqs/shengq1/vickers/vickers_201412_smallRNA_3018_batch5_liver/raw/3018-KCV-5-1_ATCACG_L005_R1_001.fastq.gz"
		],
		"3018-KCV-5-2" => [
"/gpfs21/scratch/cqs/shengq1/vickers/vickers_201412_smallRNA_3018_batch5_liver/raw/3018-KCV-5-2_CGATGT_L005_R1_001.fastq.gz"
		],
		"3018-KCV-5-3" => [
"/gpfs21/scratch/cqs/shengq1/vickers/vickers_201412_smallRNA_3018_batch5_liver/raw/3018-KCV-5-3_TTAGGC_L005_R1_001.fastq.gz"
		],
		"3018-KCV-5-4" => [
"/gpfs21/scratch/cqs/shengq1/vickers/vickers_201412_smallRNA_3018_batch5_liver/raw/3018-KCV-5-4_TGACCA_L005_R1_001.fastq.gz"
		],
		"3018-KCV-5-5" => [
"/gpfs21/scratch/cqs/shengq1/vickers/vickers_201412_smallRNA_3018_batch5_liver/raw/3018-KCV-5-5_ACAGTG_L005_R1_001.fastq.gz"
		],
		"3018-KCV-5-6" => [
"/gpfs21/scratch/cqs/shengq1/vickers/vickers_201412_smallRNA_3018_batch5_liver/raw/3018-KCV-5-6_GCCAAT_L005_R1_001.fastq.gz"
		],
		"3018-KCV-5-7" => [
"/gpfs21/scratch/cqs/shengq1/vickers/vickers_201412_smallRNA_3018_batch5_liver/raw/3018-KCV-5-7_CAGATC_L005_R1_001.fastq.gz"
		],
	},
};

my $target_dir = $def->{target_dir};
my $config     = {
	general => { "task_name" => $def->{task_name}, },
	files   => $def->{files},
	trimmer => {
		class      => "CQS::FastqTrimmer",
		perform    => 1,
		target_dir => "${target_dir}/fastq_remove_N",
		option     => "-n -z",
		extension  => "_trim.fastq.gz",
		source_ref => "files",
		cqstools   => $cqstools,
		sh_direct  => 1,
		pbs        => {
			"email"    => $email,
			"nodes"    => "1:ppn=1",
			"walltime" => "2",
			"mem"      => "10gb"
		},
	},
	fastqc_pre => {
		class      => "QC::FastQC",
		perform    => 1,
		target_dir => "${target_dir}/fastqc_pre_trim",
		option     => "",
		source_ref => "trimmer",
		pbs        => {
			"email"    => $email,
			"nodes"    => "1:ppn=1",
			"walltime" => "2",
			"mem"      => "10gb"
		},
	},
	cutadapt => {
		class      => "Cutadapt",
		perform    => 1,
		target_dir => "${target_dir}/cutadapt",
		option     => "-O 10 -m 16",
		source_ref => "trimmer",
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
	fastqc_post => {
		class      => "QC::FastQC",
		perform    => 1,
		target_dir => "${target_dir}/fastqc_post_trim",
		option     => "",
		source_ref => [ "cutadapt", ".fastq.gz" ],
		pbs        => {
			"email"    => $email,
			"nodes"    => "1:ppn=1",
			"walltime" => "2",
			"mem"      => "10gb"
		},
	},
	fastqlen => {
		class      => "FastqLen",
		perform    => 1,
		target_dir => "${target_dir}/fastq_len",
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
	identical => {
		class      => "FastqIdentical",
		perform    => 1,
		target_dir => "${target_dir}/identical",
		option     => "",
		source_ref => [ "cutadapt", ".fastq.gz" ],
		cqstools   => $cqstools,
		extension  => "_clipped_identical.fastq.gz",
		sh_direct  => 1,
		pbs        => {
			"email"    => $email,
			"nodes"    => "1:ppn=1",
			"walltime" => "24",
			"mem"      => "20gb"
		},
	},
	identical_NTA => {
		class        => "CQS::FastqMirna",
		perform      => 1,
		target_dir   => "${target_dir}/identical_NTA",
		option       => "",
		source_ref   => [ "identical", ".fastq.gz\$" ],
		seqcount_ref => [ "identical", ".dupcount\$" ],
		cqstools     => $cqstools,
		extension    => "_clipped_identical_NTA.fastq.gz",
		sh_direct    => 1,
		pbs          => {
			"email"    => $email,
			"nodes"    => "1:ppn=1",
			"walltime" => "24",
			"mem"      => "20gb"
		},
	},

	#not identical, for IGV
	bowtie1_genome_cutadapt_topN_1mm_notidentical => {
		class   => "Bowtie1",
		perform => 1,
		target_dir =>
		  "${target_dir}/topN_bowtie1_genome_cutadapt_1mm_notidentical",
		option        => $bowtie1_option_1mm,
		source_ref    => [ "cutadapt", ".fastq.gz\$" ],
		bowtie1_index => $def->{bowtie1_index},
		samonly       => 0,
		sh_direct     => 0,
		pbs           => {
			"email"    => $email,
			"nodes"    => "1:ppn=6",
			"walltime" => "72",
			"mem"      => "40gb"
		},
	},

	#1 mismatch search, NTA
	bowtie1_genome_cutadapt_topN_1mm_NTA => {
		class         => "Bowtie1",
		perform       => 1,
		target_dir    => "${target_dir}/topN_bowtie1_genome_cutadapt_1mm_NTA",
		option        => $bowtie1_option_1mm,
		source_ref    => [ "identical_NTA", ".fastq.gz\$" ],
		bowtie1_index => $def->{bowtie1_index},
		samonly       => 0,
		sh_direct     => 1,
		pbs           => {
			"email"    => $email,
			"nodes"    => "1:ppn=6",
			"walltime" => "72",
			"mem"      => "40gb"
		},
	},
	mirna_1mm_count_NTA => {
		class   => "MirnaCount",
		perform => 1,
		target_dir =>
		  "${target_dir}/topN_bowtie1_genome_cutadapt_1mm_NTA_mirna_count",
		option          => $mirnacount_option,
		source_ref      => "bowtie1_genome_cutadapt_topN_1mm_NTA",
		fastq_files_ref => "identical_NTA",
		seqcount_ref    => [ "identical_NTA", ".dupcount\$" ],
		cqs_tools       => $cqstools,
		gff_file        => $def->{mirna_coordinate},
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
	miRNA_1mm_table_NTA => {
		class   => "CQS::CQSMirnaNTATable",
		perform => 1,
		target_dir =>
		  "${target_dir}/topN_bowtie1_genome_cutadapt_1mm_NTA_mirna_table",
		option     => "",
		source_ref => [ "mirna_1mm_count_NTA", ".mapped.xml" ],
		cqs_tools  => $cqstools,
		prefix     => "miRNA_1mm_NTA_",
		sh_direct  => 1,
		pbs        => {
			"email"    => $email,
			"nodes"    => "1:ppn=1",
			"walltime" => "10",
			"mem"      => "10gb"
		},
	},

	#1 mismatch search
	bowtie1_genome_cutadapt_topN_1mm => {
		class         => "Bowtie1",
		perform       => 1,
		target_dir    => "${target_dir}/topN_bowtie1_genome_cutadapt_1mm",
		option        => $bowtie1_option_1mm,
		source_ref    => [ "identical", ".fastq.gz\$" ],
		bowtie1_index => $def->{bowtie1_index},
		samonly       => 0,
		sh_direct     => 1,
		pbs           => {
			"email"    => $email,
			"nodes"    => "1:ppn=6",
			"walltime" => "72",
			"mem"      => "40gb"
		},
	},
	miRNA_1mm_count_overlap => {
		class   => "CQSMappedCount",
		perform => 1,
		target_dir =>
		  "${target_dir}/topN_bowtie1_genome_cutadapt_1mm_count_miRNA_overlap",
		option          => $mirna_overlap_count_option,
		source_ref      => "bowtie1_genome_cutadapt_topN_1mm",
		fastq_files_ref => "identical",
		seqcount_ref    => [ "identical", ".dupcount\$" ],
		cqs_tools       => $cqstools,
		gff_file        => $def->{mirna_coordinate},
		fasta_file      => $mirna_fasta,
		samtools        => $samtools,
		sh_direct       => 1,
		pbs             => {
			"email"    => $email,
			"nodes"    => "1:ppn=1",
			"walltime" => "72",
			"mem"      => "20gb"
		},
	},
	miRNA_1mm_overlap_position => {
		class   => "CQSMappedPosition",
		perform => 1,
		target_dir =>
"${target_dir}/topN_bowtie1_genome_cutadapt_1mm_count_miRNA_overlap_position",
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
		class   => "CQSMappedCount",
		perform => 1,
		target_dir =>
		  "${target_dir}/topN_bowtie1_genome_cutadapt_1mm_count_tRNA",
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
			"mem"      => "20gb"
		},
	},
	tRNA_1mm_table => {
		class   => "CQSMappedTable",
		perform => 1,
		target_dir =>
		  "${target_dir}/topN_bowtie1_genome_cutadapt_1mm_count_tRNA_table",
		option     => "",
		source_ref => [ "tRNA_1mm_count", ".xml" ],
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
		class   => "CQSMappedPosition",
		perform => 1,
		target_dir =>
		  "${target_dir}/topN_bowtie1_genome_cutadapt_1mm_count_tRNA_position",
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
		class   => "CQSMappedCount",
		perform => 1,
		target_dir =>
		  "${target_dir}/topN_bowtie1_genome_cutadapt_1mm_count_smallRNA",
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
			"mem"      => "20gb"
		},
	},
	smallRNA_1mm_table => {
		class   => "CQSMappedTable",
		perform => 1,
		target_dir =>
		  "${target_dir}/topN_bowtie1_genome_cutadapt_1mm_count_smallRNA_table",
		option     => "",
		source_ref => [ "smallRNA_1mm_count", ".xml" ],
		cqs_tools  => $cqstools,
		prefix     => "smallRNA_1mm_",
		sh_direct  => 1,
		pbs        => {
			"email"    => $email,
			"nodes"    => "1:ppn=1",
			"walltime" => "10",
			"mem"      => "10gb"
		},
	},
	smallRNA_1mm_NTA_category => {
		class   => "CQSSmallRNACategory",
		perform => 1,
		target_dir =>
"${target_dir}/topN_bowtie1_genome_cutadapt_1mm_NTA_smallRNA_category",
		option          => "",
		source_ref      => [ "smallRNA_1mm_count", ".mapped.xml\$" ],
		mirna_count_ref => [ "mirna_1mm_count_NTA", ".mapped.xml\$" ],
		cqs_tools       => $cqstools,
		sh_direct       => 1,
		pbs             => {
			"email"    => $email,
			"nodes"    => "1:ppn=1",
			"walltime" => "72",
			"mem"      => "40gb"
		},
	},

	#2 perfect match search to mirbase only
	bowtie1_genome_cutadapt_topN_genome_pmnames => {
		class      => "Samtools::PerfectMappedReadNames",
		perform    => 1,
		target_dir => "${target_dir}/topN_bowtie1_genome_cutadapt_pmnames",
		option     => "",
		source_ref => "bowtie1_genome_cutadapt_topN_1mm",
		sh_direct  => 1,
		pbs        => {
			"email"    => $email,
			"nodes"    => "1:ppn=6",
			"walltime" => "72",
			"mem"      => "40gb"
		},
	},
	bowtie1_genome_cutadapt_topN_miRbase_pm => {
		class      => "Alignment::Bowtie1",
		perform    => 1,
		target_dir => "${target_dir}/topN_bowtie1_genome_cutadapt_miRbase_pm",
		option     => $bowtie1_option_pm,
		source_ref => [ "identical", ".fastq.gz\$" ],
		bowtie1_index =>
		  "/data/cqs/shengq1/reference/miRBase21/bowtie_index_1.0.1/mature.dna",
		samonly   => 0,
		sh_direct => 1,
		pbs       => {
			"email"    => $email,
			"nodes"    => "1:ppn=6",
			"walltime" => "72",
			"mem"      => "40gb"
		},
	},
	chromosome_count => {
		class   => "CQS::CQSChromosomeCount",
		perform => 1,
		target_dir =>
		  "${target_dir}/topN_bowtie1_genome_cutadapt_miRbase_pm_count",
		option       => "",
		source_ref   => "bowtie1_genome_cutadapt_topN_miRbase_pm",
		seqcount_ref => [ "identical", ".dupcount\$" ],
		perfect_mapped_name_ref =>
		  "bowtie1_genome_cutadapt_topN_genome_pmnames",
		cqs_tools => $cqstools,
		samtools  => $samtools,
		sh_direct => 1,
		pbs       => {
			"email"    => $email,
			"nodes"    => "1:ppn=6",
			"walltime" => "72",
			"mem"      => "40gb"
		},
	},
	chromosome_count_table => {
		class   => "CQSChromosomeTable",
		perform => 1,
		target_dir =>
		  "${target_dir}/topN_bowtie1_genome_cutadapt_miRbase_pm_table",
		option     => "",
		source_ref => [ "chromosome_count", ".xml" ],
		cqs_tools  => $cqstools,
		prefix     => "miRBase_pm_",
		sh_direct  => 1,
		pbs        => {
			"email"    => $email,
			"nodes"    => "1:ppn=1",
			"walltime" => "10",
			"mem"      => "10gb"
		},
	},
	sequencetask => {
		class      => "CQS::SequenceTask",
		perform    => 1,
		target_dir => "${target_dir}/sequencetask",
		option     => "",
		source     => {
			individual => [

				#data preparation
				"trimmer", "fastqc_pre", "cutadapt", "fastqc_post", "fastqlen",
				"identical", "identical_NTA",

				#NTA data analysis
				"bowtie1_genome_cutadapt_topN_1mm_NTA", "mirna_1mm_count_NTA",

				#non-NTA count
				"bowtie1_genome_cutadapt_topN_1mm", "miRNA_1mm_count_overlap",
				"tRNA_1mm_count",                   "smallRNA_1mm_count",

				#miRBase analysis
				"bowtie1_genome_cutadapt_topN_genome_pmnames",
				"bowtie1_genome_cutadapt_topN_miRbase_pm", "chromosome_count",

				#for IGV
				"bowtie1_genome_cutadapt_topN_1mm_notidentical",
			],
			summary => [

				#NTA table
				"miRNA_1mm_table_NTA",

				#non-NTA table and graph
				"miRNA_1mm_overlap_position", "tRNA_1mm_table",
				"smallRNA_1mm_table",         "smallRNA_1mm_NTA_category",
				"tRNA_1mm_position",

				#miRBase
				"chromosome_count_table"
			],
		},
		sh_direct => 0,
		pbs       => {
			"email"    => $email,
			"nodes"    => "1:ppn=6",
			"walltime" => "72",
			"mem"      => "40gb"
		},
	},
};

performConfig($config);

1;

