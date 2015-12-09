#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;
use CQS::FileUtils;

my $target_dir = create_directory_or_die("/scratch/cqs/shengq1/chipseq/20151208_gse53998");

my $fasta_file           = "/data/cqs/guoy1/reference/hg19/bwa_index_0.7.4/hg19_chr.fa";
my $transcript_gtf       = "/scratch/cqs/shengq1/references/hg19/Homo_sapiens.GRCh37.73.gtf";
my $transcript_gtf_index = "/scratch/cqs/shengq1/gtfindex/hg19_GRCh37_73";

my $hg19_gff = "/scratch/cqs/shengq1/references/hg19/dexseq_gff/Homo_sapiens.GRCh37.73.dexseq.gff";
my $hg19_map = "/scratch/cqs/shengq1/references/hg19/Homo_sapiens.GRCh37.73.map";

my $bowtie2_index = "/data/cqs/guoy1/reference/hg19/bowtie2_index/hg19";

my $annovar_param = "-protocol refGene,snp137,cosmic64,esp6500si_all,1000g2012apr_all -operation g,f,f,f,f --remove --otherinfo";
my $annovar_db    = "/scratch/cqs/shengq1/references/annovar/humandb/";

my $cqstools = "/home/shengq1/cqstools/CQS.Tools.exe";

my $email = "quanhu.sheng\@vanderbilt.edu";
my $task  = "gse53998";

my $config = {
	general => { task_name => $task },
	files   => {
		"SRR1106515" => [ "/gpfs21/scratch/cqs/shengq1/chipseq/20151208_gse53998/sra/SRR1106515.sra" ],
		"SRR1106516" => [ "/gpfs21/scratch/cqs/shengq1/chipseq/20151208_gse53998/sra/SRR1106516.sra" ],
		"SRR1106517" => [ "/gpfs21/scratch/cqs/shengq1/chipseq/20151208_gse53998/sra/SRR1106517.sra" ],
		"SRR1106518" => [ "/gpfs21/scratch/cqs/shengq1/chipseq/20151208_gse53998/sra/SRR1106518.sra" ],
		"SRR1106519" => [ "/gpfs21/scratch/cqs/shengq1/chipseq/20151208_gse53998/sra/SRR1106519.sra" ],
		"SRR1106520" => [ "/gpfs21/scratch/cqs/shengq1/chipseq/20151208_gse53998/sra/SRR1106520.sra" ],
		"SRR1106521" => [ "/gpfs21/scratch/cqs/shengq1/chipseq/20151208_gse53998/sra/SRR1106521.sra" ],
		"SRR1106522" => [ "/gpfs21/scratch/cqs/shengq1/chipseq/20151208_gse53998/sra/SRR1106522.sra" ],
		"SRR1106523" => [ "/gpfs21/scratch/cqs/shengq1/chipseq/20151208_gse53998/sra/SRR1106523.sra" ],
		"SRR1106524" => [ "/gpfs21/scratch/cqs/shengq1/chipseq/20151208_gse53998/sra/SRR1106524.sra" ],
		"SRR1106525" => [ "/gpfs21/scratch/cqs/shengq1/chipseq/20151208_gse53998/sra/SRR1106525.sra" ],
		"SRR1106526" => [ "/gpfs21/scratch/cqs/shengq1/chipseq/20151208_gse53998/sra/SRR1106526.sra" ],
		"SRR1106527" => [ "/gpfs21/scratch/cqs/shengq1/chipseq/20151208_gse53998/sra/SRR1106527.sra" ],
		"SRR1106528" => [ "/gpfs21/scratch/cqs/shengq1/chipseq/20151208_gse53998/sra/SRR1106528.sra" ],
		"SRR1106529" => [ "/gpfs21/scratch/cqs/shengq1/chipseq/20151208_gse53998/sra/SRR1106529.sra" ],
		"SRR1106530" => [ "/gpfs21/scratch/cqs/shengq1/chipseq/20151208_gse53998/sra/SRR1106530.sra" ],
		"SRR1106531" => [ "/gpfs21/scratch/cqs/shengq1/chipseq/20151208_gse53998/sra/SRR1106531.sra" ],
		"SRR1106532" => [ "/gpfs21/scratch/cqs/shengq1/chipseq/20151208_gse53998/sra/SRR1106532.sra" ],
		"SRR1106533" => [ "/gpfs21/scratch/cqs/shengq1/chipseq/20151208_gse53998/sra/SRR1106533.sra" ],
		"SRR1106534" => [ "/gpfs21/scratch/cqs/shengq1/chipseq/20151208_gse53998/sra/SRR1106534.sra" ],
		"SRR1106535" => [ "/gpfs21/scratch/cqs/shengq1/chipseq/20151208_gse53998/sra/SRR1106535.sra" ],
	},
	sra2fastq => {
		class      => "SRA::FastqDump",
		perform    => 1,
		ispaired   => 1,
		target_dir => "${target_dir}/FastqDump",
		option     => "",
		source_ref => "files",
		sh_direct  => 0,
		pbs        => {
			"email"    => $email,
			"nodes"    => "1:ppn=1",
			"walltime" => "10",
			"mem"      => "10gb"
		},
	},
	fastqc => {
		class      => "QC::FastQC",
		perform    => 1,
		target_dir => "${target_dir}/fastqc",
		option     => "",
		source_ref => "sra2fastq",
		sh_direct  => 0,
		pbs        => {
			"email"    => $email,
			"nodes"    => "1:ppn=2",
			"walltime" => "2",
			"mem"      => "10gb"
		},
	},
};

performConfig($config);

1;
