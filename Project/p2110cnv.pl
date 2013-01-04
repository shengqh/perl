#!/usr/bin/perl
use strict;
use warnings;

use CQS::CNV;
use CQS::FileUtils;
use CQS::SystemUtils;

my $target_dir = create_directory_or_die("/scratch/cqs/shengq1/dnaseq/2110");

my $email = "quanhu.sheng\@vanderbilt.edu";

my $config = {
	general => {
		path_file => "/home/shengq1/bin/path.txt",
		task_name => "2110"
	},
	bamfiles => {
		"2110-JP-01" => ["/scratch/cqs/shengq1/dnaseq/2110/bwa_markdup/2110-JP-1_realigned_recal_rmdup.bam"],
		"2110-JP-02" => ["/scratch/cqs/shengq1/dnaseq/2110/bwa_markdup/2110-JP-2_realigned_recal_rmdup.bam"],
		"2110-JP-03" => ["/scratch/cqs/shengq1/dnaseq/2110/bwa_markdup/2110-JP-3_realigned_recal_rmdup.bam"],
		"2110-JP-04" => ["/scratch/cqs/shengq1/dnaseq/2110/bwa_markdup/2110-JP-4_realigned_recal_rmdup.bam"],
		"2110-JP-05" => ["/scratch/cqs/shengq1/dnaseq/2110/bwa_markdup/2110-JP-5_realigned_recal_rmdup.bam"],
		"2110-JP-06" => ["/scratch/cqs/shengq1/dnaseq/2110/bwa_markdup/2110-JP-6_realigned_recal_rmdup.bam"],
		"2110-JP-07" => ["/scratch/cqs/shengq1/dnaseq/2110/bwa_markdup/2110-JP-7_realigned_recal_rmdup.bam"],
		"2110-JP-08" => ["/scratch/cqs/shengq1/dnaseq/2110/bwa_markdup/2110-JP-8_realigned_recal_rmdup.bam"],
		"2110-JP-09" => ["/scratch/cqs/shengq1/dnaseq/2110/bwa_markdup/2110-JP-9_realigned_recal_rmdup.bam"],
		"2110-JP-10" => ["/scratch/cqs/shengq1/dnaseq/2110/bwa_markdup/2110-JP-10_realigned_recal_rmdup.bam"],
		"2110-JP-11" => ["/scratch/cqs/shengq1/dnaseq/2110/bwa_markdup/2110-JP-11_realigned_recal_rmdup.bam"],
		"2110-JP-12" => ["/scratch/cqs/shengq1/dnaseq/2110/bwa_markdup/2110-JP-12_realigned_recal_rmdup.bam"],
		"2110-JP-13" => ["/scratch/cqs/shengq1/dnaseq/2110/bwa_markdup/2110-JP-13_realigned_recal_rmdup.bam"],
		"2110-JP-14" => ["/scratch/cqs/shengq1/dnaseq/2110/bwa_markdup/2110-JP-14_realigned_recal_rmdup.bam"],
		"2110-JP-15" => ["/scratch/cqs/shengq1/dnaseq/2110/bwa_markdup/2110-JP-15_realigned_recal_rmdup.bam"],
		"2110-JP-16" => ["/scratch/cqs/shengq1/dnaseq/2110/bwa_markdup/2110-JP-16_realigned_recal_rmdup.bam"],
	},
    cnvnator050 => {
        target_dir => "${target_dir}/cnvnator050",
        option     => "",
        source_ref => "bamfiles",
        binsize    => 50,
        pbs        => {
            "email"    => $email,
            "nodes"    => "1:ppn=1",
            "walltime" => "72",
            "mem"      => "10gb"
        },
    },
	cnvnator100 => {
		target_dir => "${target_dir}/cnvnator100",
		option     => "",
		source_ref => "bamfiles",
        binsize    => 100,
		pbs        => {
			"email"    => $email,
			"nodes"    => "1:ppn=1",
			"walltime" => "72",
			"mem"      => "10gb"
		},
	},
    cnvnator200 => {
        target_dir => "${target_dir}/cnvnator200",
        option     => "",
        source_ref => "bamfiles",
        binsize    => 200,
        pbs        => {
            "email"    => $email,
            "nodes"    => "1:ppn=1",
            "walltime" => "72",
            "mem"      => "10gb"
        },
    },
    cnvnator300 => {
        target_dir => "${target_dir}/cnvnator300",
        option     => "",
        source_ref => "bamfiles",
        binsize    => 300,
        pbs        => {
            "email"    => $email,
            "nodes"    => "1:ppn=1",
            "walltime" => "72",
            "mem"      => "10gb"
        },
    },
};

#cnvnator($config, "cnvnator050");
cnvnator($config, "cnvnator100");
#cnvnator($config, "cnvnator200");
#cnvnator($config, "cnvnator300");

1;
