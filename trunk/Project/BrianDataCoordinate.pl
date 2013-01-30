#!/usr/bin/perl
use strict;
use warnings;

use CQS::QC;
use CQS::RNASeq;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::TCGA;

my $target_dir = "/scratch/cqs/shengq1/briandata";

my $email = "quanhu.sheng\@vanderbilt.edu";

my $config = {
    general => {
        path_file            => "/home/shengq1/bin/path.txt",
        task_name            => "tcga"
    },
    tcga => {
        target_dir => "${target_dir}/tcga",
        option     => "",
        idfile =>"${target_dir}/tcga.txt",
        tcgaidindex =>1,
        analysisidindex => 5,
        coordinateindex => 13,
        pbs        => {
            "email"    => $email,
            "nodes"    => "1:ppn=1",
            "walltime" => "48",
            "mem"      => "10gb"
        },
    },
};

tcga_download($config, "tcga");

1;
