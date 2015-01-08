#!/usr/bin/perl
use strict;
use warnings;

use CQS::SmallRNA;

my $def = {
	#General options
	task_name  => "2795",
	email      => "quanhu.sheng\@vanderbilt.edu",
	target_dir => "/scratch/cqs/shengq1/vickers/201501_smallRNA_2795_rat/",
	max_thread => 8,

	#Data
	files => {
	    "2795-KCV-1" => ["/gpfs21/scratch/vantage_repo/Vickers/2795/2795-KCV-1_1_sequence.txt.gz"],
	},
};

performSmallRNARat($def);

1;

