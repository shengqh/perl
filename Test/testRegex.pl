#!/usr/bin/perl
use strict;
use warnings;

use Test::Simple tests => 2;

my $file1 = "/scratch/cqs/shengq1/smallRNA/20150310_guoyan_2687_human_mouse_batch2/human/identical_NTA/result/2687-GG-104_clipped_identical_NTA.fastq.gz";
my $file2 = "/scratch/cqs/shengq1/smallRNA/20150310_guoyan_2687_human_mouse_batch2/human/identical_NTA/result/2687-GG-004_clipped_identical_NTA.fastq.gz";
my $regex = "2687-GG-1.*.fastq.gz\$";
ok($file1 =~ m/$regex/, "accept file1");
ok($file2 !~ m/$regex/, "not accept file2");