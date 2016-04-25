#!/usr/bin/perl
use strict;
use warnings;

use Test::More;

my $file1 = "/scratch/cqs/shengq1/smallRNA/20150310_guoyan_2687_human_mouse_batch2/human/identical_NTA/result/2687-GG-104_clipped_identical_NTA.fastq.gz";
my $file2 = "/scratch/cqs/shengq1/smallRNA/20150310_guoyan_2687_human_mouse_batch2/human/identical_NTA/result/2687-GG-004_clipped_identical_NTA.fastq.gz";
my $regex = "2687-GG-1.*.fastq.gz\$";
ok( $file1 =~ m/$regex/, "accept file1" );
ok( $file2 !~ m/$regex/, "not accept file2" );

my $values = [
  "/scratch/cqs/shengq1/vickers/20160402_smallRNA_3018-KCV-59_62_63_human/identical_sequence_count_table/pbs/KCV-59_62_63_sequence.filelist",
  "/scratch/cqs/shengq1/vickers/20160402_smallRNA_3018-KCV-59_62_63_human/bowtie1_genome_1mm_NTA_smallRNA_table/pbs/smallRNA_1mm_KCV-59_62_63.filelist",
  "/scratch/cqs/shengq1/vickers/20160402_smallRNA_3018-KCV-59_62_63_human/bowtie1_genome_1mm_NTA_smallRNA_table/result/smallRNA_1mm_KCV-59_62_63.txt"
];

my @newvalues = grep { !/\/pbs\// } @$values;

ok( scalar(@newvalues) == 1, "should only 1" );

done_testing();