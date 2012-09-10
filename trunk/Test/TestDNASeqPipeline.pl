#!/usr/bin/perl
use strict;
use warnings;

use CQS::DNASeq;

my $fastqFile1 = "/data/cqs/guoy1/SNPtest2/rawdata/10009_1.fq";
my $fastqFile2 = "/data/cqs/guoy1/SNPtest2/rawdata/10009_2.fq";
my $fastaFile  = "/data/cqs/guoy1/reference/hg19/hg19_chr.fa";
my $rootDir    = "/scratch/cqs/shengq1/ReferenceGenomeComparison";

bwa_by_pbs_single( $fastaFile, $fastqFile1, "$rootDir/test1" );
bwa_by_pbs_double( $fastaFile, $fastqFile1, $fastqFile2, "10009",
  "$rootDir/test2" );
