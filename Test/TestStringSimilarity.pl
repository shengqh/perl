#!/usr/bin/perl
use strict;
use warnings;
use String::Similarity;
use Test::More;
use Math::Round;

my $s1 = "CAATCAGCAAGTATACTGCCCT";
my $s2 = "AATCAGCAAGTATACTGCCCT";
my $similarity12 = similarity ($s1, $s2);
my $mismatch12 = round(length($s1) * (1-$similarity12) * 2);
ok($mismatch12 eq 1);

my $s3 = "AATCAGCAAGTATACTGCCCTA";
my $similarity13 = similarity ($s1, $s3);
my $mismatch13 = round(length($s1) * (1-$similarity13) * 2);
ok($mismatch13 eq 2);

done_testing()