#!/usr/bin/perl
use strict;
use warnings;
use Test::Simple;

my $number = 54;
my $numlen = length($number);
my $actual = sprintf("%0${numlen}d", 1);

ok($actual eq "01");
