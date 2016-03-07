#!/usr/bin/perl
use strict;
use warnings;
use Test::More;
use File::Basename;

my ( $filename, $dirs, $suffix ) = fileparse("/dkjkd/addd/a.bed.c.bed", ".bed\$");

ok($filename eq "a.bed.c");
ok($suffix eq ".bed");

done_testing();