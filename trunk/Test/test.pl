#!/usr/bin/perl

use strict;
use warnings;

my $idfile = "myfile.txt";

open( DAT, $idfile ) || die("Could not open file $idfile!");
my $line     = <DAT>;
my @raw_data = <DAT>;
close(DAT);

foreach $line (@raw_data) {
  chomp($line);
  my @parts = split( '\t', $line );
  my $outputfile = $parts[0] . ".txt";
  `plink $parts[0] $parts[1] > $outputfile`;
}
