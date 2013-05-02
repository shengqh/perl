#!/usr/bin/perl 
use strict;
use warnings;
use Getopt::Long;
use CQS::FileUtils;

my $usage = "

Synopsis:

perl filtervcf.pl -i inputfile -o outputfile

Options:

  -i|--input {file}           Specify input filename
  -o|--output {file}          Specify output filename

  -h|--help                    This page.

";

Getopt::Long::Configure('bundling');

my $inputfile;
my $outputfile;
my $help;

GetOptions( 'h|help' => \$help, 'i|input=s' => \$inputfile, 'o|output=s' => \$outputfile);

if ( defined $help ) {
  print $usage;
  exit(1);
}

if ( !defined($inputfile)) {
  print $usage;
  exit(1);
}

if ( !defined($outputfile) ) {
  $outputfile = change_extension($inputfile, ".filtered.vcf");
}

`cat $inputfile | awk '(\$1 ~ \"#\") || (\$7==\"PASS\")' > $outputfile`;

1;