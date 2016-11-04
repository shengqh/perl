#!/usr/bin/perl
use strict;
use warnings;

use strict;
use warnings;
use Getopt::Long;
use CQS::FileUtils;

my $usage = "

Synopsis:

addChrToVcf.pl -i input_vcf -o output_vcf

Options:

  -i|--input {FILE}          Input vcf file without 'chr'
  -o|--output{FILE}          Output vcf file with 'chr'

  -h|--help                    This page.

";

Getopt::Long::Configure('bundling');

my $input;
my $output;
my $help;

GetOptions( 'h|help' => \$help, 'i|input=s' => \$input, 'o|output=s' => \$output);

if ( defined $help ) {
  print $usage;
  exit(1);
}

if ( !defined($input) || !defined($output) ) {
  print $usage;
  exit(1);
}

open(my $out, ">$output");
open(my $in, $input);
while(<$in>){
  if($_=~m/^#/){
    print $out $_;
  }
  else{
    print $out "chr$_";
  }
}
close($in);
close($out);
