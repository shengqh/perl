#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my $usage = "

Synopsis:

qcsummary -f fastqc_dir -r RNASeQC_dir -o outputfile

Options:

  -f|--fastqc {directory}           Fastqc directory
  -r|--rnaseqc {directory}          RNASeQC directory
  -o|--out {file}                   Output filename

  -h|--help                         This page.

";

Getopt::Long::Configure('bundling');

my $fastqcdir;
my $rnaseqcdir;
my $outfile;
my $help;

GetOptions( 'h|help' => \$help, 'f|fastqc=s' => \$fastqcdir, 'r|rnaseqc=s' => \$rnaseqcdir, 'o|out=s' => \$outfile );

if ( defined $help ) {
  print $usage;
  exit(1);
}

if ( !defined($fastqcdir) or !defined($outfile) ) {
  print $usage;
  exit(1);
}

print "fastqcdir = $fastqcdir \n";
print "rnaseqcdir = $rnaseqcdir \n";
print "outfile = $outfile \n";


