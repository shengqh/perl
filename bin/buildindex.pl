#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;

sub run_command {
  my $command = shift;
  print "$command \n";
  `$command `;
}

my $usage = "

Synopsis:

buildindex -f fastaFile

Options:

  -f|--file {fastaFile}        Fasta format sequence file
  -h|--help                    This page.
";

Getopt::Long::Configure('bundling');

my $fastaFile;
my $help;

GetOptions( 'h|help' => \$help, 'f|file=s' => \$fastaFile );

if ( defined $help ) {
  print $usage;
  exit(1);
}

if ( !defined($fastaFile) ) {
  print $usage;
  exit(1);
}

my ($base, $dir, $ext) = fileparse($fastaFile);
my $basename = basename($fastaFile);

my $bowtie2 = `bowtie2 --version | grep bowtie2 | grep version | cut -d " " -f 3`;
chomp($bowtie2);
if ( !-e "bowtie2_index_${bowtie2}" ) {
  mkdir("bowtie2_index_${bowtie2}");
}
chdir("bowtie2_index_${bowtie2}");
if ( !-e $basename ) {
  run_command("ln -s ../$fastaFile $basename ");
}
run_command("bowtie2-build $basename $base ");
chdir("..");

`bwa 2> 1`;
my $bwa = `grep Version 1 | cut -d " " -f 2 | cut -d "-" -f 1`;
chomp($bwa);
`rm 1`;
if ( !-e "bwa_index_${bwa}" ) {
  mkdir("bwa_index_${bwa}");
}
chdir("bwa_index_${bwa}");
if ( !-e $basename ) {
  run_command("ln -s ../$fastaFile $basename ");
}
print "bwa index $basename \n";
run_command("bwa index $basename");

chdir("..");

1;
