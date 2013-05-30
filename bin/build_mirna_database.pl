#!/usr/bin/perl
use strict;
use warnings;

chdir("/data/cqs/shengq1/reference/miRBase19");

my @dbs = ( "mature", "hairpin" );
my @species = ( "hsa", "mmu", "rno" );

sub run_command {
  my $command = shift;
  print "$command \n";
  `$command `;
}

foreach my $spec (@species) {
  if ( !-s "${spec}.gff3" ) {
    run_command("wget ftp://mirbase.org/pub/mirbase/CURRENT/genomes/${spec}.gff3 ");
  }
}

my $files = {};
foreach my $db (@dbs) {
  if ( !-s "${db}.fa" ) {
    run_command("wget ftp://mirbase.org/pub/mirbase/CURRENT/${db}.fa.gz; gunzip -f ${db}.fa.gz ");
  }
  if ( !-s "${db}.dna.fa" ) {
    run_command("cat ${db}.fa | perl -lane 'unless(/^>/){s/U/T/g;} print' >${db}.dna.fa ");
  }
  $files->{"${db}.dna.fa"} = "${db}.dna";

  foreach my $spec (@species) {
    run_command("cat ${db}.dna.fa | perl -lane 'if(/^>/) { \$flag=0 if \$_ !~ />${spec}/; \$flag=1 if /^>${spec}/;} print if \$flag' >${spec}.${db}.dna.fa ");
    $files->{"${spec}.${db}.dna.fa"} = "${spec}.${db}.dna";
  }
}
my %filemap = %{$files};

my $bowtie2 = `bowtie2 --version | grep bowtie2 | grep version | cut -d " " -f 3`;
chomp($bowtie2);
if ( !-e "bowtie2_index_${bowtie2}" ) {
  mkdir("bowtie2_index_${bowtie2}");
}
chdir("bowtie2_index_${bowtie2}");
for my $file ( sort keys %filemap ) {
  my $name = $filemap{$file};
  if ( !-e $file ) {
    print "ln -s ../${file} $file \n";
    run_command("ln -s ../${file} $file ");
  }
  run_command("bowtie2-build $file $name ");
}

`bwa 2> 1`;
my $bwa = `grep Version 1 | cut -d " " -f 2 | cut -d "-" -f 1`;
chomp($bwa);
`rm 1`;
if ( !-e "bwa_index_${bwa}" ) {
  mkdir("bwa_index_${bwa}");
}
chdir("bwa_index_${bwa}");
for my $file ( keys %filemap ) {
  my $name = $filemap{$file};
  if ( !-e $file ) {
    print "ln -s ../${file} $file \n";
    `ln -s ../${file} $file `;
  }
  print "bwa index $file \n";
  `bwa index $file `;
}

1;
