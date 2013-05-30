#!/usr/bin/perl
use strict;
use warnings;

chdir("/data/cqs/shengq1/reference/miRBase19");

my @dbs = ( "mature", "hairpin" );
my @species = ( "hsa", "mmu", "rno" );

foreach my $spec (@species) {
  if ( !-s "${spec}.gff3" ) {
    `wget ftp://mirbase.org/pub/mirbase/CURRENT/genomes/${spec}.gff3`;
  }
}

my $bowtie2 = `bowtie2 --version | grep bowtie2 | grep version | cut -d " " -f 3`;
chomp($bowtie2);

my $files = {};
foreach my $db (@dbs) {
  if ( !-s "${db}.fa" ) {
    `wget ftp://mirbase.org/pub/mirbase/CURRENT/${db}.fa.gz; gunzip -f ${db}.fa.gz `;
  }
  if ( !-s "${db}.dna.fa" ) {
    `cat ${db}.fa | perl -lane 'unless(/^>/){s/U/T/g;} print' >${db}.dna.fa `;
  }
  $files->{"${db}.dna.fa"} = "${db}.dna";

  foreach my $spec (@species) {
    `cat ${db}.dna.fa | perl -lane 'if(/^>/) { \$flag=0 if \$_ !~ />${spec}/; \$flag=1 if /^>${spec}/;} print if \$flag' >${spec}.${db}.dna.fa `;
    $files->{"${spec}.${db}.dna.fa"} = "${spec}.${db}.dna";
  }
}

if ( !-e "bowtie2_index_${bowtie2}" ) {
  mkdir("bowtie2_index_${bowtie2}");
}
chdir("bowtie2_index_${bowtie2}");
my %filemap = %{$files};
for my $file ( keys %filemap ) {
  my $name = $filemap{$file};
  if ( !-e $file ) {
    print "ln -s ../${file} $file \n";
    `ln -s ../${file} $file `;
  }
  `bowtie2-build $file $name `;
}

1;
