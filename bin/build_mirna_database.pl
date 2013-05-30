#!/usr/bin/perl 
use strict;
use warnings;

chdir("/data/cqs/shengq1/reference/miRBase19");

`wget ftp://mirbase.org/pub/mirbase/CURRENT/genomes/hsa.gff3`;
`wget ftp://mirbase.org/pub/mirbase/CURRENT/genomes/mmu.gff3`;
`wget ftp://mirbase.org/pub/mirbase/CURRENT/genomes/rno.gff3`;

my $bowtie2 = `bowtie2 --version | grep bowtie2 | grep version | cut -d " " -f 3`;
chomp($bowtie2);

my @dbs = ( "mature", "hairpin" );
my @species = ( "hsa", "mmu", "rno" );
my %files = {};
foreach my $db (@dbs) {
  `wget ftp://mirbase.org/pub/mirbase/CURRENT/${db}.fa.gz; gunzip -f ${db}.fa.gz `;
  `cat ${db}.fa | perl -lane 'unless(/^>/){s/U/T/g;} print' >${db}.dna.fa `;
  $files{"${db}.dna.fa"} = "${db}.dna";

  foreach my $spec (@species) {
    `cat ${db}.dna.fa | perl -lane 'if(/^>/) { \$flag=0 if \$_ !~ />${spec}/; \$flag=1 if /^>${spec}/;} print if \$flag' >${spec}.${db}.dna.fa `;
    $files{"${spec}.${db}.dna.fa"} = "${spec}.${db}.dna";
  }
}

`mkdir bowtie2_index_${bowtie2}; cd bowtie2_index_${bowtie2} `;
for my $file ( keys %files ) {
  my $name = $files{$file};
  `ln -s ../${file} $file `;
  `bowtie2-build $file $name `;
}

1;
