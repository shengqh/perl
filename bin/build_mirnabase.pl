#!/usr/bin/perl
use strict;
use warnings;

`wget ftp://mirbase.org/pub/mirbase/CURRENT/mature.fa.gz; gunzip -f mature.fa.gz`;
`wget ftp://mirbase.org/pub/mirbase/CURRENT/genomes/hsa.gff3`;
`wget ftp://mirbase.org/pub/mirbase/CURRENT/genomes/mmu.gff3`;
`wget ftp://mirbase.org/pub/mirbase/CURRENT/genomes/rno.gff3`;

`sed '/^[^>]/ y/uU/tT/' mature.fa > mature.dna.fa`;

`bwa 2> 1`;
my $bwa = `grep Version 1 | cut -d " " -f 2 | cut -d "-" -f 1`;
chomp($bwa);
`rm 1`;

my $bowtie2 = `bowtie2 --version | grep bowtie2 | grep version | cut -d " " -f 3`;
chomp($bowtie2);
