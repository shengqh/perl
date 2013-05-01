#!/usr/bin/perl 
use strict;
use warnings;

my @dbs={"mature","hairpin"};
foreach my $db (@dbs){
  `cat ${db}.fa | perl -lane 'unless(/^>/){s/U/T/g;} print' >${db}.dna.fa `;
  `bwa index ${db}.dna.fa `;
  `novoindex -m ${db}.dna.nix ${db}.dna.fa `;
}

my @species={"hsa","mmu","rno"};
foreach my $spec (@species){
  `cat mature.dna.fa | perl -lane 'if(/^>/) { \$flag=0 if \$_ !~ />${spec}/; \$flag=1 if /^>${spec}/;} print if \$flag' >${spec}.mature.dna.fa `;
  `bwa index ${spec}.mature.dna.fa `;
  `novoindex -m ${spec}.mature.dna.nix ${spec}.mature.dna.fa `;

  `cat mature.dna.fa | perl -lane 'if(/^>/) { \$flag=1 if \$_ !~ />${spec}/; \$flag=0 if /^>${spec}/;} print if \$flag' >not_${spec}.mature.dna.fa `;
  `bwa index not_${spec}.mature.dna.fa `;
  `novoindex -m not_${spec}.mature.dna.nix not_${spec}.mature.dna.fa `;
}

1;