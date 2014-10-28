#!/bin/bash

chrs=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 Pltd MT)

cd /scratch/cqs/shengq1/references/soybean

if [[ ! -s gma_ref_V1.1.fa ]]; then
  for i in "${chrs[@]}"
    do
      wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Glycine_max/Assembled_chromosomes/seq/gma_ref_V1.1_chr{i}.fa.gz
      gunzip gma_ref_V1.1_chr{i}.fa.gz
      awk '{print $1; exit}' gma_ref_V1.1_chr{i}.fa | cut -f4 -d '|' | awk '{print ">"$1}' >> gma_ref_V1.1.fa
      awk 'NR>1' gma_ref_V1.1_chr{i}.fa >> gma_ref_V1.1.fa
      rm gma_ref_V1.1_chr{i}.fa
    done
fi

wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Glycine_max/GFF/ref_V1.1_scaffolds.gff3.gz
gunzip ref_V1.1_scaffolds.gff3.gz

perl ~/program/perl/bin/buildindex.pl -f gma_ref_V1.1.fa

