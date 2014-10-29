#!/bin/bash

chrs=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 32 W Z MT LGE22C19W28_E50C23 LGE64)

cd /scratch/cqs/shengq1/references/chicken

if [[ ! -s gga_ref_Gallus_gallus-4.0.fa ]]; then
  for i in "${chrs[@]}"
    do
      wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Gallus_gallus/Assembled_chromosomes/seq/gga_ref_Gallus_gallus-4.0_chr${i}.fa.gz
      gunzip gga_ref_Gallus_gallus-4.0_chr${i}.fa.gz
      awk '{print $1; exit}' gga_ref_Gallus_gallus-4.0_chr${i}.fa | cut -f4 -d '|' | awk '{print ">"$1}' >> gga_ref_Gallus_gallus-4.0.fa
      awk 'NR>1' gga_ref_Gallus_gallus-4.0_chr${i}.fa >> gga_ref_Gallus_gallus-4.0.fa
      rm gga_ref_Gallus_gallus-4.0_chr${i}.fa
    done
fi

wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Gallus_gallus/GFF/ref_Gallus_gallus-4.0_scaffolds.gff3.gz
gunzip ref_Gallus_gallus-4.0_scaffolds.gff3.gz

perl ~/program/perl/bin/buildindex.pl -f gga_ref_Gallus_gallus-4.0.fa

