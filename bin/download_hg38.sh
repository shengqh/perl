#!/bin/bash

chrs=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y)

cd /scratch/cqs/shengq1/references/hg38_MT

if [[ ! -s hg38_MT.fa ]]; then
  rm hg38_MT.fa.tmp
  for i in "${chrs[@]}"
    do
      wget ftp://ftp.ncbi.nlm.nih.gov/genbank/genomes/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh38/Primary_Assembly/assembled_chromosomes/FASTA/chr${i}.fa.gz
      gunzip chr${i}.fa.gz
      cat ">" ${i} "\n" >> hg38_MT.fa.tmp
      awk 'NR>1' chr${i}.fa >> hg38_MT.fa.tmp
      rm chr${i}.fa
    done
  cat rCRS.fa >> hg38_MT.fa.tmp
  mv hg38_MT.fa.tmp hg38_MT.fa
fi

perl ~/program/perl/bin/buildindex.pl -f hg38_MT.fa


