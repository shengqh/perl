#!/bin/bash

chrs=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y)

cd /scratch/cqs/shengq1/references/hg20

if [[ ! -s hg20.fa ]]; then
  for i in "${chrs[@]}"
    do
      wget ftp://ftp.ncbi.nlm.nih.gov/genbank/genomes/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh38/Primary_Assembly/assembled_chromosomes/FASTA/chr${i}.fa.gz
      gunzip chr${i}.fa.gz
      cat chr${i}.fa >> hg20.fa
      rm chr${i}.fa
    done
fi

perl ~/program/perl/bin/buildindex.pl -f hg20.fa

