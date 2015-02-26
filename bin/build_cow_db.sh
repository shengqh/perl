#!/bin/bash

chrs=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 X)

cd /scratch/cqs/shengq1/references/cow

if [[ ! -s Bos_taurus_UMD_3.1.fa ]]; then
  rm Bos_taurus_UMD_3.1.fa.tmp
  for i in "${chrs[@]}"
    do
      wget ftp://ftp.ncbi.nih.gov/genbank/genomes/Eukaryotes/vertebrates_mammals/Bos_taurus/Bos_taurus_UMD_3.1/Primary_Assembly/assembled_chromosomes/FASTA/chr${i}.fa.gz
      gunzip chr${i}.fa.gz
      echo ">chr"${i} >> Bos_taurus_UMD_3.1.fa.tmp
      awk 'NR>1' chr${i}.fa >> Bos_taurus_UMD_3.1.fa.tmp
      rm chr${i}.fa
    done
  mv Bos_taurus_UMD_3.1.fa.tmp Bos_taurus_UMD_3.1.fa
fi

wget ftp://ftp.ensembl.org/pub/release-78/gtf/bos_taurus/Bos_taurus.UMD3.1.78.gtf.gz
gunzip Bos_taurus.UMD3.1.78.gtf.gz

perl ~/program/perl/bin/buildindex.pl -f Bos_taurus_UMD_3.1.fa

