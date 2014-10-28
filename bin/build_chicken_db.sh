#!/bin/bash

cd /scratch/cqs/shengq1/references/chicken

if [[ ! -s gga_ref_Gallus_gallus-4.0.fa ]]; then
  wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Gallus_gallus/Assembled_chromosomes/seq/gga_ref_Gallus_gallus-4.0_chr1.fa.gz 
  wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Gallus_gallus/Assembled_chromosomes/seq/gga_ref_Gallus_gallus-4.0_chr2.fa.gz 
  wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Gallus_gallus/Assembled_chromosomes/seq/gga_ref_Gallus_gallus-4.0_chr3.fa.gz 
  wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Gallus_gallus/Assembled_chromosomes/seq/gga_ref_Gallus_gallus-4.0_chr4.fa.gz 
  wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Gallus_gallus/Assembled_chromosomes/seq/gga_ref_Gallus_gallus-4.0_chr5.fa.gz 
  wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Gallus_gallus/Assembled_chromosomes/seq/gga_ref_Gallus_gallus-4.0_chr6.fa.gz 
  wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Gallus_gallus/Assembled_chromosomes/seq/gga_ref_Gallus_gallus-4.0_chr7.fa.gz 
  wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Gallus_gallus/Assembled_chromosomes/seq/gga_ref_Gallus_gallus-4.0_chr8.fa.gz 
  wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Gallus_gallus/Assembled_chromosomes/seq/gga_ref_Gallus_gallus-4.0_chr9.fa.gz 
  wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Gallus_gallus/Assembled_chromosomes/seq/gga_ref_Gallus_gallus-4.0_chr10.fa.gz 
  wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Gallus_gallus/Assembled_chromosomes/seq/gga_ref_Gallus_gallus-4.0_chr11.fa.gz 
  wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Gallus_gallus/Assembled_chromosomes/seq/gga_ref_Gallus_gallus-4.0_chr12.fa.gz 
  wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Gallus_gallus/Assembled_chromosomes/seq/gga_ref_Gallus_gallus-4.0_chr13.fa.gz 
  wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Gallus_gallus/Assembled_chromosomes/seq/gga_ref_Gallus_gallus-4.0_chr14.fa.gz 
  wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Gallus_gallus/Assembled_chromosomes/seq/gga_ref_Gallus_gallus-4.0_chr15.fa.gz 
  wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Gallus_gallus/Assembled_chromosomes/seq/gga_ref_Gallus_gallus-4.0_chr16.fa.gz 
  wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Gallus_gallus/Assembled_chromosomes/seq/gga_ref_Gallus_gallus-4.0_chr17.fa.gz 
  wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Gallus_gallus/Assembled_chromosomes/seq/gga_ref_Gallus_gallus-4.0_chr18.fa.gz 
  wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Gallus_gallus/Assembled_chromosomes/seq/gga_ref_Gallus_gallus-4.0_chr19.fa.gz 
  wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Gallus_gallus/Assembled_chromosomes/seq/gga_ref_Gallus_gallus-4.0_chr20.fa.gz 
  wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Gallus_gallus/Assembled_chromosomes/seq/gga_ref_Gallus_gallus-4.0_chr21.fa.gz 
  wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Gallus_gallus/Assembled_chromosomes/seq/gga_ref_Gallus_gallus-4.0_chr22.fa.gz 
  wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Gallus_gallus/Assembled_chromosomes/seq/gga_ref_Gallus_gallus-4.0_chr23.fa.gz 
  wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Gallus_gallus/Assembled_chromosomes/seq/gga_ref_Gallus_gallus-4.0_chr24.fa.gz 
  wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Gallus_gallus/Assembled_chromosomes/seq/gga_ref_Gallus_gallus-4.0_chr25.fa.gz 
  wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Gallus_gallus/Assembled_chromosomes/seq/gga_ref_Gallus_gallus-4.0_chr26.fa.gz 
  wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Gallus_gallus/Assembled_chromosomes/seq/gga_ref_Gallus_gallus-4.0_chr27.fa.gz 
  wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Gallus_gallus/Assembled_chromosomes/seq/gga_ref_Gallus_gallus-4.0_chr28.fa.gz 
  wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Gallus_gallus/Assembled_chromosomes/seq/gga_ref_Gallus_gallus-4.0_chr32.fa.gz 
  wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Gallus_gallus/Assembled_chromosomes/seq/gga_ref_Gallus_gallus-4.0_chrW.fa.gz
  wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Gallus_gallus/Assembled_chromosomes/seq/gga_ref_Gallus_gallus-4.0_chrZ.fa.gz
  wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Gallus_gallus/Assembled_chromosomes/seq/gga_ref_Gallus_gallus-4.0_chrMT.fa.gz
  wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Gallus_gallus/Assembled_chromosomes/seq/gga_ref_Gallus_gallus-4.0_chrLGE22C19W28_E50C23.fa.gz
  wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Gallus_gallus/Assembled_chromosomes/seq/gga_ref_Gallus_gallus-4.0_chrLGE64.fa.gz
  
  gunzip *.fa.gz 
fi  
