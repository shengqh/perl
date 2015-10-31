#!/bin/bash

cd /scratch/cqs/shengq1/references/dbsnp/mm10 

if [[ ! -s mouse_GRCm38_v142_MT.vcf ]]; then
  wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/mouse_10090/VCF/vcf_chr_1.vcf.gz 
  wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/mouse_10090/VCF/vcf_chr_2.vcf.gz 
  wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/mouse_10090/VCF/vcf_chr_3.vcf.gz 
  wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/mouse_10090/VCF/vcf_chr_4.vcf.gz 
  wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/mouse_10090/VCF/vcf_chr_5.vcf.gz 
  wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/mouse_10090/VCF/vcf_chr_6.vcf.gz 
  wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/mouse_10090/VCF/vcf_chr_7.vcf.gz 
  wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/mouse_10090/VCF/vcf_chr_8.vcf.gz 
  wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/mouse_10090/VCF/vcf_chr_9.vcf.gz 
  wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/mouse_10090/VCF/vcf_chr_10.vcf.gz 
  wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/mouse_10090/VCF/vcf_chr_11.vcf.gz 
  wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/mouse_10090/VCF/vcf_chr_12.vcf.gz 
  wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/mouse_10090/VCF/vcf_chr_13.vcf.gz 
  wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/mouse_10090/VCF/vcf_chr_14.vcf.gz 
  wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/mouse_10090/VCF/vcf_chr_15.vcf.gz 
  wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/mouse_10090/VCF/vcf_chr_16.vcf.gz 
  wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/mouse_10090/VCF/vcf_chr_17.vcf.gz 
  wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/mouse_10090/VCF/vcf_chr_18.vcf.gz 
  wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/mouse_10090/VCF/vcf_chr_19.vcf.gz 
  wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/mouse_10090/VCF/vcf_chr_20.vcf.gz 
  wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/mouse_10090/VCF/vcf_chr_X.vcf.gz 
  wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/mouse_10090/VCF/vcf_chr_Y.vcf.gz 
  wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/mouse_10090/VCF/vcf_chr_MT.vcf.gz 
  gunzip *.vcf.gz 
  grep "^#" vcf_chr_1.vcf > mouse_GRCm38_v142_MT.vcf
  cat vcf_chr_1.vcf vcf_chr_2.vcf vcf_chr_3.vcf vcf_chr_4.vcf vcf_chr_5.vcf vcf_chr_6.vcf vcf_chr_7.vcf vcf_chr_8.vcf vcf_chr_9.vcf vcf_chr_10.vcf vcf_chr_11.vcf vcf_chr_12.vcf vcf_chr_13.vcf vcf_chr_14.vcf vcf_chr_15.vcf vcf_chr_16.vcf vcf_chr_17.vcf vcf_chr_18.vcf vcf_chr_19.vcf vcf_chr_X.vcf vcf_chr_Y.vcf vcf_chr_MT.vcf | grep -v "^#" >> mouse_GRCm38_v142_MT.vcf 
  rm vcf_chr*.vcf   
fi

if [[ ! -s mouse_GRCm38_v142_M.vcf ]]; then
  cat mouse_GRCm38_v142_MT.vcf | awk 'BEGIN {OFS="\t"} {if($1=="MT")$1="M"; print }' > mouse_GRCm38_v142_M.vcf
fi 
