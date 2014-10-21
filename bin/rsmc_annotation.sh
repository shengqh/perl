#!/bin/bash

for i in `ls -d1 /scratch/cqs/shengq1/somaticmutation_comparison/TCGA_rsmc_positionInRead_RNA/result/*/*NT.tsv`
do
  echo $i
  mono-sgen /home/shengq1/rsmc/rsmc.exe annotation \
    -i $i \
    --annovar --annovar_buildver hg19 \
    --rnaediting --rnaediting_db /data/cqs/shengq1/reference/rnaediting/hg19.txt \
    --distance --distance_exon_gtf /data/cqs/shengq1/reference/ensembl_gtf/Homo_sapiens.GRCh37.75.M.gtf
done


for i in `ls -d1 /scratch/cqs/shengq1/somaticmutation_comparison/TCGA_rsmc_positionInRead_DNA/result/*/*NT.tsv`
do
  echo $i
  mono-sgen /home/shengq1/rsmc/rsmc.exe annotation \
    -i $i \
    --annovar --annovar_buildver hg19 \
    --rnaediting --rnaediting_db /data/cqs/shengq1/reference/rnaediting/hg19.txt \
    --distance --distance_exon_gtf /data/cqs/shengq1/reference/ensembl_gtf/Homo_sapiens.GRCh37.75.M.gtf
done
