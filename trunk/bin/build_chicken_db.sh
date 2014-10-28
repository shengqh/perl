#!/bin/bash

chrs=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 32 W Z MT LGE22C19W28_E50C23 LGE64)

rm gga_ref_Gallus_gallus-4.0.fa

for i in "${chrs[@]}"
    do
        wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Gallus_gallus/Assembled_chromosomes/seq/gga_ref_Gallus_gallus-4.0_chr${i}.fa.gz
        gunzip gga_ref_Gallus_gallus-4.0_chr${i}.fa.gz
        echo ">chr${i}" >> gga_ref_Gallus_gallus-4.0.fa
        awk 'NR>1' gga_ref_Gallus_gallus-4.0_chr${i}.fa >> gga_ref_Gallus_gallus-4.0.fa
        rm gga_ref_Gallus_gallus-4.0_chr${i}.fa
    done
