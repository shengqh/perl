cd /data/cqs/shengq1/reference/cosmic
version="v69"

grep "^#" CosmicCodingMuts.vcf > VCF_Header 
grep -v "^#" CosmicCodingMuts.vcf > Coding.clean 
grep -v "^#" CosmicNonCodingVariants.vcf > NonCoding.clean 
cat Coding.clean NonCoding.clean | sort -gk 2,2 | perl /scratch/cqs/shengq1/local/bin/gatk/public/perl/sortByRef.pl --k 1 - /scratch/cqs/shengq1/somaticmutation_comparison/tcga_fasta/hg19_TCGA_rCRS.fa.fai > Cosmic.hg19.16569MT
cat VCF_Header Cosmic.hg19.16569MT > cosmic_${version}_hg19_16569_MT.vcf
cat cosmic_v69_hg19_16569_MT.vcf | awk '{if($1=="MT")$1="M"; print }' > cosmic_v69_hg19_16569_M.vcf

rm Cosmic.hg19.16569MT Coding.clean NonCoding.clean VCF_Header
