cd /data/cqs/shengq1/reference/cosmic
version="v69"

grep "^#" CosmicCodingMuts.vcf > VCF_Header 
grep -v "^#" CosmicCodingMuts.vcf > Coding.clean 
grep -v "^#" CosmicNonCodingVariants.vcf > NonCoding.clean 
cat Coding.clean NonCoding.clean | sort -gk 2,2 | perl /scratch/cqs/shengq1/local/bin/gatk/public/perl/sortByRef.pl --k 1 - /data/cqs/guoy1/reference/hg19/hg19_chr.fa.fai > Cosmic.hg19.16571
cat VCF_Header Cosmic.hg19.16571 > cosmic_${version}.hg19.16571.vcf
cat Coding.clean NonCoding.clean | sort -gk 2,2 | perl /scratch/cqs/shengq1/local/bin/gatk/public/perl/sortByRef.pl --k 1 - /scratch/cqs/shengq1/somaticmutation_comparison/tcga_fasta/hg19_TCGA_rCRS.fa.fai > Cosmic.hg19.16569MT
cat VCF_Header Cosmic.hg19.16569MT > cosmic_${version}.hg19.16569MT.vcf
cat Coding.clean NonCoding.clean | sort -gk 2,2 | perl /scratch/cqs/shengq1/local/bin/gatk/public/perl/sortByRef.pl --k 1 - /data/cqs/guoy1/reference/hg19/hg19_rCRS/hg19_rCRS.fa.fai > Cosmic.hg19.16569M
cat VCF_Header Cosmic.hg19.16569M > cosmic_${version}.hg19.16569M.vcf
rm Coding.clean Cosmic.hg19.* NonCoding.clean VCF_Header

