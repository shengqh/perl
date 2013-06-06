version="v65_28052013"
wget ftp://ngs.sanger.ac.uk/production/cosmic/CosmicNonCodingVariants_${version}_noLimit.vcf.gz
wget ftp://ngs.sanger.ac.uk/production/cosmic/CosmicCodingMuts_${version}_noLimit.vcf.gz
gunzip Cosmic*.gz
grep "^#" CosmicCodingMuts_${version}_noLimit.vcf > VCF_Header
grep -v "^#" CosmicCodingMuts_${version}_noLimit.vcf > Coding.clean
grep -v "^#" CosmicNonCodingVariants_${version}_noLimit.vcf > NonCoding.clean
cat Coding.clean NonCoding.clean | sort -gk 2,2 | perl /scratch/cqs/shengq1/local/bin/gatk/public/perl/sortByRef.pl --k 1 - /data/cqs/guoy1/reference/hg19/hg19_chr.fa.fai > Cosmic.hg19.16571
cat VCF_Header Cosmic.hg19.16571 > cosmic_${version}.hg19.16571.vcf
cat Coding.clean NonCoding.clean | sort -gk 2,2 | perl /scratch/cqs/shengq1/local/bin/gatk/public/perl/sortByRef.pl --k 1 - /data/cqs/guoy1/reference/hg19/hg19_rCRS/hg19_rCRS.fa.fai > Cosmic.hg19.16569
cat VCF_Header Cosmic.hg19.16569 > cosmic_${version}.hg19.16569.vcf
rm Coding.clean Cosmic.hg19.* NonCoding.clean VCF_Header CosmicNonCodingVariants_${version}_noLimit.vcf CosmicCodingMuts_${version}_noLimit.vcf
