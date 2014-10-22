cd /scratch/cqs/shengq1/references/dbsnp 

wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b142_GRCh37p13/VCF/All.vcf.gz 
gunzip All.vcf.gz  
mv All.vcf human_GRCh37_v142_16569_MT.vcf 
cat human_GRCh37_v142_16569_MT.vcf | awk 'BEGIN {OFS="\t"} {if($1=="MT")$1="M"; print }' > human_GRCh37_v142_16569_M.vcf 

wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b142_GRCh37p13/BED/bed_chr_1.bed.gz 
wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b142_GRCh37p13/BED/bed_chr_2.bed.gz 
wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b142_GRCh37p13/BED/bed_chr_3.bed.gz 
wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b142_GRCh37p13/BED/bed_chr_4.bed.gz 
wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b142_GRCh37p13/BED/bed_chr_5.bed.gz 
wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b142_GRCh37p13/BED/bed_chr_6.bed.gz 
wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b142_GRCh37p13/BED/bed_chr_7.bed.gz 
wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b142_GRCh37p13/BED/bed_chr_8.bed.gz 
wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b142_GRCh37p13/BED/bed_chr_9.bed.gz 
wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b142_GRCh37p13/BED/bed_chr_10.bed.gz 
wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b142_GRCh37p13/BED/bed_chr_11.bed.gz 
wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b142_GRCh37p13/BED/bed_chr_12.bed.gz 
wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b142_GRCh37p13/BED/bed_chr_13.bed.gz 
wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b142_GRCh37p13/BED/bed_chr_14.bed.gz 
wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b142_GRCh37p13/BED/bed_chr_15.bed.gz 
wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b142_GRCh37p13/BED/bed_chr_16.bed.gz 
wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b142_GRCh37p13/BED/bed_chr_17.bed.gz 
wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b142_GRCh37p13/BED/bed_chr_18.bed.gz 
wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b142_GRCh37p13/BED/bed_chr_19.bed.gz 
wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b142_GRCh37p13/BED/bed_chr_20.bed.gz 
wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b142_GRCh37p13/BED/bed_chr_21.bed.gz 
wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b142_GRCh37p13/BED/bed_chr_22.bed.gz 
wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b142_GRCh37p13/BED/bed_chr_X.bed.gz 
wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b142_GRCh37p13/BED/bed_chr_Y.bed.gz 
wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b142_GRCh37p13/BED/bed_chr_MT.bed.gz 
gunzip *.bed.gz 
cat bed_chr_1.bed bed_chr_2.bed bed_chr_3.bed bed_chr_4.bed bed_chr_5.bed bed_chr_6.bed bed_chr_7.bed bed_chr_8.bed bed_chr_9.bed bed_chr_10.bed bed_chr_11.bed bed_chr_12.bed bed_chr_13.bed bed_chr_14.bed bed_chr_15.bed bed_chr_16.bed bed_chr_17.bed bed_chr_18.bed bed_chr_19.bed bed_chr_20.bed bed_chr_21.bed bed_chr_22.bed bed_chr_X.bed bed_chr_Y.bed bed_chr_MT.bed >human_GRCh37_v142_16569_MT.bed 
rm bed_chr*.bed   
cat human_GRCh37_v142_16569_MT.bed | awk 'BEGIN {OFS="\t"} {if($1=="chrMT")$1="chrM"; print }' > human_GRCh38_v142_16569_M.bed 

wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/All.vcf.gz 
gunzip All.vcf.gz 
mv All.vcf human_GRCh38_v142_16569_MT.vcf 
cat human_GRCh38_v142_16569_MT.vcf | awk 'BEGIN {OFS="\t"} {if($1=="MT")$1="M"; print }' > human_GRCh38_v142_16569_M.vcf 
