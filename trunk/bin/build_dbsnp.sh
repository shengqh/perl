cd /data/cqs/shengq1/reference/dbsnp

wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b141_GRCh37p13/VCF/All.vcf.gz
gunzip All.vcf.gz
mv All.vcf human_GRCh37_v141_16569_MT.vcf
cat human_GRCh37_v141_16569_MT.vcf | awk '{if($1=="MT")$1="M"; print }' > human_GRCh37_v141_16569_M.vcf

wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/All.vcf.gz
gunzip All.vcf.gz
mv All.vcf human_GRCh38_v141_16569_MT.vcf
cat human_GRCh38_v141_16569_MT.vcf | awk '{if($1=="MT")$1="M"; print }' > human_GRCh38_v141_16569_M.vcf
