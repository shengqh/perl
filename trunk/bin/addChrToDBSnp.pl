open(REF, "/data/cqs/guoy1/reference/dbsnp138/00-All.vcf");
open(ILL, ">/data/cqs/shengq1/reference/hg19_illumina/illumina-dbsnp138.vcf");
while(<REF>){
    if($_=~m/^#/){
	print ILL $_;
    }
    else{
	print ILL "chr$_";
    }
}
close(ILL);
