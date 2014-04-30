$oldvcf="/data/cqs/guoy1/reference/dbsnp138/00-All.vcf";
$newvcf="/data/cqs/shengq1/reference/hg19_illumina/illumina-dbsnp138.vcf";

open(ILL, ">$newvcf");

open(REF, $oldvcf);
while(<REF>){
  if($_=~m/^#/){
    print ILL $_;
  }
  else{
    last;
  }
}
close(REF);

open(REF, $oldvcf);
while(<REF>){
  if($_=~m/^M/){
    print ILL "chr$_";
  }
}
close(REF);

open(REF, $oldvcf);
while(<REF>){
  if($_=~m/^#/){
    next;
  }
  else{
    if($_=~m/^MT/){
      next;
    }
    else{
      print ILL "chr$_";
    }
  }
}
close(REF);

close(ILL);
