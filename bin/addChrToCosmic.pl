$oldvcf="/data/cqs/shengq1/reference/cosmic/cosmic_v67_20131024.hg19.16571.vcf";
$newvcf="/data/cqs/shengq1/reference/hg19_illumina/illumina-cosmicv67.vcf";

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
  if($_=~m/^MT/){
    $data = substr $_, 2;
    print ILL "chrM$data";
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
