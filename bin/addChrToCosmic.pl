open(REF, "/data/cqs/shengq1/reference/cosmic/cosmic_v67_20131024.hg19.16571.vcf");
open(ILL, ">/data/cqs/shengq1/reference/hg19_illumina/illumina-cosmicv67.vcf");
while(<REF>){
	if($_=~m/^#/){
		print ILL $_;
    	}
    	else{
		if($_=~m/^MT/){
			$data = substr $_, 2;
			print ILL "chrM$data";
    		}
   		else{
			print ILL "chr$_";
		}
	}
}
close(ILL);
