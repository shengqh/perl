#!/usr/bin/perl
use strict;
use warnings;

my $rootdir = "/scratch/cqs/shengq1/groseq";

my $exp    = "mtg8";
my @ranges = (
  "chr1:20,850,679-21,106,250", "chr1:55,098,491-55,162,388",
  "chr2:80,092,416-80,563,519", "chr3:121,886,990-122,714,378"
);
my $bamfile = $rootdir . "/" . "2083-SH-uniq-merged.sorted.bam";

foreach my $range (@ranges) {
  my $name = $range;
  $name =~ s/[:,-]//g;
  $name = "$rootdir/${exp}_${name}";

  my $bamfileP = "${name}.p.bam";
  my $bamfileR = "${name}.r.bam";

  system(
"samtools view -u -F 16 $bamfile $range | samtools mpileup - > ${name}.p.basecount"
  );
  system(
"samtools view -u -f 16 $bamfile $range | samtools mpileup - > ${name}.r.basecount"
  );

  my $rfile = "${name}.R";
  print "$rfile\n";
  open( OUT, ">$rfile" ) or die $!;
  print OUT "setwd(\"${rootdir}\"\);\n";
  print OUT "p<-read.table(\"${name}.p.basecount\",header=F,sep=\"\t\");\n";
  print OUT "r<-read.table(\"${name}.r.basecount\",header=F,sep=\"\t\");\n";
  print OUT
"png(filename=\"$name.png\",width=4000,height=3000,units=\"px\",res=300);\n";
  print OUT
"plot(p[,2],p[,4],type=\"h\",col=\"red\",xlab=\"Base\",ylab=\"Read Count\",ylim=c(-50,50));\n";
  print OUT "lines(r[,2],-r[,4],type=\"h\", col=\"blue\");\n";
  print OUT "dev.off();\n";
  close OUT;

  system("R --no-save < $rfile");
}

print "Finished";
