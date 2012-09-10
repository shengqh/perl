use strict;
use warnings;

use CQS::FileUtils;
use File::Spec::Functions qw(all);

my $dir  = "D:/sqh/projects/BreastCancer/Dataset";
my @ebis = ("E-TABM-158/E-TABM-158.sample.txt");

foreach my $ebi (@ebis) {
  my ( $volume, $directories, $file ) = splitpath($dir . "/" . $ebi);
  $ebi =~ /(.+)\//;
  my $dataset = $1;
  my $samplefile = $volume . $directories . "/${dataset}.sample";
  
  
}

