#!/usr/bin/perl
use strict;
use warnings;

use LWP::UserAgent;
use LWP::Simple;
use File::Basename;

use CQS::SystemUtils;
sub run_command {
  my $command = shift;
  print "$command \n";
  `$command `;
}

if ( is_linux() ) {
  chdir("/data/cqs/shengq1/reference/trna");
  run_command("wget http://gtrnadb.ucsc.edu/download/tRNAs/GtRNAdb-all-tRNAs.fa.gz; gunzip -f GtRNAdb-all-tRNAs.fa.gz ");
}
else{
  chdir("e:/test");
}

my $species = {
  "hg19" => "http://gtrnadb.ucsc.edu/Hsapi19/Hsapi19-by-locus-txt.html",
  "mm10" => "http://gtrnadb.ucsc.edu/Mmusc10/Mmusc10-by-locus-txt.html",
  "rn4"  => "http://gtrnadb.ucsc.edu/Rnorv/Rnorv-by-locus-txt.html"
};

my $files = {};
foreach my $db ( sort keys %{$species} ) {
  my $dbFile = "${db}_tRNA.bed";
  my $url   = $species->{$db};
  my $filename = basename($url);
  
  run_command("wget $url ");
  
  open(DATA, "<${filename}");
   
#  open(DB, ">$dbFile");
  while (<DATA>) {
    if ($_ =~ /^chr/) {
      print "$_";
      if ($_ =~ /(\S+(chr\S+)\.\S+)\s+\((\d+)-(\d+)/){
        print $1, "\t", $2, "\t", $3, "\t", $4, "\n";
      }
    };
  }
  close(DATA);
}

1;
