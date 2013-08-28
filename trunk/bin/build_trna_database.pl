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
}
else {
  chdir("e:/test");
}

my $species = {
  "hg19" => "http://gtrnadb.ucsc.edu/Hsapi19/Hsapi19-by-locus-txt.html",
  "mm10" => "http://gtrnadb.ucsc.edu/Mmusc10/Mmusc10-by-locus-txt.html",
  "rn4"  => "http://gtrnadb.ucsc.edu/Rnorv/Rnorv-by-locus-txt.html"
};

my $files = {};
foreach my $db ( sort keys %{$species} ) {
  my $dbFile   = "${db}_tRNA.bed";
  my $url      = $species->{$db};
  my $filename = basename($url);

  if(! -e $filename){
    run_command("wget $url ");
  }

  open( DATA, "<${filename}" );

  open(DB, ">$dbFile");
  while (<DATA>) {
    if ( $_ =~ /^chr/ ) {
      print "$_";
      my @parts  = split(/\s+/,  $_);
      my $start  = $parts[2];
      my $end    = $parts[3];
      my $strand = '+';
      if ( $start > $end ) {
        $strand = '-';
        $start  = $parts[3];
        $end    = $parts[2];
      }
      $start = $start -1;
      print DB $parts[0], "\t", $start, "\t", $end, "\t", $parts[0], ".tRNA", $parts[1], "-", $parts[4], $parts[5], "\t1000\t", $strand, "\n";
    }
  }
  close(DATA);
  close(DB);
}

1;
