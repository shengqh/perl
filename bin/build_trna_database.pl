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

  run_command("wget $url ");

  open( DATA, "<${filename}" );

  #  open(DB, ">$dbFile");
  while (<DATA>) {
    if ( $_ =~ /^chr/ ) {
      print "$_";
      my @parts  = split(/\s+/,  $_);
      my $start  = $parts[3];
      my $end    = $parts[4];
      my $strand = '+';
      if ( $start > $end ) {
        $strand = '-';
        $start  = $parts[4];
        $end    = $parts[3];
      }
      print $parts[1], "\t", $start, "\t", $end, "\t", $parts[1], ".tRNA", $parts[2], "-", $parts[5], $parts[6], "\t1000\t", $strand, "\n";
    }
  }
  close(DATA);
}

1;
