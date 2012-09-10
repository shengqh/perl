#!/usr/local/bin/perl

use strict;

use LWP::Simple;
use HTML::Entities;

sub mychomp {
  my ($result) = @_;
  chomp($result);
  return ($result);
}

sub return_unique {
  my @array = @_;
  my %count;
  map { $count{$_} = 1 } @array;
  return sort keys(%count);
}

my $id = "10033900";

print $id . "\n";
my $url =
  "http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=" . $id;
my $content = get($url);
if ( defined($content) ) {
  if ( $content =~ m/Discordant Genotypes:.+?(<TR.+?)<\/TABLE>/sg ) {
    my $table = $1;
    print $table;
  }
}

print "Finished!\n";

1;
