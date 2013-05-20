#!/usr/bin/perl 
use strict;
use warnings;
use Getopt::Long;
use CQS::FileUtils;

my $usage = "

Synopsis:

extract_count_file.pl -r rootdir -o outputfile

Options:

  -r|--root {dirroot}          Root of task directory
  -o|--out {file}              Specify output filename

  -h|--help                    This page.

";

Getopt::Long::Configure('bundling');

my $dirroot;
my $outfile;
my $help;

GetOptions( 'h|help' => \$help, 'r|root=s' => \$dirroot, 'o|out=s' => \$outfile );

if ( defined $help ) {
  print $usage;
  exit(1);
}

if ( !defined($dirroot) || !defined($outfile) ) {
  print $usage;
  exit(1);
}

print "dirroot = $dirroot \n";
print "outfile = $outfile \n";

my @subdirs = list_directories($dirroot);
my %ids;
my %samples;
my %data;
foreach my $subdir (@subdirs) {
  my $path = "${dirroot}/${subdir}";
  my $f    = "${path}/${subdir}.sam.count";

  if ( !-s $f ) {
    next;
  }

  print "Reading $f...\n";
  open( IN, $f ) or die $!;
  $samples{$subdir} = 1;
  while (<IN>) {
    s/\r|\n//g;
    my ( $g, $count ) = split "\t";
    $data{$g}->{$subdir} = $count;
  }
}

open( OUT, ">$outfile" ) || die $!;
print OUT join ",", ( "miRNA", sort( keys %samples ) );
print OUT "\n";
foreach my $g ( sort keys %data ) {
  print OUT $g;
  foreach my $s ( sort ( keys %samples ) ) {
    if ( exists( $data{$g}->{$s} ) ) {
      print OUT ",", $data{$g}->{$s};
    }
    else {
      print OUT ",0";
    }
  }
  print OUT "\n";
}
close OUT;

print "done\n";
