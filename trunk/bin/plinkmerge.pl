#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my $usage = "

Synopsis:

plinkmerge -f firstfile_root -m listfile -o outputfile_root

Options:

  -f|--file {fileroot}         Specify .ped and .map files
  -m|--merge-list {list file}  Merge multiple standard and/or binary filesets
  -o|--out {fileroot}          Specify output root filename

  -h|--help                    This page.

";

Getopt::Long::Configure('bundling');

my $firstfile;
my $listfile;
my $outfile;
my $help;

GetOptions( 'h|help' => \$help, 'f|file=s' => \$firstfile, 'm|merge-list=s' => \$listfile, 'o|out=s' => \$outfile );

if ( defined $help ) {
  print $usage;
  exit(1);
}

if ( !defined($firstfile) || !defined($listfile) || !defined($outfile) ) {
  print $usage;
  exit(1);
}

print "firstfile = $firstfile \n";
print "listfile = $listfile \n";
print "outfile = $outfile \n";

my $firstped = $firstfile . ".ped";
my $firstmap = $firstfile . ".map";
my $outped   = $outfile . ".ped";
my $outmap   = $outfile . ".map";

open LIST_FILE, $listfile or die $!;

my @beds = ();
my @maps = ();

push( @maps, $firstmap );
while (<LIST_FILE>) {
  my @parts = split(/\s+/);
  if ( scalar(@parts) > 1 ) {
    push( @beds, $parts[0] );
    push( @maps, $parts[1] );
    print "$parts[0] $parts[1] \n";
  }
}

print "merging bed files to $outped ... \n";

my @bedfiles = ();
open FIRST_FILE, $firstped || die $!;
my $filecount = scalar(@beds);
for ( my $index = 0 ; $index < $filecount ; $index++ ) {
  local *FILE;
  open FILE, $beds[$index] || die $!;
  push( @bedfiles, *FILE );
}

open( OUT_FILE, '>$outped' ) || die $!;
while ( my $line = <FIRST_FILE> ) {
  chomp($line);
  if ( $line =~ m/^\S+\s\S+\s\S+\s\S+\s\S+\s\S+\s.+$/ ) {
    print OUT_FILE $line;
    foreach my $file (@bedfiles) {
      my $otherline = <$file>;
      chomp($otherline);
      $otherline =~ m/^\S+\s\S+\s\S+\s\S+\s\S+\s\S+\s(.+)$/;
      print OUT_FILE ' $1';
    }
    print OUT_FILE "\n";
  }
}
close(OUT_FILE);
close(FIRST_FILE);
foreach my $file (@bedfiles) {
  close($file);
}

print "merging map files to $outmap ... \n";

open( OUT_FILE, '>$outmap' ) || die $!;
foreach my $mapfile (@maps) {
  open FILE, $mapfile || die $!;
  while ( my $line = <FILE> ) {
    print OUT_FILE $line;
  }
  close(FILE);
}
close(OUT_FILE);

print "done\n";