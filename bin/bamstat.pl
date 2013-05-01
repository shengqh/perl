#!/usr/bin/perl 
use strict;
use warnings;
use Getopt::Long;
use CQS::FileUtils;

my $usage = "

Synopsis:

bamstat -r rootdir -o outputfile_root

Options:

  -r|--root {dirroot}          Root of task directory
  -o|--out {fileroot}          Specify output root filename

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
open( OUT_FILE, ">$outfile" ) || die $!;
print OUT_FILE "sample\ttotal_reads\tmapped_reads\tpercentage\n";
foreach my $subdir (@subdirs){
  my $path = $dirroot . "/" . $subdir;
  my $total = `cat $path/*.stat|grep "in total ("| cut -d ' ' -f1`;
  chomp($total);
  my $mapped = `cat $path/*.stat|grep "mapped ("| cut -d ' ' -f1`;
  chomp($mapped);
  print OUT_FILE $subdir . "\t" . $total . "\t" . $mapped . "\t" . ($mapped / $total) . "\n";
}
close (OUT_FILE);

print "done\n";