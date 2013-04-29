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

foreach my $subdir (@subdirs){
  print $subdir;
  my $path = $dirroot . "/" . $subdir;
  print `cat $path/*.stat|grep "in total ("| cut -d ' ' -f1`;
  print `cat $path/*.stat|grep "mapped ("| cut -d ' ' -f1`;
  print "\n";
}

print "done\n";