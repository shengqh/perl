#!/usr/bin/perl 
use strict;
use warnings;
use Getopt::Long;
use Text::CSV;

my $usage = "

Synopsis:

epadata -i input_file -o output_file

Options:

  -i|--input {file}            Input EPA file
  -o|--out {file}              Preprocessed EPA file
  -h|--help                    This page.

";

Getopt::Long::Configure('bundling');

my $input_file;
my $output_file;
my $help;

GetOptions( 'h|help' => \$help, 'i|input=s' => \$input_file, 'o|out=s' => \$output_file );

if ( defined $help ) {
  print $usage;
  exit(1);
}

if ( !defined($input_file) || !defined($output_file) ) {
  print $usage;
  exit(1);
}

print "input_file = $input_file \n";
print "output_file = $output_file \n";

my @rows;
my $csv = Text::CSV->new( { binary => 1 } )    # should set binary attribute.
  or die "Cannot use CSV: " . Text::CSV->error_diag();

open my $in,  "<:encoding(utf8)", $input_file  or die "$input_file: $!";
open my $out, ">:encoding(utf8)", $output_file or die "$output_file: $!";
while ( my $row = $csv->getline($in) ) {
  $row->[0] !~ m/^END OF FILE/ or last;
  $csv->print( $out, $row );
}
close $in;
close $out;

print "done\n";
