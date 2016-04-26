#!/usr/bin/perl 
use strict;
use warnings;
use Getopt::Long;

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

open my $in,  "<$input_file"  or die "$input_file: $!";
open my $out, ">$output_file" or die "$output_file: $!";

my $header = <$in>;
$header =~ s/[\r\n]+//;

print $out $header, "\n";

my @headers = split( ',', $header );
my $header_count = scalar(@headers);

sub is_correct_line {
  my $line = shift;
  if ( $line !~ /\"[^,]+\"$/ ) {
    return 0;
  }

  my @parts = split( ',', $line );
  my $part_count = scalar(@parts);
  
  if($part_count == $header_count && $parts[0] =~ /[\d.]+/){
    return 1;
  }elsif($part_count > $header_count){
    die "Long line: " . $line, "\n";
  }else{
    return 0;
  }
}

is_correct_line('38.895572,-76.958072,"NAD83",10,"11","001","0041","44201",1,"Ozone","2006-03-21","22:00","2006-03-22","03:00",2006,81,.025,"Parts per million","1 HOUR","",.005,,"","Equivalent","INSTRUMENTAL-ULTRA VIOLET"') or die;

sub is_insert {
  my $line = shift;
  return $line =~ /[^,\s"]+/;
}

my $lastline;
while ( my $row = <$in> ) {
  $row =~ s/[\r\n]+//;
  if ( $row =~ /END OF FILE/ ) {
    last;
  }

  my $printout = 1;
  while ( !is_correct_line($row) ) {
    #print $row, "\n";
    
    my $nextline1 = <$in>;
    $nextline1 =~ s/[\r\n]+//;
    if ( length($nextline1) == 0 ) {
      if ( is_insert($row) && is_correct_line($lastline) ) {
        $printout = 0;
        last;
      }
      else {
        die "lastline=" . $lastline . "\ncurrentline=" . $row . "\nnextline is empty\n";
      }
    }

    if ( !is_insert($nextline1) ) {
      die "lastline=" . $lastline . "\ncurrentline=" . $row . "\nnextline=" . $nextline1 . "\n";
    }

    my $nextline2 = <$in>;
    $nextline2 =~ s/[\r\n]+//;

    $row = $row . $nextline2;
    
    if($row =~ /^[^.]{4-6}\./){
      die "lastline=" . $lastline . "\ncurrentline=" . $row . "\nnextline1=" . $nextline1. "\nnextline2=" . $nextline2 . "\n";
    }
  }

  if ($printout) {
    print $out $row, "\n";
  }
  $lastline = $row;
}
close $in;
close $out;

print "done\n";
