#!/usr/bin/perl 
use strict;
use warnings;
use File::Copy;
use Getopt::Long;
use CQS::FileUtils;

my $usage = "

Synopsis:

perl extract_comparison.pl -i input_dir -o output_dir

Options:

  -i|--input {dir}             Specify input directory including all cuffdiff result directories 
  -o|--output {dir}            Specify output directory
  -g|--geneonly                Gene only (default false)
  -f|--foldchange {double}     Filter by absolute value of log2(fold_change)
  -s|--significant             Filter by significant and foldchange (only work with -f)
  -d|--different {double}      Filter by difference of value between two samples (only work with -f and only check with the one has inf or -inf fold change)

  -h|--help                    This page.

";

Getopt::Long::Configure('bundling');

my $input_root;
my $output_dir;
my $help;
my $gene_only;
my $fold_change;
my $significant;
my $different;

GetOptions(
  'h|help'         => \$help,
  'g|geneonly'     => \$gene_only,
  'i|input=s'      => \$input_root,
  'o|output=s'     => \$output_dir,
  'f|foldchange=s' => \$fold_change,
  's|significant'  => \$significant,
  's|different=s'  => \$different,
);

if ( defined $help ) {
  print $usage;
  exit(1);
}

if ( !defined($input_root) || !defined($output_dir) ) {
  print $usage;
  exit(1);
}

$gene_only = 0 if !defined($gene_only);

die "directory $input_root is not exists" if !-d $input_root;

create_directory_or_die($output_dir) if !-d $output_dir;

my @subdirs = list_directories($input_root);
if ( 0 == scalar(@subdirs) ) {
  die "$input_root has no sub CuffDiff directories";
}

my @filenames;
if ($gene_only) {
  @filenames = ("gene_exp.diff");
  $gene_only = "(\$3 != \"-\") && ";
}
else {
  @filenames = ( "gene_exp.diff", "splicing.diff" );
  $gene_only = "";
}

my $filter;
if ( defined($fold_change) ) {
  $filter = "(\$10 != \"inf\" && \$10 != \"-inf\" && (\$10 <= -$fold_change || \$10 >= $fold_change))";

  if ( defined($different) ) {
    $filter = "($filter || (\$10 == \"inf\" &&  (\$9-\$8) >= $different) || (\$10 == \"-inf\" &&  (\$8-\$9) >= $different))";
  }

  if ( defined($significant) ) {
    $filter = $filter . " && (\$14==\"yes\")";
  }
}
else {
  $filter = "(\$14==\"yes\")";
}

for my $subdir (@subdirs) {
  foreach my $filename (@filenames) {
    my $file = "${input_root}/${subdir}/${filename}";
    if ( -s $file ) {
      print "Processing " . $file . "\n";
      open IN, "<$file" or die "Cannot open file $file";
      my $line = <IN>;
      $line = <IN>;
      close(IN);

      if ( defined($line) ) {
        my @parts      = split( /\t/, $line );
        my $partcount  = scalar(@parts);
        my $targetname = "${output_dir}/${subdir}.${filename}";

        copy( $file, $targetname ) or die "copy failed : $!";

        my $target_sign_name = $targetname . ".sig";
        my $cmd              = "cat $targetname | awk '$gene_only ((\$14==\"significant\") || ($filter))' > $target_sign_name";

        print $cmd . "\n";
        `$cmd`;
      }
    }
  }
}

1;
