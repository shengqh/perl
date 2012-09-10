use strict;
use warnings;

use CQS::FileUtils;

require "SampleGeneInfo.pl";

my $root = "D:/sqh/projects/BreastCancer/Dataset";

my @subdirs = ListDirectories($root);

my $summary = $root . "/summary.txt";

open (my $summaryfile, ">$summary") or die "Cannot create file $summary: $!"; 
print $summaryfile "Dataset\tFileCount\tSampleCount\tGeneCount\tGeneSymbol\n"; 

foreach my $subdir (@subdirs) {
  my $realdir  = $root . "/" . $subdir;
  my @txtfiles = ListFiles(
    $realdir,
    sub {
      my $link = shift;
      return $link =~ /\.txt$/;
    }
  );
  my $totalfile = @txtfiles;
  
  if ($totalfile == 0){
    next;
  }

  my $keypattern;
  my $ignorepattern;

  if ( $subdir =~ /^E-\S{4}-\d+/ ) {    #EBI data
    $keypattern    = "Scan\\sREF";
    $ignorepattern = "\\sREF";
  }
  else {
    if ( ( $subdir =~ "^GSE" ) || ( $subdir =~ "^MDA" ) ) {    #NCBI data
      $keypattern    = "ID_REF";
      $ignorepattern = "^!";
    }
    else {
      die "Don't know how to parse data of $subdir";
    }
  }

  my @totalsamples;
  my @totalgenes;
  my @genesymbols;

  foreach my $file (@txtfiles) {
    my $realfile = $realdir . "/" . $file;
    my %ret = GetSampleGeneInfo($realfile, $keypattern, $ignorepattern);
    
    my $sampleref = $ret{'sample'};
    my $generef = $ret{'gene'};
    
    @totalsamples = (@totalsamples, @$sampleref);
    @totalgenes = (@totalgenes,  @$generef);
    push(@genesymbols, $ret{'symbol'});
  }
  
  my $totalsample = return_unique_count(@totalsamples); 
  my $totalgenes = return_unique_count(@totalgenes);
  my $totalsymbols = join('/', @genesymbols);
  
  print $summaryfile "$subdir\t$totalfile\t$totalsample\t$totalgenes\t$totalsymbols\n"; 
}

close $summaryfile;

print "\nFinished!";

1;
