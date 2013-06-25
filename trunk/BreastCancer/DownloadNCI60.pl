
use strict;
use warnings;

require "DownloadGEO.pl";

my @datasets = (
  "GSE5846"
);

my $targetdir   = "H:/shengquanhu/projects/miRNA/NCI60/microarray_";

DownloadGeoDatasets($targetdir, @datasets);

print "Finished!\n";
