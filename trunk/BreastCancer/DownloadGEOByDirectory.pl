use strict;
use warnings;

use CQS::FileUtils;

require "DownloadGEO.pl";

my $targetdir   = "D:/projects/BreastCancer/Test/";
my @datasets = grep(/GSE/, list_directories($targetdir));

DownloadGeoDatasets($targetdir, @datasets);

print "Finished!\n";
