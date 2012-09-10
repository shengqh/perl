use strict;
use warnings;

use CQS::FileUtils;

require "DownloadGEO.pl";

my $targetdir   = "I:/projects/BreastCancer/Dataset/";
my @datasets = grep(/GSE/, list_directories($targetdir));

DownloadGeoDatasets(@datasets);

print "Finished!\n";
