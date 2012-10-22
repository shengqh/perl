use strict;
use warnings;

use CQS::FileUtils;
use CQS::EBI;

require "DownloadGEO.pl";

my $targetdir   = "D:/projects/BreastCancer/Cellline/";

download_ebi_dataset( $targetdir, "E-MTAB-783" );

my @datasets = grep(/GSE/, list_directories($targetdir));

DownloadGeoDatasets($targetdir, @datasets);

print "Finished!\n";
