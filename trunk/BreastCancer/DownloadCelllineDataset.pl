use strict;
use warnings;

use CQS::FileUtils;
use CQS::EBI;

require "DownloadGEO.pl";

my $targetdir   = "I:/projects/BreastCancer/Cellline/";

my @datasets = ( "E-MTAB-7","E-TABM-244","E-TABM-157","E-MEXP-440","E-MEXP-232","E-MTAB-783" );
foreach my $dataset (@datasets) {
  download_ebi_dataset( $targetdir, $dataset );
}

@datasets = grep(/GSE/, list_directories($targetdir));

DownloadGeoDatasets($targetdir, @datasets);

print "Finished!\n";
