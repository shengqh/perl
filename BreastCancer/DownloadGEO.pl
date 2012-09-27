use strict;
use warnings;

use CQS::GEO;
use CQS::FileUtils;

my $excludeddir = "I:/projects/BreastCancer/Excluded/";

sub DownloadGeoDatasets {
  my ($targetdir, @datasets) = @_;

  foreach my $dataset (@datasets) {
    my $localdir       = $targetdir . $dataset;
    my $curexcludeddir = $excludeddir . $dataset;

    if ( !( -d $localdir ) && ( -d $curexcludeddir ) ) {
      next;
    }

    unless ( -d $localdir ) {
      mkdir $localdir or die "failed to create directory $localdir";
    }

    if (
      !has_file(
        $localdir,
        sub {
          my $file = shift;
          return $file =~ /\.cel.gz$/i;
        }
      )
      )
    {
      download_geo_supplementary( $dataset, $localdir );
    }

    download_geo_matrix( $dataset, $localdir );
  }
}

1;
