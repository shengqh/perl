use strict;
use warnings;

use CQS::GEO;
use CQS::FileUtils;

my @training = ( "GSE5851", "GSE18088", "GSE13067", "GSE13294", "GSE14333", "GSE26682", "GSE26906", "GSE28702", "GSE33113", "GSE37892", "GSE39582" );

my @testing = (
      "GSE17536",
      #"GSE17537",
      #"GSE17538"
);

my @celllines = ( "GSE8332", "GSE36133" );

my $targetdir = "H:/shengquanhu/projects/chenxi/20130724-chenxi_colon_cancer";

#download_geo_ids( "${targetdir}/training/",  "", @training );
download_geo_ids( "${targetdir}/testing/",   "", @testing );
#download_geo_ids( "${targetdir}/celllines/", "", @celllines );

print "Finished!\n";
