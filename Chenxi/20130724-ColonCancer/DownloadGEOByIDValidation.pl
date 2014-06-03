use strict;
use warnings;

use CQS::GEO;
use CQS::FileUtils;

my @training = ( "GSE13067", "GSE13294", "GSE14333", "GSE17536", "GSE17537", "GSE18088", "GSE26682", "GSE33113" );

my $targetdir = "H:/shengquanhu/projects/DupChecker/";

download_geo_ids( $targetdir,   "", @training );

print "Finished!\n";
