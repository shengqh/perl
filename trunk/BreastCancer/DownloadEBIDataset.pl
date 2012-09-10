use strict;
use warnings;

use LWP::Simple;

use CQS::Download;
use CQS::HashUtils;

my $filter = sub {
  my @url_array = @_;
  my @myresult  = ();

  foreach my $link (@url_array) {
    if ( $link =~ /(.+\.raw\.\d+\.zip)$/ ) {
      push( @myresult, "http://www.ebi.ac.uk" . $1 );
    }
    elsif ( $link =~ /(.+\.idf.txt)$/ ) {
      push( @myresult, "http://www.ebi.ac.uk" . $1 );
    }
    elsif ( $link =~ /(.+\.sdrf.txt)$/ ) {
      push( @myresult, "http://www.ebi.ac.uk" . $1 );
    }
  }

  return (@myresult);
};

my $makefilename = sub {
  my $oldurl = shift;
  $oldurl =~ /\/([^\/]+)$/;
  return ($1);
};

sub download_ebi_dataset {
  my ( $rootdir, $dataset ) = @_;
  my $url = "http://www.ebi.ac.uk/arrayexpress/experiments/$dataset";
  my $dir = "$rootdir/$dataset";

  if ( !( -d $dir ) ) {
    mkdir $dir;
  }

  download_files( $url, $dir, $filter, $makefilename );
}

#my @datasets = ( "E-MTAB-365", "E-NCMF-3", "E-TABM-158", "E-UCON-1","E-TABM-43" );
my @datasets = ( "E-MTAB-365", "E-TABM-158");
foreach my $dataset (@datasets) {
  download_ebi_dataset( "I:/projects/BreastCancer/Dataset", $dataset )
    ;
}

#my @celllines = ("E-MTAB-7","E-TABM-244","E-TABM-157","E-MEXP-440","E-MEXP-232");
#foreach my $dataset (@celllines) {
#  download_ebi_dataset( "I:/projects/BreastCancer/Cellline", $dataset )
#    ;
#}

1;
