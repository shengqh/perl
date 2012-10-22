use strict;
use warnings;

use CQS::EBI;

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
