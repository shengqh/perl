#!/usr/local/bin/perl

use CQS::Download;

my $rnaurl = 'https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/brca/cgcc/unc.edu/illuminahiseq_rnaseq/rnaseq/';
my $rnadir = 'I:/projects/TNBC/data/rnaseqv1/';

download_dirs($rnaurl, $rnadir);

print "\ndownload rnaseqv1 finished.\n";

1;