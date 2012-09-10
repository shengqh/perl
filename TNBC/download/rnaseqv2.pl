#!/usr/local/bin/perl

use CQS::Download;

my $rnaurl = 'https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/brca/cgcc/unc.edu/illuminahiseq_rnaseqv2/rnaseqv2/';
my $rnadir = 'I:/projects/TNBC/data/rnaseqv2/';


download_dirs($rnaurl, $rnadir);

print "download rnaseqv2 finished.";

1;