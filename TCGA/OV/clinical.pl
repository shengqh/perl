#!/usr/local/bin/perl

use CQS::Download;

my $clinicalurl = 'https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/ov/bcr/biotab/clin/';
my $clinicaldir = 'I:/projects/OV/data/clinical/';

download_files($clinicalurl, $clinicaldir);

1;