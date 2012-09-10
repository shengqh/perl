#!/usr/local/bin/perl

use CQS::Download;

my $clinicalurl = 'https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/lusc/bcr/biotab/clin/';
my $clinicaldir = 'I:/projects/LUSC/data/clinical/';

download_files($clinicalurl, $clinicaldir);

1;