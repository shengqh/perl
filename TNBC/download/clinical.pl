#!/usr/local/bin/perl

require "./common/common.pl";

my $clinicalurl = 'https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/brca/bcr/nationwidechildrens.org/biotab/clin/';
my $clinicaldir = 'D:/sqh/projects/TNBC/data/clinical/';

download_files($clinicalurl, $clinicaldir);

1;