#!/usr/local/bin/perl

require "./common/common.pl";

my $proteinurl = 'https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/brca/cgcc/mdanderson.org/mda_rppa_core/protein_exp/mdanderson.org_BRCA.MDA_RPPA_Core.Level_3.1.0.0/';
my $proteindir = 'D:/sqh/projects/TNBC/data/protein/';

download_files($proteinurl, $proteindir);

1;