#!/usr/local/bin/perl

use strict;
use CQS::Download;

my $level2url = 'https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/laml/cgcc/genome.wustl.edu/hg-u133_plus_2/transcriptome/genome.wustl.edu_LAML.HG-U133_Plus_2.Level_2.1.2.0/';
my $level2dir = 'I:/projects/LAML/data/microarray/genome.wustl.edu_LAML.HG-U133_Plus_2.Level_2.1.2.0/';

download_files($level2url, $level2dir);

my $celfilter = sub{
  my @url_array = @_;
  
  my @myresult=();
  
  foreach my $link(@url_array)
  {
    if (!($link =~ /.CEL$/)){
      next;
    }
    push(@myresult, $link);
  }
  
  return (@myresult);  
};


my $level1url = 'https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/laml/cgcc/genome.wustl.edu/hg-u133_plus_2/transcriptome/genome.wustl.edu_LAML.HG-U133_Plus_2.Level_1.1.2.0/';
my $level1dir = 'I:/projects/LAML/data/microarray/genome.wustl.edu_LAML.HG-U133_Plus_2.Level_1.1.2.0/';

download_files($level1url, $level1dir, $celfilter);

1;