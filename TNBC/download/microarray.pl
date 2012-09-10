#!/usr/local/bin/perl

use CQS::Download;

my $microarrayurl = 'https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/brca/cgcc/unc.edu/agilentg4502a_07_3/transcriptome/';
my $microarraydir = 'I:/projects/TNBC/data/microarray/';

my $microarrayfilter = sub{
  my @url_array = @_;
  
  my @myresult=();
  
  foreach my $link(@url_array)
  {
    if (!($link =~ /\/$/)){
      next;
    }
    
    if ($link =~ /Level_[12]/) 
    {
      next;
    }
    
    push(@myresult, $link);
  }
  
  return (@myresult);  
};

download_dirs($microarrayurl, $microarraydir, $microarrayfilter);

1;