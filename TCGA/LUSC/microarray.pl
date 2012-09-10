#!/usr/local/bin/perl

use CQS::Download;

my $url = 'https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/lusc/cgcc/broad.mit.edu/ht_hg-u133a/transcriptome/';
my $dir = 'I:/projects/LUSC/data/microarray/';

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

download_dirs($url, $dir, $microarrayfilter);

print "\ndownload finished.\n";

1;