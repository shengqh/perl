#!/usr/local/bin/perl

use CQS::Download;

my $dirfilter = sub{
  my @url_array = @_;
  
  my @myresult=();
  
  foreach my $link(@url_array)
  {
    if ($link =~ /Level_[12]/){
      next;
    }
    push(@myresult, $link);
  }
  
  return (@myresult);  
};

my $filefilter = sub{
  my @url_array = @_;
  
  my @myresult=();
  
  foreach my $link(@url_array)
  {
    if ($link =~ /gene.quantification.txt$/){
      push(@myresult, $link);
      next;
    }

    if ($link =~ /sdrf.txt$/){
      push(@myresult, $link);
      next;
    }

    if ($link =~ /genes.results$/){
      push(@myresult, $link);
      next;
    }
  }
  
  return (@myresult);  
};


my $rnav1url = 'https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/ov/cgcc/bcgsc.ca/illuminahiseq_rnaseq/rnaseq/';
my $rnav1dir = 'I:/projects/OV/data/rnaseqv1/';
download_dirs($rnav1url, $rnav1dir, $dirfilter, $filefilter);

print "\ndownload rnaseq finished.\n";

1;