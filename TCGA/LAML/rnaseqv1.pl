#!/usr/local/bin/perl

use CQS::Download;

my $rnaurl = 'https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/laml/cgcc/bcgsc.ca/illuminaga_rnaseq/rnaseq/bcgsc.ca_LAML.IlluminaGA_RNASeq.Level_3.1.5.0/';
my $rnadir = 'I:/projects/laml/data/rnaseqv1/bcgsc.ca_LAML.IlluminaGA_RNASeq.Level_3.1.5.0/';

my $rnaseqfilter = sub{
  my @url_array = @_;
  
  my @myresult=();
  
  foreach my $link(@url_array)
  {
    if (!($link =~ /gene.quantification.txt$/)){
      next;
    }
    push(@myresult, $link);
  }
  
  return (@myresult);  
};


download_files($rnaurl, $rnadir, $rnaseqfilter);

print "\ndownload rnaseqv1 finished.\n";

1;