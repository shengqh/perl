#!/usr/local/bin/perl

require "./common/common.pl";

my $snpurl = 'https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/brca/cgcc/broad.mit.edu/genome_wide_snp_6/snp/';
my $snpdir = 'D:/sqh/projects/TNBC/data/snp/';

my $snpfilter = sub{
  my @url_array = @_;
  
  my %result = ();
  
  foreach my $link(@url_array)
  {
    if (!($link =~ /\/$/)){
      next;
    }
    
    if ($link =~ /(.+)[0-9]{4}\.[0-9]\/$/) 
    {
      my $keylink = $1;
      if (exists($result{$keylink})){
        my $oldlink = $result{$keylink};
        if ($oldlink lt $link){
          print "$oldlink => $link \n";
          $result{$1} = $link;
        }
        else{
          print "kept $oldlink than $link \n";
        }
      }
      else{
        $result{$1} = $link;
      }
    }
    else{
      next;
    }
  }
  
  my @myresult=();
  
  foreach $value (values %result){
    push(@myresult, $value);  
  }
  
  return (@myresult);  
};

download_dirs($snpurl,$snpdir,$snpfilter);

1;