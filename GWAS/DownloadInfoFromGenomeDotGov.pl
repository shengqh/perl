#!/usr/local/bin/perl

use strict;

use LWP::Simple;
use HTML::Entities;
use LWP::UserAgent;

my $ua = new LWP::UserAgent;
$ua->agent( "AgentName/0.1 " . $ua->agent );

sub mychomp {
  my ($result) = @_;
  chomp($result);
  return ($result);
}

sub return_unique {
  my @array = @_;
  my %count;
  map { $count{$_} = 1 } @array;
  return sort keys(%count);
}

sub return_clean_text {
  my ($result) = @_;
  $result =~ s/<br \/>/ /g;
  $result =~ s/\R/ /g;
  $result =~ s/<.+?>//g;
  $result =~ s/^\s+|\s+$//g;
  $result =~ s/\s\s+/ /g;
  $result = decode_entities($result);
  if ( $result eq "\"" ) {
    $result = "";
  }
  return ($result);
}

sub download_data {
  my ($id) = @_;
  chomp($id);

  my $post = 'snp=' . $id;

  print $post . "\n";

  # Create a request
  my $req = new HTTP::Request POST =>
    'https://www.genome.gov/page.cfm?pageid=26525384#searchForm';
  $req->content_type('application/x-www-form-urlencoded');
  $req->content($post);

  # Pass request to the user agent and get a response back
  my $res = $ua->request($req);

  if ( $res->is_success ) {
    return ( 1, $res->content );
  }
  else {
    return ( 0, "" );
  }
}

open( MYFILE, 'I:/projects/guoyan/gwas/rsids.txt' );
my @array = <MYFILE>;
close(MYFILE);

open( SAVEFILE, '>I:/projects/guoyan/gwas/trait.genome.gov.txt' );
my $isfirst = 1;

my @ids = return_unique(@array);

foreach (@ids) {
  my $id = $_;
  #my $id = "rs10033900";

  chomp($id);

  my ( $succeed, $content ) = download_data($id);

  if ($succeed) {

    #print $content;

    if ($isfirst) {
      if ( $content =~ m/(<TR VALIGN="TOP" CLASS="tableheader">.+?<\/TR>)/sg ) {
        print SAVEFILE "id";
        my $header = $1;
        while ( $header =~ m/<TH\s+ID.+?>(.+?)<\/TH>/sg ) {
          print SAVEFILE "\t" . return_clean_text($1);
        }
        print SAVEFILE "\n";
        $isfirst = 0;
      }
    }

    while ( $content =~ m/(<TR bgcolor.+?<\/TR>)/sg ) {
      my $table = $1;
      print SAVEFILE $id;

      while ( $table =~ m/<TD\s+HEADERS.+?>(.*?)<\/TD>/sg ) {
        my $value = return_clean_text($1);
        print SAVEFILE "\t" . $value;
      }
      print SAVEFILE "\n";
    }

    sleep(1);

    #exit;
  }
}
close(SAVEFILE);

print "Finished!\n";

1;
