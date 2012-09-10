#!/usr/local/bin/perl

use strict;

use LWP::Simple;
use HTML::Entities;

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

sub guess_trait_by_odds_and_title {
  my ( $trait, $title, $odds ) = @_;
  if ( $trait eq "" ) {
    if ( ( $odds =~ /meters/ ) or ( $odds =~ /\scm\s/ ) ) {
      $trait = "Height";
    }
    elsif ( $odds =~ /kg\/m2/ ) {
      $trait = "BMI";
    }
    elsif ( $odds =~ /\skg\s/ ) {
      $trait = "Weight";
    }
    elsif ( $odds =~ /\smmHg\s/ ) {
      $trait = "Blood Pressure";
    }
  }

  if ( $trait eq "" ) {
    if ( ( $title =~ /body mass index/ )
      or ( $title =~ /length\/height/ ) )
    {
      $trait = "BMI";
    }
    elsif ( $title =~ /\sheight\s/ ) {
      $trait = "P-Height";
    }
  }
  return ($trait);
}

open( MYFILE, 'I:/projects/guoyan/gwas/rsids.txt' );
my @array = <MYFILE>;
close(MYFILE);

open( SAVEFILE, '>I:/projects/guoyan/gwas/trait.txt' );
print SAVEFILE "ID\tPMID\tTrait\tTitle\tRisk Allele\tPval\tOdds Ratio\n";

my @ids = return_unique(@array);

foreach (@ids) {
  my $id = $_;
  chomp($id);
  print $id . "\n";
  my $url     = "http://www.snpedia.com/index.php/" . $id;
  my $content = get($url);
  if ( defined($content) ) {
    while ( $content =~ m/(<table style.+?<\/table>)/sg ) {
      my $table = $1;
      if ( $table =~
/PMID\s(\d+).+?Trait.+?<td>(.+?)<\/td>.+?Title.+?<td>(.+?)<\/td>.+?Risk Allele.+?<td>(.+?)<\/td>.+?P-val.+?<td>(.+?)<\/td>.+?Odds Ratio.+?<td>(.+?)<\/td>/s
        )
      {
        my $pmid  = mychomp($1);
        my $trait = mychomp($2);
        my $title = mychomp($3);
        my $risk  = mychomp($4);
        my $pval  = mychomp($5);
        my $odds  = decode_entities( mychomp($6) );
        
        if ( $trait =~ /<a href=.+>(.+)<\/a>$/sg ) {
          $trait = $1;
        }

        $trait = guess_trait_by_odds_and_title( $trait, $title, $odds );

        print SAVEFILE $id . "\t" 
          . $pmid . "\t" 
          . $trait . "\t" 
          . $title . "\t"
          . $risk . "\t"
          . $pval . "\t"
          . $odds . "\n";
      }
    }
    sleep(1);
  }
}
close(SAVEFILE);

print "Finished!\n";

1;
