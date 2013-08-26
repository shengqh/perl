#!/usr/bin/perl
use strict;
use warnings;
use threads;
use threads::shared;

my $current : shared;

sub addCurrent {
  my ($num) = @_;
  while (1) {
    {
      lock($current);
      $current++;
      print $num, "\t", $current, "\n";
      if ( $current > 10 ) {
        return;
      }
    }
    sleep(1);
  }
}

my $t1 = threads->create( 'addCurrent', '1' );
my $t2 = threads->create( 'addCurrent', "2" );
my $t3 = threads->create( 'addCurrent', "3" );

$t1->join();
$t2->join();
$t3->join();

