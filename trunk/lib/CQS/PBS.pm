#!/usr/bin/perl
package CQS::PBS;

use strict;
use warnings;
use File::Basename;

require Exporter;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(init_dir get_pbs_desc)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

sub init_dir {
  my $rootDir = shift;

  #defined several folders
  my $pbsDir    = "$rootDir/pbs";
  my $resultDir = "$rootDir/result";
  my $logDir    = "$rootDir/log";

  if ( !( -d $rootDir ) ) {
    system "mkdir $rootDir";
  }
  if ( !( -d $pbsDir ) ) {
    system "mkdir $pbsDir";
  }
  if ( !( -d $resultDir ) ) {
    system "mkdir $resultDir";
  }
  if ( !( -d $logDir ) ) {
    system "mkdir $logDir";
  }

  return ( $logDir, $pbsDir, $resultDir );
}

sub get_pbs_desc {
	my ($hour) = @_;
	if (not defined $hour ){
		$hour = "24";
	}
  my $pbsDesc = <<PBS;
#!/bin/bash
#Beginning of PBS bash script
#PBS -M quanhu.sheng\@vanderbilt.edu
#Status/Progress Emails to be sent
#PBS -m bae
#Email generated at b)eginning, a)bort, and e)nd of jobs
#PBS -l mem=15000mb
#Total job memory required (specify how many megabytes)
#PBS -l walltime=${hour}:00:00
#You must specify Wall Clock time (hh:mm:ss) [Maximum allowed 30 days = 720:00:00]
#PBS -q batch
PBS

  return ($pbsDesc);
}

1;

