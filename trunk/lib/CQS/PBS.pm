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
	my $hour   = "48";
	my $email  = "quanhu.sheng\@vanderbilt.edu";
	my $memory = "15000mb";
	my $nodes  = "8";

	my ($refHash) = @_;
	if ( defined $refHash ) {
		my %hash = %{$refHash};
		foreach ( keys %hash ) {
			if ( $_ eq "hour" ) {
				$hour = $hash{$_};
			}
			elsif ( $_ eq "email" ) {
				$email = $hash{$_};
			}
			elsif ( $_ eq "memory" ) {
				$memory = $hash{$_};
			}
			elsif ( $_ eq "nodes" ) {
				$nodes = $hash{$_};
			}
		}
	}
	my $pbsDesc = <<PBS;
#!/bin/bash
#Beginning of PBS bash script
#PBS -M $email
#Status/Progress Emails to be sent
#PBS -m bae
#Email generated at b)eginning, a)bort, and e)nd of jobs
#PBS -l nodes=$nodes
#Processors needed
#PBS -l mem=$memory
#Total job memory required (specify how many megabytes)
#PBS -l walltime=${hour}:00:00
#You must specify Wall Clock time (hh:mm:ss) [Maximum allowed 30 days = 720:00:00]
#PBS -q batch
PBS

	return ($pbsDesc);
}

1;

