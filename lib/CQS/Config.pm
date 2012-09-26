#!/usr/bin/perl
package CQS::Config;

use strict;
use warnings;
use File::Basename;
use CQS::FileUtils;
use CQS::PBS;

require Exporter;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(get_parameter get_param_file get_raw_files)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

sub get_parameter {
	my ( $config, $section ) = @_;

    my $task_name = $config->{general}{task_name} or die "define general::task_name first";
	my $path_file = get_param_file( $config->{general}{path_file}, "path_file", 0 );
	my $refPbs     = $config->{$section}{pbs}        or die "define ${section}::pbs parameters first";
	my $target_dir = $config->{$section}{target_dir} or die "define ${section}::target_dir parameters first";
	my ( $logDir, $pbsDir, $resultDir ) = init_dir($target_dir);
	my ($pbsDesc) = get_pbs_desc($refPbs);

	die "define ${section}::option first" if ( !defined $config->{$section}{option} );
	my $option = $config->{$section}{option};

	return ( $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option );
}

#get parameter which indicates a file. If required, not defined or not exists, die. If defined but not exists, die.
#returned file either undefined or exists.
sub get_param_file {
	my ( $file, $name, $required ) = @_;

	my $result = $file;

	if ($required) {
		if ( !defined $file ) {
			die "$name was not defined!";
		}

		if ( !-e $file ) {
			die "$name $file defined but not exists!";
		}
	}
	else {
		if ( defined($file) ) {
			if ( $file eq "" ) {
				undef($result);
			}
			elsif ( !-e $file ) {
				die "$name $file defined but not exists!";
			}
		}
	}
	return ($result);
}


sub get_raw_files {
    my ( $config, $section ) = @_;

    if ( defined $config->{$section}{source} ) {
        return ( $config->{$section}{source} );
    }

    if ( defined $config->{$section}{source_ref} ) {
        my $result = $config->{ $config->{$section}{source_ref} };
        if ( !defined $result ) {
            die "section ${result} was not defined!";
        }
        return ($result);
    }

    die "define ${section}::source or ${section}::source_ref first!";
}

1;

