#!/usr/bin/perl
package CQS::Config;

use strict;
use warnings;
use File::Basename;
use CQS::FileUtils;

require Exporter;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(get_parameter get_param_file)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

sub get_parameter {
    my ( $config, $section ) = @_;

    my $path_file = get_param_file( $config->{general}{path_file}, "path_file", 0 );
    my $refPbs     = $config->{$section}{pbs}        or die "define ${section}::pbs parameters first";
    my $target_dir = $config->{$section}{target_dir} or die "define ${section}::target_dir parameters first";
    my ( $logDir, $pbsDir, $resultDir ) = init_dir($target_dir);
    my ($pbsDesc) = get_pbs_desc($refPbs);
    my $option = $config->{$section}{option} or die "define ${section}::option first";

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

1;


