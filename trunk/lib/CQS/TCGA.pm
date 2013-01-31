#!/usr/bin/perl
package CQS::TCGA;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;

require Exporter;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(tcga_download tcga_get_coordinate)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

use Cwd;

sub output_header {
	my ( $pbsFile, $pbsDesc, $path_file, $log ) = @_;
	open( OUT, ">$pbsFile" ) or die $!;
	print OUT $pbsDesc;
	print OUT "#PBS -o $log\n";
	print OUT "#PBS -j oe\n\n";
	if ( defined $path_file ) {
		print OUT "source $path_file\n";
	}
}

sub output_footer() {
	print OUT "echo finished=`date`\n";
	close OUT;
}

my $bamfilter = sub {
	my $filename = shift;

	return ( $filename =~ /.bam$/ );
};

sub tcga_download {
	my ( $config, $section ) = @_;

	my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option ) = get_parameter( $config, $section );

	my $idfile = get_param_file( $config->{$section}{idfile}, "analysis id file", 1 );
	my $analysisidindex = $config->{$section}{analysisidindex} or die "Define analysisidindex at section $section";

	open( DAT, $idfile ) || die("Could not open file $idfile!");
	my $line     = <DAT>;
	my @raw_data = <DAT>;
	close(DAT);

	my $pbsFile = $pbsDir . "/${task_name}_download.pbs";
	my $log     = $logDir . "/${task_name}_download.log";

	output_header( $pbsFile, $pbsDesc, $path_file, $log );

	my $rawdir        = create_directory_or_die( $resultDir . "/raw" );

	print OUT "echo download=`date` \n";
    print OUT "cd $rawdir \n";

	foreach $line (@raw_data) {
		chomp($line);
		my @parts      = split( '\t', $line );
		my $partSize   = @parts;
        my $tcga       = $parts[$tcgaidindex];
		my $analysisid = $parts[$analysisidindex];

		my $url = "https://cghub.ucsc.edu/cghub/data/analysis/download/" . $analysisid;
        print OUT "echo $tcga `date` \n";
		print OUT "GeneTorrent -v -c ~/.ssh/mykey.pem -C ~/pylibs/share/GeneTorrent -d $url \n";
	}
	output_footer();
	print "$pbsFile created\n";
}

sub tcga_get_coordinate {
    my ( $config, $section ) = @_;

    my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option ) = get_parameter( $config, $section );

    my $idfile = get_param_file( $config->{$section}{idfile}, "analysis id file", 1 );
    my $tcgaidindex     = $config->{$section}{tcgaidindex}     or die "Define tcgaidindex at section $section";
    my $analysisidindex = $config->{$section}{analysisidindex} or die "Define analysisidindex at section $section";
    my $coordinateindex = $config->{$section}{coordinateindex} or die "Define coordinateindex at section $section";

    open( DAT, $idfile ) || die("Could not open file $idfile!");
    my $line     = <DAT>;
    my @raw_data = <DAT>;
    close(DAT);

    my $pbsFile = $pbsDir . "/${task_name}_coordidate.pbs";
    my $log     = $logDir . "/${task_name}_coordidate.log";

    output_header( $pbsFile, $pbsDesc, $path_file, $log );

    my $rawdir        = create_directory_or_die( $resultDir . "/raw" );
    my $coordinatedir = create_directory_or_die( $resultDir . "/coordindates" );

    print OUT "echo download=`date` \n";

    foreach $line (@raw_data) {
        chomp($line);
        my @parts      = split( '\t', $line );
        my $partSize   = @parts;
        my $tcga       = $parts[$tcgaidindex];
        my $analysisid = $parts[$analysisidindex];
        my $coordinate = $parts[$coordinateindex];
        
        if(!defined($coordinate) ){
        	next;
        }
        
        if(length($coordinate) == 0){
            next;
        }
        
        print $coordinate . "\n";

        my $subdir = $rawdir . '/' . $analysisid;

        my $targetfile = $coordinatedir . '/' . $tcga . ".coordinates.sam";

        my @bamfiles     = list_files( $subdir, $bamfilter );
        my $bamfile      = $subdir . "/" . $bamfiles[0];
        my $bamfileindex = $bamfile . ".bai";

        if ( !-e $bamfileindex ) {
            my $cmd = "samtools index " . $bamfile . " ";
            print OUT $cmd . "\n";
        }

        my $cmd = "samtools view " . $bamfile . " " . $coordinate . " > " . $targetfile . " ";
        print OUT $cmd . "\n";
    }
    output_footer();
    print "$pbsFile created\n";
}

1;
