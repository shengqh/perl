#!/usr/bin/perl  -w
############################################################################
#   File:   findtext.pl
#   Brief:  Match text content in a directory
#
#   Usage:  findtext.pl  DIR PATTERN
#           - DIR      work directory
#           - PATTERN  pattern for text search
#
#   Author: Thinkhy
#   Date:   2011-11-8
#   Changes:
#
############################################################################


use strict;
use Encode;
use File::Glob ':glob';

my $cnt = @ARGV;
if ($cnt <= 0)
{
    print " Usage:  findtext.pl  DIR PATTERN
             - DIR      work directory
             - PATTERN  pattern for text search
          ";
    exit -1;
}

my $dir = decode("gb2312", $ARGV[0]);
my $pattern0 = $ARGV[1];
my $pattern = decode("gb2312", $ARGV[1]);

print $pattern;

chomp($dir);
chomp($pattern);

$dir =~ s/\\/\//g;

print "Search content in directory: ".$dir."\n";
print "===================================================\n";
print "\n";

Find($dir);

my $deep = 0;
sub Find
{
    my $path = shift;

    $deep++;

    my @pathes = listpath($path);
    my @files = listfile($path);

    foreach(@files)
    {
        my $file = $_;
        if (isPatternOnFile($file, $pattern))
        {
            print $file."\n"
        }
    }

    if ($pathes[0])
    {
        foreach(@pathes)
        {
            Find($_);
        }
    }

    $deep--;
}

sub listpath
{
    my $path = shift;
    my @list = bsd_glob "$path/*";
    #my @list = <$path/*>;
    my @pathes = grep {  -d  } @list;
    return @pathes;
}

sub listfile
{
    my $path = shift;
    my @list = bsd_glob "$path/*";
    #my @list = <$path/*.*>;
    my @files = grep {  -f } @list;
    return @files;
}

sub isPatternOnFile
{
    my ($file, $pattern) = @_;
    my $temp;

    #read the entire file
    open my $fh, $file or die $!;
    {
        local $/;
        $temp = <$fh>;
    }
    return 1 if ($temp =~ /$pattern/i) or return 0;
}
