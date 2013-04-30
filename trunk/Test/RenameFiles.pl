#!/usr/bin/perl
use strict;
use warnings;
use CQS::FileUtils;

my $dir = "D:/The Big Bang Theory/mp4/MP4/";
my @files = list_files($dir);
for my $file (@files) {
	print $file . "\n";
	my $newfile = $file;
#	if ( $newfile =~ /(S\d+E\d+)/ ) {
	if ( $newfile =~ /([12]\d)/ ) {
		$newfile = "TheBigBangTheoryS05E" . $1 . ".mp4";
		print $newfile . "\n";
		if(!($file eq $newfile)){
			print $file . "->" . $newfile . "\n";
			if (!rename ($dir . $file , $dir . $newfile)){
				print $!, "\n";
			} 
		}
	}
}
