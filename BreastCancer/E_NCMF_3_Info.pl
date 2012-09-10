use strict;
use warnings;

require "SampleGeneInfo.pl";

my $uncompressedfile="D:/sqh/projects/BreastCancer/Dataset/E-NCMF-3/E-NCMF-3-processed-data-2516405675.txt";
my ($samplecount, $genecount, $genesymbol) = GetSampleGeneInfo($uncompressedfile, "Scan\\sREF", "Reporter\\sREF");
print "\n$samplecount\t$genecount\t$genesymbol";

