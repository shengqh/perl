use strict;
use warnings;

require "SampleGeneInfo.pl";

my $uncompressedfile="D:/sqh/projects/BreastCancer/Dataset/E-TABM-158/E-TABM-158.txt";
my ($samplecount, $genecount, $genesymbol) = GetSampleGeneInfo($uncompressedfile, "Scan\\sREF", "\\sREF");
print "\n$samplecount\t$genecount\t$genesymbol";

