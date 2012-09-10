use strict;
use warnings;

require "SampleGeneInfo.pl";

my $uncompressedfile="D:/sqh/projects/BreastCancer/Dataset/MDA133/MDA133PredictorTrainAndValidation.txt";
my ($samplecount, $genecount, $genesymbol) = GetSampleGeneInfo($uncompressedfile, "ID_REF");
print "\n$samplecount\t$genecount\t$genesymbol";

