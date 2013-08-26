#!/usr/bin/perl
package CQS::DESeq;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::Task;
use CQS::NGSCommon;

our @ISA = qw(CQS::Task);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name} = "DESeq";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $pairs = get_raw_files( $config, $section );

  my $groups = get_raw_files( $config, $section, "groups" );

  my $countfile = parse_param_file( $config, $section, "countfile", 1 );

  my $rfile = $resultDir . "/${task_name}.r";
  open( RF, ">$rfile" ) or die "Cannot create $rfile";
  print RF "
library(\"DESeq\")

setwd(\"$resultDir\")  
  
data<-read.table(\"$countfile\",row.names=1, header=T, check.names=F)
if(is.numeric(data[1,1])){
  countTable<-data
}
else{
  countTable<-data[,c(2:ncol(data))]
}

groups=list(
";
  my $first = 1;
  for my $groupName ( sort keys %{$groups} ) {
    my @samples = @{ $groups->{$groupName} };
    for my $sample (sort @samples){
      if (! $first){
        print RF ",";
      }
      print RF "\"$sample\"=\"$groupName\"\n";
      $first = 0;
    }
  }
  print RF ")

condition=factor(unlist(groups[colnames(countTable)]))

cds = newCountDataSet(countTable, condition)
cds = estimateSizeFactors(cds)
cds = estimateDispersions(cds)

rs=rowSums(counts(cds))
theta=0.4
use=(rs > quantile(rs, probs=theta))
cdsFilt = cds[use,]

";

  for my $pairName ( sort keys %{$pairs} ) {
    my @groupNames = @{ $pairs->{$pairName} };
    
    if(scalar(@groupNames) != 2){
      die "Comparison in pair $pairName should contains and only contains two groups!";
    }
    
    print RF "tb=nbinomTest(cdsFilt, \"" . $groupNames[0], "\", \"", $groupNames[1], "\")
tbb<-tb[order(tb\$padj),]
write.csv(tbb, paste0($pairName, \".csv\"))

";
  }
   print "!!!R file $rfile created, you can run this R file to do calculation.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $pairs = get_raw_files( $config, $section, "pairs" );

  my $result = {};
  for my $pairName ( sort keys %{$pairs} ) {
    my $curDir      = $resultDir . "/$pairName";
    my @resultFiles = ();
    push( @resultFiles, $curDir . "/gene_exp.diff" );
    push( @resultFiles, $curDir . "/genes.read_group_tracking" );
    push( @resultFiles, $curDir . "/splicing.diff" );
    push( @resultFiles, $resultDir . "/${task_name}_group_sample.map" );

    $result->{$pairName} = filter_array( \@resultFiles, $pattern );
  }
  return $result;
}

1;
