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
library(\"DESeq2\")
library(\"heatmap.plus\")

hmcols <- colorRampPalette(c(\"green\", \"black\", \"red\"))(256)

setwd(\"$resultDir\")  
  
data<-read.table(\"$countfile\",row.names=1, header=T, check.names=F)

hasname <- (! is.numeric(data[1,1]))
if(hasname){
  countData<-data[,c(2:ncol(data))]
}else{
  countData<-data
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

condition=factor(unlist(groups[colnames(countData)]))

colData=data.frame(condition=condition)
";

  for my $pairName ( sort keys %{$pairs} ) {
    my @groupNames = @{ $pairs->{$pairName} };
    
    if(scalar(@groupNames) != 2){
      die "Comparison in pair $pairName should contains and only contains two groups!";
    }
    
    my $g1 = $groupNames[0];
    my $g2 = $groupNames[1];
    
    print RF "
#$pairName
    
pairCountData=cbind(countData[,colData\$condition==\"$g1\"], 
                    countData[,colData\$condition==\"$g2\"])
                    
pairColData=data.frame( condition=factor(unlist(groups[colnames(pairCountData)]), 
                        levels=c(\"$g2\", \"$g1\")))

pairColorDef=list(\"$g1\"=\"RED\", \"$g2\"=\"BLUE\")
pairColors<-unlist(pairColorDef[pairColData\$condition])

dds=DESeqDataSetFromMatrix(countData = pairCountData,
                           colData = pairColData,
                           design = ~ condition)

dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds, cooksCutoff=FALSE)

res<-results(dds)

if(hasname){
  res\$name<-data[,1]
  res<-res[,c(ncol(res), 1:(ncol(res)-1))]
}

tbb<-res[order(res\$padj),]

write.csv(as.data.frame(tbb),\"${pairName}.csv\"))

select<- (!is.na(res\$padj)) & (res\$padj<0.05)

vsd<-varianceStabilizingTransformation(dds,blind=TRUE)
vsdmatrix<-as.matrix(assay(vsd))
vsdselect<-vsdmatrix[select,]

png(filename=\"${pairName}.png\", width=4000, height=3000, res=300)

clab<-matrix(c(rep(\"white\", ncol(vsdselect)), pairColors), ncol=2, byrow=FALSE)
colnames(clab)<-c(\"\", \"Group\")
par(mar=c(12, 10, 10, 10))
heatmap.plus(vsdselect, col = hmcols, ColSideColors = clab, margins=c(10,15))

dev.off()
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
