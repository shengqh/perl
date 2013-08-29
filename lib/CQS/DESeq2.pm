#!/usr/bin/perl
package CQS::DESeq2;

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
  $self->{_name} = "DESeq2";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $pairs = get_raw_files( $config, $section );

  my $groups = get_raw_files( $config, $section, "groups" );

  my $countfile = parse_param_file( $config, $section, "countfile", 1 );

  my %tpgroups = ();
  for my $groupName ( sort keys %{$groups} ) {
    my @samples = @{ $groups->{$groupName} };
    $tpgroups{$groupName} = "\"" . join( "\",\"", @samples ) . "\"";
  }

  my $rfile = $resultDir . "/${task_name}.r";
  open( RF, ">$rfile" ) or die "Cannot create $rfile";
  print RF "
library(\"DESeq2\")
library(\"heatmap.plus\")
library(\"gplots\")

hmcols <- colorRampPalette(c(\"green\", \"black\", \"red\"))(256)

setwd(\"$resultDir\")  
  
data<-read.table(\"$countfile\",row.names=1, header=T, check.names=F)

hasname <- (! is.numeric(data[1,1]))
if(hasname){
  countData<-data[,c(2:ncol(data))]
}else{
  countData<-data
}

pairs=list(
";
  my $first = 1;
  for my $pairName ( sort keys %{$pairs} ) {
    my @groupNames = @{ $pairs->{$pairName} };

    if ( scalar(@groupNames) != 2 ) {
      die "Comparison in pair $pairName should contains and only contains two groups!";
    }

    my $g1 = $groupNames[0];
    my $g2 = $groupNames[1];
    my $s1 = $tpgroups{$g1};
    my $s2 = $tpgroups{$g2};
    if ( !$first ) {
      print RF ",";
    }
    print RF "\"$pairName\" = list(\"$g1\" = c($s1), \"$g2\" = c($s2)) \n";
    $first = 0;
  }

  print RF ")
pairnames=names(pairs)

for(pairname in pairnames){
  str(pairname)
  gs=pairs[[pairname]]
  gnames=names(gs)
  g1name=gnames[1]
  g2name=gnames[2]
  g1=gs[[g1name]]
  g2=gs[[g2name]]
  c1=countData[,colnames(countData) %in% g1]
  c2=countData[,colnames(countData) %in% g2]
  
  if(ncol(c1) != length(g1)){
    warning(paste0(\"There are only \", ncol(c1), \" samples in group \", g1name, \" but \", length(g1), \" required!\"))
    next
  }
  
  if(ncol(c2) != length(g2)){
    warning(paste0(\"There are only \", ncol(c2), \" samples in group \", g2name, \" but \", length(g2), \" required!\"))
    next
  }
  
  pairCountData=cbind(c1, c2)
  
  pairColData=data.frame(condition=factor(c(rep(g1name, ncol(c1)), rep(g2name, ncol(c2)))))
  pairColors<-c(rep(\"RED\", ncol(c1)), rep(\"BLUE\", ncol(c2)))
  
  dds=DESeqDataSetFromMatrix(countData = pairCountData,
                             colData = pairColData,
                             design = ~ condition)
  
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersions(dds)
  dds <- nbinomWaldTest(dds)
  
  res<-results(dds)
  
  if(hasname){
    res\$name<-data[,1]
    res<-res[,c(ncol(res), 1:(ncol(res)-1))]
  }
  
  tbb<-res[order(res\$padj),]
  write.csv(as.data.frame(tbb),paste0(pairname, \"_DESeq2.csv\"))
  
  select<- (!is.na(res\$padj)) & (res\$padj<0.05) & ((res\$log2FoldChange >= 1) | (res\$log2FoldChange <= -1))
  
  vsd<-varianceStabilizingTransformation(dds,blind=TRUE)
  vsdmatrix<-as.matrix(assay(vsd))
  vsdselect<-vsdmatrix[select,]
  colnames(vsdselect)<-colnames(pairCountData)
  
  png(filename=paste0(pairname, \".png\"), width=4000, height =3000, res=300)
  
  clab<-matrix(c(rep(\"white\", ncol(vsdselect)), pairColors), ncol=2, byrow=FALSE)
  colnames(clab)<-c(\"\", \"Group\")
  par(mar=c(12, 10, 10, 10))
  heatmap.plus(vsdselect, col = hmcols, ColSideColors = clab, margins=c(10,15))
  
  grid.text(g1name, x = unit(0.02, \"npc\"), y = unit(0.90, \"npc\"), just = \"left\", gp=gpar(fontsize=20, col=\"RED\"))
  grid.text(g2name, x = unit(0.02, \"npc\"), y = unit(0.85, \"npc\"), just = \"left\", gp=gpar(fontsize=20, col=\"BLUE\"))
  
  dev.off()
}
";
  print "!!!R file $rfile created, you can run this R file to do calculation.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $pairs = get_raw_files( $config, $section );

  my $result = {};
  for my $pairName ( sort keys %{$pairs} ) {
    my @resultFiles = ();
    push( @resultFiles, $resultDir . "/${pairName}.csv" );
    push( @resultFiles, $resultDir . "/${pairName}.png" );
    $result->{$pairName} = filter_array( \@resultFiles, $pattern );
  }
  return $result;
}

1;
