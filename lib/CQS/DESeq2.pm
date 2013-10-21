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
setwd(\"$resultDir\")  
  
data<-read.table(\"$countfile\",row.names=1, header=T, check.names=F)

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
library(\"DESeq2\")
library(\"heatmap3\")
library(\"lattice\")
library(\"reshape\")
library(\"ggplot2\")

hmcols <- colorRampPalette(c(\"green\", \"black\", \"red\"))(256)

hasname <- (! is.numeric(data[1,1]))
if(hasname){
  countData<-data[,c(2:ncol(data))]
}else{
  countData<-data
}

pairnames=names(pairs)
pairname=pairnames[1]

for(pairname in pairnames){
  str(pairname)
  gs=pairs[[pairname]]
  gnames=names(gs)
  g1name=gnames[1]
  g2name=gnames[2]
  g1=gs[[g1name]]
  g2=gs[[g2name]]
  c1=countData[,colnames(countData) %in% g1,drop=F]
  c2=countData[,colnames(countData) %in% g2,drop=F]
  
  if(ncol(c1) != length(g1)){
    warning(paste0(\"There are only \", ncol(c1), \" samples in group \", g1name, \" but \", length(g1), \" required!\"))
    next
  }
  
  if(ncol(c2) != length(g2)){
    warning(paste0(\"There are only \", ncol(c2), \" samples in group \", g2name, \" but \", length(g2), \" required!\"))
    next
  }
  
  pairCountData=cbind(c1, c2)
  
  pairColData=data.frame(condition=factor(c(rep(g1name, ncol(c1)), rep(g2name, ncol(c2))), levels=gnames))
  rownames(pairColData)<-colnames(pairCountData)
  pairColors<-as.matrix(data.frame(Group=c(rep(\"red\", ncol(c1)), rep(\"blue\", ncol(c2)))))
  
  #different expression analysis
  dds=DESeqDataSetFromMatrix(countData = pairCountData,
                             colData = pairColData,
                             design = ~ condition)
  
  dds <- DESeq(dds)
  res<-results(dds,cooksCutoff=FALSE)
  
  select<- (!is.na(res\$padj)) & (res\$padj<0.05) & ((res\$log2FoldChange >= 1) | (res\$log2FoldChange <= -1))
  
  if(hasname){
    tbb<-cbind(data[,1,drop=F], pairCountData, res)
  }else{
    tbb<-cbind(pairCountData, res)
  }
  tbbselect<-tbb[select,,drop=F]
  
  tbb<-tbb[order(tbb\$padj),,drop=F]
  write.csv(as.data.frame(tbb),paste0(pairname, \"_DESeq2.csv\"))
  
  tbbselect<-tbbselect[order(tbbselect\$padj),,drop=F]
  write.csv(as.data.frame(tbbselect),paste0(pairname, \"_DESeq2_sig.csv\"))

  #transform the count data
  rld<-rlogTransformation(dds, blind=TRUE)
  rldmatrix=as.matrix(assay(rld))
  
  #draw density graph
  rsdata<-melt(rldmatrix)
  colnames(rsdata)<-c(\"Gene\", \"Sample\", \"logCount\")
  png(filename=paste0(pairname, \".logdensity.png\"), width=4000, height=3000, res=300)
  g<-ggplot(rsdata, aes(x=logCount, colour=Sample)) + geom_density() + xlab(\"log transformed count\")
  print(g)
  dev.off()
    
  #draw heatmap
  rldselect<-rldmatrix[select,,drop=F]
  if(nrow(rldselect) > 2){
    png(filename=paste0(pairname, \".heatmap.png\"), width=3000, height =3000, res=300)
    heatmap3(rldselect, col = hmcols, ColSideColors = pairColors, margins=c(12,5), scale=\"r\", dist=dist, labRow=\"\",
             legendfun=function() showLegend(legend=paste0(\"Group \", gnames),col=c(\"red\",\"blue\"),cex=1.5,x=\"center\"))
    dev.off()
  }
}
";

  my $shfile = $pbsDir . "/${task_name}_de2.sh";
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  print SH "cd $resultDir
R --vanilla -f $rfile
";
  close(SH);

  print "!!!shell file $shfile created, you can run this shell file to do calculation.\n";
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
