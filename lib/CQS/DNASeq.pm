#!/usr/bin/perl
package CQS::DNASeq;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;

require Exporter;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(bwa_by_pbs_single bwa_by_pbs_double samtools_index get_sorted_bam refine_bam_file)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

use Cwd;

sub get_sorted_bam {
  my $bamFile = shift;
  my ( $name, $path, $suffix ) = fileparse( $bamFile, qr/\Q.bam\E/ );
  my $bamSorted     = $path . $name . "_sorted";
  my $bamSortedFile = $bamSorted . ".bam";
  return ( $bamSortedFile, $bamSorted );
}

sub bwa_by_pbs_single {
	my ( $config, $section ) = @_;

	my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option ) = get_parameter( $config, $section );

	my $faFile = get_param_file( $config->{$section}{fasta_file}, "fasta_file", 1 );

	die "define ${section}::option_samse first" if ( !defined $config->{$section}{option_samse} );
	my $option_samse = $config->{$section}{option_samse};

	my %rawFiles = %{ get_raw_files( $config, $section ) };

	my $shfile = $pbsDir . "/${task_name}.sh";
	open( SH, ">$shfile" ) or die "Cannot create $shfile";
	print SH "type -P qsub &>/dev/null && export MYCMD=\"qsub\" || export MYCMD=\"bash\" \n";

	for my $sampleName ( sort keys %rawFiles ) {
		my @sampleFiles = @{ $rawFiles{$sampleName} };

		my $sampleFile1 = $sampleFiles[0];

		my ( $sampleName1, $directories1, $suffix1 ) = fileparse($sampleFile1);
		my $saiFile1      = $sampleName1 . ".sai";
		my $samFile       = $sampleName . ".sam";
		my $bamFile       = $sampleName . ".bam";
    my $sortedBamPrefix = $sampleName . "_sorted";
    my $sortedBamFile = $sortedBamPrefix . ".bam";

		my $pbsName = "${sampleName}_bwa.pbs";
		my $pbsFile = "${pbsDir}/$pbsName";

		print SH "\$MYCMD ./$pbsName \n";

		my $log = "${logDir}/${sampleName}_bwa.log";

		open( OUT, ">$pbsFile" ) or die $!;
		print OUT $pbsDesc;
		print OUT "#PBS -o $log\n";
		print OUT "#PBS -j oe\n\n";

		if ( -e $path_file ) {
			print OUT "source $path_file\n";
		}
		print OUT "echo bwa=`date`\n";

		my $curDir = create_directory_or_die( $resultDir . "/$sampleName" );

		#my $tag="'\@RG\tID:$sample\tLB:$sample\tSM:$sample\tPL:ILLUMINA'";
		print OUT "cd $curDir\n\n";

    print OUT "if [ ! -e $saiFile1 ]; then\n";
		print OUT "  echo sai1=`date` \n";
		print OUT "  bwa aln $option $faFile $sampleFile1 > $saiFile1 \n";
		print OUT "fi\n\n";
		
		print OUT "echo aln=`date` \n";
		print OUT "bwa samse -r '\@RG\tID:${sampleName}\tLB:${sampleName}\tSM:${sampleName}\tPL:ILLUMINA' $option_samse $faFile $saiFile1 $sampleFile1 > $samFile \n\n";
		
		print OUT "echo sam2bam=`date`\n";
		print OUT "samtools view -b -S $samFile -o $bamFile \n\n";
		
		print OUT "echo sortbam=`date`\n";
		print OUT "samtools sort $bamFile $sortedBamPrefix \n\n";
    
    print OUT "echo indexbam=`date`\n";
    print OUT "samtools index $sortedBamFile \n\n";
		
		print OUT "echo bamstat=`date`\n";
		print OUT "samtools flagstat $sortedBamFile > ${sortedBamFile}.stat \n\n";
		
		print OUT "echo finished=`date` \n";
		close OUT;

		print "$pbsFile created \n";
	}
	close(SH);

	if ( is_linux() ) {
		chmod 0755, $shfile;
	}

	print "!!!shell file $shfile created, you can run this shell file to submit all bwa tasks.\n";

	#`qsub $pbsFile`;
}

sub bwa_by_pbs_double {
	my ( $config, $section ) = @_;

	my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option ) = get_parameter( $config, $section );

	my $faFile = get_param_file( $config->{$section}{fasta_file}, "fasta_file", 1 );
	my $inserts = $config->{$section}{estimate_insert};

	die "define ${section}::option_sampe first" if ( !defined $config->{$section}{option_sampe} );
	my $option_sampe = $config->{$section}{option_sampe};

	my %rawFiles = %{ get_raw_files( $config, $section ) };

	my $shfile = $pbsDir . "/${task_name}.sh";
	open( SH, ">$shfile" ) or die "Cannot create $shfile";
	print SH "type -P qsub &>/dev/null && export MYCMD=\"qsub\" || export MYCMD=\"bash\" \n";

	for my $sampleName ( sort keys %rawFiles ) {
		my @sampleFiles = @{ $rawFiles{$sampleName} };

		my $sampleFile1 = $sampleFiles[0];
		my $sampleFile2 = $sampleFiles[1];

		my ( $sampleName1, $directories1, $suffix1 ) = fileparse($sampleFile1);
		my $saiFile1 = $sampleName1 . ".sai";
		my ( $sampleName2, $directories2, $suffix2 ) = fileparse($sampleFile2);
		my $saiFile2      = $sampleName2 . ".sai";
		my $samFile       = $sampleName . ".sam";
		my $bamFile       = $sampleName . ".bam";
    my $sortedBamPrefix = $sampleName . "_sort";
		my $sortedBamFile = $sortedBamPrefix . ".bam";

		my $pbsName = "${sampleName}_bwa.pbs";
		my $pbsFile = "${pbsDir}/$pbsName";

		print SH "\$MYCMD ./$pbsName \n";

		my $log = "${logDir}/${sampleName}_bwa.log";

		open( OUT, ">$pbsFile" ) or die $!;
		print OUT $pbsDesc;
		print OUT "#PBS -o $log\n";
		print OUT "#PBS -j oe\n\n";

		if ( -e $path_file ) {
			print OUT "source $path_file\n";
		}
		print OUT "echo bwa=`date`\n";

		my $curDir = create_directory_or_die( $resultDir . "/$sampleName" );

		my $tag="'\@RG\tID:$sampleName\tLB:$sampleName\tSM:$sampleName\tPL:ILLUMINA'";
		print OUT "cd $curDir\n\n";

		print OUT "if [ -e $sortedBamFile ]; then\n";
		print OUT "  echo job has already been done. if you want to do again, delete $sortedBamFile and submit job again.\n";
		print OUT "else\n";
    print OUT "  if [ ! -e $saiFile1 ]; then\n";
		print OUT "    echo sai1=`date` \n";
		print OUT "    bwa aln $option $faFile $sampleFile1 >$saiFile1 \n";
    print OUT "  fi\n\n";
    
    print OUT "  if [ ! -e $saiFile2 ]; then\n";
		print OUT "    echo sai2=`date` \n";
		print OUT "    bwa aln $option $faFile $sampleFile2 >$saiFile2 \n";
    print OUT "  fi\n\n";
		
		print OUT "  echo aln=`date` \n";
		print OUT "  bwa sampe -r $tag $option_sampe $faFile $saiFile1 $saiFile2 $sampleFile1 $sampleFile2 > $samFile \n\n";
		
		print OUT "  echo sam2bam=`date`\n";
		print OUT "  samtools view -b -S $samFile -o $bamFile \n\n";
		
		print OUT "  echo sortbam=`date`\n";
		print OUT "  samtools sort $bamFile $sortedBamPrefix \n\n";
		
    print OUT "  echo indexbam=`date`\n";
    print OUT "  samtools index $sortedBamFile \n\n";
    
		print OUT "  echo bamstat=`date`\n";
		print OUT "  samtools flagstat $sortedBamFile > ${sortedBamFile}.stat \n\n";

		if ($inserts) {
			print OUT "  echo insertsize=`date`\n";
			print OUT "  samtools view $sortedBamFile | awk 'and (\$2, 0x0002) && and (\$2, 0x0040)' | cut -f 9 | sed 's/^-//' > ${sortedBamPrefix}.len \n";
			print OUT "  sort -n ${sortedBamPrefix}.len | awk ' { x[NR]=\$1; s+=\$1; } END {mean=s/NR; for (i in x){ss+=(x[i]-mean)^2}; sd=sqrt(ss/NR); if(NR %2) {median=x[(NR+1)/2];}else{median=(x[(NR/2)]+x[(NR/2)+1])/2.0;} print \"mean=\"mean \"; stdev=\"sd \"; median=\"median }' > ${sortedBamPrefix}.inserts \n";
		}
		print OUT "fi\n\n";

		print OUT "echo finished=`date`\n";
		close OUT;

		print "$pbsFile created\n";
	}
	close(SH);

	if ( is_linux() ) {
		chmod 0755, $shfile;
	}

	print "!!!shell file $shfile created, you can run this shell file to submit all bwa tasks.\n";

	#`qsub $pbsFile`;
}

sub samtools_index {
	my ( $config, $section ) = @_;

	my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option ) = get_parameter( $config, $section );

	my %rawFiles = %{ get_raw_files( $config, $section ) };

	my $isbamsorted = $config->{$section}{isbamsorted};
	if ( !defined($isbamsorted) ) {
		$isbamsorted = 0;
	}

	my $shfile = $pbsDir . "/${task_name}.sh";
	open( SH, ">$shfile" ) or die "Cannot create $shfile";
	print SH "type -P qsub &>/dev/null && export MYCMD=\"qsub\" || export MYCMD=\"bash\" \n";

	for my $sampleName ( sort keys %rawFiles ) {
		my @sampleFiles = @{ $rawFiles{$sampleName} };

		my $pbsName = "${sampleName}_index.pbs";
		my $pbsFile = "${pbsDir}/$pbsName";

		print SH "\$MYCMD ./$pbsName \n";

		my $log = "${logDir}/${sampleName}_index.log";

		open( OUT, ">$pbsFile" ) or die $!;
		print OUT $pbsDesc;
		print OUT "#PBS -o $log\n";
		print OUT "#PBS -j oe\n\n";

		if ( -e $path_file ) {
			print OUT "source $path_file\n";
		}
		print OUT "echo index=`date`\n";

		my $bamFile = $sampleFiles[0];

		my $bamSortedFile;
		if ($isbamsorted) {
			$bamSortedFile = $bamFile;
		}
		else {
			( $bamSortedFile, my $bamSorted ) = get_sorted_bam($bamFile);
			print OUT "if [ ! -s $bamSortedFile ]; then\n";
			print OUT "  echo samtools_sort=`date`\n";
			print OUT "  samtools sort $bamFile $bamSorted \n";
			print OUT "fi\n";
		}

		my $bamIndexFile = $bamSortedFile . ".bai";
		print OUT "if [ ! -s $bamIndexFile ]; then\n";
		print OUT "  echo samtools_index=`date`\n";
		print OUT "  samtools index $bamSortedFile \n";
		print OUT "fi\n";
		print OUT "echo finished=`date`\n";
		close OUT;

		print "$pbsFile created\n";
	}
	close(SH);

	if ( is_linux() ) {
		chmod 0755, $shfile;
	}

	print "!!!shell file $shfile created, you can run this shell file to submit all samtools index tasks.\n";

	#`qsub $pbsFile`;
}

sub refine_bam_file {
  my ( $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option ) = get_parameter( $config, $section );

  my $faFile             = get_param_file( $config->{$section}{fasta_file},         "fasta_file",         1 );
  my @vcfFiles           = $config->{$section}{vcf_files};
  my $gatk_jar           = get_param_file( $config->{$section}{gatk_jar},           "gatk_jar",           1 );
  my $markDuplicates_jar = get_param_file( $config->{$section}{markDuplicates_jar}, "markDuplicates_jar", 1 );
  my $thread_count       = $config->{$section}{thread_count};
  if ( !defined($thread_count) ) {
    $thread_count = 1;
  }

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $shfile = $pbsDir . "/${task_name}.sh";
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  print SH "type -P qsub &>/dev/null && export MYCMD=\"qsub\" || export MYCMD=\"bash\" \n";

  for my $sampleName ( sort keys %rawFiles ) {
    my $curDir = create_directory_or_die( $resultDir . "/$sampleName" );

    my @sampleFiles = @{ $rawFiles{$sampleName} };

    my $sampleFile   = $sampleFiles[0];
    
    my $sFile = $curDir . "/" . basename($sampleFile);
    if($sFile eq $sampleFile){
      $sFile = basename($sampleFile);
    }
    
    my $intervalFile  = $sFile . ".intervals";
    my $realignedFile = change_extension( $sFile, ".realigned.bam" );
    my $grpFile       = $realignedFile . ".grp";
    my $recalFile     = change_extension( $realignedFile, ".recal.bam" );
    my $rmdupFile     = change_extension( $recalFile, ".rmdup.bam" );

    my $pbsName = "${sampleName}_refine.pbs";
    my $pbsFile = "${pbsDir}/$pbsName";

    print SH "\$MYCMD ./$pbsName \n";

    my $log    = "${logDir}/${sampleName}_refine.log";

    open( OUT, ">$pbsFile" ) or die $!;

    print OUT "
$pbsDesc
#PBS -o $log
#PBS -j oe
";

    if ( -e $path_file ) {
      print OUT "source $path_file\n";
    }

  my $knownvcf = "";
  my $knownsitesvcf="";
  
  foreach my $vcf(@vcfFiles){
    $knownvcf = $knownvcf . " -known $vcf";
    $knownsitesvcf = $knownsitesvcf . " -knownSites $vcf";
  }
  
    print OUT "
echo bwa=`date`
cd $curDir

if [ ! -s tmpdir ]; then
  mkdir tmpdir
fi

if [ ! -e $intervalFile ]; then
  echo RealignerTargetCreator=`date` 
  java $option -jar $gatk_jar -T RealignerTargetCreator -I $sFile -R $faFile $knownvcf -nt $thread_count -o $intervalFile
fi

if [[ -e $intervalFile && ! -e $realignedFile ]]; then
  echo IndelRealigner=`date` 
  java $option -Djava.io.tmpdir=tmpdir -jar $gatk_jar -T IndelRealigner -I $sFile -R $faFile -targetIntervals $intervalFile $knownvcf --consensusDeterminationModel KNOWNS_ONLY -LOD 0.4 -o $realignedFile 
fi

if [[ -e $realignedFile && ! -e $grpFile ]]; then
  echo BaseRecalibrator=`date` 
  java $option -jar $gatk_jar -T BaseRecalibrator -R $faFile -I $realignedFile $knownsitesvcf -o $grpFile -plots ${grpFile}.pdf
fi

if [[ -e $grpFile && ! -e $recalFile ]]; then
  echo PrintReads=`date`
  java $option -jar $gatk_jar -T PrintReads -R $faFile -I $realignedFile -BQSR $grpFile -o $recalFile 
fi

if [[ -e $recalFile && ! -e $rmdupFile ]]; then
  echo RemoveDuplicate=`date` 
  java $option -jar $markDuplicates_jar I=$recalFile O=$rmdupFile M=${rmdupFile}.matrix VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true REMOVE_DUPLICATES=true
fi

if [[ -e $rmdupFile && ! -e ${rmdupFile}.bai ]]; then
  echo BamIndex=`date` 
  samtools index $rmdupFile
fi

echo finished=`date`
";

    close OUT;

    print "$pbsFile created\n";
  }
  close(SH);

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all bwa tasks.\n";

  #`qsub $pbsFile`;
}

1;
