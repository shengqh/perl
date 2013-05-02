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

our %EXPORT_TAGS = ( 'all' => [qw(bwa_by_pbs_single bwa_by_pbs_double samtools_index get_sorted_bam refine_bam_file gatk_snpindel)] );

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
	die "define ${section}::option_samse first" if !defined( $config->{$section}{option_samse} );
	my $option_samse = $config->{$section}{option_samse};
	my $sort_bam;
	$sort_bam = $config->{$section}{sort_bam} or $sort_bam = 0;

	my %rawFiles = %{ get_raw_files( $config, $section ) };

	my $shfile = $pbsDir . "/${task_name}_bwa.sh";
	open( SH, ">$shfile" ) or die "Cannot create $shfile";
	print SH "type -P qsub &>/dev/null && export MYCMD=\"qsub\" || export MYCMD=\"bash\" \n";

	for my $sampleName ( sort keys %rawFiles ) {
		my @sampleFiles = @{ $rawFiles{$sampleName} };

		my $sampleFile1 = $sampleFiles[0];

		my ( $sampleName1, $directories1, $suffix1 ) = fileparse($sampleFile1);
		my $saiFile1        = $sampleName1 . ".sai";
		my $samFile         = $sampleName . ".sam";
		my $bamFile         = $sampleName . ".bam";
		my $sortedBamPrefix = $sampleName . "_sorted";
		my $sortedBamFile   = $sortedBamPrefix . ".bam";

		my $pbsName = "${sampleName}_bwa.pbs";
		my $pbsFile = "${pbsDir}/$pbsName";

		print SH "\$MYCMD ./$pbsName \n";

		my $log    = "${logDir}/${sampleName}_bwa.log";
		my $curDir = create_directory_or_die( $resultDir . "/$sampleName" );
		my $tag    = "'\@RG\tID:$sampleName\tLB:$sampleName\tSM:$sampleName\tPL:ILLUMINA'";

		open( OUT, ">$pbsFile" ) or die $!;
		print OUT "$pbsDesc
#PBS -o $log
#PBS -j oe

$path_file

echo bwa=`date`

cd $curDir

if [ ! -s $saiFile1 ]; then
  echo sai1=`date` 
  bwa aln $option $faFile $sampleFile1 > $saiFile1 
fi

echo aln=`date` 
bwa samse -r $tag $option_samse $faFile $saiFile1 $sampleFile1 > $samFile 

echo sam2bam=`date`
samtools view -b -S $samFile -o $bamFile 

";

		if ($sort_bam) {
			print OUT "
echo sortbam=`date`
samtools sort $bamFile $sortedBamPrefix 

echo indexbam=`date`
samtools index $sortedBamFile 

echo bamstat=`date`
samtools flagstat $sortedBamFile > ${sortedBamFile}.stat 

echo finished=`date` 
";
		}
		else {
			print OUT "
echo indexbam=`date`
samtools index $bamFile 

echo bamstat=`date`
samtools flagstat $bamFile > ${bamFile}.stat 

echo finished=`date` 
";
		}
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

	die "define ${section}::option_sampe first" if !defined( $config->{$section}{option_sampe} );
	my $option_sampe = $config->{$section}{option_sampe};
	my $sort_bam;
	$sort_bam = $config->{$section}{sort_bam} or $sort_bam = 0;

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
		my $saiFile2        = $sampleName2 . ".sai";
		my $samFile         = $sampleName . ".sam";
		my $bamFile         = $sampleName . ".bam";
		my $sortedBamPrefix = $sampleName . "_sort";
		my $sortedBamFile   = $sortedBamPrefix . ".bam";

		my $pbsName = "${sampleName}_bwa.pbs";
		my $pbsFile = "${pbsDir}/$pbsName";
		my $curDir  = create_directory_or_die( $resultDir . "/$sampleName" );
		my $tag     = "'\@RG\tID:$sampleName\tLB:$sampleName\tSM:$sampleName\tPL:ILLUMINA'";

		my $inserts_str = "";
		if ($inserts) {
			$inserts_str = "  echo insertsize=`date`
  samtools view $sortedBamFile | awk 'and (\$2, 0x0002) && and (\$2, 0x0040)' | cut -f 9 | sed 's/^-//' > ${sortedBamPrefix}.len 
  sort -n ${sortedBamPrefix}.len | awk ' { x[NR]=\$1; s+=\$1; } END {mean=s/NR; for (i in x){ss+=(x[i]-mean)^2}; sd=sqrt(ss/NR); if(NR %2) {median=x[(NR+1)/2];}else{median=(x[(NR/2)]+x[(NR/2)+1])/2.0;} print \"mean=\"mean \"; stdev=\"sd \"; median=\"median }' > ${sortedBamPrefix}.inserts 
  ";
		}

		print SH "\$MYCMD ./$pbsName \n";

		my $log = "${logDir}/${sampleName}_bwa.log";

		open( OUT, ">$pbsFile" ) or die $!;
		print OUT "$pbsDesc
#PBS -o $log
#PBS -j oe

$path_file

echo bwa=`date`

cd $curDir

if [ -s $sortedBamFile ]; then
  echo job has already been done. if you want to do again, delete $sortedBamFile and submit job again.
else
  if [ ! -s $saiFile1 ]; then
    echo sai1=`date` 
    bwa aln $option $faFile $sampleFile1 >$saiFile1 
  fi

  if [ ! -s $saiFile2 ]; then
    echo sai2=`date` 
    bwa aln $option $faFile $sampleFile2 >$saiFile2 
  fi

  echo aln=`date` 
  bwa sampe -r $tag $option_sampe $faFile $saiFile1 $saiFile2 $sampleFile1 $sampleFile2 > $samFile 

  echo sam2bam=`date`
  samtools view -b -S $samFile -o $bamFile 
";

		if ($sort_bam) {
			print OUT "
  echo sortbam=`date`
  samtools sort $bamFile $sortedBamPrefix 

  echo indexbam=`date`
  samtools index $sortedBamFile 

  echo bamstat=`date`
  samtools flagstat $sortedBamFile > ${sortedBamFile}.stat 
  
  $inserts_str
fi

echo finished=`date`
";
		}
		else {
			print OUT "
  echo indexbam=`date`
  samtools index $bamFile 

  echo bamstat=`date`
  samtools flagstat $bamFile > ${bamFile}.stat 
  
  $inserts_str
fi

echo finished=`date`
";
		}
		close OUT;

		print "$pbsFile created\n";
	}
	close(SH);

	if ( is_linux() ) {
		chmod 0755, $shfile;
	}

	print "!!!shell file $shfile created, you can run this shell file to submit all bwa tasks.\n";
}

sub samtools_index {
	my ( $config, $section ) = @_;

	my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option ) = get_parameter( $config, $section );

	my %rawFiles = %{ get_raw_files( $config, $section ) };

	my $isbamsorted = $config->{$section}{isbamsorted};
	if ( !defined($isbamsorted) ) {
		$isbamsorted = 0;
	}

	my $shfile = $pbsDir . "/${task_name}_index.sh";
	open( SH, ">$shfile" ) or die "Cannot create $shfile";
	print SH "type -P qsub &>/dev/null && export MYCMD=\"qsub\" || export MYCMD=\"bash\" \n";

	for my $sampleName ( sort keys %rawFiles ) {
		my @sampleFiles = @{ $rawFiles{$sampleName} };

		my $pbsName = "${sampleName}_index.pbs";
		my $pbsFile = "${pbsDir}/$pbsName";

		print SH "\$MYCMD ./$pbsName \n";

		my $log = "${logDir}/${sampleName}_index.log";

		open( OUT, ">$pbsFile" ) or die $!;
		print OUT "$pbsDesc
#PBS -o $log
#PBS -j oe

$path_file

echo index=`date`
";

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
		print OUT "if [ ! -s $bamIndexFile ]; then
  echo samtools_index=`date`
  samtools index $bamSortedFile 
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

	print "!!!shell file $shfile created, you can run this shell file to submit all samtools index tasks.\n";
}

sub refine_bam_file {
	my ( $config, $section ) = @_;

	my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option ) = get_parameter( $config, $section );

	my $faFile             = get_param_file( $config->{$section}{fasta_file}, "fasta_file", 1 );
	my @vcfFiles           = @{ $config->{$section}{vcf_files} };
	my $gatk_jar           = get_param_file( $config->{$section}{gatk_jar}, "gatk_jar", 1 );
	my $markDuplicates_jar = get_param_file( $config->{$section}{markDuplicates_jar}, "markDuplicates_jar", 1 );
	my $thread_count       = $config->{$section}{thread_count};
	if ( !defined($thread_count) ) {
		$thread_count = 1;
	}

	my %rawFiles = %{ get_raw_files( $config, $section ) };

	my $shfile = $pbsDir . "/${task_name}_refine.sh";
	open( SH, ">$shfile" ) or die "Cannot create $shfile";
	print SH "type -P qsub &>/dev/null && export MYCMD=\"qsub\" || export MYCMD=\"bash\" \n";

	for my $sampleName ( sort keys %rawFiles ) {
		my $curDir = create_directory_or_die( $resultDir . "/$sampleName" );

		my @sampleFiles = @{ $rawFiles{$sampleName} };

		my $sampleFile = $sampleFiles[0];

		my $sFile = $curDir . "/" . basename($sampleFile);
		if ( $sFile eq $sampleFile ) {
			$sFile = basename($sampleFile);
		}

		my $intervalFile  = $sFile . ".intervals";
		my $realignedFile = change_extension( $sFile, ".realigned.bam" );
		my $grpFile       = $realignedFile . ".grp";
		my $recalFile     = change_extension( $realignedFile, ".recal.bam" );
		my $rmdupFile     = change_extension( $recalFile, ".rmdup.bam" );
		my $sortedPrefix     = change_extension( $recalFile, ".sorted" );
		my $sortedFile = $sortedPrefix + ".bam";

		my $pbsName = "${sampleName}_refine.pbs";
		my $pbsFile = "${pbsDir}/$pbsName";

		print SH "\$MYCMD ./$pbsName \n";

		my $log           = "${logDir}/${sampleName}_refine.log";
		my $knownvcf      = "";
		my $knownsitesvcf = "";

		foreach my $vcf (@vcfFiles) {
			$knownvcf      = $knownvcf . " -known $vcf";
			$knownsitesvcf = $knownsitesvcf . " -knownSites $vcf";
		}

		open( OUT, ">$pbsFile" ) or die $!;

		print OUT "
$pbsDesc
#PBS -o $log
#PBS -j oe

$path_file

echo bwa=`date`
cd $curDir

if [ ! -e tmpdir ]; then
  mkdir tmpdir
fi

if [ ! -s $intervalFile ]; then
  echo RealignerTargetCreator=`date` 
  java $option -jar $gatk_jar -T RealignerTargetCreator -I $sFile -R $faFile $knownvcf -nt $thread_count -o $intervalFile
fi

if [[ -s $intervalFile && ! -s $realignedFile ]]; then
  echo IndelRealigner=`date` 
  java $option -Djava.io.tmpdir=tmpdir -jar $gatk_jar -T IndelRealigner -I $sFile -R $faFile -targetIntervals $intervalFile $knownvcf --consensusDeterminationModel KNOWNS_ONLY -LOD 0.4 -o $realignedFile 
fi

if [[ -s $realignedFile && ! -s $grpFile ]]; then
  echo BaseRecalibrator=`date` 
  java $option -jar $gatk_jar -T BaseRecalibrator -rf BadCigar -R $faFile -I $realignedFile $knownsitesvcf -o $grpFile -plots ${grpFile}.pdf
fi

if [[ -s $grpFile && ! -s $recalFile ]]; then
  echo PrintReads=`date`
  java $option -jar $gatk_jar -T PrintReads -rf BadCigar -R $faFile -I $realignedFile -BQSR $grpFile -o $recalFile 
fi

if [[ -s $recalFile && ! -s $rmdupFile ]]; then
  echo RemoveDuplicate=`date` 
  java $option -jar $markDuplicates_jar I=$recalFile O=$rmdupFile M=${rmdupFile}.matrix VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true REMOVE_DUPLICATES=true
fi

if [[ -s $rmdupFile && ! -s $sortedFile ]]; then
  echo BamSort=`date` 
  samtools sort $rmdupFile $sortedPrefix 
fi

if [[ -s $sortedFile && ! -s ${sortedFile}.bai ]]; then
  echo BamIndex=`date` 
  samtools index $sortedFile
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
}

sub gatk_snpindel {
	my ( $config, $section ) = @_;

	my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option ) = get_parameter( $config, $section );

	my $faFile = get_param_file( $config->{$section}{fasta_file}, "fasta_file", 1 );

	my @vcfFiles = @{ $config->{$section}{vcf_files} };
	my $knownvcf = "";
	foreach my $vcf (@vcfFiles) {
		if ( $knownvcf eq "" ) {
			$knownvcf = "-D $vcf";
		}
		else {
			$knownvcf = $knownvcf . " -comp $vcf";
		}
	}

	my $gatk_jar = get_param_file( $config->{$section}{gatk_jar}, "gatk_jar", 1 );
	my $gatk_option = $config->{$section}{gatk_option} or die "define ${section}::gatk_option first";

	my %rawFiles = %{ get_raw_files( $config, $section ) };

	my $shfile = $pbsDir . "/${task_name}_snpindel.sh";
	open( SH, ">$shfile" ) or die "Cannot create $shfile";
	print SH "type -P qsub &>/dev/null && export MYCMD=\"qsub\" || export MYCMD=\"bash\" \n";

	for my $sampleName ( sort keys %rawFiles ) {
		my $curDir       = create_directory_or_die( $resultDir . "/$sampleName" );
		my $listfilename = "${sampleName}.list";
		my $listfile     = $curDir . "/$listfilename";
		open( LIST, ">$listfile" ) or die "Cannot create $listfile";
		my @sampleFiles = @{ $rawFiles{$sampleName} };
		foreach my $sampleFile (@sampleFiles) {
			print LIST $sampleFile . "\n";
		}
		close(LIST);

		my $snpOut  = $sampleName . "_snp.vcf";
		my $snpStat = $sampleName . "_snp.stat";

		my $pbsName = "${sampleName}_snp.pbs";

		print SH "\$MYCMD ./$pbsName \n";

		my $log = "${logDir}/${sampleName}_snp.log";

		my $pbsFile = "${pbsDir}/$pbsName";
		open( OUT, ">$pbsFile" ) or die $!;
		print OUT "
$pbsDesc
#PBS -o $log
#PBS -j oe

$path_file

cd $curDir

echo SNP=`date` 
java -jar $option $gatk_jar -T UnifiedGenotyper -R $faFile -I $listfilename $knownvcf --out $snpOut -metrics $snpStat -glm SNP $gatk_option

echo finished=`date`
";
		close OUT;
		print "$pbsFile created\n";

		my $indelOut  = $sampleName . "_indel.vcf";
		my $indelStat = $sampleName . "_indel.stat";

		$pbsName = "${sampleName}_indel.pbs";

		print SH "\$MYCMD ./$pbsName \n";

		$log = "${logDir}/${sampleName}_indel.log";

		$pbsFile = "${pbsDir}/$pbsName";
		open( OUT, ">$pbsFile" ) or die $!;
		print OUT "
$pbsDesc
#PBS -o $log
#PBS -j oe

$path_file

cd $curDir

echo InDel=`date` 
java -jar $option $gatk_jar -T UnifiedGenotyper -R $faFile -I $listfilename $knownvcf --out $indelOut -metrics $indelStat -glm INDEL $gatk_option

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
}

1;
