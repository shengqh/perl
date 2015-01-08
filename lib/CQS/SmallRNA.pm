#!/usr/bin/perl
package CQS::SmallRNA;

use strict;
use warnings;
use Pipeline::SmallRNA;

require Exporter;
our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(performSmallRNAHuman performSmallRNAMouse performSmallRNAHumanTask performSmallRNAMouseTask performSmallRNARat performSmallRNARatTask)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

my $human_genome = {

	#genome database
	mirbase_count_option => "-p hsa",
	mirna_coordinate     => "/data/cqs/shengq1/reference/miRBase20/hsa.gff3",
	trna_coordinate      => "/data/cqs/guoy1/reference/smallrna/hg19_tRNA_ucsc_ensembl.bed",
	trna_fasta           => "/data/cqs/guoy1/reference/smallrna/hg19_tRNA_ucsc_ensembl.bed.fa",
	smallrna_coordinate  => "/data/cqs/guoy1/reference/smallrna/hg19_smallRNA_ucsc_ensembl.bed",
	bowtie1_index        => "/data/cqs/guoy1/reference/hg19/bowtie_index_hg19_rCRS_1.0.0/hg19_rCRS",
};

my $mouse_genome = {

	#genome database
	mirbase_count_option => "-p mmu",
	mirna_coordinate     => "/data/cqs/shengq1/reference/miRBase20/mmu.gff3",
	trna_coordinate     => "/data/cqs/guoy1/reference/smallrna/mm10_tRNA_ucsc_ensembl.bed",
	trna_fasta          => "/data/cqs/guoy1/reference/smallrna/mm10_tRNA_ucsc_ensembl.bed.fa",
	smallrna_coordinate => "/data/cqs/guoy1/reference/smallrna/mm10_smallRNA_ucsc_ensembl.bed",
	bowtie1_index       => "/data/cqs/shengq1/reference/mm10/bowtie_index/mm10",
};

my $rat_genome = {

	#genome database
	mirbase_count_option => "-p rno",
	mirna_coordinate     => "/data/cqs/shengq1/reference/miRBase20/rno.gff3",
	trna_coordinate      => "/data/cqs/guoy1/reference/smallrna/rn4_tRNA_ucsc_ensembl.bed",
	trna_fasta           => "/data/cqs/guoy1/reference/smallrna/rn4_tRNA_ucsc_ensembl.bed.fa",
	smallrna_coordinate  => "/data/cqs/guoy1/reference/smallrna/rn4_smallRNA_ucsc_ensembl.bed",
	bowtie1_index        => "/data/cqs/shengq1/reference/rn4/bowtie1_index/rn4",
};

sub getDefinition {
	my ( $userdef, $genome ) = @_;

	my $cluster="torque";
	
	if(defined $userdef->{cluster}){
		$cluster = $userdef->{cluster};
	}
	
	my $def = {

		#General options
		task_name  => $userdef->{task_name},
		email      => $userdef->{email},
		target_dir => $userdef->{target_dir},
		max_thread => $userdef->{max_thread},
		cluster    => $cluster,
		
		#Data
		files => $userdef->{files},

		#Default software parameter (don't change it except you really know it)
		bowtie1_option_1mm         => "-a -m 100 --best --strata -v 1 -p 8",
		bowtie1_option_pm          => "-a -m 100 --best --strata -v 0 -p 8",
		mirnacount_option          => "-s",                                         #ignore score
		smallrnacount_option       => "-s --min_overlap 0.5 --length --sequence",
		mirna_overlap_count_option => "-s --min_overlap 0.5 --gtf_key miRNA",
		min_read_length            => 16,

		#Software and miRBase database options
		samtools              => "/home/shengq1/local/bin/samtools/samtools",
		cqstools              => "/home/shengq1/cqstools/CQS.Tools.exe",
		mirna_fasta           => "/data/cqs/shengq1/reference/miRBase20/mature.dna.fa",
		bowtie1_miRBase_index => "/data/cqs/shengq1/reference/miRBase21/bowtie_index_1.0.1/mature.dna",

		#genome database
		mirbase_count_option => $genome->{mirbase_count_option},
		mirna_coordinate     => $genome->{mirna_coordinate},
		trna_coordinate      => $genome->{trna_coordinate},
		trna_fasta           => $genome->{trna_fasta},
		smallrna_coordinate  => $genome->{smallrna_coordinate},
		bowtie1_index        => $genome->{bowtie1_index},
	};
	
	return $def;
}

sub performSmallRNAHuman {
	my ($userdef) = @_;
	my $def = getDefinition($userdef, $human_genome);

	performSmallRNA($def);
}

sub performSmallRNAHumanTask {
	my ($userdef, $task) = @_;
	my $def = getDefinition($userdef, $human_genome);

	performSmallRNATask($def, $task);
}

sub performSmallRNAMouse {
	my ($userdef) = @_;
	my $def = getDefinition($userdef, $mouse_genome);

	performSmallRNA($def);
}

sub performSmallRNAMouseTask {
	my ($userdef, $task) = @_;
	my $def = getDefinition($userdef, $mouse_genome);

	performSmallRNATask($def, $task);
}

sub performSmallRNARat {
	my ($userdef) = @_;
	my $def = getDefinition($userdef, $rat_genome);

	performSmallRNA($def);
}

sub performSmallRNARatTask {
	my ($userdef, $task) = @_;
	my $def = getDefinition($userdef, $rat_genome);

	performSmallRNATask($def, $task);
}

1;