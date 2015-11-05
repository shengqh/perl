#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;
use CQS::PerformSmallRNA;
use CQS::FileUtils;

my $def = {

  #General options
  task_name  => "2797_mouse",
  email      => "quanhu.sheng\@vanderbilt.edu",
  target_dir => create_directory_or_die("/scratch/cqs/shengq1/vickers/20150911_parclip_gsnap_2797-mouse/"),
  max_thread => 8,
  cluster    => "slurm",

  cqstools => "/home/shengq1/cqstools/CQS.Tools.exe",

  gsnap_index_directory => "/scratch/cqs/shengq1/references/mm10/gsnap_index_k14_2015-06-23/",
  gsnap_index_name      => "mm10",

  #Data
  identicals => {
    "RPI47" => [ "/gpfs21/scratch/cqs/shengq1/vickers/20150401_parclip_2797_rat_mouse_human/mm10/identical/result/RPI47_clipped_identical.fastq.gz", ],
    "RPI48" => [ "/gpfs21/scratch/cqs/shengq1/vickers/20150401_parclip_2797_rat_mouse_human/mm10/identical/result/RPI48_clipped_identical.fastq.gz", ],
  },
  identicals_count => {
    "RPI47" => [ "/gpfs21/scratch/cqs/shengq1/vickers/20150401_parclip_2797_rat_mouse_human/mm10/identical/result/RPI47_clipped_identical.fastq.dupcount", ],
    "RPI48" => [ "/gpfs21/scratch/cqs/shengq1/vickers/20150401_parclip_2797_rat_mouse_human/mm10/identical/result/RPI48_clipped_identical.fastq.dupcount", ],
  },
  files => {
    "RPI47" => [ "/gpfs21/scratch/cqs/shengq1/vickers/20150401_parclip_2797_rat_mouse_human/mm10/identical_NTA/result/RPI47_clipped_identical_NTA.fastq.gz", ],
    "RPI48" => [ "/gpfs21/scratch/cqs/shengq1/vickers/20150401_parclip_2797_rat_mouse_human/mm10/identical_NTA/result/RPI48_clipped_identical_NTA.fastq.gz", ],
  },
  counts => {
    "RPI47" => ["/gpfs21/scratch/cqs/shengq1/vickers/20150401_parclip_2797_rat_mouse_human/mm10/identical_NTA/result/RPI47_clipped_identical_NTA.fastq.dupcount"],
    "RPI48" => ["/gpfs21/scratch/cqs/shengq1/vickers/20150401_parclip_2797_rat_mouse_human/mm10/identical_NTA/result/RPI48_clipped_identical_NTA.fastq.dupcount"],
  },
};

my $config = {
  general => { "task_name" => $def->{task_name}, },
  gsnap   => {
    class                 => "Alignment::Gsnap",
    perform               => 0,
    target_dir            => $def->{target_dir} . "/gsnap",
    option                => "-y 0 -z 0 -Y 0 -Z 0 -m 1 -Q --trim-mismatch-score 0 --trim-indel-score 0 --mode ttoc-nonstranded --gunzip",
    gsnap_index_directory => $def->{gsnap_index_directory},
    gsnap_index_name      => $def->{gsnap_index_name},
    source                => $def->{files},
    sh_direct             => 1,
    pbs                   => {
      "email"    => $def->{email},
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    }
  },
  gsnap_smallRNA_count => {
    class           => 'CQS::SmallRNACount',
    perform         => 0,
    target_dir      => $def->{target_dir} . "/gsnap_smallRNA_count",
    option          => '-s -e 4',
    cluster         => 'slurm',
    sh_direct       => 1,
    fastq_files     => $def->{files},
    coordinate_file => "/scratch/cqs/shengq1/references/smallrna/mm10_miRBase21_ucsc-tRNA_ensembl78.bed",
    fasta_file      => "/scratch/cqs/shengq1/references/smallrna/mm10_miRBase21_ucsc-tRNA_ensembl78.bed.fa",
    source_ref      => 'gsnap',
    cqs_tools       => $def->{cqstools},
    seqcount        => $def->{counts},
    pbs             => {
      'email'    => $def->{email},
      'walltime' => '72',
      'mem'      => '40gb',
      'nodes'    => '1:ppn=1'
    },
  },
  'gsnap_smallRNA_t2c' => {
    'pbs' => {
      'email'    => 'quanhu.sheng@vanderbilt.edu',
      'walltime' => '72',
      'mem'      => '20gb',
      'nodes'    => '1:ppn=1'
    },
    'sh_direct'  => 1,
    'perform'    => 1,
    'target_dir' => $def->{target_dir} . '/gsnap_smallRNA_t2c',
    'source_ref' => [ 'gsnap_smallRNA_count', '.mapped.xml$' ],
    'cqs_tools'  => '/home/shengq1/cqstools/CQS.Tools.exe',
    'option'     => '-p 0.05 -e 0.013',
    'class'      => 'CQS::ParclipT2CFinder'
  },
  gsnap_smallRNA_t2c_table => {
    class      => 'SmallRNA::T2CSummary',
    perform    => 1,
    target_dir => $def->{target_dir} . "/gsnap_smallRNA_t2c_table",
    option     => '',
    cluster    => 'slurm',
    sh_direct  => 1,
    source_ref => [ 'gsnap_smallRNA_count', '.mapped.xml$' ],
    cqs_tools  => $def->{cqstools},
    pbs        => {
      'email'    => 'quanhu.sheng@vanderbilt.edu',
      'walltime' => '72',
      'mem'      => '40gb',
      'nodes'    => '1:ppn=1'
    },
  },
  'gsnap_smallRNA_category' => {
    'pbs' => {
      'email'    => $def->{email},
      'walltime' => '72',
      'mem'      => '40gb',
      'nodes'    => '1:ppn=1'
    },
    'cluster'    => 'slurm',
    'sh_direct'  => 1,
    'perform'    => 1,
    'target_dir' => $def->{target_dir} . '/gsnap_smallRNA_category',
    'source_ref' => [ 'gsnap_smallRNA_count', '.info$' ],
    'cqs_tools'  => $def->{cqstools},
    'class'      => 'CQS::SmallRNACategory',
    'option'     => ''
  },
  'unmappedReads' => {
    'source2_ref' => [ 'gsnap_smallRNA_count', '.mapped.xml' ],
    'pbs'         => {
      'email'    => 'quanhu.sheng@vanderbilt.edu',
      'walltime' => '1',
      'mem'      => '10gb',
      'nodes'    => '1:ppn=1'
    },
    'perlFile'   => 'unmappedReadsToFastq.pl',
    'sh_direct'  => 1,
    'perform'    => 1,
    'target_dir' => $def->{target_dir} . '/unmappedReads',
    'source'     => $def->{identicals},
    'output_ext' => '_clipped_identical.unmapped.fastq.gz',
    'class'      => 'CQS::Perl'
  },
  'unmappedReads_bowtie1_genome_1mm' => {
    'pbs' => {
      'email'    => 'quanhu.sheng@vanderbilt.edu',
      'walltime' => '72',
      'mem'      => '40gb',
      'nodes'    => '1:ppn=8'
    },
    'cluster'       => 'slurm',
    'sh_direct'     => 1,
    'perform'       => 1,
    'target_dir'    => $def->{target_dir} . '/unmappedReads_bowtie1_genome_1mm',
    'samonly'       => 0,
    'source_ref'    => [ 'unmappedReads', '.fastq.gz$' ],
    'bowtie1_index' => '/scratch/cqs/shengq1/references/hg19_16569_MT/bowtie_index_1.1.2/hg19_16569_MT',
    'option'        => '-a -m 100 --best --strata -v 1 -p 8',
    'class'         => 'Bowtie1',
    'mappedonly'    => 1
  },
  'unmappedReads_bowtie1_genome_1mm_3utr_count' => {
    'pbs' => {
      'email'    => 'quanhu.sheng@vanderbilt.edu',
      'walltime' => '72',
      'mem'      => '20gb',
      'nodes'    => '1:ppn=1'
    },
    'sh_direct'       => 1,
    'perform'         => 1,
    'target_dir'      => $def->{target_dir} . '/unmappedReads_bowtie1_genome_1mm_3utr_count',
    'fastq_files_ref' => [ 'unmappedReads', '.fastq.gz$' ],
    'coordinate_file' => '/data/cqs/shengq1/reference/utr3/20140612_ucsc_hg19_3UTR.txt',
    'source_ref'      => [ 'unmappedReads_bowtie1_genome_1mm', '.bam$' ],
    'cqs_tools'       => '/home/shengq1/cqstools/CQS.Tools.exe',
    'seqcount'        => $def->{identicals_count},
    'option'          => '',
    'class'           => 'CQS::SmallRNACount'
  },

  'unmappedReads_bowtie1_genome_1mm_3utr_count_target' => {
    'pbs' => {
      'email'    => 'quanhu.sheng@vanderbilt.edu',
      'walltime' => '72',
      'mem'      => '20gb',
      'nodes'    => '1:ppn=1'
    },
    'fasta_file'   => '/gpfs21/scratch/cqs/shengq1/references/hg19_16569_M/hg19_16569_M.fa',
    'sh_direct'    => 1,
    'perform'      => 1,
    'target_dir'   => $def->{target_dir} . '/unmappedReads_bowtie1_genome_1mm_3utr_count_target',
    'source_ref'   => [ 'gsnap_smallRNA_t2c', '.xml$' ],
    'target_ref'   => [ 'unmappedReads_bowtie1_genome_1mm_3utr_count', '.xml$' ],
    'cqs_tools'    => '/home/shengq1/cqstools/CQS.Tools.exe',
    'option'       => '',
    'class'        => 'CQS::ParclipTarget',
    'refgene_file' => '/gpfs21/scratch/cqs/shengq1/references/hg19_16569_M/hg19_refgene.tsv'
  },
  'sequencetask' => {
    'pbs' => {
      'email'    => 'quanhu.sheng@vanderbilt.edu',
      'walltime' => '72',
      'mem'      => '40gb',
      'nodes'    => '1:ppn=8'
    },
    'source' =>
      { 'step1' => [ 'gsnap_smallRNA_t2c', 'unmappedReads', 'unmappedReads_bowtie1_genome_1mm', 'unmappedReads_bowtie1_genome_1mm_3utr_count', 'unmappedReads_bowtie1_genome_1mm_3utr_count_target' ] }
    ,
    'cluster'    => 'slurm',
    'sh_direct'  => 1,
    'perform'    => 1,
    'target_dir' => $def->{target_dir} . '/sequencetask',
    'class'      => 'CQS::SequenceTask',
    'option'     => ''
  },
};

performConfig($config);

1;

