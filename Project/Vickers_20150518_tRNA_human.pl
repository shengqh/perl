#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;

my $config = {
  'general' => {
    'cluster'   => 'slurm',
    'task_name' => '20150518_tRNA'
  },
  'files' => {
    'KCVH08' => ['/gpfs21/scratch/cqs/shengq1/vickers/data/20150515_tRNA/KCVH8_S14_R1_001.fastq.gz'],
    'KCVH02' => ['/gpfs21/scratch/cqs/shengq1/vickers/data/20150515_tRNA/KCVH2_S8_R1_001.fastq.gz'],
    'KCVH03' => ['/gpfs21/scratch/cqs/shengq1/vickers/data/20150515_tRNA/KCVH3_S9_R1_001.fastq.gz'],
    'KCVH04' => ['/gpfs21/scratch/cqs/shengq1/vickers/data/20150515_tRNA/KCVH4_S10_R1_001.fastq.gz'],
    'KCVH09' => ['/gpfs21/scratch/cqs/shengq1/vickers/data/20150515_tRNA/KCVH9_S15_R1_001.fastq.gz'],
    'KCVH06' => ['/gpfs21/scratch/cqs/shengq1/vickers/data/20150515_tRNA/KCVH6_S12_R1_001.fastq.gz'],
    'KCVH07' => ['/gpfs21/scratch/cqs/shengq1/vickers/data/20150515_tRNA/KCVH7_S13_R1_001.fastq.gz'],
    'KCVH10' => ['/gpfs21/scratch/cqs/shengq1/vickers/data/20150515_tRNA/KCVH10_S7_R1_001.fastq.gz'],
    'KCVH05' => ['/gpfs21/scratch/cqs/shengq1/vickers/data/20150515_tRNA/KCVH5_S11_R1_001.fastq.gz'],
    'KCVH01' => ['/gpfs21/scratch/cqs/shengq1/vickers/data/20150515_tRNA/KCVH1_S6_R1_001.fastq.gz']
  },
  'fastqc_pre_trim' => {
    'pbs' => {
      'email'    => 'quanhu.sheng@vanderbilt.edu',
      'walltime' => '2',
      'mem'      => '10gb',
      'nodes'    => '1:ppn=1'
    },
    'source_ref' => 'files',
    'cluster'    => 'slurm',
    'perform'    => 1,
    'class'      => 'QC::FastQC',
    'target_dir' => '/scratch/cqs/shengq1/vickers/20150518_tRNA_human/fastqc_pre_trim',
    'option'     => ''
  },
  'fastqc_pre_trim_summary' => {
    'pbs' => {
      'email'    => 'quanhu.sheng@vanderbilt.edu',
      'walltime' => '2',
      'mem'      => '10gb',
      'nodes'    => '1:ppn=1'
    },
    'cluster'    => 'slurm',
    'perform'    => 1,
    'cqstools'   => '/home/shengq1/cqstools/CQS.Tools.exe',
    'class'      => 'QC::FastQCSummary',
    'target_dir' => '/scratch/cqs/shengq1/vickers/20150518_tRNA_human/fastqc_pre_trim',
    'option'     => ''
  },
  'cutadapt' => {
    'pbs' => {
      'email'    => 'quanhu.sheng@vanderbilt.edu',
      'walltime' => '24',
      'mem'      => '20gb',
      'nodes'    => '1:ppn=1'
    },
    'extension'  => '_clipped.fastq',
    'cluster'    => 'slurm',
    'sh_direct'  => 1,
    'perform'    => 1,
    'target_dir' => '/scratch/cqs/shengq1/vickers/20150518_tRNA_human/cutadapt',
    'source_ref' => 'files',
    'adapter'    => 'AGATCGGAAGAG',
    'option'     => '-O 10 -m 16',
    'class'      => 'Cutadapt'
  },
  'fastqc_post_trim' => {
    'pbs' => {
      'email'    => 'quanhu.sheng@vanderbilt.edu',
      'walltime' => '2',
      'mem'      => '10gb',
      'nodes'    => '1:ppn=1'
    },
    'source_ref' => [ 'cutadapt', '.fastq.gz' ],
    'cluster'    => 'slurm',
    'perform'    => 1,
    'class'      => 'QC::FastQC',
    'target_dir' => '/scratch/cqs/shengq1/vickers/20150518_tRNA_human/fastqc_post_trim',
    'option'     => ''
  },
  'fastqc_post_trim_summary' => {
    'pbs' => {
      'email'    => 'quanhu.sheng@vanderbilt.edu',
      'walltime' => '2',
      'mem'      => '10gb',
      'nodes'    => '1:ppn=1'
    },
    'cluster'    => 'slurm',
    'perform'    => 1,
    'cqstools'   => '/home/shengq1/cqstools/CQS.Tools.exe',
    'class'      => 'QC::FastQCSummary',
    'target_dir' => '/scratch/cqs/shengq1/vickers/20150518_tRNA_human/fastqc_post_trim',
    'option'     => ''
  },
  'identical' => {
    'pbs' => {
      'email'    => 'quanhu.sheng@vanderbilt.edu',
      'walltime' => '24',
      'mem'      => '20gb',
      'nodes'    => '1:ppn=1'
    },
    'cluster'    => 'slurm',
    'extension'  => '_clipped_identical.fastq.gz',
    'sh_direct'  => 1,
    'perform'    => 1,
    'target_dir' => '/scratch/cqs/shengq1/vickers/20150518_tRNA_human/identical',
    'source_ref' => [ 'cutadapt', '.fastq.gz' ],
    'cqstools'   => '/home/shengq1/cqstools/CQS.Tools.exe',
    'class'      => 'FastqIdentical',
    'option'     => ''
  },
  star => {
    class      => "Alignment::STAR",
    perform    => 1,
    target_dir => "/scratch/cqs/shengq1/vickers/20150518_tRNA_human/star",
    option     => "",
    source_ref => "identical",
    genome_dir => "/scratch/cqs/shengq1/references/hg19_16569_M/STAR_index_v37.75_2.4.0j_sjdb75",
    sh_direct  => 1,
    pbs        => {
      "email"    => 'quanhu.sheng@vanderbilt.edu',
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "30gb"
    },
  },
  'fastq_len' => {
    'pbs' => {
      'email'    => 'quanhu.sheng@vanderbilt.edu',
      'walltime' => '24',
      'mem'      => '20gb',
      'nodes'    => '1:ppn=1'
    },
    'cluster'    => 'slurm',
    'sh_direct'  => 1,
    'perform'    => 1,
    'target_dir' => '/scratch/cqs/shengq1/vickers/20150518_tRNA_human/fastq_len',
    'source_ref' => 'cutadapt',
    'cqstools'   => '/home/shengq1/cqstools/CQS.Tools.exe',
    'class'      => 'FastqLen',
    'option'     => ''
  },
  'identical_sequence_table' => {
    'pbs' => {
      'email'    => 'quanhu.sheng@vanderbilt.edu',
      'walltime' => '10',
      'mem'      => '10gb',
      'nodes'    => '1:ppn=1'
    },
    'cluster'    => 'slurm',
    'sh_direct'  => 1,
    'perform'    => 1,
    'target_dir' => '/scratch/cqs/shengq1/vickers/20150518_tRNA_human/identical_sequence_count_table',
    'source_ref' => [ 'identical', '.dupcount$' ],
    'cqs_tools'  => '/home/shengq1/cqstools/CQS.Tools.exe',
    'suffix'     => '_sequence',
    'class'      => 'CQS::SmallRNASequenceCountTable',
    'option'     => ''
  },
  'bowtie1_genome_3mm' => {
    'pbs' => {
      'email'    => 'quanhu.sheng@vanderbilt.edu',
      'walltime' => '72',
      'mem'      => '40gb',
      'nodes'    => '1:ppn=8'
    },
    'cluster'       => 'slurm',
    'sh_direct'     => 1,
    'perform'       => 1,
    'target_dir'    => '/scratch/cqs/shengq1/vickers/20150518_tRNA_human/bowtie1_genome_3mm',
    'samonly'       => 0,
    'source_ref'    => [ 'identical', '.fastq.gz$' ],
    'class'         => 'Bowtie1',
    'option'        => '-a -m 100 --best --strata -v 3 -p 8',
    'bowtie1_index' => '/data/cqs/guoy1/reference/hg19/bowtie_index_hg19_rCRS_1.0.0/hg19_rCRS'
  },
  'bowtie1_genome_3mm_smallRNA_count' => {
    'pbs' => {
      'email'    => 'quanhu.sheng@vanderbilt.edu',
      'walltime' => '72',
      'mem'      => '40gb',
      'nodes'    => '1:ppn=1'
    },
    'fasta_file'      => '/scratch/cqs/shengq1/references/smallrna/hg19_miRBase20_ucsc-tRNA_ensembl75.bed.fa',
    'cluster'         => 'slurm',
    'sh_direct'       => 1,
    'perform'         => 1,
    'target_dir'      => '/scratch/cqs/shengq1/vickers/20150518_tRNA_human/bowtie1_genome_3mm_smallRNA_count',
    'fastq_files_ref' => 'identical',
    'coordinate_file' => '/scratch/cqs/shengq1/references/smallrna/hg19_miRBase20_ucsc-tRNA_ensembl75.bed',
    'source_ref'      => 'bowtie1_genome_3mm',
    'cqs_tools'       => '/home/shengq1/cqstools/CQS.Tools.exe',
    'seqcount_ref'    => [ 'identical', '.dupcount$' ],
    'option'          => '-s -m 3',
    'class'           => 'CQS::SmallRNACount',
    'samtools'        => '/scratch/cqs/shengq1/local/bin/samtools'
  },

  'bowtie1_genome_3mm_smallRNA_category' => {
    'pbs' => {
      'email'    => 'quanhu.sheng@vanderbilt.edu',
      'walltime' => '72',
      'mem'      => '40gb',
      'nodes'    => '1:ppn=1'
    },
    'cluster'    => 'slurm',
    'sh_direct'  => 1,
    'perform'    => 1,
    'target_dir' => '/scratch/cqs/shengq1/vickers/20150518_tRNA_human/bowtie1_genome_3mm_smallRNA_category',
    'source_ref' => [ 'bowtie1_genome_3mm_smallRNA_count', '.info$' ],
    'cqs_tools'  => '/home/shengq1/cqstools/CQS.Tools.exe',
    'option'     => '',
    'class'      => 'CQS::SmallRNACategory'
  },
  'bowtie1_genome_3mm_smallRNA_table' => {
    'pbs' => {
      'email'    => 'quanhu.sheng@vanderbilt.edu',
      'walltime' => '10',
      'mem'      => '10gb',
      'nodes'    => '1:ppn=1'
    },
    'cluster'    => 'slurm',
    'sh_direct'  => 1,
    'perform'    => 1,
    'target_dir' => '/scratch/cqs/shengq1/vickers/20150518_tRNA_human/bowtie1_genome_3mm_smallRNA_table',
    'source_ref' => [ 'bowtie1_genome_3mm_smallRNA_count', '.mapped.xml' ],
    'cqs_tools'  => '/home/shengq1/cqstools/CQS.Tools.exe',
    'option'     => '',
    'class'      => 'CQS::SmallRNATable',
    'prefix'     => 'smallRNA_3mm_'
  },
  'bowtie1_genome_3mm_notidentical' => {
    'pbs' => {
      'email'    => 'quanhu.sheng@vanderbilt.edu',
      'walltime' => '72',
      'mem'      => '40gb',
      'nodes'    => '1:ppn=8'
    },
    'cluster'       => 'slurm',
    'sh_direct'     => 0,
    'perform'       => 1,
    'target_dir'    => '/scratch/cqs/shengq1/vickers/20150518_tRNA_human/bowtie1_genome_3mm_notidentical',
    'samonly'       => 0,
    'source_ref'    => [ 'cutadapt', '.fastq.gz' ],
    'class'         => 'Bowtie1',
    'option'        => '-a -m 100 --best --strata -v 3 -p 8',
    'bowtie1_index' => '/data/cqs/guoy1/reference/hg19/bowtie_index_hg19_rCRS_1.0.0/hg19_rCRS'
  },
  'sequencetask' => {
    'source' => {
      'step2' => [ 'fastqc_pre_trim_summary', 'fastqc_post_trim_summary', 'identical_sequence_table', 'bowtie1_genome_3mm_smallRNA_table', 'bowtie1_genome_3mm_smallRNA_category' ],
      'step1' => [ 'fastqc_pre_trim', 'cutadapt', 'fastqc_post_trim', 'fastq_len', 'identical', 'bowtie1_genome_3mm', 'bowtie1_genome_3mm_notidentical', 'bowtie1_genome_3mm_smallRNA_count' ]
    },
    'pbs' => {
      'email'    => 'quanhu.sheng@vanderbilt.edu',
      'walltime' => '72',
      'mem'      => '40gb',
      'nodes'    => '1:ppn=8'
    },
    'cluster'    => 'slurm',
    'sh_direct'  => 0,
    'perform'    => 1,
    'target_dir' => '/scratch/cqs/shengq1/vickers/20150518_tRNA_human/sequencetask',
    'class'      => 'CQS::SequenceTask',
    'option'     => ''
  },
};

performConfig($config);
1;

