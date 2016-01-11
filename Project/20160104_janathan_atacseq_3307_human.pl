#!/usr/bin/perl
use strict;
use warnings;

use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;
use Hash::Merge qw( merge );

my $target_dir = create_directory_or_die("/scratch/cqs/shengq1/atacseq/20160104_janathan_atacseq_3307_human");
my $email      = "quanhu.sheng\@vanderbilt.edu";

my $cqstools   = "/home/shengq1/cqstools/CQS.Tools.exe";
my $picard_jar = "/scratch/cqs/shengq1/local/bin/picard/picard.jar";
my $gatk_jar   = "/home/shengq1/local/bin/GATK/GenomeAnalysisTK.jar";
my $dbsnp      = "/scratch/cqs/shengq1/references/dbsnp/human_GRCh37_v142_16569_M.vcf";

my $macs2call_option = "-f BED -g mm -B -q 0.01 --nomodel";

my $bwa_fasta     = "/scratch/cqs/shengq1/references/hg19_tmp/bwa_index_0.7.12/hg19_16569_M.fa";
my $bowtie2_index = "/scratch/cqs/shengq1/references/hg19_16569_M/bowtie2_index_2.2.4/hg19_16569_M";

my $config = {
  general => { task_name => "3307" },
  files   => {
    "JDB1_1K_NoTNF"  => [ "/data/cqs/shengq1/data/3307/PE-150/3307-JDB-1_1_sequence.txt.gz", "/data/cqs/shengq1/data/3307/PE-150/3307-JDB-1_2_sequence.txt.gz" ],
    "JDB2_50K_NoTNF" => [ "/data/cqs/shengq1/data/3307/PE-150/3307-JDB-2_1_sequence.txt.gz", "/data/cqs/shengq1/data/3307/PE-150/3307-JDB-2_2_sequence.txt.gz" ],
    "JDB3_1K_TNF"    => [ "/data/cqs/shengq1/data/3307/PE-150/3307-JDB-3_1_sequence.txt.gz", "/data/cqs/shengq1/data/3307/PE-150/3307-JDB-3_2_sequence.txt.gz" ],
    "JDB4_50K_TNF"   => [ "/data/cqs/shengq1/data/3307/PE-150/3307-JDB-4_1_sequence.txt.gz", "/data/cqs/shengq1/data/3307/PE-150/3307-JDB-4_2_sequence.txt.gz" ],
    "JDB6_1K_NoTNF"  => [ "/data/cqs/shengq1/data/3307/PE-150/3307-JDB-6_1_sequence.txt.gz", "/data/cqs/shengq1/data/3307/PE-150/3307-JDB-6_2_sequence.txt.gz" ],
    "JDB7_1K_TNF"    => [ "/data/cqs/shengq1/data/3307/PE-150/3307-JDB-7_1_sequence.txt.gz", "/data/cqs/shengq1/data/3307/PE-150/3307-JDB-7_2_sequence.txt.gz" ],
    "JDB8_50K_NoTNF" => [ "/data/cqs/shengq1/data/3307/PE-150/3307-JDB-8_1_sequence.txt.gz", "/data/cqs/shengq1/data/3307/PE-150/3307-JDB-8_2_sequence.txt.gz" ],
    "JDB9_50K_TNF"   => [ "/data/cqs/shengq1/data/3307/PE-150/3307-JDB-9_1_sequence.txt.gz", "/data/cqs/shengq1/data/3307/PE-150/3307-JDB-9_2_sequence.txt.gz" ],
  },
  individual_comparison => {
    "1K_NoTNF_vs_1K_TNF_1"   => [ "JDB1_1K_NoTNF",  "JDB3_1K_TNF" ],
    "50K_NoTNF_vs_50K_TNF_1" => [ "JDB2_50K_NoTNF", "JDB4_50K_TNF" ],
    "1K_NoTNF_vs_1K_TNF_2"   => [ "JDB6_1K_NoTNF",  "JDB7_1K_TNF" ],
    "50K_NoTNF_vs_50K_TNF_2" => [ "JDB8_50K_NoTNF", "JDB9_50K_TNF" ],
  },
  fastqc => {
    class      => "QC::FastQC",
    perform    => 0,
    target_dir => "${target_dir}/fastqc",
    option     => "",
    source_ref => "files",
    sh_direct  => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=2",
      "walltime" => "2",
      "mem"      => "40gb"
    },
  },
  fastqc_summary => {
    class      => "QC::FastQCSummary",
    perform    => 0,
    target_dir => "${target_dir}/fastqc",
    option     => "",
    cqstools   => $cqstools,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "2",
      "mem"      => "10gb"
    },
  },
  bwa_pretrim => {
    class      => "Alignment::BWA",
    perform    => 0,
    target_dir => "${target_dir}/bwa_pretrim",
    option     => "",
    bwa_index  => $bwa_fasta,
    picard_jar => $picard_jar,
    source_ref => "files",
    sh_direct  => 0,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  bwa_pretrim_cleanbam => {
    class             => "ATACseq::CleanBam",
    perform           => 1,
    target_dir        => "${target_dir}/bwa_pretrim_cleanbam",
    option            => "",
    source_ref        => "bwa_pretrim",
    picard_jar        => $picard_jar,
    remove_chromosome => "M",
    sh_direct         => 0,
    pbs               => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "240",
      "mem"      => "40gb"
    },
  },
  bwa_pretrim_bam2bed => {
    class          => "ATACseq::BamToBed",
    perform        => 1,
    target_dir     => "${target_dir}/bwa_pretrim_bam2bed",
    option         => "",
    source_ref     => "bwa_pretrim_cleanbam",
    blacklist_file => "/scratch/cqs/shengq1/references/mappable_region/hg19/wgEncodeDacMapabilityConsensusExcludable.bed",
    sh_direct      => 0,
    pbs            => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  bwa_pretrim_macs2callpeak_individual_nomodel => {
    class      => "Chipseq::MACS2Callpeak",
    perform    => 1,
    target_dir => "${target_dir}/bwa_pretrim_macs2callpeak_individual_nomodel",
    option     => $macs2call_option,
    source_ref => "bwa_pretrim_bam2bed",
    sh_direct  => 0,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  bwa_pretrim_macs2bdgdiff_individual_nomodel => {
    class      => "Chipseq::MACS2Bdgdiff",
    perform    => 1,
    target_dir => "${target_dir}/bwa_pretrim_macs2bdgdiff_individual_nomodel",
    option     => "",
    source_ref => "bwa_pretrim_macs2callpeak_individual_nomodel",
    groups_ref => "individual_comparison",
    sh_direct  => 0,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  bowtie2_pretrim => {
    class         => "Alignment::Bowtie2",
    perform       => 1,
    target_dir    => "${target_dir}/bowtie2_pretrim",
    option        => "-k 1",
    bowtie2_index => $bowtie2_index,
    picard_jar    => $picard_jar,
    source_ref    => "files",
    sh_direct     => 0,
    pbs           => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  bowtie2_pretrim_cleanbam => {
    class             => "ATACseq::CleanBam",
    perform           => 1,
    target_dir        => "${target_dir}/bowtie2_pretrim_cleanbam",
    option            => "",
    source_ref        => "bowtie2_pretrim",
    picard_jar        => $picard_jar,
    remove_chromosome => "M",
    sh_direct         => 0,
    pbs               => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "240",
      "mem"      => "40gb"
    },
  },
  bowtie2_pretrim_bam2bed => {
    class          => "ATACseq::BamToBed",
    perform        => 1,
    target_dir     => "${target_dir}/bowtie2_pretrim_bam2bed",
    option         => "",
    source_ref     => "bowtie2_pretrim_cleanbam",
    blacklist_file => "/scratch/cqs/shengq1/references/mappable_region/hg19/wgEncodeDacMapabilityConsensusExcludable.bed",
    sh_direct      => 0,
    pbs            => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  bowtie2_pretrim_macs2callpeak_individual_nomodel => {
    class      => "Chipseq::MACS2Callpeak",
    perform    => 1,
    target_dir => "${target_dir}/bowtie2_pretrim_macs2callpeak_individual_nomodel",
    option     => $macs2call_option,
    source_ref => "bowtie2_pretrim_bam2bed",
    sh_direct  => 0,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  bowtie2_pretrim_macs2bdgdiff_individual_nomodel => {
    class      => "Chipseq::MACS2Bdgdiff",
    perform    => 1,
    target_dir => "${target_dir}/bowtie2_pretrim_macs2bdgdiff_individual_nomodel",
    option     => "",
    source_ref => "bowtie2_pretrim_macs2callpeak_individual_nomodel",
    groups_ref => "individual_comparison",
    sh_direct  => 0,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  cutadapt => {
    class      => "Trimmer::Cutadapt",
    perform    => 0,
    target_dir => "${target_dir}/cutadapt",
    option     => "-O 10 -m 50",
    source_ref => "files",
    adapter    => "CTGTCTCTTATA",
    extension  => "_clipped.fastq.gz",
    sh_direct  => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "24",
      "mem"      => "20gb"
    },
  },
  fastqlen => {
    class      => "CQS::FastqLen",
    perform    => 0,
    target_dir => "${target_dir}/fastqlen",
    option     => "",
    source_ref => "cutadapt",
    cqstools   => $cqstools,
    sh_direct  => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "24",
      "mem"      => "20gb"
    },
  },
  bowtie2_posttrim => {
    class         => "Alignment::Bowtie2",
    perform       => 1,
    target_dir    => "${target_dir}/bowtie2_posttrim",
    option        => "-k 1",
    bowtie2_index => $bowtie2_index,
    picard_jar    => $picard_jar,
    source_ref    => "cutadapt",
    sh_direct     => 0,
    pbs           => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  bowtie2_posttrim_cleanbam => {
    class             => "ATACseq::CleanBam",
    perform           => 1,
    target_dir        => "${target_dir}/bowtie2_posttrim_cleanbam",
    option            => "",
    source_ref        => "bowtie2_posttrim",
    picard_jar        => $picard_jar,
    remove_chromosome => "M",
    sh_direct         => 0,
    pbs               => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "240",
      "mem"      => "40gb"
    },
  },
  bowtie2_posttrim_bam2bed => {
    class          => "ATACseq::BamToBed",
    perform        => 1,
    target_dir     => "${target_dir}/bowtie2_posttrim_bam2bed",
    option         => "",
    source_ref     => "bowtie2_posttrim_cleanbam",
    blacklist_file => "/scratch/cqs/shengq1/references/mappable_region/hg19/wgEncodeDacMapabilityConsensusExcludable.bed",
    sh_direct      => 0,
    pbs            => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  bowtie2_posttrim_macs2callpeak_individual_nomodel => {
    class      => "Chipseq::MACS2Callpeak",
    perform    => 1,
    target_dir => "${target_dir}/bowtie2_posttrim_macs2callpeak_individual_nomodel",
    option     => $macs2call_option,
    source_ref => "bowtie2_posttrim_bam2bed",
    sh_direct  => 0,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  bowtie2_posttrim_macs2bdgdiff_individual_nomodel => {
    class      => "Chipseq::MACS2Bdgdiff",
    perform    => 1,
    target_dir => "${target_dir}/bowtie2_posttrim_macs2bdgdiff_individual_nomodel",
    option     => "",
    source_ref => "bowtie2_posttrim_macs2callpeak_individual_nomodel",
    groups_ref => "individual_comparison",
    sh_direct  => 0,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },

  sequencetask => {
    class      => "CQS::SequenceTask",
    perform    => 1,
    target_dir => "${target_dir}/sequencetask",
    option     => "",
    source     => {
      T1 => [ "bwa_pretrim",      "bwa_pretrim_cleanbam",      "bwa_pretrim_bam2bed",      "bwa_pretrim_macs2callpeak_individual_nomodel" ],
      T2 => [ "bowtie2_pretrim",  "bowtie2_pretrim_cleanbam",  "bowtie2_pretrim_bam2bed",  "bowtie2_pretrim_macs2callpeak_individual_nomodel" ],
      T3 => [ "bowtie2_posttrim", "bowtie2_posttrim_cleanbam", "bowtie2_posttrim_bam2bed", "bowtie2_posttrim_macs2callpeak_individual_nomodel" ],

    },
    sh_direct => 0,
    pbs       => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
    }

};

performConfig($config);

#performTask( $config, "bowtie2_pretrim" );

1;

