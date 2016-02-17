#!/usr/bin/perl
use strict;
use warnings;

use CQS::PerformSmallRNA;
use CQS::ClassFactory;


#!/usr/bin/perl
use strict;
use warnings;

use CQS::PerformSmallRNA;

my $def = {

        #General options
        task_name => "3018-KCV-52_53_54",
        email     => "quanhu.sheng\@vanderbilt.edu",
        target_dir =>
          "/scratch/cqs/zhaos/vickers/20160211_smallRNA_3018-KCV-52_53_54_sheep",
        max_thread => 8,
        cqstools   => "/home/shengq1/cqstools/CQS.Tools.exe",

        #Default software parameter (don't change it except you really know it)
        fastq_remove_N => 0,

        remove_sequences      => "'CCACGTTCCCGTGG;ACAGTCCGACGATC'",
        search_unmapped_reads => 1,
        fastq_remove_random=>4,
          blast_unmapped_reads => 1,

        #Data
        files => {
  "SheepHDL" => ["/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-52/3018-KCV-52-i9_S9_R1_001.fastq.gz"],
  "SheepLDL" => ["/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-53/3018-KCV-53-i25_S9_R1_001.fastq.gz"],
  "SheepVLDL" => ["/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-54/3018-KCV-54-i40_S9_R1_001.fastq.gz"],
        },
};

my $config = performSmallRNA_oar3($def, 0);
$config->{bowtie1_unmapped_reads}{target_dir} = "/scratch/cqs/shengq1/vickers/20160211_smallRNA_3018-KCV-52_53_54_sheep_blastn/bowtie1_unmapped_reads";
$config->{bowtie1_unmapped_sequence_count_table}{target_dir} = "/scratch/cqs/shengq1/vickers/20160211_smallRNA_3018-KCV-52_53_54_sheep_blastn/bowtie1_unmapped_sequence_count_table";
$config->{bowtie1_unmapped_sequence_blast}{target_dir} = "/scratch/cqs/shengq1/vickers/20160211_smallRNA_3018-KCV-52_53_54_sheep_blastn/bowtie1_unmapped_sequence_blast";

performTask($config, "bowtie1_unmapped_reads");
performTask($config, "bowtie1_unmapped_sequence_count_table");
performTask($config, "bowtie1_unmapped_sequence_blast");

1;

