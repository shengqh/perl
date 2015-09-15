#!/usr/bin/perl
use strict;
use warnings;

use CQS::PerformSmallRNA;
use CQS::ClassFactory;

my $config = {
  'general' => {
    'cluster'   => 'slurm',
    'task_name' => 'wael'
  },
  'files' =>{
    "T169"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/B/D/T169.fq.gz"],
    "T173"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/B/D/T173.fq.gz"],
    "T176"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/B/D/T176.fq.gz"],
    "T178"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/B/D/T178.fq.gz"],
    "T179"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/B/D/T179.fq.gz"],
    "T180"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/B/D/T180.fq.gz"],
    "T193"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/B/F/T193.fq.gz"],
  },
  'valid_fastq' => {
    'pbs' => {
      'email'    => 'quanhu.sheng@vanderbilt.edu',
      'walltime' => '24',
      'mem'      => '20gb',
      'nodes'    => '1:ppn=1'
    },
    'cluster'    => 'slurm',
    'sh_direct'  => 1,
    'perform'    => 1,
    'target_dir' => '/scratch/cqs/shengq1/smallRNA/20150910_guoyan_smallRNA_human/valid_fastq',
    'source_ref' => 'files',
    'option'     => '',
    'class'      => 'CQS::ValidFastqExtractor'
  },
};

performConfig($config);

my $def = {

  #General options
  task_name  => "wael",
  email      => "quanhu.sheng\@vanderbilt.edu",
  target_dir => "/scratch/cqs/shengq1/smallRNA/20150910_guoyan_smallRNA_human",
  max_thread => 8,

  #Default software parameter (don't change it except you really know it)
  fastq_remove_N => 0,

  #Data
  files => {
    "1T4"      => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/B/D/1T4.fq.gz"],
    "1T40"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/B/D/1T40.fq.gz"],
    "1T42"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/B/D/1T42.fq.gz"],
    "N103"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/C/N103.fq.gz"],
    "N111"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/C/N111.fq.gz"],
    "N112"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/C/N112.fq.gz"],
    "N114"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/C/N114.fq.gz"],
    "N119"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/C/N119.fq.gz"],
    "N125"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/C/N125.fq.gz"],
    "N128"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/E/N128.fq.gz"],
    "N132"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/E/N132.fq.gz"],
    "N134"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/E/N134.fq.gz"],
    "N145"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/E/N145.fq.gz"],
    "N148"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/E/N148.fq.gz"],
    "N151"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/E/N151.fq.gz"],
    "N153"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/E/N153.fq.gz"],
    "N156"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/E/N156.fq.gz"],
    "N166"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/E/N166.fq.gz"],
    "N169"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/C/N169.fq.gz"],
    "N170"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/C/N170.fq.gz"],
    "N172"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/C/N172.fq.gz"],
    "N173"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/C/N173.fq.gz"],
    "N174"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/C/N174.fq.gz"],
    "N175"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/C/N175.fq.gz"],
    "N177"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/C/N177.fq.gz"],
    "N178"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/C/N178.fq.gz"],
    "N179"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/C/N179.fq.gz"],
    "N180"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/C/N180.fq.gz"],
    "N181"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/C/N181.fq.gz"],
    "N187"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/C/N187.fq.gz"],
    "N188"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/C/N188.fq.gz"],
    "N190"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/E/N190.fq.gz"],
    "N191"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/E/N191.fq.gz"],
    "N192"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/E/N192.fq.gz"],
    "N193"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/E/N193.fq.gz"],
    "N194"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/E/N194.fq.gz"],
    "N197"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/E/N197.fq.gz"],
    "N198"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/E/N198.fq.gz"],
    "N199"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/E/N199.fq.gz"],
    "N200"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/E/N200.fq.gz"],
    "N201"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/E/N201.fq.gz"],
    "normal01" => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/A/normal1.fq.gz"],
    "normal02" => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/A/normal2.fq.gz"],
    "normal03" => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/A/normal3.fq.gz"],
    "normal04" => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/A/normal4.fq.gz"],
    "normal05" => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/A/normal5.fq.gz"],
    "normal06" => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/A/normal6.fq.gz"],
    "normal08" => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/A/normal8.fq.gz"],
    "normal09" => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/A/normal9.fq.gz"],
    "normal10" => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/A/normal10.fq.gz"],
    "T029"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/B/D/T29.fq.gz"],
    "T042"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/B/D/T42.fq.gz"],
    "T103"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/B/D/T103.fq.gz"],
    "T111"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/B/D/T111.fq.gz"],
    "T112"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/B/D/T112.fq.gz"],
    "T114"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/B/D/T114.fq.gz"],
    "T119"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/B/D/T119.fq.gz"],
    "T120"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/B/D/T120.fq.gz"],
    "T121"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/B/D/T121.fq.gz"],
    "T125"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/B/D/T125.fq.gz"],
    "T128"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/B/F/T128.fq.gz"],
    "T129"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/B/F/T129.fq.gz"],
    "T130"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/B/F/T130.fq.gz"],
    "T131"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/B/F/T131.fq.gz"],
    "T132"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/B/F/T132.fq.gz"],
    "T134"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/B/F/T134.fq.gz"],
    "T135"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/B/F/T135.fq.gz"],
    "T136"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/B/F/T136.fq.gz"],
    "T138"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/B/F/T138.fq.gz"],
    "T142"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/B/F/T142.fq.gz"],
    "T143"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/B/F/T143.fq.gz"],
    "T144"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/B/F/T144.fq.gz"],
    "T145"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/B/F/T145.fq.gz"],
    "T146"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/B/F/T146.fq.gz"],
    "T147"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/B/F/T147.fq.gz"],
    "T148"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/B/F/T148.fq.gz"],
    "T149"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/B/F/T149.fq.gz"],
    "T150"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/B/F/T150.fq.gz"],
    "T151"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/B/F/T151.fq.gz"],
    "T152"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/B/F/T152.fq.gz"],
    "T153"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/B/F/T153.fq.gz"],
    "T154"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/B/F/T154.fq.gz"],
    "T155"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/B/F/T155.fq.gz"],
    "T156"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/B/F/T156.fq.gz"],
    "T157"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/B/F/T157.fq.gz"],
    "T159"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/B/F/T159.fq.gz"],
    "T162"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/B/F/T162.fq.gz"],
    "T166"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/B/F/T166.fq.gz"],

    #"T169"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/B/D/T169.fq.gz"],
    "T172" => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/B/D/T172.fq.gz"],

    #"T173"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/B/D/T173.fq.gz"],
    "T174" => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/B/D/T174.fq.gz"],
    "T175" => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/B/D/T175.fq.gz"],

    #"T176"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/B/D/T176.fq.gz"],
    "T177" => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/B/D/T177.fq.gz"],

    #"T178"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/B/D/T178.fq.gz"],
    #"T179"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/B/D/T179.fq.gz"],
    #"T180"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/B/D/T180.fq.gz"],
    "T181" => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/B/D/T181.fq.gz"],
    "T182" => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/B/D/T182.fq.gz"],
    "T183" => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/B/D/T183.fq.gz"],
    "T184" => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/B/D/T184.fq.gz"],
    "T185" => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/B/D/T185.fq.gz"],
    "T186" => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/B/D/T186.fq.gz"],
    "T187" => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/B/D/T187.fq.gz"],
    "T188" => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/B/D/T188.fq.gz"],
    "T189" => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/B/D/T189.fq.gz"],
    "T190" => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/B/F/T190.fq.gz"],
    "T191" => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/B/F/T191.fq.gz"],
    "T192" => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/B/F/T192.fq.gz"],

    #"T193"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/B/F/T193.fq.gz"],
    #"T194"     => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/B/F/T194.fq.gz"],
    "T196" => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/B/F/T196.fq.gz"],
    "T197" => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/B/F/T197.fq.gz"],
    "T198" => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/B/F/T198.fq.gz"],
    "T199" => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/B/F/T199.fq.gz"],
    "T200" => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/B/F/T200.fq.gz"],
    "T201" => ["/gpfs21/scratch/cqs/guoy1/wael/miRNA/B/F/T201.fq.gz"],
  },
};

#performSmallRNA_hg19($def);

1;

