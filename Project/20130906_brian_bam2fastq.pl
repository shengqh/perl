#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;

my $target_dir = create_directory_or_die("/scratch/shengq1/20130906_brian_bam2fastq");

my $email    = "quanhu.sheng\@vanderbilt.edu";
my $cqstools = "/home/shengq1/cqstools/CQS.Tools.exe";
my $samtools = "/home/shengq1/local/bin/samtools/samtools";

my $config = {
  general  => { task_name => "bam2fastq" },
  rnafiles => {
    bamfiles => {
      "UNCID_1092669.31757dda-12d6-4f35-87f1-d00464f2815b" =>
        ["/gpfs20/data/lehmanbd/VirusSeq/MMTV_RNAseq/UNCID_1092669.31757dda-12d6-4f35-87f1-d00464f2815b.sorted_genome_alignments.resort.fixed.bam"],
      "UNCID_1115718.20486352-7d5f-4065-ab10-7a3f792ec1e1" =>
        ["/gpfs20/data/lehmanbd/VirusSeq/MMTV_RNAseq/UNCID_1115718.20486352-7d5f-4065-ab10-7a3f792ec1e1.sorted_genome_alignments.resort.fixed.bam"],
      "UNCID_1118137.a74b20ca-06f9-47ec-8888-cb01e0ba3257" =>
        ["/gpfs20/data/lehmanbd/VirusSeq/MMTV_RNAseq/UNCID_1118137.a74b20ca-06f9-47ec-8888-cb01e0ba3257.sorted_genome_alignments.resort.fixed.bam"],
      "UNCID_1118269.ab3e7f05-9fb2-4c19-a6a9-3ffe1e4c1525" =>
        ["/gpfs20/data/lehmanbd/VirusSeq/MMTV_RNAseq/UNCID_1118269.ab3e7f05-9fb2-4c19-a6a9-3ffe1e4c1525.sorted_genome_alignments.resort.fixed.bam"],
      "UNCID_1118498.526c9170-0f25-4374-9287-1fc40acb6555" =>
        ["/gpfs20/data/lehmanbd/VirusSeq/MMTV_RNAseq/UNCID_1118498.526c9170-0f25-4374-9287-1fc40acb6555.sorted_genome_alignments.resort.fixed.bam"],
      "UNCID_1125970.5213da80-35ac-4172-b2e8-de686e73a498" =>
        ["/gpfs20/data/lehmanbd/VirusSeq/MMTV_RNAseq/UNCID_1125970.5213da80-35ac-4172-b2e8-de686e73a498.sorted_genome_alignments.resort.fixed.bam"],
      "UNCID_1126558.8131b824-c3db-4c83-99b9-aa0c60747c5d" =>
        ["/gpfs20/data/lehmanbd/VirusSeq/MMTV_RNAseq/UNCID_1126558.8131b824-c3db-4c83-99b9-aa0c60747c5d.sorted_genome_alignments.resort.fixed.bam"],
      "UNCID_1126711.75b2227e-0279-4dfb-a1a1-e64c53b48d7a" =>
        ["/gpfs20/data/lehmanbd/VirusSeq/MMTV_RNAseq/UNCID_1126711.75b2227e-0279-4dfb-a1a1-e64c53b48d7a.sorted_genome_alignments.resort.fixed.bam"],
      "UNCID_1127283.d31b476c-c047-4e5f-a0cb-eff1065f0121" =>
        ["/gpfs20/data/lehmanbd/VirusSeq/MMTV_RNAseq/UNCID_1127283.d31b476c-c047-4e5f-a0cb-eff1065f0121.sorted_genome_alignments.resort.fixed.bam"],
      "UNCID_1127508.6d066a72-f59f-45a8-ab90-216000b36da4" =>
        ["/gpfs20/data/lehmanbd/VirusSeq/MMTV_RNAseq/UNCID_1127508.6d066a72-f59f-45a8-ab90-216000b36da4.sorted_genome_alignments.resort.fixed.bam"],
      "UNCID_1128643.fb03a31f-6611-4343-816d-37a332147eac" =>
        ["/gpfs20/data/lehmanbd/VirusSeq/MMTV_RNAseq/UNCID_1128643.fb03a31f-6611-4343-816d-37a332147eac.sorted_genome_alignments.resort.fixed.bam"],
      "UNCID_1128735.e15ff180-38d3-4614-88ee-3645afaab54b" =>
        ["/gpfs20/data/lehmanbd/VirusSeq/MMTV_RNAseq/UNCID_1128735.e15ff180-38d3-4614-88ee-3645afaab54b.sorted_genome_alignments.resort.fixed.bam"],
      "UNCID_1129979.2aa34033-a6f5-45a8-833b-b433352e58ba" =>
        ["/gpfs20/data/lehmanbd/VirusSeq/MMTV_RNAseq/UNCID_1129979.2aa34033-a6f5-45a8-833b-b433352e58ba.sorted_genome_alignments.resort.fixed.bam"],
      "UNCID_1131613.b05828bc-eb75-4220-9972-6c6d62fa5303" =>
        ["/gpfs20/data/lehmanbd/VirusSeq/MMTV_RNAseq/UNCID_1131613.b05828bc-eb75-4220-9972-6c6d62fa5303.sorted_genome_alignments.resort.fixed.bam"],
      "UNCID_1138621.ddb62dc3-293b-4021-905a-83ae10aa0f8e" =>
        ["/gpfs20/data/lehmanbd/VirusSeq/MMTV_RNAseq/UNCID_1138621.ddb62dc3-293b-4021-905a-83ae10aa0f8e.sorted_genome_alignments.resort.fixed.bam"],
      "UNCID_1141399.72c34d2c-f255-492b-9249-e9e9c506ece2" =>
        ["/gpfs20/data/lehmanbd/VirusSeq/MMTV_RNAseq/UNCID_1141399.72c34d2c-f255-492b-9249-e9e9c506ece2.sorted_genome_alignments.resort.fixed.bam"],
      "UNCID_1148831.f091d050-4847-4d4d-b94c-1efacecdef6b" =>
        ["/gpfs20/data/lehmanbd/VirusSeq/MMTV_RNAseq/UNCID_1148831.f091d050-4847-4d4d-b94c-1efacecdef6b.sorted_genome_alignments.resort.fixed.bam"],
      "UNCID_1148929.d73174ef-182c-48b9-af79-dfb457d5b7be" =>
        ["/gpfs20/data/lehmanbd/VirusSeq/MMTV_RNAseq/UNCID_1148929.d73174ef-182c-48b9-af79-dfb457d5b7be.sorted_genome_alignments.resort.fixed.bam"],
      "UNCID_1148937.7c7b9232-aadf-476f-b9c7-2fdd67efe1e7" =>
        ["/gpfs20/data/lehmanbd/VirusSeq/MMTV_RNAseq/UNCID_1148937.7c7b9232-aadf-476f-b9c7-2fdd67efe1e7.sorted_genome_alignments.resort.fixed.bam"],
      "UNCID_1149453.ae592879-5689-4818-8ac3-23a00c5d3934" =>
        ["/gpfs20/data/lehmanbd/VirusSeq/MMTV_RNAseq/UNCID_1149453.ae592879-5689-4818-8ac3-23a00c5d3934.sorted_genome_alignments.resort.fixed.bam"],
      "UNCID_1151896.01370d42-f75c-4532-9b9c-24ff7302b033" =>
        ["/gpfs20/data/lehmanbd/VirusSeq/MMTV_RNAseq/UNCID_1151896.01370d42-f75c-4532-9b9c-24ff7302b033.sorted_genome_alignments.resort.fixed.bam"],
      "UNCID_1155522.95b5d7ca-4078-4af9-ac6b-4ed1c3cd9c60" =>
        ["/gpfs20/data/lehmanbd/VirusSeq/MMTV_RNAseq/UNCID_1155522.95b5d7ca-4078-4af9-ac6b-4ed1c3cd9c60.sorted_genome_alignments.resort.fixed.bam"],
      "UNCID_1155915.0a9215ea-d5d2-41ad-a00d-bd7722f41229" =>
        ["/gpfs20/data/lehmanbd/VirusSeq/MMTV_RNAseq/UNCID_1155915.0a9215ea-d5d2-41ad-a00d-bd7722f41229.sorted_genome_alignments.resort.fixed.bam"],
      "UNCID_1157461.619ffb2b-a591-4f3b-aa48-b0c44f4ee36a" =>
        ["/gpfs20/data/lehmanbd/VirusSeq/MMTV_RNAseq/UNCID_1157461.619ffb2b-a591-4f3b-aa48-b0c44f4ee36a.sorted_genome_alignments.resort.fixed.bam"],
      "UNCID_1161272.acd59617-770c-4653-891f-71ea87d67c41" =>
        ["/gpfs20/data/lehmanbd/VirusSeq/MMTV_RNAseq/UNCID_1161272.acd59617-770c-4653-891f-71ea87d67c41.sorted_genome_alignments.resort.fixed.bam"],
      "UNCID_1167393.a2405d64-34eb-4915-abf7-8530151d5cb0" =>
        ["/gpfs20/data/lehmanbd/VirusSeq/MMTV_RNAseq/UNCID_1167393.a2405d64-34eb-4915-abf7-8530151d5cb0.sorted_genome_alignments.resort.fixed.bam"],
      "UNCID_1169322.97168a1f-abf8-414b-af10-b63f5daa7023" =>
        ["/gpfs20/data/lehmanbd/VirusSeq/MMTV_RNAseq/UNCID_1169322.97168a1f-abf8-414b-af10-b63f5daa7023.sorted_genome_alignments.resort.fixed.bam"],
      "UNCID_1338259.42ba8e87-8d37-4cbd-82f0-bca0cbe5dcd3" =>
        ["/gpfs20/data/lehmanbd/VirusSeq/MMTV_RNAseq/UNCID_1338259.42ba8e87-8d37-4cbd-82f0-bca0cbe5dcd3.sorted_genome_alignments.resort.fixed.bam"],
      "UNCID_1673778.09d1b9e5-a828-4e12-8a66-2a90c18d2f3e" =>
        ["/gpfs20/data/lehmanbd/VirusSeq/MMTV_RNAseq/UNCID_1673778.09d1b9e5-a828-4e12-8a66-2a90c18d2f3e.sorted_genome_alignments.resort.fixed.bam"],
      "UNCID_1806145.9a9d7ee5-e3fa-491f-ab5b-0ab69652b8a5" =>
        ["/gpfs20/data/lehmanbd/VirusSeq/MMTV_RNAseq/UNCID_1806145.9a9d7ee5-e3fa-491f-ab5b-0ab69652b8a5.sorted_genome_alignments.resort.fixed.bam"],
    },
  },
  bam2fastq => {
    class      => "Bam2Fastq",
    perform    => 1,
    target_dir => "${target_dir}/bam2fastq",
    option     => "",
    source_ref => "bamfiles",
    cqstools   => $cqstools,
    samtools   => $samtools,
    ispaired   => 1,
    sh_direct  => 0,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "20gb"
    },
  },
};

performConfig($config);

1;
