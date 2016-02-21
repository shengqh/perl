#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use File::Spec;
use Getopt::Long;

sub run_command {
  my $command = shift;
  print "$command \n";
  `$command `;
}

my $usage = "

Synopsis:

buildindex -f fastaFile [Options]

Options:

  -f|--file {fastaFile}        Fasta format sequence file
  -b|--bowtie                  Build bowtie index
  -w|--bwa                     Build bwa index
  -s|--star                    Build STAR index
  --sjdbGTFfile                Option sjdbGTFfile for STAR index
  --sjdbGTFfileVersion         Option sjdbGTFfile version for STAR index
  --sjdbOverhang               Option sjdbOverhang for STAR index
  --thread                     Option thread for STAR index
  -B|--bowtie2                 Build bowtie2 index
  -g|--gsnap                   Build gsnap index
  -h|--help                    This page.
";

Getopt::Long::Configure('bundling');

my $fastaFile;
my $dogsnap;
my $dobowtie;
my $dobowtie2;
my $dobwa;
my $dostar;
my $sjdbGTFfile;
my $sjdbOverhang;
my $thread;
my $help;
my $sjdbGTFfileVersion;

GetOptions(
  'h|help'               => \$help,
  'f|file=s'             => \$fastaFile,
  'g|gsnap'              => \$dogsnap,
  'b|bowtie'             => \$dobowtie,
  'B|bowtie2'            => \$dobowtie2,
  'w|bwa'                => \$dobwa,
  's|star'               => \$dostar,
  'sjdbGTFfile=s'        => \$sjdbGTFfile,
  'sjdbGTFfileVersion=s' => \$sjdbGTFfileVersion,
  'sjdbOverhang=i'       => \$sjdbOverhang,
  'thread=i'             => \$thread,
);

my $pass = 1;

if ( defined $help ) {
  print $usage;
  exit(1);
}

if ( !defined($fastaFile) ) {
  print STDERR "file is required for build index." . "\n";
  $pass = 0;
}
elsif ( !-e $sjdbGTFfile ) {
  print STDERR "sjdbGTFfile is not exist " . $sjdbGTFfile . "\n";
  $pass = 0;
}

if ( defined $dostar ) {
  if ( !defined $sjdbGTFfile ) {
    print STDERR "sjdbGTFfile is required for STAR index." . "\n";
    $pass = 0;
  }
  elsif ( !-e $sjdbGTFfile ) {
    print STDERR "sjdbGTFfile is not exist " . $sjdbGTFfile . "\n";
    $pass = 0;
  }

  if ( !defined $sjdbGTFfileVersion ) {
    print STDERR "sjdbGTFfileVersion is required for STAR index." . "\n";
    $pass = 0;
  }

  if ( !defined $sjdbOverhang ) {
    $sjdbOverhang = 99;
  }

  if ( !defined $thread ) {
    $thread = 8;
  }
}

if ( !$pass ) {
  print $usage;
  exit(1);
}

my $basename = basename($fastaFile);
( my $base = $basename ) =~ s/\.[^.]+$//;

# index fasta file
if ( !-e "${basename}.fai" ) {
  run_command("samtools faidx $fastaFile ");
}

if ( !-e "${base}.dict" ) {
  run_command("java -jar /scratch/cqs/shengq1/local/bin/picard/picard.jar CreateSequenceDictionary R=$fastaFile O=${base}.dict");
}

if ( !-e "${base}.len" ) {
  run_command("perl /home/shengq1/local/bin/get_fasta_lengths.pl $fastaFile");
  run_command("mv res_${basename} ${base}.len ");
}

my $absolute_dir = File::Spec->rel2abs(".");

if ( defined $dobowtie ) {

  # bowtie
  my $bowtie = `bowtie --version | grep bowtie | grep version | cut -d " " -f 3`;
  chomp($bowtie);
  if ( !-e "bowtie_index_${bowtie}" ) {
    mkdir("bowtie_index_${bowtie}");
    chdir("bowtie_index_${bowtie}");
    if ( !-e $basename ) {
      run_command("ln -s ../$fastaFile $basename ");
      run_command("ln -s ../${base}.dict ${base}.dict ");
      run_command("ln -s ../${basename}.fai ${basename}.fai ");
      run_command("ln -s ../${base}.len ${base}.len ");
    }
    run_command("bowtie-build $basename $base ");
    chdir($absolute_dir);
  }
}

if ( defined $dobwa ) {

  #bwa
  `bwa 2> 1`;
  my $bwa = `grep Version 1 | cut -d " " -f 2 | cut -d "-" -f 1`;
  chomp($bwa);
  `rm 1`;
  if ( !-e "bwa_index_${bwa}" ) {
    mkdir("bwa_index_${bwa}");
    chdir("bwa_index_${bwa}");
    if ( !-e $basename ) {
      run_command("ln -s ../$fastaFile $basename ");
      run_command("ln -s ../${base}.dict ${base}.dict ");
      run_command("ln -s ../${basename}.fai ${basename}.fai ");
      run_command("ln -s ../${base}.len ${base}.len ");
    }
    print "bwa index $basename \n";
    run_command("bwa index $basename");

    chdir($absolute_dir);
  }
}

if ( defined $dostar ) {

  # STAR
  my $star = `STAR --version`;
  chomp($star);
  my $star_dir = "STAR_index_${star}_${sjdbGTFfileVersion}_sjdb${sjdbOverhang}";
  
  my $absolute_gtf = File::Spec->rel2abs($sjdbGTFfile);

  if ( !-e $star_dir ) {
    mkdir($star_dir);
    chdir($star_dir);
    if ( !-e $basename ) {
      run_command("ln -s ../$fastaFile $basename ");
      run_command("ln -s ../${base}.dict ${base}.dict ");
      run_command("ln -s ../${basename}.fai ${basename}.fai ");
      run_command("ln -s ../${base}.len ${base}.len ");
    }
    run_command("STAR --runThreadN $thread --runMode genomeGenerate --genomeDir . --genomeFastaFiles $basename --sjdbGTFfile $absolute_gtf --sjdbOverhang $sjdbOverhang");
    chdir($absolute_dir);
  }
}

if ( defined $dobowtie2 ) {

  # bowtie2
  my $bowtie2 = `bowtie2 --version | grep bowtie2 | grep version | cut -d " " -f 3`;
  chomp($bowtie2);
  if ( !-e "bowtie2_index_${bowtie2}" ) {
    mkdir("bowtie2_index_${bowtie2}");
    chdir("bowtie2_index_${bowtie2}");
    if ( !-e $basename ) {
      run_command("ln -s ../$fastaFile $basename ");
      run_command("ln -s ../${base}.dict ${base}.dict ");
      run_command("ln -s ../${basename}.fai ${basename}.fai ");
      run_command("ln -s ../${base}.len ${base}.len ");
    }
    run_command("bowtie2-build $basename $base ");
    chdir("..");
  }
}

if ( defined $dogsnap ) {

  # gsnap
  `gsnap 2> 1`;
  my $gsnap = `grep version 1 | cut -d " " -f 3`;
  chomp($gsnap);
  `rm 1`;
  if ( !-e "gsnap_index_k14_${gsnap}" ) {
    mkdir("gsnap_index_k14_${gsnap}");
    chdir("gsnap_index_k14_${gsnap}");
    if ( !-e $basename ) {
      run_command("ln -s ../$fastaFile $basename ");
      run_command("ln -s ../${base}.dict ${base}.dict ");
      run_command("ln -s ../${basename}.fai ${basename}.fai ");
      run_command("ln -s ../${base}.len ${base}.len ");
    }

    #min_read_length = kmers + interval - 1
    #in order to control the min_read_length = 16, we have to smaller the kmers from 15 to 14 when keep the "sampling interval for genome" equals 3
    run_command("gmap_build -D . -d $base -k 14 -s none $basename");
    run_command("atoiindex -F . -d $base");
    chdir($absolute_dir);
  }
}

1;
