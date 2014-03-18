#svn propset svn:eol-style LF bin/install_mRNAseq_TCGA.sh

cd /scratch/cqs/shengq1/local/bin/TCGA

wget -r -nc --no-parent --reject "index.html*,MapSplice_multi_threads_2.0.1.9.zip,hg19_M_rCRS.fa.tgz" -X "/public/mRNAseq_TCGA/hg19_M_rCRS/chromosomes/togo,/public/mRNAseq_TCGA/rsem_ref/bak" -e robots=off https://webshare.bioinf.unc.edu/public/mRNAseq_TCGA/

ln -s webshare.bioinf.unc.edu/public/mRNAseq_TCGA/ mRNAseq_TCGA

cd mRNAseq_TCGA

tar -xzvf MapSplice_multithreads_2.0.1.9.tar.gz
tar -xzvf MapSplice_multithreads_12_07.tar.gz

#ubu tool and required libraries

wget -nc https://github.com/downloads/mozack/ubu/ubu-1.0.jar
wget -nc http://repo1.maven.org/maven2/net/sf/jopt-simple/jopt-simple/4.6/jopt-simple-4.6.jar
wget -nc https://raw.github.com/mozack/ubu/master/src/perl/sort_bam_by_reference_and_name.pl
wget -nc http://apache.mirrors.hoobly.com//commons/collections/binaries/commons-collections-3.2.1-bin.tar.gz
tar -zxvf commons-collections-3.2.1-bin.tar.gz
mv commons-collections-3.2.1/commons-collections-3.2.1.jar .
rm -rf commons-collections-3.2.1 

#rsem

wget -nc http://deweylab.biostat.wisc.edu/rsem/src/rsem-1.1.13.tar.gz
tar -xzvf rsem-1.1.13.tar.gz
cd rsem-1.1.13
make

#samtools, the version of the one that mapsplice attached missing a function cat

wget -nc https://downloads.sourceforge.net/project/samtools/samtools/0.1.19/samtools-0.1.19.tar.bz2
bzip2 -cd samtools-0.1.19.tar.bz2 | tar xf -
cd samtools-0.1.19
make



