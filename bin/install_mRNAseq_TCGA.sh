#svn propset svn:eol-style LF bin/install_mRNAseq_TCGA.sh

cd /scratch/cqs/shengq1/local/bin/TCGA

wget -r -nc --no-parent --reject "index.html*,MapSplice_multi_threads_2.0.1.9.zip,hg19_M_rCRS.fa.tgz" -X "/public/mRNAseq_TCGA/hg19_M_rCRS/chromosomes/togo,/public/mRNAseq_TCGA/rsem_ref/bak,/public/mRNAseq_TCGA/rsem_ref/rsem" -e robots=off https://webshare.bioinf.unc.edu/public/mRNAseq_TCGA/

ln -s webshare.bioinf.unc.edu/public/mRNAseq_TCGA/ mRNAseq_TCGA

cd mRNAseq_TCGA
wget -nc https://github.com/downloads/mozack/ubu/ubu-1.0.jar

wget -nc http://deweylab.biostat.wisc.edu/rsem/src/rsem-1.1.13.tar.gz
tar -xzvf rsem-1.1.13.tar.gz
cd rsem-1.1.13
make

tar -xzvf MapSplice_multithreads_2.0.1.9.tar.gz
tar -xzvf MapSplice_multithreads_12_07.tar.gz

cd MapSplice_multi_threads_2.0.1.9/samtools-0.1.9/
make



