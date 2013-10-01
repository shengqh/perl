cd /data/cqs/shengq1/reference/trna
mono-sgen /home/shengq1/cqstools/CQS.Tools.exe bed2fasta -i hg19_tRNA.bed -f /data/cqs/guoy1/reference/hg19/hg19_chr.fa -o hg19_tRNA.bed.fa
mono-sgen /home/shengq1/cqstools/CQS.Tools.exe bed2fasta -i mm10_tRNA.bed -f /data/cqs/guoy1/reference/mm10/mm10.fa -o mm10_tRNA.bed.fa --keepChrInName
mono-sgen /home/shengq1/cqstools/CQS.Tools.exe bed2fasta -i rn4_tRNA.bed -f /data/cqs/shengq1/reference/rn4/rn4.fa -o rn4_tRNA.bed.fa --keepChrInName

