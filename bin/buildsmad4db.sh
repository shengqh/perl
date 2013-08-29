cd /data/cqs/shengq1/reference/mm10/smad4 

grep Smad4 /data/cqs/guoy1/reference/annotation2/mm10/Mus_musculus.GRCm38.68.gtf > smad4.gtf 

python /home/shengq1/pylibs/bin/dexseq_prepare_annotation.py smad4.gtf smad4.gff 