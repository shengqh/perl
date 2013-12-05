cd /data/cqs/shengq1/reference/miRBase20
awk 'BEGIN{name="undefined"} {if($1 ~ "^>"){name=substr($1,2)}else{print name","$1}} ' hsa.mature.dna.fa > hsa.mature.dna.db
awk 'BEGIN{name="undefined"} {if($1 ~ "^>"){name=substr($1,2)}else{print name","$1}} ' mmu.mature.dna.fa > mmu.mature.dna.db
awk 'BEGIN{name="undefined"} {if($1 ~ "^>"){name=substr($1,2)}else{print name","$1}} ' rno.mature.dna.fa > rno.mature.dna.db

