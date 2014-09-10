for i in $(ls -d */); do
        for j in $(ls ${i}*.tmp); do
                grep "ref:" ${j} >> ${i}/${i%%/}.gen
                grep "snp:" ${j} >> ${i}/${i%%/}.gen
        done
done
