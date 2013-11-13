l *.log| awk '{if($7 == "5"){print $NF}}' |xargs rm
