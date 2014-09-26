#!/bin/bash

echo "Usage: $0 xml_directory fasta_file ";

if [[ $# -eq 0 || $# -gt 2 ]]
then
  echo "No/wrong ($#) arguments detected "
  exit 1 #exit shell script
fi

for i in `ls -d1 $1/*.pep.xml`
do
  echo $i
  if grep -v -Fxq "alternative_protein" $i
  then
    RefreshParser $i $2
  fi
done
