#!/bin/bash

BaseDir=$1

FastqExtension=".fastq.gz"
FastqDir=${BaseDir}"/Exp"
OutDir=${BaseDir}"/Exp/FastQC"
Nthreads=4

if [[ ! -d $OutDir ]]; then
  mkdir -p $OutDir
fi

cd $FastqDir

for Ind1 in *${FastqExtension}
do
  if (( i % Nthreads == 0 )); then
    wait
  fi
  ((i++))
  fastqc $Ind1 &> /dev/null &
done

wait

mv *_fastqc*  $OutDir


