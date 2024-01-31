#!/bin/bash

BaseDir=$1

DataExt=".srt.bam"
DataDir=${BaseDir}"/MapRaw"
OutDir=${DataDir}"/FastQC"
Nthreads=4

if [[ ! -d $OutDir ]]; then
  mkdir -p $OutDir
fi

cd $DataDir

for Ind1 in *${DataExt}
do
  if (( i % Nthreads == 0 )); then
    wait
  fi
  ((i++))
  fastqc $Ind1 &> /dev/null &
done

wait

mv *_fastqc*  $OutDir

DataExt=".rmd.bam"
DataDir=${BaseDir}"/MapDdp"
OutDir=$DataDir"/FastQC"

if [[ ! -d $OutDir ]]; then
    mkdir -p $OutDir
fi

cd $DataDir

for Ind1 in *${DataExt}
do
  if (( i % Nthreads == 0 )); then
    wait
  fi
  ((i++))
  fastqc $Ind1 &> /dev/null &
done

wait

mv *_fastqc*  $OutDir


