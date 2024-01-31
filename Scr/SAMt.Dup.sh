#!/bin/bash

BaseDir=$1

OutDir="MapDdp"
OutPath=$BaseDir"/"$OutDir
Nthreads=16

cd $BaseDir

if [[ ! -d $OutPath ]]; then
  mkdir -p $OutPath
fi

for Ind1 in MapRaw/*bam
do
  if (( i % Nthreads == 0 )); then
    wait
  fi
  ((i++))
  Pref=$(echo $Ind1 | cut -d "/" -f 2 | cut -d "." -f 1,2)
  (
  samtools rmdup $Ind1 $OutPath/$Pref.srt.rmd.bam
  samtools index $OutPath/$Pref.srt.rmd.bam
  ) &
done

wait


