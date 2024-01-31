#!/bin/bash

# preparing data for GC content normalization

Ref1=$1
Ref2=$2
MainDir=$3

Ref="$Ref1 $Ref2"
MyChr="chrI chrII chrIII chrIV chrIX chrV chrVI chrVII chrVIII chrX chrXI chrXII chrXIII chrXIV chrXV chrXVI"
WorkDir=$MainDir"/Ref/Mod"

cd $WorkDir

for Ind1 in $Ref
do
  PathRefDir=$MainDir"/CNV/GCdata/"$Ind1
  if [[ ! -d $$PathRefDir"/ChrRef" ]]; then
    mkdir -p $PathRefDir"/ChrRef"
  fi
  # prepara il file LenChr.txt
  grep -v "chrMT" $Ind1.genome.chrref.fa.fai | awk 'BEGIN {FS="\t"} {OFS="\t"} {print NR, $1, $2}' > $PathRefDir"/LenChr.txt"

  for Ind2 in $MyChr
  do
    # prepara chrFiles
    samtools faidx $Ind1".genome.chrref.fa" $Ind2"_"$Ind1 > $PathRefDir"/ChrRef/"$Ind2"_"$Ind1".fa"
  done
done


