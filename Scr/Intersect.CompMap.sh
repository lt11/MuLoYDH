#!/bin/bash

Ref1Label=$1
Ref2Label=$2
BaseDir=$3

Nthreads=8

# segmenti in eterozigosi
CallsDir=${BaseDir}"/VariantCalls/Ploidy1"

OutPath=$CallsDir"/Intersect"
if [[ -d $OutPath ]]; then
  rm -rf $OutPath
fi
mkdir -p $OutPath

# retrieve sample name from FreeB folder
# Ind1 is FreeB vcf
# Ind2 is SAMt vcf
# Ind3 the path to Intersect output
for Ind1 in $(ls ${CallsDir}/FreeB/*vcf.gz) 
do
  if (( i % Nthreads == 0 )); then
    wait
  fi
  ((i++))
  (
  Ind2=$(echo $Ind1 | sed "s/FreeB/SAMt/")
  Ind3=$(echo $Ind1 | sed "s/FreeB/Intersect/")
  vcf-isec -n +2 $Ind2 $Ind1 | bgzip > $Ind3  
  ) &
done

# segmenti in omozigosi
CallsDir=${BaseDir}"/VariantCalls/Ploidy2"

OutPath=$CallsDir"/Intersect"
if [[ -d $OutPath ]]; then
  rm -rf $OutPath
fi
mkdir -p $OutPath

# retrieve sample name from FreeB folder
# Ind1 is FreeB vcf
# Ind2 is SAMt vcf
# Ind3 the path to Intersect output
for Ind1 in $(ls ${CallsDir}/FreeB/*vcf.gz) 
do
  if (( i % Nthreads == 0 )); then
    wait
  fi
  ((i++))
  (
  Ind2=$(echo $Ind1 | sed "s/FreeB/SAMt/")
  Ind3=$(echo $Ind1 | sed "s/FreeB/Intersect/")
  vcf-isec -n +2 $Ind2 $Ind1 | bgzip > $Ind3  
  ) &
done

# varianti unphased
CallsDir=${BaseDir}"/VariantCalls/Unphased"
BothRef=$Ref1Label" "$Ref2Label

OutPath=$CallsDir"/Intersect"
if [[ -d $OutPath ]]; then
  rm -rf $OutPath
fi
mkdir -p $OutPath

# for each reference
# retrieve sample name from FreeB folder
# Ind1 is FreeB vcf
# Ind2 is SAMt vcf
# Ind3 the path to Intersect output
for IndR in $BothRef
do
  for Ind1 in $(ls ${CallsDir}/FreeB/*vcf.gz | grep "$IndR.vcf.gz")
  do
    if (( i % Nthreads == 0 )); then
      wait
    fi
    ((i++))
    (
    Ind2=$(echo $Ind1 | sed "s/FreeB/SAMt/")
    Ind3=$(echo $Ind1 | sed "s/FreeB/Intersect/")
    vcf-isec -n +2 $Ind2 $Ind1 | bgzip > $Ind3  
    ) &
  done
  
  wait
done

wait



