#!/bin/bash

Ref1Label=$1
Ref2Label=$2
ParLabel=$3
BaseDir=$4

Nthreads=4

# segmenti in eterozigosi
OutDir=$BaseDir"/VariantCalls/Ploidy1/ParSubtracted"
InputDir=$BaseDir"/VariantCalls/Ploidy1/Intersect"

# per evitare che si pianti nei re-run
rm -f $InputDir/*tbi

if [[ ! -d $OutDir ]]; then
  mkdir -p $OutDir
else
  rm -rf $OutDir
  mkdir -p $OutDir
fi

if [[ $ParLabel == "-" ]]; then
  cp $InputDir/*.vcf.gz $OutDir
else
  cd $InputDir
  # espandi il file Par, poi leggi i vcf non espansi
  gunzip $ParLabel.$Ref1Label+$Ref2Label.vcf.gz
  # leggi i vcf non espansi
  for File1 in $(ls *vcf.gz)
  do
    if (( i % Nthreads == 0 )); then
      wait
    fi
    ((i++))
    (
    File1ID=$(echo $File1 | cut -d "." -f 1)
    sed "s|$ParLabel|$File1ID|g" < $ParLabel.$Ref1Label+$Ref2Label.vcf > temp.4.$File1ID.vcf
    bgzip temp.4.$File1ID.vcf
    tabix temp.4.$File1ID.vcf.gz
    tabix -f $File1
    vcf-isec -c $File1 temp.4.$File1ID.vcf.gz | bgzip > $OutDir"/"$File1
    rm -f temp.4.$File1ID*
    ) &
  done
  
  wait
  bgzip $ParLabel.$Ref1Label+$Ref2Label.vcf
  cp $ParLabel.$Ref1Label+$Ref2Label.vcf.gz $OutDir
fi

# segmenti in omozigosi
OutDir=$BaseDir"/VariantCalls/Ploidy2/ParSubtracted"
InputDir=$BaseDir"/VariantCalls/Ploidy2/Intersect"

# per evitare che si pianti nei re-run
rm -f $InputDir/*tbi

if [[ ! -d $OutDir ]]; then
  mkdir -p $OutDir
else
  rm -rf $OutDir
  mkdir $OutDir
fi

if [[ $ParLabel == "-" ]]; then
  cp $InputDir/*.vcf.gz $OutDir
else
  cd $InputDir
  # espandi il file Par, poi leggi i vcf non espansi
  gunzip $ParLabel.$Ref1Label+$Ref2Label.vcf.gz
  # leggi i vcf non espansi
  for File1 in $(ls *vcf.gz)
  do
    if (( i % Nthreads == 0 )); then
      wait
    fi
    ((i++))
    (
    File1ID=$(echo $File1 | cut -d "." -f 1)
    sed "s|$ParLabel|$File1ID|g" < $ParLabel.$Ref1Label+$Ref2Label.vcf > temp.4.$File1ID.vcf
    bgzip temp.4.$File1ID.vcf
    tabix temp.4.$File1ID.vcf.gz
    tabix -f $File1
    vcf-isec -c $File1 temp.4.$File1ID.vcf.gz | bgzip > $OutDir"/"$File1
    rm -f temp.4.$File1ID*
    ) &
  done
  
  wait
  bgzip $ParLabel.$Ref1Label+$Ref2Label.vcf
  cp $ParLabel.$Ref1Label+$Ref2Label.vcf.gz $OutDir
fi

# varianti unphased
OutDir=$BaseDir"/VariantCalls/Unphased/ParSubtracted"
InputDir=$BaseDir"/VariantCalls/Unphased/Intersect"
BothRef=$Ref1Label" "$Ref2Label

# per evitare che si pianti nei re-run
rm -f $InputDir/*tbi

if [[ ! -d $OutDir ]]; then
  mkdir -p $OutDir
else
  rm -rf $OutDir
  mkdir $OutDir
fi

for IndR in $BothRef
do
  if [[ $ParLabel == "-" ]]; then
    cp $InputDir/*.vcf.gz $OutDir
  else
    cd $InputDir
    # espandi il file Par, poi leggi i vcf non espansi
    gunzip $ParLabel.$IndR.vcf.gz
    # leggi i vcf non espansi
    for File1 in $(ls *vcf.gz | grep $IndR.vcf.gz)
    do
      if (( i % Nthreads == 0 )); then
        wait
      fi
      ((i++))
      (
      File1ID=$(echo $File1 | cut -d "." -f 1)
      sed "s|$ParLabel|$File1ID|g" < $ParLabel.$IndR.vcf > temp.4.$File1ID.vcf
      bgzip temp.4.$File1ID.vcf
      tabix -f temp.4.$File1ID.vcf.gz
      tabix -f $File1
      vcf-isec -c $File1 temp.4.$File1ID.vcf.gz | bgzip > $OutDir"/"$File1
      rm -f temp.4.$File1ID*
      ) &
    done
    
    wait
    bgzip $ParLabel.$IndR.vcf
    cp $ParLabel.$IndR.vcf.gz $OutDir
  fi
done



