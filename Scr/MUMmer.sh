#!/bin/bash

Ref1Label=$1
Ref2Label=$2
GenosStruct=$3
BaseDir=$4

RefDir=${BaseDir}"/Ref"

# reference 1 paths
RefPair12=$Ref1Label"_"$Ref2Label
OutDir12=${BaseDir}"/MUMmer/$RefPair12"
# reference 2 paths
RefPair21=$Ref2Label"_"$Ref1Label
OutDir21=${BaseDir}"/MUMmer/$RefPair21"

if [[ $GenosStruct == "collinear" ]]; then
  # estrae i cromosomi SOLO dal reference 1
  MyChromosomes=$(grep "^>" $RefDir/$Ref1Label.genome.fa | sed 's|.||' )
  
  # reference 1 calculations
  
  if [[ ! -d $OutDir12 ]]; then
    mkdir -p $OutDir12
  fi
  
  cd $OutDir12
  
  if [[ -e $RefPair12.prt.all ]]; then
      rm -f $RefPair12.prt.all
  fi
  
  Flag=0
  for IndC in $MyChromosomes
  do
    samtools faidx $RefDir/$Ref1Label.genome.fa $IndC > $RefDir/$Ref1Label.genome.$IndC.fa
    samtools faidx $RefDir/$Ref2Label.genome.fa $IndC > $RefDir/$Ref2Label.genome.$IndC.fa
    nucmer --prefix=$RefPair12.$IndC $RefDir/$Ref1Label.genome.$IndC.fa $RefDir/$Ref2Label.genome.$IndC.fa
    if [[ $Flag == 0 ]]; then
      show-snps -ClrT $RefPair12.$IndC.delta >> $RefPair12.prt.all
      Flag=1
    else
      show-snps -ClrT $RefPair12.$IndC.delta | sed '1,4d' >> $RefPair12.prt.all
    fi
  done
  
  # separate SNM from indels
  # questi marker (RefPair12.prt.all) li uso per filtrare le varianti simulate
  head -4 $RefPair12.prt.all > $RefPair12.prt.snps
  sed '1,4d' $RefPair12.prt.all | awk 'BEGIN{FS="\t"}{if ($2 != "." && $3 != ".") print $0}' >> $RefPair12.prt.snps
  head -4 $RefPair12.prt.all > $RefPair12.prt.indel
  sed '1,4d' $RefPair12.prt.all | awk 'BEGIN{FS="\t"}{if ($2 == "." || $3 == ".") print $0}' >> $RefPair12.prt.indel
  
  # make header file for Marker.Intersect.R (used for appending results)
  head -4 $RefPair12.prt.all > $RefPair12.intersect.snps
  sed 's|NUCMER|NUCMER intersection|' $RefPair12.intersect.snps > temp.snps
  mv temp.snps $RefPair12.intersect.snps
  
  # reference 2 calculations
  
  if [[ ! -d $OutDir21 ]]; then
    mkdir -p $OutDir21
  fi
  
  cd $OutDir21
  
  if [[ -e $RefPair21.prt.all ]]; then
    rm -f $RefPair21.prt.all
  fi
  
  Flag=0
  for IndC in $MyChromosomes
  do
    nucmer --prefix=$RefPair21.$IndC $RefDir/$Ref2Label.genome.$IndC.fa $RefDir/$Ref1Label.genome.$IndC.fa
    if [[ $Flag == 0 ]]; then
      show-snps -ClrT $RefPair21.$IndC.delta >> $RefPair21.prt.all
      Flag=1
    else
      show-snps -ClrT $RefPair21.$IndC.delta | sed '1,4d' >> $RefPair21.prt.all
    fi
  done
  
  # separate SNM from indels
  # queste indels non servono a una sega
  head -4 $RefPair21.prt.all > $RefPair21.prt.snps
  sed '1,4d' $RefPair21.prt.all | awk 'BEGIN{FS="\t"}{if ($2 != "." && $3 != ".") print $0}' >> $RefPair21.prt.snps
  head -4 $RefPair21.prt.all > $RefPair21.prt.indel
  sed '1,4d' $RefPair21.prt.all | awk 'BEGIN{FS="\t"}{if ($2 == "." || $3 == ".") print $0}' >> $RefPair21.prt.indel

elif [[ $GenosStruct == "rearranged" ]]; then
  # reference 1 calculations
  if [[ ! -d $OutDir12 ]]; then
    mkdir -p $OutDir12
  fi
  
  cd $OutDir12
  nucmer --prefix=$RefPair12 $RefDir/$Ref1Label.genome.fa $RefDir/$Ref2Label.genome.fa
  show-snps -ClrT $RefPair12.delta > $RefPair12.prt.all
  
  # separate SNM from indels
  # questi marker (RefPair12.prt.all) li uso per filtrare le varianti simulate
  head -4 $RefPair12.prt.all > $RefPair12.prt.snps
  sed '1,4d' $RefPair12.prt.all | awk 'BEGIN{FS="\t"}{if ($2 != "." && $3 != ".") print $0}' >> $RefPair12.prt.snps
  head -4 $RefPair12.prt.all > $RefPair12.prt.indel
  sed '1,4d' $RefPair12.prt.all | awk 'BEGIN{FS="\t"}{if ($2 == "." || $3 == ".") print $0}' >> $RefPair12.prt.indel
  
  # make header file for Marker.Intersect.R (used for appending results)
  head -4 $RefPair12.prt.all > $RefPair12.intersect.snps
  sed 's|NUCMER|NUCMER intersection|' $RefPair12.intersect.snps > temp.snps
  mv temp.snps $RefPair12.intersect.snps
  
  # reference 2 calculations
  if [[ ! -d $OutDir21 ]]; then
    mkdir -p $OutDir21
  fi
  
  cd $OutDir21
  nucmer --prefix=$RefPair21 $RefDir/$Ref2Label.genome.fa $RefDir/$Ref1Label.genome.fa
  show-snps -ClrT $RefPair21.delta > $RefPair21.prt.all
  
  # separate SNM from indels
  # queste indels non servono a una sega
  head -4 $RefPair21.prt.all > $RefPair21.prt.snps
  sed '1,4d' $RefPair21.prt.all | awk 'BEGIN{FS="\t"}{if ($2 != "." && $3 != ".") print $0}' >> $RefPair21.prt.snps
  head -4 $RefPair21.prt.all > $RefPair21.prt.indel
  sed '1,4d' $RefPair21.prt.all | awk 'BEGIN{FS="\t"}{if ($2 == "." || $3 == ".") print $0}' >> $RefPair21.prt.indel
fi

cd $BaseDir"/Scr"
Rscript Marker.Intersect.R $Ref1Label $Ref2Label $BaseDir


