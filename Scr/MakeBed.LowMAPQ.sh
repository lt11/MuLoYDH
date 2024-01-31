#!/bin/bash

# interseca le regioni a bassa MAPQ con quelle calcolate da MakeBed.HomoHeteroReg.R
# le regioni prodotte verranno chiamate da mapping standard
# poi, corregge le regioni calcolate da MakeBed.HomoHeteroReg.R
# sottraendo quelle a bassa MAPQ

Ref1Label=$1
Ref2Label=$2
BaseDir=$3

WorkDir=$BaseDir"/VariantCalls/HomoHeteroReg"

cd $WorkDir

for IndS in $(ls Hetero*bed | cut -d"." -f2)
do
  # calcolo regioni da chiamare standard (occhio, il file LowMAPQ.$IndS.bed contiene le regioni di ENTRAMBI i reference)
  # calcolo regioni non in LOH che chiamiamo randomly sul reference 1
  bedtools intersect -a <(grep "_$Ref1Label" Hetero.$IndS.bed) -b RegionsLowMAPQ.bed > LowMAPQ.$IndS.bed
  # faccio un evento fake sul reference 2 per evitare che FreeBayes si pianti (cosa che accade se non ci sono regioni)
  echo -e "chrMT_$Ref2Label\t0\t1" >> LowMAPQ.$IndS.bed
  # calcolo regioni da chiamare standard in LOH (che chiamo sul reference in LOH)
  bedtools intersect -a Homo.$IndS.bed -b RegionsLowMAPQ.bed >> LowMAPQ.$IndS.bed
  # correzione regioni hetero/homo da chiamare competitive
  bedtools subtract -a Hetero.$IndS.bed -b RegionsLowMAPQ.bed > HeteroCorr.$IndS.bed
  bedtools subtract -a Homo.$IndS.bed -b RegionsLowMAPQ.bed > HomoCorr.$IndS.bed
  
  sort -k1,1 -k2,2n LowMAPQ.$IndS.bed > temp
  mv -f temp LowMAPQ.$IndS.bed
done



