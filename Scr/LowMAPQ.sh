#!/bin/bash

# prende il primo bam controllo
# o parental o uno a caso 
# toglie le read a MAPQ < 5
# calcola il coverage
# estrae le regioni a coverage < 5

HybridZero=$1
ParentalZero=$2
BaseDir=$3

OutDir=$BaseDir"/VariantCalls/HomoHeteroReg"

if [[ ! -d $OutDir ]]; then
  mkdir -p $OutDir
fi

cd $BaseDir/MapDdp

if [[ $HybridZero !=  "-" ]]; then
  H0Label=$(echo $HybridZero | cut -d" " -f1)
  H0Sample=$(ls $H0Label*bam | grep "+")
  samtools view -b -q 5 $H0Sample > Appoggio.bam
  bedtools genomecov -ibam Appoggio.bam -bga | awk 'BEGIN {FS="\t"} {OFS="\t"} {if ($4 < 5) print $1,$2,$3}' > $OutDir"/RegionsLowMAPQ.bed"
elif [[ $ParentalZero !=  "-"  ]]; then
  P0Sample=$(ls $ParentalZero*bam | grep "+")
  samtools view -b -q 5 $P0Sample > Appoggio.bam
  bedtools genomecov -ibam Appoggio.bam -bga | awk 'BEGIN {FS="\t"} {OFS="\t"} {if ($4 < 5) print $1,$2,$3}' > $OutDir"/RegionsLowMAPQ.bed"
else
  Sample=$(ls *bam | grep "+" | head -1)
  samtools view -b -q 5 $Sample > Appoggio.bam
  bedtools genomecov -ibam Appoggio.bam -bga | awk 'BEGIN {FS="\t"} {OFS="\t"} {if ($4 < 5) print $1,$2,$3}' > $OutDir"/RegionsLowMAPQ.bed"
fi

rm -f Appoggio.bam


