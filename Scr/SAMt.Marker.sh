#!/bin/bash

Ref1=$1
Ref2=$2
BaseDir=$3

Nthreads=8
BamDir=${BaseDir}"/MapDdp"
DataExt=".rmd.bam"

RefPair=$Ref1"_"$Ref2
PosDir=${BaseDir}"/MUMmer/$RefPair"
PosPath="$PosDir/$RefPair.intersect.snps"

OutDir=${BaseDir}"/Markers"

# controlla che Ref1 e Ref2 matchino con le prima riga di $PosPath

if [[ ! -d $OutDir ]]; then
  mkdir -p $OutDir
fi

# Format site files
sed '1,4d' $PosPath | awk -v Aref1=$Ref1 'BEGIN{FS="\t"} {print $11"_"Aref1,$1}' > $OutDir/Sites.$Ref1.bed
sed '1,4d' $PosPath | awk -v Aref2=$Ref2 'BEGIN{FS="\t"} {print $12"_"Aref2,$4}' > $OutDir/Sites.$Ref2.bed

cd $BamDir

for Ind1 in $(ls *${DataExt} | grep -v "+")
do
  if (( i % Nthreads == 0 )); then
    wait
  fi
  ((i++))
  (
  RefPath=$(samtools view -H $Ind1 | grep "^@PG" | head -1 | cut -f 5 | cut -d " " -f 6)
  RefName=$(basename $RefPath | cut -d "." -f 1)
  SampleName=$(echo $Ind1 | cut -d "." -f 1)
  bcftools mpileup --regions-file <(awk 'BEGIN {OFS="\t"} {print $1,$2}' $OutDir/Sites.$RefName.bed) -min-MQ5 -a FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP --skip-indels -D -f $RefPath $Ind1 | bcftools call -c -Oz > $OutDir"/"$SampleName"."$RefName".vcf.gz"  
  ) &
done

wait


