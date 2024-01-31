#!/bin/bash

Ref1Label=$1
Ref2Label=$2
BaseDir=$3

Nthreads=8
BamDir=${BaseDir}"/MapDdp"
DataExt=".rmd.bam"
PathReg=${BaseDir}"/VariantCalls/HomoHeteroReg"

cd $BamDir

# segmenti in eterozigosi (corretti togliendo le regioni a bassa MAPQ), $PathReg"/HeteroCorr."$SampleName".bed"
OutDir=${BaseDir}"/VariantCalls/Ploidy1"
OutSAMt=$OutDir"/SAMt"
if [[ ! -d $OutSAMt ]]; then
  mkdir -p $OutSAMt
else
  rm -rf $OutSAMt
  mkdir -p $OutSAMt
fi

for Ind1 in $(ls *${DataExt} | grep "+") 
do
  if (( i % Nthreads == 0 )); then
    wait
  fi
  ((i++))
  (
  RefPath=$(samtools view -H $Ind1 | grep "^@PG" | head -1 | cut -f 5 | cut -d " " -f 6)
  RefName=$(basename $RefPath | cut -d "." -f 1)
  SampleName=$(echo $Ind1 | cut -d "." -f 1)
  # se non ci sono regioni nel bed non fare le chiamate
  if [[ -s $PathReg"/HeteroCorr."$SampleName".bed" ]]; then
    bcftools mpileup --regions-file $PathReg"/HeteroCorr."$SampleName".bed" -min-MQ5 -a FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP -D -f $RefPath $Ind1 | bcftools call --ploidy 1 -mv -Oz > $OutSAMt"/"$SampleName"."$RefName".vcf.gz"
    tabix -p vcf $OutSAMt"/"$SampleName"."$RefName".vcf.gz"
  fi
  ) &
done

# segmenti in omozigosi (corretti togliendo le regioni a bassa MAPQ), $PathReg"/HomoCorr."$SampleName".bed"
OutDir=${BaseDir}"/VariantCalls/Ploidy2"
OutSAMt=$OutDir"/SAMt"
if [[ ! -d $OutSAMt ]]; then
  mkdir -p $OutSAMt
else
  rm -rf $OutSAMt
  mkdir -p $OutSAMt
fi

for Ind1 in $(ls *${DataExt} | grep "+") 
do
  if (( i % Nthreads == 0 )); then
    wait
  fi
  ((i++))
  (
  RefPath=$(samtools view -H $Ind1 | grep "^@PG" | head -1 | cut -f 5 | cut -d " " -f 6)
  RefName=$(basename $RefPath | cut -d "." -f 1)
  SampleName=$(echo $Ind1 | cut -d "." -f 1)
  # se non ci sono regioni nel bed non fare le chiamate
  if [[ -s $PathReg"/HomoCorr."$SampleName".bed" ]]; then   
    bcftools mpileup --regions-file $PathReg"/HomoCorr."$SampleName".bed" -min-MQ5 -a FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP -D -f $RefPath $Ind1 | bcftools call -mv -Oz > $OutSAMt"/"$SampleName"."$RefName".vcf.gz"
    tabix -p vcf $OutSAMt"/"$SampleName"."$RefName".vcf.gz"
  fi
  ) &
done

# chiamate unphased
OutDir=${BaseDir}"/VariantCalls/Unphased"
OutSAMt=$OutDir"/SAMt"
if [[ ! -d $OutSAMt ]]; then
  mkdir -p $OutSAMt
else
  rm -rf $OutSAMt
  mkdir -p $OutSAMt
fi

# occhio, qui l'indice corre sugli ID dei campioni, non sui file bam; infatti SampleName=$IndS
for IndS in $(ls *${DataExt} | grep "+" | cut -d "." -f1) 
do
  if (( i % Nthreads == 0 )); then
    wait
  fi
  ((i++))
  (
  # chiamate unphased sul Ref1
  BamName=$IndS.$Ref1Label.srt.rmd.bam
  RefPath=$(samtools view -H $BamName | grep "^@PG" | head -1 | cut -f 5 | cut -d " " -f 6)
  RefName=$(basename $RefPath | cut -d "." -f 1)
  SampleName=$IndS
  # se non ci sono regioni nel bed non fare le chiamate
  if [[ -s $PathReg"/LowMAPQ."$SampleName".bed" ]]; then
    grep "_$Ref1Label" $PathReg"/LowMAPQ."$SampleName".bed" > $PathReg"/LowMAPQ."$SampleName"."$Ref1Label".bed"
    bcftools mpileup --regions-file $PathReg"/LowMAPQ."$SampleName"."$Ref1Label".bed" -min-MQ5 -a FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP -D -f $RefPath $BamName | bcftools call -mv -Oz > $OutSAMt"/"$SampleName"."$RefName".vcf.gz"
    tabix -p vcf $OutSAMt"/"$SampleName"."$RefName".vcf.gz"
  fi
  ) &
done

# occhio, qui l'indice corre sugli ID dei campioni, non sui file bam; infatti SampleName=$IndS
for IndS in $(ls *${DataExt} | grep "+" | cut -d "." -f1) 
do
  if (( i % Nthreads == 0 )); then
    wait
  fi
  ((i++))
  (
  # chiamate unphased sul Ref2
  BamName=$IndS.$Ref2Label.srt.rmd.bam
  RefPath=$(samtools view -H $BamName | grep "^@PG" | head -1 | cut -f 5 | cut -d " " -f 6)
  RefName=$(basename $RefPath | cut -d "." -f 1)
  SampleName=$IndS
  # se non ci sono regioni nel bed non fare le chiamate
  if [[ -s $PathReg"/LowMAPQ."$SampleName".bed" ]]; then
    grep "_$Ref2Label" $PathReg"/LowMAPQ."$SampleName".bed" > $PathReg"/LowMAPQ."$SampleName"."$Ref2Label".bed"
    bcftools mpileup --regions-file $PathReg"/LowMAPQ."$SampleName"."$Ref2Label".bed" -min-MQ5 -a FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP -D -f $RefPath $BamName | bcftools call -mv -Oz > $OutSAMt"/"$SampleName"."$RefName".vcf.gz"
    tabix -p vcf $OutSAMt"/"$SampleName"."$RefName".vcf.gz"
  fi
  ) &
done

wait


