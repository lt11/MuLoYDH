#!/bin/bash

Ref1Label=$1
Ref2Label=$2
BaseDir=$3

Nthreads=8
BamDir=${BaseDir}"/MapDdp"
DataExt=".rmd.bam"
PathReg=${BaseDir}"/VariantCalls/HomoHeteroReg"

cd $BamDir

# chiamate da competitive nei segmenti in eterozigosi (corretti togliendo le regioni a bassa MAPQ), $PathReg"/HeteroCorr."$SampleName".bed"
OutDir=${BaseDir}"/VariantCalls/Ploidy1"
OutFreeB=$OutDir"/FreeB"
if [[ ! -d $OutFreeB ]]; then
  mkdir -p $OutFreeB
else 
  rm -rf $OutFreeB
  mkdir -p $OutFreeB
fi

for Ind1 in $(ls *${DataExt} | grep "+") 
do
  if (( i % Nthreads == 0 )); then
    wait
  fi
  ((i++))
  (
  RefPath=$(samtools view -H $Ind1 | grep "^@PG" | head -1 | cut -f 5 | cut -d " " -f 6)
  RefName=$(basename $RefPath | cut -d"." -f 1)
  SampleName=$(echo $Ind1 | cut -d "." -f 1)
  # se non ci sono regioni nel bed FreeBayes si pianta
  if [[ -s  $PathReg"/HeteroCorr."$SampleName".bed" ]]; then
    freebayes -p 1 -t $PathReg"/HeteroCorr."$SampleName".bed" --fasta-reference $RefPath $Ind1 > $OutFreeB"/"$SampleName"."$RefName".vcf"
    # fix sample name
    sed "s|unknown|$Ind1|" $OutFreeB"/"$SampleName"."$RefName".vcf" > $OutFreeB"/"$SampleName"."$RefName".temp.vcf"
    mv -f $OutFreeB"/"$SampleName"."$RefName".temp.vcf" $OutFreeB"/"$SampleName"."$RefName".vcf"
    bgzip $OutFreeB"/"$SampleName"."$RefName".vcf"
    tabix -p vcf $OutFreeB"/"$SampleName"."$RefName".vcf.gz"
  fi
  ) &
done

# chiamate da competitive nei segmenti in omozigosi (corretti togliendo le regioni a bassa MAPQ), $PathReg"/HomoCorr."$SampleName".bed"
OutDir=${BaseDir}"/VariantCalls/Ploidy2"
OutFreeB=$OutDir"/FreeB"
if [[ ! -d $OutFreeB ]]; then
  mkdir -p $OutFreeB
else 
  rm -rf $OutFreeB
  mkdir -p $OutFreeB
fi

for Ind1 in $(ls *${DataExt} | grep "+") 
do
  if (( i % Nthreads == 0 )); then
    wait
  fi
  ((i++))
  (
  RefPath=$(samtools view -H $Ind1 | grep "^@PG" | head -1 | cut -f 5 | cut -d " " -f 6)
  RefName=$(basename $RefPath | cut -d"." -f 1)
  SampleName=$(echo $Ind1 | cut -d "." -f 1)
  # se non ci sono regioni nel bed FreeBayes si pianta
  if [[ -s $PathReg"/HomoCorr."$SampleName".bed" ]]; then
    freebayes -p 2 -t $PathReg"/HomoCorr."$SampleName".bed" --fasta-reference $RefPath $Ind1 > $OutFreeB"/"$SampleName"."$RefName".vcf"
    # fix sample name
    sed "s|unknown|$Ind1|" $OutFreeB"/"$SampleName"."$RefName".vcf" > $OutFreeB"/"$SampleName"."$RefName".temp.vcf"
    mv -f $OutFreeB"/"$SampleName"."$RefName".temp.vcf" $OutFreeB"/"$SampleName"."$RefName".vcf"
    bgzip $OutFreeB"/"$SampleName"."$RefName".vcf"
    tabix -p vcf $OutFreeB"/"$SampleName"."$RefName".vcf.gz"
  fi
  ) &
done

# chiamate unphased
OutDir=${BaseDir}"/VariantCalls/Unphased"
OutFreeB=$OutDir"/FreeB"
if [[ ! -d $OutFreeB ]]; then
  mkdir -p $OutFreeB
else 
  rm -rf $OutFreeB
  mkdir -p $OutFreeB
fi

# occhio, qui l'indice corre sugli ID dei campioni, non sui file bam; infatti SampleName=$IndS
for IndS in $(ls *${DataExt} | grep "+" | cut -d"." -f1) 
do
  if (( i % Nthreads == 0 )); then
    wait
  fi
  ((i++))
  (
  # chiamate unphased sul Ref1
  BamName=$IndS.$Ref1Label.srt.rmd.bam
  RefPath=$(samtools view -H $BamName | grep "^@PG" | head -1 | cut -f 5 | cut -d " " -f 6)
  RefName=$(basename $RefPath | cut -d"." -f 1)
  SampleName=$IndS
  # se non ci sono regioni nel bed FreeBayes si pianta
  if [[ -s $PathReg"/LowMAPQ."$SampleName".bed" ]]; then
    freebayes -p 2 -t <(grep "_$Ref1Label" $PathReg"/LowMAPQ."$SampleName".bed") --fasta-reference $RefPath $BamName > $OutFreeB"/"$SampleName"."$RefName".vcf"
    # fix sample name
    sed "s|unknown|$BamName|" $OutFreeB"/"$SampleName"."$RefName".vcf" > $OutFreeB"/"$SampleName"."$RefName".temp.vcf"
    mv -f $OutFreeB"/"$SampleName"."$RefName".temp.vcf" $OutFreeB"/"$SampleName"."$RefName".vcf"
    bgzip $OutFreeB"/"$SampleName"."$RefName".vcf"
    tabix -p vcf $OutFreeB"/"$SampleName"."$RefName".vcf.gz"
  fi
  ) &
done

# occhio, qui l'indice corre sugli ID dei campioni, non sui file bam; infatti SampleName=$IndS
for IndS in $(ls *${DataExt} | grep "+" | cut -d"." -f1) 
do
  if (( i % Nthreads == 0 )); then
    wait
  fi
  ((i++))
  (
  # chiamate unphased sul Ref2
  BamName=$IndS.$Ref2Label.srt.rmd.bam
  RefPath=$(samtools view -H $BamName | grep "^@PG" | head -1 | cut -f 5 | cut -d " " -f 6)
  RefName=$(basename $RefPath | cut -d"." -f 1)
  SampleName=$IndS
  # se non ci sono regioni nel bed FreeBayes si pianta
  if [[ -s $PathReg"/LowMAPQ."$SampleName".bed" ]]; then
    freebayes -p 2 -t <(grep "_$Ref2Label" $PathReg"/LowMAPQ."$SampleName".bed") --fasta-reference $RefPath $BamName > $OutFreeB"/"$SampleName"."$RefName".vcf"
    # fix sample name
    sed "s|unknown|$BamName|" $OutFreeB"/"$SampleName"."$RefName".vcf" > $OutFreeB"/"$SampleName"."$RefName".temp.vcf"
    mv -f $OutFreeB"/"$SampleName"."$RefName".temp.vcf" $OutFreeB"/"$SampleName"."$RefName".vcf"
    bgzip $OutFreeB"/"$SampleName"."$RefName".vcf"
    tabix -p vcf $OutFreeB"/"$SampleName"."$RefName".vcf.gz"
  fi
  ) &
done

wait


