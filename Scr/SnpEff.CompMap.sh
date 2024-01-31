#!/bin/bash

################################################
# annota le varianti nella cartella 
# "Region" con SnpEff
# poi fa lo stesso con "Intersect" ma 
# l'output è una sottocartella
# (meno visibile)

# occhio: crossref!
################################################

Ref1Label=$1
Ref2Label=$2
Ref1=$3
Ref2=$4
SnpEff=$5
BaseDir=$6

Ref1Version="1"
Ref2Version="1"
SnpEffConfig=$(echo $SnpEff | sed 's|...$|config|')

# segmenti in eterozigosi e omozigosi, varianti filtrate per regione e quality

WorkDir=$BaseDir"/VariantCalls/Merged/Regions"
OutDir=$BaseDir"/VariantCalls/Merged/Effect"

if [[ ! -d $OutDir ]]; then
  mkdir -p $OutDir
fi

cd $WorkDir
# occhio, questi vcf NON sono compressi
for IndSamp in $(ls *.$Ref1.vcf | grep -v "Ann")
do
  SampName=$(echo $IndSamp | cut -d"." -f1)
  # ricerco nel file config il nome del genoma e la versione da dare in pasto a SnpEff
  Ref1Name=$(grep $Ref1Label $SnpEffConfig | grep -v ^# | cut -d" " -f1 | grep "\.$Ref1Version\.genome$" | sed 's|\.genome||')
  java -Xmx4g -jar $SnpEff -v $Ref1Name -s $SampName.$Ref1.Summary.html $IndSamp > $SampName.$Ref1.Ann.vcf
  mv $SampName.$Ref1.Ann.vcf $SampName.$Ref1.Summary* $OutDir
done &
# occhio, questi vcf NON sono compressi
for IndSamp in $(ls *.$Ref2.vcf | grep -v "Ann")
do
  SampName=$(echo $IndSamp | cut -d"." -f1)
  # ricerco nel file config il nome del genoma e la versione da dare in pasto a SnpEff
  Ref2Name=$(grep $Ref2Label $SnpEffConfig | grep -v ^# | cut -d" " -f1 | grep "\.$Ref2Version\.genome$" | sed 's|\.genome||')
  java -Xmx4g -jar $SnpEff -v $Ref2Name -s $SampName.$Ref2.Summary.html $IndSamp > $SampName.$Ref2.Ann.vcf
  mv $SampName.$Ref2.Ann.vcf $SampName.$Ref2.Summary* $OutDir
done &

# segmenti in eterozigosi, varianti NON filtrate per regione e quality

WorkDir=$BaseDir"/VariantCalls/Ploidy1/Intersect"
OutDir=$BaseDir"/VariantCalls/Ploidy1/Intersect/Effect"

if [[ ! -d $OutDir ]]; then
  mkdir -p $OutDir
fi

cd $WorkDir
# occhio, questi vcf sono compressi e contengono le varianti su entrambi i reference
for IndSamp in $(ls *.vcf.gz | grep -v "Ann")
do
  SampName=$(echo $IndSamp | cut -d"." -f1)
  # tolgo dal vcf le varianti e le stringhe contig nell'header dell'altro reference col grep -v
  # tolgo _RefNLabel dalle entry del campo CHROM (e dall'header) con sed (sennò il DB di SnpEff si incazza)
  zcat $IndSamp | grep -v "_$Ref2Label" | sed "s|_${Ref1Label}||g" > $SampName.Ref1.vcf.tmp
  # ricerco nel file config il nome del genoma e la versione da dare in pasto a SnpEff
  Ref1Name=$(grep $Ref1Label $SnpEffConfig | grep -v ^# | cut -d" " -f1 | grep "\.$Ref1Version\.genome$" | sed 's|\.genome||')
  java -Xmx4g -jar $SnpEff -v $Ref1Name -s $SampName.$Ref1.Summary.html $SampName.Ref1.vcf.tmp > $SampName.$Ref1.Ann.vcf
  mv $SampName.$Ref1.Ann.vcf $SampName.$Ref1.Summary* $OutDir
  rm $SampName.Ref1.vcf.tmp
done &
# occhio, questi vcf sono compressi e contengono le varianti su entrambi i reference
for IndSamp in $(ls *.vcf.gz | grep -v "Ann")
do
  SampName=$(echo $IndSamp | cut -d"." -f1)
  # tolgo dal vcf le varianti e le stringhe contig nell'header dell'altro reference col grep -v
  # tolgo _RefNLabel dalle entry del campo CHROM (e dall'header) con sed (sennò il DB di SnpEff si incazza)
  zcat $IndSamp | grep -v "_$Ref1Label" | sed "s|_${Ref2Label}||g" > $SampName.Ref2.vcf.tmp
  # ricerco nel file config il nome del genoma e la versione da dare in pasto a SnpEff
  Ref2Name=$(grep $Ref2Label $SnpEffConfig | grep -v ^# | cut -d" " -f1 | grep "\.$Ref2Version\.genome$" | sed 's|\.genome||')
  java -Xmx4g -jar $SnpEff -v $Ref2Name -s $SampName.$Ref2.Summary.html $SampName.Ref2.vcf.tmp > $SampName.$Ref2.Ann.vcf
  mv $SampName.$Ref2.Ann.vcf $SampName.$Ref2.Summary* $OutDir
  rm $SampName.Ref2.vcf.tmp
done &

# segmenti in omozigosi, varianti NON filtrate per regione e quality

WorkDir=$BaseDir"/VariantCalls/Ploidy2/Intersect"
OutDir=$BaseDir"/VariantCalls/Ploidy2/Intersect/Effect"

if [[ ! -d $OutDir ]]; then
  mkdir -p $OutDir
fi

cd $WorkDir
# occhio, questi vcf sono compressi e contengono le varianti su entrambi i reference
for IndSamp in $(ls *.vcf.gz | grep -v "Ann")
do
  SampName=$(echo $IndSamp | cut -d"." -f1)
  # tolgo dal vcf le varianti e le stringhe contig nell'header dell'altro reference col grep -v
  # tolgo _RefNLabel dalle entry del campo CHROM (e dall'header) con sed (sennò il DB di SnpEff si incazza)
  zcat $IndSamp | grep -v "_$Ref2Label" | sed "s|_${Ref1Label}||g" > $SampName.Ref1.vcf.tmp
  # ricerco nel file config il nome del genoma e la versione da dare in pasto a SnpEff
  Ref1Name=$(grep $Ref1Label $SnpEffConfig | grep -v ^# | cut -d" " -f1 | grep "\.$Ref1Version\.genome$" | sed 's|\.genome||')
  java -Xmx4g -jar $SnpEff -v $Ref1Name -s $SampName.$Ref1.Summary.html $SampName.Ref1.vcf.tmp > $SampName.$Ref1.Ann.vcf
  mv $SampName.$Ref1.Ann.vcf $SampName.$Ref1.Summary* $OutDir
  rm $SampName.Ref1.vcf.tmp
done &
# occhio, questi vcf sono compressi e contengono le varianti su entrambi i reference
for IndSamp in $(ls *.vcf.gz | grep -v "Ann")
do
  SampName=$(echo $IndSamp | cut -d"." -f1)
  # tolgo dal vcf le varianti e le stringhe contig nell'header dell'altro reference col grep -v
  # tolgo _RefNLabel dalle entry del campo CHROM (e dall'header) con sed (sennò il DB di SnpEff si incazza)
  zcat $IndSamp | grep -v "_$Ref1Label" | sed "s|_${Ref2Label}||g" > $SampName.Ref2.vcf.tmp
  # ricerco nel file config il nome del genoma e la versione da dare in pasto a SnpEff
  Ref2Name=$(grep $Ref2Label $SnpEffConfig | grep -v ^# | cut -d" " -f1 | grep "\.$Ref2Version\.genome$" | sed 's|\.genome||')
  java -Xmx4g -jar $SnpEff -v $Ref2Name -s $SampName.$Ref2.Summary.html $SampName.Ref2.vcf.tmp > $SampName.$Ref2.Ann.vcf
  mv $SampName.$Ref2.Ann.vcf $SampName.$Ref2.Summary* $OutDir
  rm $SampName.Ref2.vcf.tmp
done &

wait


