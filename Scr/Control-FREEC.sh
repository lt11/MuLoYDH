#!/bin/bash

# works with the config file template in Lib

MainDir=$1

TemplateFile=$MainDir"/Lib/Control-FREEC.Config.Template.txt"
Path2CFlib=$MainDir"/Lib"

Nthreads=4
BamDir=$MainDir"/MapDdp"
DataExt=".rmd.bam"

# retrieve ploidy from template for plotting
Ploidy=$(grep "ploidy =" $TemplateFile | cut -d "=" -f 2 | sed 's|^ ||')

if [[ ! -s $TemplateFile ]]; then
  echo "Missing template file"
  exit 1
fi

cd $BamDir

for IndS in $(ls *${DataExt} | grep -v "+")
do
  if (( i % Nthreads == 0 )); then
    wait
  fi
  ((i++))
  (
  RefPath=$(samtools view -H $IndS | grep "^@PG" | head -1 | cut -f 5 | cut -d " " -f 6)
  Reference=$(basename $RefPath | cut -d "." -f 1)
  SampleName=$(echo $IndS | cut -d "." -f 1)
  
  OutDir=$MainDir"/CNV/Results/"$Reference"/"$SampleName
  if [[ ! -d $OutDir ]]; then
    mkdir -p $OutDir
  fi
  
  Config=$MainDir"/CNV/Config/"$Reference"/"$SampleName
  if [[ ! -d $Config ]]; then
    mkdir -p $Config
  fi
  
  # crea il file config_GC.txt specifico per il campione/reference
  sed "s|STRINGchrLenFile|$MainDir/CNV/GCdata/$Reference/LenChr.txt|" < $TemplateFile | sed "s|STRINGchrFiles|$MainDir/CNV/GCdata/$Reference/ChrRef|" | sed "s|STRINGmateFile|$MainDir/MapDdp/$SampleName.$Reference.srt.rmd.bam|" | sed "s|STRINGoutputDir|$OutDir|" | sed "s|STRINGgemMappabilityFile|$MainDir/CNV/Mappability/$Reference/$Reference.mappability|" > $Config"/CF."$SampleName".Input.GC.txt"
  
  # running Control-FREEC
  freec -conf $Config"/CF."$SampleName".Input.GC.txt" &> $OutDir"/"$SampleName".log"
  
  # calculate significance
  cat $Path2CFlib"/AssessSignificance.R" | R --slave --args $OutDir"/"$SampleName"."$Reference".srt.rmd.bam_CNVs" $OutDir"/"$SampleName"."$Reference".srt.rmd.bam_ratio.txt"
  
  # plot: in primis, trovo i nomi dei cromosomi, poi li passo allo script R
  AllChr=$(cut -f2 $MainDir"/CNV/GCdata/"$Reference"/LenChr.txt" | sed 's|chr||')
  Nchrom=$(echo ${AllChr} | wc -w)
  cat $Path2CFlib"/MakeGraph.R" | R --slave --args $Ploidy $OutDir"/"$SampleName"."$Reference".srt.rmd.bam_ratio.txt" $Nchrom $AllChr
  ) &
done

wait


