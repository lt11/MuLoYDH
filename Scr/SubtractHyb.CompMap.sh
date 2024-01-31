#!/bin/bash

Ref1Label=$1
Ref2Label=$2
H0Label=$3
BaseDir=$4

Nthreads=4

# segmenti in eterozigosi
OutDir=$BaseDir"/VariantCalls/Ploidy1/HybSubtracted"
InputDir=$BaseDir"/VariantCalls/Ploidy1/ParSubtracted"

if [[ ! -d $OutDir ]]; then
  mkdir -p $OutDir
else
  rm -rf $OutDir
  mkdir -p $OutDir
fi
cd $InputDir
# per evitare che si pianti nei re-run
rm -f $InputDir/*tbi $InputDir/Merged.txt.gz
if [[ $H0Label == "-" ]]; then
  # index vcf per bcftools merge
  # serve un for perché tabix
  # è scarso e non legge multi input
  for IndS in $(ls | grep -v Merged.txt.gz); do tabix -f $IndS; done
  AllFiles=$(ls *vcf.gz | grep -v ParHyb)
  bcftools merge --output-type v $AllFiles | grep -v "^##" | bgzip > Merged.txt.gz
  Rscript $BaseDir"/Scr/FiltSharedVar.R" "Ploidy1" $BaseDir
  cd $OutDir
  # anche bgzip è molto scarso
  for i in *; do bgzip $i; done
elif [[ $H0Label == "skip" ]]; then
  cp $InputDir/*vcf.gz $OutDir
  break
else
  # apro tutti i gz dei controlli per correggere il nome del campione
  for IndH0 in $H0Label
  do
    gunzip $IndH0.$Ref1Label+$Ref2Label.vcf.gz
  done

  for File1 in $(ls *vcf.gz)
  do
    if (( i % Nthreads == 0 )); then
      wait
    fi
    ((i++))
    (
    File1ID=$(echo $File1 | cut -d "." -f 1)
    # correzione nome campione
    for IndH0 in $H0Label
    do
      sed "s|$IndH0|$File1ID|g" < $IndH0.$Ref1Label+$Ref2Label.vcf > $IndH0.temp.4.$File1ID.vcf
    done
    bgzip *temp.4.$File1ID.vcf
    tabix *temp.4.$File1ID.vcf.gz
    tabix -f $File1
    # sottrazione varianti
    vcf-isec -c $File1 *temp.4.$File1ID.vcf.gz | bgzip > $OutDir"/"$File1
    rm -f *temp.4.$File1ID*
    ) &
  done
  
  wait
  
  for IndH0 in $H0Label
  do
    bgzip $IndH0.$Ref1Label+$Ref2Label.vcf
    cp $IndH0.$Ref1Label+$Ref2Label.vcf.gz $OutDir
  done
fi

# segmenti in omozigosi
OutDir=$BaseDir"/VariantCalls/Ploidy2/HybSubtracted"
InputDir=$BaseDir"/VariantCalls/Ploidy2/ParSubtracted"

if [[ ! -d $OutDir ]]; then
  mkdir -p $OutDir
else
  rm -rf $OutDir
  mkdir -p $OutDir
fi
cd $InputDir
# per evitare che si pianti nei re-run
rm -f $InputDir/*tbi $InputDir/Merged.txt.gz
if [[ $H0Label == "-" ]]; then
  # index vcf per bcftools merge
  # serve un for perché tabix
  # è scarso e non legge multi input
  for IndS in $(ls | grep -v Merged.txt.gz); do tabix -f $IndS; done
  AllFiles=$(ls *vcf.gz | grep -v ParHyb)
  bcftools merge --output-type v $AllFiles | grep -v "^##" | bgzip > Merged.txt.gz
  Rscript $BaseDir"/Scr/FiltSharedVar.R" "Ploidy2" $BaseDir
  cd $OutDir
  # anche bgzip è molto scarso
  for i in *; do bgzip $i; done
elif [[ $H0Label == "skip" ]]; then
  cp $InputDir/*vcf.gz $OutDir
  break
else
  # apro tutti i gz dei controlli per correggere il nome del campione
  for IndH0 in $H0Label
  do
    gunzip $IndH0.$Ref1Label+$Ref2Label.vcf.gz
  done

  for File1 in $(ls *vcf.gz)
  do
    if (( i % Nthreads == 0 )); then
      wait
    fi
    ((i++))
    (
    File1ID=$(echo $File1 | cut -d "." -f 1)
    # correzione nome campione
    for IndH0 in $H0Label
    do
      sed "s|$IndH0|$File1ID|g" < $IndH0.$Ref1Label+$Ref2Label.vcf > $IndH0.temp.4.$File1ID.vcf
    done
    bgzip *temp.4.$File1ID.vcf
    tabix *temp.4.$File1ID.vcf.gz
    tabix -f $File1
    # sottrazione varianti
    vcf-isec -c $File1 *temp.4.$File1ID.vcf.gz | bgzip > $OutDir"/"$File1
    rm -f *temp.4.$File1ID*
    ) &
  done
  
  wait
  
  for IndH0 in $H0Label
  do
    bgzip $IndH0.$Ref1Label+$Ref2Label.vcf
    cp $IndH0.$Ref1Label+$Ref2Label.vcf.gz $OutDir
  done
fi

# varianti unphased
OutDir=$BaseDir"/VariantCalls/Unphased/HybSubtracted"
InputDir=$BaseDir"/VariantCalls/Unphased/ParSubtracted"
BothRef=$Ref1Label" "$Ref2Label

if [[ ! -d $OutDir ]]; then
  mkdir -p $OutDir
else
  rm -rf $OutDir
  mkdir -p $OutDir
fi
for IndR in $BothRef
do
  cd $InputDir
  # per evitare che si pianti nei re-run (e in questo caso anche per evitare casini tra Ref1 e Ref2)
  rm -f $InputDir/*tbi $InputDir/Merged.txt.gz
  if [[ $H0Label == "-" ]]; then
    # index vcf per bcftools merge
    # serve un for perché tabix
    # è scarso e non legge multi input
    for IndS in $(ls *$IndR.vcf.gz | grep -v Merged.txt.gz); do tabix -f $IndS; done
    AllFiles=$(ls *$IndR.vcf.gz | grep -v ParHyb)
    bcftools merge --output-type v $AllFiles | grep -v "^##" | bgzip > Merged.txt.gz
    Rscript $BaseDir"/Scr/FiltSharedVar.R" "Unphased" $BaseDir
    cd $OutDir
    # anche bgzip è molto scarso
    for i in *vcf; do bgzip $i; done
  elif [[ $H0Label == "skip" ]]; then
    cp $InputDir/*vcf.gz $OutDir
    break
  else
    # apro tutti i gz dei controlli per correggere il nome del campione
    for IndH0 in $H0Label
    do
      gunzip $IndH0.$IndR.vcf.gz
    done
  
    for File1 in $(ls *vcf.gz | grep $IndR.vcf.gz)
    do
      if (( i % Nthreads == 0 )); then
        wait
      fi
      ((i++))
      (
      File1ID=$(echo $File1 | cut -d "." -f 1)
      # correzione nome campione
      for IndH0 in $H0Label
      do
        sed "s|$IndH0|$File1ID|g" < $IndH0.$IndR.vcf > $IndH0.temp.4.$File1ID.$IndR.vcf
      done
      bgzip *temp.4.$File1ID.$IndR.vcf
      tabix *temp.4.$File1ID.$IndR.vcf.gz
      tabix -f $File1
      # sottrazione varianti
      vcf-isec -c $File1 *temp.4.$File1ID.$IndR.vcf.gz | bgzip > $OutDir"/"$File1
      rm -f *temp.4.$File1ID*
      ) &
    done
    
    wait
    
    for IndH0 in $H0Label
    do
      bgzip $IndH0.$IndR.vcf
      cp $IndH0.$IndR.vcf.gz $OutDir
    done
  fi
done



