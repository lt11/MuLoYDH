#!/bin/bash

Ref1=$1
Ref2=$2
ReadLen=$3
BaseDir=$4

Nthreads=8
CopyNumDir=${BaseDir}"/CNV"
MapDir=${BaseDir}"/CNV/Mappability"
ModDir=${BaseDir}"/Ref/Mod"

AllRef="${Ref1} ${Ref2}"

if [[ ! -d ${MapDir} ]]; then
  mkdir -p ${MapDir}
fi

for IndRef in ${AllRef}
do
  WorkDir=$MapDir"/"$IndRef
  if [[ ! -d ${WorkDir} ]]; then
    mkdir -p ${WorkDir}
  fi
  
  cp $ModDir"/"$IndRef".genome.chrref.fa" $WorkDir
  cd $WorkDir
  
  gem-indexer -T $Nthreads -c dna -i $IndRef".genome.chrref.fa" -o $IndRef
  gem-mappability -T $Nthreads -I  $IndRef".gem" -l $ReadLen -o  $IndRef
  # output $IndRef.mappability
  rm -f $IndRef".genome.chrref.fa"
done


