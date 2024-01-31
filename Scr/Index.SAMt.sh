#!/bin/bash

WorkDir=$1"/Ref/Mod"

cd $WorkDir

for Ind1 in *.fa
do 
  samtools faidx $Ind1
done


