#!/bin/bash

RefDir=$1"/Ref/Mod"

for Ind in $RefDir/*fa
do
  bwa index $Ind &
done 

wait


