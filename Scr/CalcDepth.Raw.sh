#!/bin/bash

# Calculate depth of coverage statistics and check REF name consistency
# between file name and bam file header

BaseDir=$1

BamDir=${BaseDir}"/MapRaw"
Depth=${BaseDir}"/DepthRaw"
Nthreads=1

if [[ ! -d $Depth ]]; then
  mkdir -p $Depth
fi

cd $BamDir

for Ind1 in *bam
do
  if (( i % Nthreads == 0 )); then
    wait
  fi
  ((i++))
  (
  AlnName=$(echo $Ind1 | cut -d "." -f 1,2)
  samtools depth $Ind1 > $Depth/$AlnName.depth.txt
  RefLen=$(samtools view -H $Ind1 | grep "^@SQ" | cut -f 3 | cut -d ":" -f 2 | awk '{sum += $1} END {print sum}')
  awk -v RL=$RefLen 'BEGIN{FS="\t"}{sum+=$3; sumsq+=$3*$3} END {print "N covered sites=" NR; print "% reference covered=" 100*NR/RL; print "Reference length=" RL; print "Average=",sum/NR; print "Stdev=",sqrt(sumsq/NR - (sum/NR)**2)}' $Depth/$AlnName.depth.txt > $Depth/$AlnName.mean.txt
  rm -f $Depth/$AlnName.depth.txt
  
  # Check reference consistency: file name VS bam header
  RefHeader=$(samtools view -H $Ind1 | grep "^@PG" | head -1 | cut -f 5 | cut -d " " -f 6 | cut -d "/" -f 7 | cut -d "." -f 1)
  RefFileName=$(echo $Ind1 | cut -d "." -f 2)
  if [[ "$RefHeader" != "$RefFileName"  ]]; then
    echo "ERROR: REF file name did not match REF bam header, $AlnName"
  fi
  ) &
  
done

wait


