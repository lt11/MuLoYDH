#!/bin/bash

# genera gli esperimenti 
# con le read usate per gli assembly
# dei parent

# occhio!
# l'ID di S288C è sbagliato!!!
# seqtk sample /home/Share/jxyue/12Strains_Illumina_FASTQ/S288c.R1.fastq.gz 0.12 > PAR$ParID2.R1.fastq &
# seqtk sample /home/Share/jxyue/12Strains_Illumina_FASTQ/S288c.R2.fastq.gz 0.12 > PAR$ParID2.R2.fastq &

# variables
ParID1="UWOPS034614"
ShortID1="MA"
ParID2="YPS128"
ShortID2="NA"

BaseDir=$(pwd | sed 's|....$||')

if [[ ! -d $BaseDir/Exp ]]; then
  mkdir -p $BaseDir/Exp
fi

cd $BaseDir/Exp

seqtk sample /home/Share/jxyue/12Strains_Illumina_FASTQ/$ParID1.R1.fastq.gz 0.17 > PAR$ParID1.R1.fastq &
seqtk sample /home/Share/jxyue/12Strains_Illumina_FASTQ/$ParID1.R2.fastq.gz 0.17 > PAR$ParID1.R2.fastq &

seqtk sample /home/Share/jxyue/12Strains_Illumina_FASTQ/$ParID2.R1.fastq.gz 0.08 > PAR$ParID2.R1.fastq &
seqtk sample /home/Share/jxyue/12Strains_Illumina_FASTQ/$ParID2.R2.fastq.gz 0.08 > PAR$ParID2.R2.fastq &

wait

cat PAR$ParID1.R1.fastq PAR$ParID2.R1.fastq > ParHyb$ShortID1$ShortID2.R1.fastq &
cat PAR$ParID1.R2.fastq PAR$ParID2.R2.fastq > ParHyb$ShortID1$ShortID2.R2.fastq &

wait

rm -f PAR$ParID1* PAR$ParID2*

gzip *fastq


