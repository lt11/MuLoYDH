#!/bin/bash

#####################
### user settings ###
#####################
 
# assemblies labels (must match prefix of assemblies stored in "Rep" folder)
Ref1Label="YPS128"
Ref2Label="DBVPG6765"
# short labels (used to name file and plotting)
Ref1="NA"
Ref2="WE"
# prefix of optional parental short read data, if it's not the variable must be set to "-"
ParentalZero="-"
# prefix of optional control short read data
# if it's not the variable can be set to "-" (to filter out variants shared by at least 95% of the samples)
# or "skip" (to bypass this step)
HybridZero="A887R46"
# label of optional chromosome bearing selection marker (will be filtered before calculation of LOH statistics)
# if no selection marker has been used the variable must be set to "-"
MarkerChrom="-"
# "collinear" or "rearranged" parental genomes
GenomesStructure="collinear"
# reciprocal overlap for LOH filtering
OverlapThreshold=0.5
# read length
ReadLength="150"

#####################
### settings' end ###
#####################

### part 0

# base folder
BaseDir=$(dirname "$(pwd)")
# path to SnpEff.jar
SnpEff="${BaseDir}"/Lib/snpEff/snpEff.jar
# check Logs folder
if [[ ! -d Logs ]]; then mkdir Logs; fi

### part I: references & annotations initialization

mkdir -p "${BaseDir}"/Ref/Ann
cp "${BaseDir}"/Rep/Asm/"${Ref1Label}".genome.fa "${BaseDir}"/Ref &
cp "${BaseDir}"/Rep/Asm/"${Ref2Label}".genome.fa "${BaseDir}"/Ref &
cp "${BaseDir}"/Rep/Ann/"${Ref1Label}".* "${BaseDir}"/Ref/Ann &
cp "${BaseDir}"/Rep/Ann/"${Ref2Label}".* "${BaseDir}"/Ref/Ann &

wait

### part II: mapping

# bash QualCheck1.sh "${BaseDir}" &> Logs/QualCheck1.log &

/usr/bin/time -v Rscript FormatRef.CompMap.R "${Ref1Label}" "${Ref2Label}" "${BaseDir}"/Ref > Logs/FormatRef.CompMap.out 2> Logs/Time.FormatRef.CompMap.err

(
/usr/bin/time -v bash Index.BWA.sh "${BaseDir}" > Logs/Index.BWA.out 2> Logs/Time.Index.BWA.err &

/usr/bin/time -v bash Index.SAMt.sh "${BaseDir}" &> Logs/Index.SAMt.log &

wait
)

/usr/bin/time -v bash BWA.Mapping.sh "${BaseDir}" > Logs/BWA.Mapping.out 2> Logs/Time.BWA.Mapping.err

/usr/bin/time -v bash SAMt.Dup.sh "${BaseDir}" > Logs/SAMt.Dup.out 2> Logs/Time.SAMt.Dup.err

/usr/bin/time -v bash LowMAPQ.sh "${HybridZero}" "${ParentalZero}" "${BaseDir}" > Logs/LowMAPQ.out 2> Logs/Time.LowMAPQ.err

(
/usr/bin/time -v bash CalcDepth.Raw.sh "${BaseDir}" > Logs/CalcDepth.Raw.out 2> Logs/Time.CalcDepth.Raw.err &

/usr/bin/time -v bash CalcDepth.Ddp.sh "${BaseDir}" > Logs/CalcDepth.Ddp.out 2> Logs/Time.CalcDepth.Ddp.err &

wait

/usr/bin/time -v Rscript DepthParser.R "${BaseDir}" > Logs/DepthParser.out 2> Logs/Time.DepthParser.err &
) &

bash QualCheck2.sh "${BaseDir}" &> Logs/QualCheck2.log &

### part III: markers calculation, genotyping & variants calling

## markers genotyping

(
/usr/bin/time -v bash MUMmer.sh "${Ref1Label}" "${Ref2Label}" "${GenomesStructure}" "${BaseDir}" > Logs/MUMmer.out 2> Logs/Time.MUMmer.err

/usr/bin/time -v Rscript TableLMDR.R "${ReadLength}" "${BaseDir}" > Logs/TableLMDR.out 2> Logs/Time.TableLMDR.err &

/usr/bin/time -v bash SAMt.Marker.sh "${Ref1Label}" "${Ref2Label}" "${BaseDir}" > Logs/SAMt.Marker.out 2> Logs/Time.SAMt.Marker.err

/usr/bin/time -v Rscript Marker.Parser.R "${BaseDir}" > Logs/Marker.Parser.out 2> Logs/Time.Marker.Parser.err
) &

## copy number variants call

(
/usr/bin/time -v bash Mappability.GEM.sh "${Ref1Label}" "${Ref2Label}" "${ReadLength}" "${BaseDir}" > Logs/Mappability.GEM.out 2> Logs/Time.Mappability.GEM.err

/usr/bin/time -v bash PreControl-FREEC.sh "${Ref1Label}" "${Ref2Label}" "${BaseDir}" > Logs/PreControl-FREEC.out 2> Logs/Time.PreControl-FREEC.err

/usr/bin/time -v bash Control-FREEC.sh "${BaseDir}" > /dev/null 2> Logs/Time.Control-FREEC.err
) &

## LOH detection & plotting

wait

/usr/bin/time -v Rscript ClrS.V8.R "${Ref1Label}" "${Ref2Label}" "${Ref1}" "${Ref2}" "${BaseDir}" > Logs/ClrS.V8.out 2> Logs/Time.ClrS.V8.err

/usr/bin/time -v Rscript Calc.Marker.R "${Ref1Label}" "${Ref2Label}" "${Ref1}" "${Ref2}" "${BaseDir}" > Logs/Calc.Marker.out 2> Logs/Time.Calc.Marker.err

/usr/bin/time -v Rscript BAFplot.R "${Ref1Label}" "${Ref2Label}" "${Ref1}" "${Ref2}" "${BaseDir}" > Logs/BAFplot.out 2> Logs/Time.BAFplot.err &

/usr/bin/time -v Rscript Marker.Dist.R "${Ref1Label}" "${Ref2Label}" "${Ref1}" "${Ref2}" "${BaseDir}" > Logs/Marker.Dist.out 2> Logs/Time.Marker.Dist.err &

/usr/bin/time -v Rscript AllSeg.R "${Ref1Label}" "${Ref2Label}" "${Ref1}" "${Ref2}" "${BaseDir}" > Logs/AllSeg.out 2> Logs/Time.AllSeg.err

/usr/bin/time -v Rscript MakeBed.HomoHeteroReg.R "${Ref1Label}" "${Ref2Label}" "${Ref1}" "${Ref2}" "${BaseDir}" > Logs/MakeBed.HomoHeteroReg.out 2> Logs/Time.MakeBed.HomoHeteroReg.err

/usr/bin/time -v bash MakeBed.LowMAPQ.sh "${Ref1Label}" "${Ref2Label}" "${BaseDir}" > Logs/MakeBed.LowMAPQ.out 2> Logs/Time.MakeBed.LowMAPQ.err

/usr/bin/time -v Rscript HeatMap.R "${Ref1Label}" "${Ref1}" "${BaseDir}" > Logs/HeatMap.out 2> Logs/Time.HeatMap.Ref1.err &

/usr/bin/time -v Rscript HeatMap.R "${Ref2Label}" "${Ref2}" "${BaseDir}" > Logs/HeatMap.out 2> Logs/Time.HeatMap.Ref2.err &

## small variants

(
/usr/bin/time -v bash SAMt.Call.CompMap.sh "${Ref1Label}" "${Ref2Label}" "${BaseDir}" > Logs/SAMt.Call.CompMap.out 2> Logs/Time.SAMt.Call.CompMap.err &

/usr/bin/time -v bash FreeBayes.Call.CompMap.sh "${Ref1Label}" "${Ref2Label}" "${BaseDir}" > Logs/FreeBayes.Call.CompMap.out 2> Logs/Time.FreeBayes.Call.CompMap.err &

wait

/usr/bin/time -v bash Intersect.CompMap.sh "${Ref1Label}" "${Ref2Label}" "${BaseDir}" > Logs/Intersect.CompMap.out 2> Logs/Time.Intersect.CompMap.err

/usr/bin/time -v bash SubtractPar.CompMap.sh "${Ref1Label}" "${Ref2Label}" "${ParentalZero}" "${BaseDir}" > Logs/SubtractPar.CompMap.out 2> Logs/Time.SubtractPar.CompMap.err

/usr/bin/time -v bash SubtractHyb.CompMap.sh "${Ref1Label}" "${Ref2Label}" "${HybridZero}" "${BaseDir}" > Logs/SubtractHyb.CompMap.out 2> Logs/Time.SubtractHyb.CompMap.err
)

/usr/bin/time -v Rscript VariantParser.R "${BaseDir}" > Logs/VariantParser.out 2> Logs/Time.VariantParser.err

### part IV: annotations, statistics

/usr/bin/time -v Rscript Annotation.SNV.CompMap.R "${Ref1Label}" "${Ref2Label}" "${Ref1}" "${Ref2}" "${BaseDir}" > Logs/Annotation.SNV.CompMap.out 2> Logs/Time.Annotation.SNV.CompMap.err

/usr/bin/time -v bash SnpEff.CompMap.sh "${Ref1Label}" "${Ref2Label}" "${Ref1}" "${Ref2}" "${SnpEff}" "${BaseDir}" > Logs/SnpEff.CompMap.out 2> Logs/Time.SnpEff.CompMap.err &

(
/usr/bin/time -v Rscript Annotation.InterTerm.V2.R "${Ref1Label}" "${Ref2Label}" "${Ref1}" "${Ref2}" "${BaseDir}" > Logs/Annotation.InterTerm.V2.out 2> Logs/Time.Annotation.InterTerm.V2.err

/usr/bin/time -v Rscript StatLOH.V1.R "${OverlapThreshold}" "${MarkerChrom}" "${BaseDir}" > Logs/StatLOH.V1.out 2> Logs/Time.StatLOH.V1.err
) &

/usr/bin/time -v Rscript Annotation.Features.V2.R "${Ref1Label}" "${Ref2Label}" "${Ref1}" "${Ref2}" "${BaseDir}" > Logs/Annotation.Features.V2.out 2> Logs/Time.Annotation.Features.V2.err

wait

# rm -rf "${BaseDir}"/MapRaw

echo "Finished!" Logs/Finished.txt
