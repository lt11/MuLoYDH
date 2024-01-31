# header ------------------------------------------------------------------

# calculate low marker density regions
# ovvero le regioni con marker distanti più di 500 bp

rm(list = ls())
options(stringsAsFactors = F)

SpecDec <- function(x, k) as.numeric(format(round(x, k), nsmall = k))

# settings ----------------------------------------------------------------

# arguments
argsVal <- commandArgs(trailingOnly = T)
baseDir <- argsVal[1]

distThreshold <- 1000
outDir <- file.path(baseDir, "MUMmer", "LowDensityRegions")
dir.create(path = outDir, showWarnings = F, recursive = T)

mumDir <- file.path(baseDir, "MUMmer")
inputDir <- list.files(mumDir, pattern = "_", full.names = T)

allChr <- c("chrI", "chrII", "chrIII", "chrIV", "chrV", "chrVI", "chrVII", "chrVIII", 
            "chrIX", "chrX", "chrXI", "chrXII", "chrXIII", "chrXIV", "chrXV", "chrXVI")
colNames <- c("posR1", "allR1", "allR2", "posR2", "buff", "dist", "lenR1", "lenR2", "frm", "tags", "chromR1", "chromR2")

inputData <- list.files(inputDir, pattern = "intersect\\.snps$", full.names = T)
bothRef <- unlist(strsplit(x = gsub(pattern = "\\.intersect\\.snps", replacement = "", basename(inputData)), split = "_"))

nucmerData <- read.table(file = inputData, header = F, sep = "\t", col.names = colNames, skip = 4)

# init out files
outFileR1 <- file.path(outDir, paste0(bothRef[1], ".LowDensityRegions.SubTel.txt"))
headerString <- c("chrom", "start", "end", "\n")
cat(x = headerString, file = outFileR1, append = F, sep = "\t")
outFileR2 <- file.path(outDir, paste0(bothRef[2], ".LowDensityRegions.SubTel.txt"))
headerString <- c("chrom", "start", "end", "\n")
cat(x = headerString, file = outFileR2, append = F, sep = "\t")

# init out summary files by chromosome
outSummaryFileR1 <- file.path(outDir, paste0(bothRef[1], ".Summary.LMDR.SubTel.txt"))
headerString <- c("chrom", "LMDR [bp]", "length [bp]", "LMDR fraction", "\n")
cat(x = headerString, file = outSummaryFileR1, append = F, sep = "\t")
outSummaryFileR2 <- file.path(outDir, paste0(bothRef[2], ".Summary.LMDR.SubTel.txt"))
headerString <- c("chrom", "LMDR [bp]", "length [bp]", "LMDR fraction", "\n")
cat(x = headerString, file = outSummaryFileR2, append = F, sep = "\t")

# remove MT markers
indNoMT <- which(nucmerData$chromR1 != "chrMT" & nucmerData$chromR2 != "chrMT")
nucmerDataFilt <- nucmerData[indNoMT, ]

# load chromosome length data
inFileChLenR1 <- file.path(baseDir, "Ref", "Mod", paste0(bothRef[1], ".genome.chrref.fa.fai"))
refR1ChromLen <- read.table(inFileChLenR1, nrows = 16)
inFileChLenR2 <- file.path(baseDir, "Ref", "Mod", paste0(bothRef[2], ".genome.chrref.fa.fai"))
refR2ChromLen <- read.table(inFileChLenR2, nrows = 16)

# load (sub)telomers lengths
inFileSubLenR1 <- file.path(baseDir, "Ref", "Ann", paste0(bothRef[1], ".subtel.txt"))
lenSubTelR1 <- read.table(inFileSubLenR1)
totSubTelR1 <- lenSubTelR1[, 1] + refR1ChromLen[, 2] - lenSubTelR1[, 2]
inFileSubLenR2 <- file.path(baseDir, "Ref", "Ann", paste0(bothRef[2], ".subtel.txt"))
lenSubTelR2 <- read.table(inFileSubLenR2)
totSubTelR2 <- lenSubTelR2[, 1] + refR2ChromLen[, 2] - lenSubTelR2[, 2]

# reference 1 -------------------------------------------------------------

counterChr <- 1
for (indC in allChr) {
  # estreai posizioni cromosoma
  posRefR1 <- nucmerDataFilt$posR1[which(nucmerDataFilt$chromR1 == indC)]
  
  # splendido esempio di posizioni su cromosoma lungo 5150
  # posRefR1 <- c(151, 1000, 512, 1020, 1029, 1571, 3011, 3012, 3030, 2010, 
  #               3040, 3200, 3232, 3244, 3333, 3500, 3700, 4000, 4049)
  # posRefR1Sorted <- sort(posRefR1)
  # posRefR1SortedExtreme <- c(1, posRefR1Sorted, 5150)
  # diffPosRefR1SortedExtreme <- diff(posRefR1SortedExtreme)
  # startReg <- c()
  # endReg <- c()
  # flagReg <- 0
  
  # sort by pos
  posRefR1Sorted <- sort(posRefR1)
  # concatena 1, posizioni, length chromosome in posRefR1
  posRefR1SortedExtreme <- c(1, posRefR1Sorted, refR1ChromLen[counterChr, 2])
  # distance
  diffPosRefR1SortedExtreme <- diff(posRefR1SortedExtreme)
  startReg <- c()
  endReg <- c()
  flagReg <- 0
  # loop agghiacciante
  for (indD in 1:length(diffPosRefR1SortedExtreme)) {
    if (diffPosRefR1SortedExtreme[indD] > distThreshold & flagReg == 0) {
      # inizio regione
      startReg <- c(startReg, posRefR1SortedExtreme[indD])
      flagReg <- 1
    } else if (diffPosRefR1SortedExtreme[indD] > distThreshold & flagReg == 1) {
      # proseguimento regione
      next()
    } else if (diffPosRefR1SortedExtreme[indD] <= distThreshold & flagReg == 1) {
      # la regione è finita 
      endReg <- c(endReg, posRefR1SortedExtreme[indD])
      flagReg <- 0
    }
  } # loop agghiacciante
  
  if (flagReg == 1) {
    endReg <- c(endReg, posRefR1SortedExtreme[length(posRefR1SortedExtreme)])
  }
  dfReg <- data.frame(startReg, endReg)
  # accumula stats by chrom
  chrVect <- rep(indC, nrow(dfReg))
  dfOut <- data.frame(chrVect, dfReg)
  # sottrai (sub)telomeri
  # indBon <- which(dfOut$startReg > lenSubTelR1[counterChr, 1] & dfOut$endReg < lenSubTelR1[counterChr, 2])
  # dfOut <- dfOut[indBon, ]
  # cat(length(indBon), nrow(dfOut) , "\n")
  # save dfOut
  write.table(x = dfOut, file = outFileR1, append = T, row.names = F, col.names = F, quote = F, sep = "\t")
  # save whole-chromosome stats
  if (is.data.frame(dfOut)) {
    sumLowDenReg <- sum(dfOut$endReg - dfOut$startReg)
    outSummaryR1 <- data.frame(indC, sumLowDenReg, refR1ChromLen[counterChr, 2], 
                               SpecDec(x = sumLowDenReg / refR1ChromLen[counterChr, 2], 5))
    write.table(x = outSummaryR1, file = outSummaryFileR1, append = T, row.names = F, col.names = F, quote = F, sep = "\t")
  } else {
    sumLowDenReg <- 0
    outSummaryR1 <- data.frame(indC, sumLowDenReg, refR1ChromLen[counterChr, 2], 
                               SpecDec(x = sumLowDenReg / refR1ChromLen[counterChr, 2], 5))
    write.table(x = outSummaryR1, file = outSummaryFileR1, append = T, row.names = F, col.names = F, quote = F, sep = "\t")
  }
  # increment chromosome
  counterChr <- counterChr + 1
} # chromosome loop

# reference 2 -------------------------------------------------------------

counterChr <- 1
for (indC in allChr) {
  # estreai posizioni cromosoma
  posRefR2 <- nucmerDataFilt$posR2[which(nucmerDataFilt$chromR2 == indC)]
  
  # splendido esempio di posizioni su cromosoma lungo 5150
  # posRefR2 <- c(151, 1000, 512, 1020, 1029, 1571, 3011, 3012, 3030, 2010, 
  #               3040, 3200, 3232, 3244, 3333, 3500, 3700, 4000, 4049)
  # posRefR2Sorted <- sort(posRefR2)
  # posRefR2SortedExtreme <- c(1, posRefR2Sorted, 5150)
  # diffPosRefR2SortedExtreme <- diff(posRefR2SortedExtreme)
  # startReg <- c()
  # endReg <- c()
  # flagReg <- 0
  
  # sort by pos
  posRefR2Sorted <- sort(posRefR2)
  # concatena 1, posizioni, length chromosome in posRefR2
  posRefR2SortedExtreme <- c(1, posRefR2Sorted, refR2ChromLen[counterChr, 2])
  # distance
  diffPosRefR2SortedExtreme <- diff(posRefR2SortedExtreme)
  startReg <- c()
  endReg <- c()
  flagReg <- 0
  # loop agghiacciante
  for (indD in 1:length(diffPosRefR2SortedExtreme)) {
    if (diffPosRefR2SortedExtreme[indD] > distThreshold & flagReg == 0) {
      # inizio regione
      startReg <- c(startReg, posRefR2SortedExtreme[indD])
      flagReg <- 1
    } else if (diffPosRefR2SortedExtreme[indD] > distThreshold & flagReg == 1) {
      # proseguimento regione
      next()
    } else if (diffPosRefR2SortedExtreme[indD] <= distThreshold & flagReg == 1) {
      # la regione è finita 
      endReg <- c(endReg, posRefR2SortedExtreme[indD])
      flagReg <- 0
    }
  } # loop agghiacciante
  
  if (flagReg == 1) {
    endReg <- c(endReg, posRefR2SortedExtreme[length(posRefR2SortedExtreme)])
  }
  dfReg <- data.frame(startReg, endReg)
  # accumula stats by chrom
  chrVect <- rep(indC, nrow(dfReg))
  dfOut <- data.frame(chrVect, dfReg)
  # sottrai (sub)telomeri
  # indBon <- which(dfOut$startReg > lenSubTelR2[counterChr, 1] & dfOut$endReg < lenSubTelR2[counterChr, 2])
  # dfOut <- dfOut[indBon, ]
  # cat(length(indBon), nrow(dfOut) , "\n")
  # save dfOut
  write.table(x = dfOut, file = outFileR2, append = T, row.names = F, col.names = F, quote = F, sep = "\t")
  # save whole-chromosome stats
  if (is.data.frame(dfOut)) {
    sumLowDenReg <- sum(dfOut$endReg - dfOut$startReg)
    outSummaryR2 <- data.frame(indC, sumLowDenReg, refR2ChromLen[counterChr, 2], 
                               SpecDec(x = sumLowDenReg / refR2ChromLen[counterChr, 2], 5))
    write.table(x = outSummaryR2, file = outSummaryFileR2, append = T, row.names = F, col.names = F, quote = F, sep = "\t")
  } else {
    sumLowDenReg <- 0
    outSummaryR2 <- data.frame(indC, sumLowDenReg, refR2ChromLen[counterChr, 2], 
                               SpecDec(x = sumLowDenReg / refR2ChromLen[counterChr, 2], 5))
    write.table(x = outSummaryR2, file = outSummaryFileR2, append = T, row.names = F, col.names = F, quote = F, sep = "\t")
  }
  # increment chromosome
  counterChr <- counterChr + 1
} # chromosome loop


