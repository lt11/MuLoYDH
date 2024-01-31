# header ------------------------------------------------------------------

# sottrae le LOH dei controlli dai campioni (se ce ne sono)
# calcola la frazione di genoma in LOH
# e tutte le altre statistiche,
# per chrII si intende il cromosoma con il marker di selezione

# dati:
# myData -> dati segmenti status 0, 1 da Annotations.SNV.IntTer.txt
# da cui genero myControls e myData
# da myControls prendo i segmenti status 0 -> myControlsStatusZero
# (notare che qui c'è anche chrII)
# da myData -> tolgo chrII (myDataFiltChrII) e poi prendo solo gli status 0 (myDataFiltChrIIStatusZero)

rm(list = ls())
options(stringsAsFactors = F)
library(ggplot2)

# function(s) -------------------------------------------------------------

ReciprocalOverlap  <-  function(intervalA, intervalB, frcOverlap) {
  # calculates reciprocal overlap between two intervals
  # arguments:
  #   intervalA: a vector of length 2
  #   intervalB: a vector of length 2
  #   frcOverlap: minimum fraction of reciprocal overlap
  # returns:
  #   a list of [[1]] a binary number, 0 for no overlap, 1 for reciprocal overlap
  #   and [[2]] the new coordinates of the overlapping segment
  
  startA <- intervalA[1]
  endA <- intervalA[2]
  startB <- intervalB[1]
  endB <- intervalB[2]
  lengthA <- endA - startA
  lengthB <- endB - startB
  binOverlap <- 0
  newCoord <- c(0, 0)
  if (startA <= startB & endA >= startB & endA <= endB) {
    lenCommon <- endA - startB
    lenCommonFractA <- lenCommon / lengthA
    lenCommonFractB <- lenCommon / lengthB
    if (lenCommonFractA > frcOverlap & lenCommonFractB > frcOverlap) {
      binOverlap <- 1
      newCoord <- c(startA, endB)
    }
  }
  if (startA <= startB & endA >= endB) {
    lenCommon <- endB - startB
    lenCommonFractA <- lenCommon / lengthA
    lenCommonFractB <- lenCommon / lengthB
    if (lenCommonFractA > frcOverlap & lenCommonFractB > frcOverlap) {
      binOverlap <- 1
      newCoord <- c(startA, endA)
    }
  }
  if (startA >= startB & startA <= endB & endA >= endB) {
    lenCommon <- endB - startA
    lenCommonFractA <- lenCommon / lengthA
    lenCommonFractB <- lenCommon / lengthB
    if (lenCommonFractA > frcOverlap & lenCommonFractB > frcOverlap) {
      binOverlap <- 1
      newCoord <- c(startB, endA)
    }
  }
  if (startA >= startB & startA <= endB & endA <= endB) {
    lenCommon <- endA - startA
    lenCommonFractA <- lenCommon / lengthA
    lenCommonFractB <- lenCommon / lengthB
    if (lenCommonFractA > frcOverlap & lenCommonFractB > frcOverlap) {
      binOverlap <- 1
      newCoord <- c(startB, endB)
    }
  }
  lstRes <- list()
  lstRes$binOverlap <- binOverlap
  lstRes$newCoord <- newCoord
  return(lstRes)
}

# settings ----------------------------------------------------------------

# significant figures
nDecimal <- 5

argsVal <- commandArgs(trailingOnly = T)
overlapThreshold <- as.numeric(argsVal[1])
markerChrom <- argsVal[2]
baseDir <- argsVal[3]
# commentami
# baseDir <- "/Users/lorenzo/Desktop/Out\ of\ Control"
# commentami
# overlapThreshold <- 0.5
# commentami
# markerChrom <- "-"

# make output dirs
dirOutTab <- file.path(baseDir, "LOH", "Summary", "Tables")
dir.create(path = dirOutTab, showWarnings = F, recursive = T)
dirOutStat <- file.path(baseDir, "LOH", "Summary", "Stat")
dir.create(path = dirOutStat, showWarnings = F, recursive = T)
dirOutPlot <- file.path(baseDir, "LOH", "Summary", "Plot")
dir.create(path = dirOutPlot, showWarnings = F, recursive = T)
# output table file
fileTable <- file.path(dirOutTab, "LOH.filtered.variants.txt")
# output RData file for TestRearranged.R
fileForTestRearranged <- file.path(dirOutTab, "DataFiltChrSelectionMarker.StatusZero.RData")
# output stat file
fileStat <- file.path(dirOutStat, paste(basename(baseDir), "LOH.stats.txt", sep = "."))
unlink(fileStat, recursive = T)

allChr <- c("chrI", "chrII", "chrIII", "chrIV", "chrV", "chrVI", "chrVII", "chrVIII", 
            "chrIX", "chrX", "chrXI", "chrXII", "chrXIII", "chrXIV", "chrXV", "chrXVI")
distThreshold <- 25000
# numero di marker delle LOH di cui calcolo
# la stat "numero di LOH with 1:5 markers"
nMarkerThreshold <- 5

indF <- list.files(path = file.path(baseDir, "LOH"),
                   pattern = "^Annotations\\.SNV\\.IntTer\\.txt",
                   full.names = T, recursive = T)
# read summary table
myData <- read.table(file = indF, header = T, sep = "\t", na.strings = "antani")

# read T0, ref1, ref2
inputFilePath <- list.files(path = file.path(baseDir, "Scr"), pattern = "^MuLoYDH.*sh$", full.names = T)
inputFile <- readLines(con = inputFilePath)

# lettura argmenti che non "entrano" nella riga di comando
# (per evitare di farla lunga un km)
controlIDs <- gsub("^\"|\"$", "", 
                   unlist(strsplit(grep(pattern = "HybridZero=", x = inputFile, value = T)[1], split = "=", fixed = T))[2])
booControlID <- controlIDs == "-"
controlIDs <- unlist(strsplit(controlIDs, split = " "))

ref1Label <- gsub("^\"|\"$", "", 
                  unlist(strsplit(grep(pattern = "Ref1Label=", x = inputFile, value = T)[1], split = "=", fixed = T))[2])
ref2Label <- gsub("^\"|\"$", "", 
                  unlist(strsplit(grep(pattern = "Ref2Label=", x = inputFile, value = T)[1], split = "=", fixed = T))[2])
ref1 <- gsub("^\"|\"$", "", 
             unlist(strsplit(grep(pattern = "Ref1=", x = inputFile, value = T)[1], split = "=", fixed = T))[2])
ref2 <- gsub("^\"|\"$", "", 
             unlist(strsplit(grep(pattern = "Ref2=", x = inputFile, value = T)[1], split = "=", fixed = T))[2])

chrLenFilesRef1 <- list.files(path = file.path(baseDir, "CNV", "GCdata", ref1Label),
                              pattern = "LenChr\\.txt$", recursive = T, full.names = T)
chrLenFilesRef2 <- list.files(path = file.path(baseDir, "CNV", "GCdata", ref2Label),
                              pattern = "LenChr\\.txt$", recursive = T, full.names = T)

# read chromosome length data
allChrLenRef1 <- read.table(file = chrLenFilesRef1, header = F, sep = "\t")
allChrLenRef2 <- read.table(file = chrLenFilesRef2, header = F, sep = "\t")

if (markerChrom == "-") {
  chrLenRef1 <- allChrLenRef1
  chrLenRef2 <- allChrLenRef2
} else {
  chrLenRef1 <- allChrLenRef1[-which(allChr == markerChrom), ]
  chrLenRef2 <- allChrLenRef2[-which(allChr == markerChrom), ]
}
allChrLen <- list(allChrLenRef1, allChrLenRef2)

# genome length without chrII
genomeLength <- sum(chrLenRef1[, 3]) + sum(chrLenRef2[, 3])

# filters -----------------------------------------------------------------

# cancellami
# fakeEvent1 <- c("A314R20", "N44", "chrX", 332735, 332744, 338370, 338380, 0, 503, 0.0894, 11, 5645, 0.0896, 11, 5626, 0, 0, 0, 0, 0, 0, "int")
# fakeEvent2 <- c("A314R21", "N44", "chrX", 332735, 332744, 332786, 332886, 0, 5, 0.0894, 11, 151, 0.0896, 11, 42, 0, 0, 0, 0, 0, 0, "int")
# fakeEvent3 <- c("A314R20", "N44", "chrXVI", 432577, 432643, 432643, 432823,  0, 1, 0.0894, 11, 151, 0.0896, 11, 42, 0, 0, 0, 0, 0, 0, "int")
# myData <- rbind(myData, fakeEvent1, fakeEvent2, fakeEvent3)

# prepara la tabella dei controlli (H0 + ParHyb) per sottrarre le LOH dai campioni
indH0 <- c()
# se il controllo letto dal wrapper NON è == "-"
if (!booControlID) {
  for (labelH0 in controlIDs) {
    indH0 <- c(indH0, which(myData$sample == labelH0))
  }
}
indPar <- grep(pattern = "^ParHyb", myData$sample)
indControls <- c(indH0, indPar)

if (length(indControls) != 0) {
  myControls <- myData[indControls, ]
  myData <- myData[-indControls, ]
  
  indControlStatusZero <- which(myControls$status == 0)
  # init myControlsStatusZero così il for di indR gira anche se non ci sono controlli a status 0
  myControlsStatusZero <- c()
  if (length(indControlStatusZero) != 0) {
    myControlsStatusZero <- myControls[indControlStatusZero, ]
  }
}

# experiment
expName <- unlist(strsplit(indF, split = "/", fixed = T))
cat("Experiments ID:", expName[length(expName) - 3], "\n", file = fileStat, append = T)

# numero di esperimenti
mySamples <- unique(myData$sample)
nSamp <- length(mySamples)
cat("# of experiments:", nSamp, "\n", file = fileStat, append = T)

# calculation without chrII and heterozygous segments
if (markerChrom == "-") {
  myDataFiltChrII <- myData
  myDataFiltChrIIStatusZero <- myDataFiltChrII[which(myDataFiltChrII$status == 0), ]
} else {
  myDataFiltChrII <- myData[which(myData$chr != markerChrom), ]
  myDataFiltChrIIStatusZero <- myDataFiltChrII[which(myDataFiltChrII$status == 0), ]
}

indOutSamp <- c()
if (length(indControls) != 0) {
  # filter LOHs in myDataFiltChrIIStatusZero from LOHs in myControlsStatusZero
  # indici delle LOH nei campioni (in myDataFiltChrIIStatusZero) da buttare in Arno
  for (indR in c(ref1, ref2)) {
    indBon <- which(myControlsStatusZero$ref == indR)
    if (length(indBon) > 0) {
      myControlsStatusZeroRef <- myControlsStatusZero[indBon, ]
      for (indC in unique(myControlsStatusZeroRef$chr)) {
        myControlsStatusZeroRefChr <-myControlsStatusZeroRef[which(myControlsStatusZeroRef$chr == indC), ]
        # indOver = indici delle LOH nei campioni
        indOver <- which(myDataFiltChrIIStatusZero$ref == indR & myDataFiltChrIIStatusZero$chr == indC & myDataFiltChrIIStatusZero$len != 1)
        if (length(indOver) != 0) {
          # indPluto è una delle LOH indOver
          for (indPluto in indOver) {
            sampleEvent <- as.numeric(c(myDataFiltChrIIStatusZero$first[indPluto], myDataFiltChrIIStatusZero$last[indPluto]))
            # qui facciamo l'overlap di una LOH alla volta dei campioni VS quelle dei controlli
            # e accumuliamo l'output [[1]], ovvero 0 se non c'è overlap
            # o 1 se c'è overlap
            outRepOv <- c()
            # indAntani è l'indice di una delle LOH nei controlli
            for (indAntani in 1:nrow(myControlsStatusZeroRefChr)) {
              controlEvent <- as.numeric(c(myControlsStatusZeroRefChr$first[indAntani], 
                                           myControlsStatusZeroRefChr$last[indAntani]))
              outRepOv <- c(outRepOv, ReciprocalOverlap(sampleEvent, controlEvent, overlapThreshold)[[1]])
            }
            # se nell'output cumulativo c'è almeno un 1
            # allora l'LOH indPluto ha overlap con un controllo
            if (length(which(outRepOv == 1)) > 0) {
              indOutSamp <- c(indOutSamp, indPluto)
            }
          }
        }
      }
    }
  }
}

# init overlappedWithControlsSegments per poi concatenare
overlappedWithControlsSegments <- c()
if (length(indOutSamp) != 0) {
  # le LOH che hanno overlap con i controlli diventano 
  # segmenti in eterozigosi
  myDataFiltChrIIStatusZero$status[indOutSamp] <- 1
  myDataFiltChrIIStatusZero$TI[indOutSamp] <- "sub"
  # li metto da parte per fare poi i plot
  overlappedWithControlsSegments <- myDataFiltChrIIStatusZero[indOutSamp, ]
  # cancella l'LOH
  myDataFiltChrIIStatusZero <- myDataFiltChrIIStatusZero[-indOutSamp, ]
}

if (length(indControls) != 0) {
  # cerca le LOH single marker
  # in un modo rapido
  # se fai match con ref, chr, first, last, len possono uscire solo le len 1
  # il resto, ovvero roba con overlap 100% è già stata filtrata
  stringSamp <- paste(myDataFiltChrIIStatusZero$ref, myDataFiltChrIIStatusZero$chr, myDataFiltChrIIStatusZero$first, myDataFiltChrIIStatusZero$last, myDataFiltChrIIStatusZero$len, sep = "_")
  stringCont <- paste(myControlsStatusZero$ref, myControlsStatusZero$chr, myControlsStatusZero$first, myControlsStatusZero$last, myControlsStatusZero$len, sep = "_")
  # indici delle LOH con un solo marker nei campioni presenti anche nei controlli
  indOutSingMark <- which(is.element(stringSamp, stringCont))
  if (length(indOutSingMark) != 0) {
    myDataFiltChrIIStatusZero$status[indOutSingMark] <- 1
    myDataFiltChrIIStatusZero$TI[indOutSingMark] <- "sub"
    overlappedWithControlsSegments <- rbind(overlappedWithControlsSegments, myDataFiltChrIIStatusZero[indOutSingMark, ])
    myDataFiltChrIIStatusZero <- myDataFiltChrIIStatusZero[-indOutSingMark, ]
  }
}

# calculations ------------------------------------------------------------

lenLOH <- myDataFiltChrIIStatusZero$distES

singSampPerc <- c()
for (indS in mySamples) {
  # % genoma in LOH - single sample
  singSampPerc <- c(singSampPerc, 
                    sum(myDataFiltChrIIStatusZero$distES[which(myDataFiltChrIIStatusZero$sample == indS)]) / genomeLength)
}

# frazione markers HOM ref1/ref2 per campione
for (indS in unique(myDataFiltChrII$sample)) {
  nHom <- sum(myDataFiltChrIIStatusZero$len[which(myDataFiltChrIIStatusZero$ref == ref1 & myDataFiltChrIIStatusZero$sample == indS)])
  nHet <- sum(myDataFiltChrII$len[which(myDataFiltChrII$ref == ref1 & myDataFiltChrII$status == 1 & myDataFiltChrII$sample == indS)])
  cat("Fraction of markers HOM", ref1, paste0(indS, ": "), file = fileStat, append = T)
  cat(signif(nHom / (nHom + nHet), nDecimal), "\n", file = fileStat, append = T)
  
  nHom <- sum(myDataFiltChrIIStatusZero$len[which(myDataFiltChrIIStatusZero$ref == ref2 & myDataFiltChrIIStatusZero$sample == indS)])
  nHet <- sum(myDataFiltChrII$len[which(myDataFiltChrII$ref == ref2 & myDataFiltChrII$status == 1 & myDataFiltChrII$sample == indS)])
  cat("Fraction of markers HOM", ref2, paste0(indS, ": "), file = fileStat, append = T)
  cat(signif(nHom / (nHom + nHet), nDecimal), "\n", file = fileStat, append = T)
}

# mean fraction of genome in LOH (per sample)
cat("##########################\n", file = fileStat, append = T)
cat("Mean statistics per sample\n", file = fileStat, append = T)
cat("##########################\n", file = fileStat, append = T)
cat("Fraction of genome in LOH: ", file = fileStat, append = T)
cat(signif(mean(singSampPerc), nDecimal), signif(sd(singSampPerc), nDecimal), signif(median(singSampPerc), nDecimal), "\n", 
    file = fileStat, append = T)
cat("mean", "sd", "median", "\n", file = fileStat, append = T)

# numero di LOH > distThreshold
indDistMajRef1 <- which(myDataFiltChrIIStatusZero$distLF > distThreshold & myDataFiltChrIIStatusZero$ref == ref1)
cat("# LOH >", distThreshold / 1000, "kb", ref1, "\n", file = fileStat, append = T)
cat(signif(length(indDistMajRef1) / nSamp, nDecimal), "\n", file = fileStat, append = T)
indDistMajRef2 <- which(myDataFiltChrIIStatusZero$distLF > distThreshold & myDataFiltChrIIStatusZero$ref == ref2)
cat("# LOH >", distThreshold / 1000, "kb", ref2, "\n", file = fileStat, append = T)
cat(signif(length(indDistMajRef2) / nSamp, nDecimal), "\n", file = fileStat, append = T)
# numero di LOH <= distThreshold
cat("# LOH <=", distThreshold / 1000, "kb", ref1, "\n", file = fileStat, append = T)
cat(signif(length(which(myDataFiltChrIIStatusZero$distLF <= distThreshold & myDataFiltChrIIStatusZero$ref == ref1)) / nSamp, nDecimal), "\n", file = fileStat, append = T)
cat("# LOH <=", distThreshold / 1000, "kb", ref2, "\n", file = fileStat, append = T)
cat(signif(length(which(myDataFiltChrIIStatusZero$distLF <= distThreshold & myDataFiltChrIIStatusZero$ref == ref2)) / nSamp, nDecimal), "\n", file = fileStat, append = T)

# numero di LOH with more than 5 markers
cat("# LOH with more than 5 markers", ref1, "\n", file = fileStat, append = T)
cat(signif(length(which(myDataFiltChrIIStatusZero$len > nMarkerThreshold & myDataFiltChrIIStatusZero$ref == ref1)) / nSamp, nDecimal), "\n", file = fileStat, append = T)
cat("# LOH with more than 5 markers", ref2, "\n", file = fileStat, append = T)
cat(signif(length(which(myDataFiltChrIIStatusZero$len > nMarkerThreshold & myDataFiltChrIIStatusZero$ref == ref2)) / nSamp, nDecimal), "\n", file = fileStat, append = T)

# numero di LOH with 5 or less markers
cat("# LOH with 5 or less markers", ref1, "\n", file = fileStat, append = T)
cat(signif(length(which(myDataFiltChrIIStatusZero$len <= nMarkerThreshold & myDataFiltChrIIStatusZero$ref == ref1)) / nSamp, nDecimal), "\n", file = fileStat, append = T)
cat("# LOH with 5 or less markers", ref2, "\n", file = fileStat, append = T)
cat(signif(length(which(myDataFiltChrIIStatusZero$len <= nMarkerThreshold & myDataFiltChrIIStatusZero$ref == ref2)) / nSamp, nDecimal), "\n", file = fileStat, append = T)

# numero di LOH with 1:5 markers
cat("# LOH with 1-5 markers & mean marker distance", ref1, "\n", file = fileStat, append = T)
for (indN in 1:nMarkerThreshold) {
  indbon <- which(myDataFiltChrIIStatusZero$len == indN & myDataFiltChrIIStatusZero$ref == ref1)
  cat(signif(length(indbon) / nSamp, nDecimal), "\t", signif(mean(myDataFiltChrIIStatusZero$evrLF[indbon]), nDecimal), "\n", file = fileStat, append = T)
}
cat("# LOH with 1-5 markers & mean marker distance", ref2, "\n", file = fileStat, append = T)
for (indN in 1:nMarkerThreshold) {
  indbon <- which(myDataFiltChrIIStatusZero$len == indN & myDataFiltChrIIStatusZero$ref == ref2)
  cat(signif(length(indbon) / nSamp, nDecimal), "\t", signif(mean(myDataFiltChrIIStatusZero$evrLF[indbon]), nDecimal), "\n", file = fileStat, append = T)
}

# terminal and interstitial
cat("Number of interstitial events: ", file = fileStat, append = T)
cat(signif(length(which(myDataFiltChrIIStatusZero$TI == "int")) / nSamp, nDecimal), "\n", file = fileStat, append = T)
cat("Number of terminal events: ", file = fileStat, append = T)
cat(signif(length(which(myDataFiltChrIIStatusZero$TI == "ter")) / nSamp, nDecimal), "\n", file = fileStat, append = T)

cat("Interstitial events in", ref1, "\n", file = fileStat, append = T)
cat(signif(length(which(myDataFiltChrIIStatusZero$TI == "int" & myDataFiltChrIIStatusZero$ref == ref1)) / nSamp, nDecimal), "\n", file = fileStat, append = T)
cat("Interstitial events <", distThreshold / 1000, "kb", ref1, "\n", file = fileStat, append = T)
cat(signif(length(which(myDataFiltChrIIStatusZero$TI == "int" & myDataFiltChrIIStatusZero$ref == ref1 & myDataFiltChrIIStatusZero$distLF <= distThreshold)) / nSamp, nDecimal), "\n", file = fileStat, append = T)
cat("Interstitial events >", distThreshold / 1000, "kb", ref1, "\n", file = fileStat, append = T)
cat(signif(length(which(myDataFiltChrIIStatusZero$TI == "int" & myDataFiltChrIIStatusZero$ref == ref1 & myDataFiltChrIIStatusZero$distLF > distThreshold)) / nSamp, nDecimal), "\n", file = fileStat, append = T)
cat("Terminal events in", ref1, "\n", file = fileStat, append = T)
cat(signif(length(which(myDataFiltChrIIStatusZero$TI == "ter" & myDataFiltChrIIStatusZero$ref == ref1)) / nSamp, nDecimal), "\n", file = fileStat, append = T)
cat("Terminal events <", distThreshold / 1000, "kb", ref1, "\n", file = fileStat, append = T)
cat(signif(length(which(myDataFiltChrIIStatusZero$TI == "ter" & myDataFiltChrIIStatusZero$ref == ref1 & myDataFiltChrIIStatusZero$distLF <= distThreshold)) / nSamp, nDecimal), "\n", file = fileStat, append = T)
cat("Terminal events >", distThreshold / 1000, "kb", ref1, "\n", file = fileStat, append = T)
cat(signif(length(which(myDataFiltChrIIStatusZero$TI == "ter" & myDataFiltChrIIStatusZero$ref == ref1 & myDataFiltChrIIStatusZero$distLF > distThreshold)) / nSamp, nDecimal), "\n", file = fileStat, append = T)

cat("Interstitial events in", ref2, "\n", file = fileStat, append = T)
cat(signif(length(which(myDataFiltChrIIStatusZero$TI == "int" & myDataFiltChrIIStatusZero$ref == ref2)) / nSamp, nDecimal), "\n", file = fileStat, append = T)
cat("Interstitial events <", distThreshold / 1000, "kb", ref2, "\n", file = fileStat, append = T)
cat(signif(length(which(myDataFiltChrIIStatusZero$TI == "int" & myDataFiltChrIIStatusZero$ref == ref2 & myDataFiltChrIIStatusZero$distLF <= distThreshold)) / nSamp, nDecimal), "\n", file = fileStat, append = T)
cat("Interstitial events >", distThreshold / 1000, "kb", ref2, "\n", file = fileStat, append = T)
cat(signif(length(which(myDataFiltChrIIStatusZero$TI == "int" & myDataFiltChrIIStatusZero$ref == ref2 & myDataFiltChrIIStatusZero$distLF > distThreshold)) / nSamp, nDecimal), "\n", file = fileStat, append = T)
cat("Terminal events in", ref2, "\n", file = fileStat, append = T)
cat(length(which(myDataFiltChrIIStatusZero$TI == "ter" & myDataFiltChrIIStatusZero$ref == ref2)) / nSamp, "\n", file = fileStat, append = T)
cat("Terminal events <", distThreshold / 1000, "kb", ref2, "\n", file = fileStat, append = T)
cat(signif(length(which(myDataFiltChrIIStatusZero$TI == "ter" & myDataFiltChrIIStatusZero$ref == ref2 & myDataFiltChrIIStatusZero$distLF <= distThreshold)) / nSamp, nDecimal), "\n", file = fileStat, append = T)
cat("Terminal events >", distThreshold / 1000, "kb", ref2, "\n", file = fileStat, append = T)
cat(signif(length(which(myDataFiltChrIIStatusZero$TI == "ter" & myDataFiltChrIIStatusZero$ref == ref2 & myDataFiltChrIIStatusZero$distLF > distThreshold)) / nSamp, nDecimal), "\n", file = fileStat, append = T)

# by chromosome
# misleading?
# a normalizzato con la lunghezza dell'altro cromosoma?
# cat("# LOH per chromosome per Mb", ref1, "\n", file = fileStat, append = T)
# for (indH in unique(myDataFiltChrIIStatusZero$chr)) {
#   numRef1 <- length(which(myDataFiltChrIIStatusZero$chr == indH & myDataFiltChrIIStatusZero$ref == ref1)) / nSamp / chrLenRef1$V3[which(chrLenRef1$V2 == paste(indH, ref1Label, sep = "_"))] * 1E6
#   cat(indH, "\t", signif(numRef1, nDecimal), "\n", file = fileStat, append = T)
# }
# cat("# LOH per chromosome per Mb", ref2, "\n", file = fileStat, append = T)
# for (indH in unique(myDataFiltChrIIStatusZero$chr)) {
#   numRef2 <- length(which(myDataFiltChrIIStatusZero$chr == indH & myDataFiltChrIIStatusZero$ref == ref2)) / nSamp / chrLenRef2$V3[which(chrLenRef2$V2 == paste(indH, ref2Label, sep = "_"))] * 1E6
#   cat(indH, "\t", signif(numRef2, nDecimal), "\n", file = fileStat, append = T)
# }

if (markerChrom == "-") {
  # save tables: rbind di le LOH (con chrII), i segmenti in eterozigosi (con chrII), le LOH filtrate perché presenti nei controlli (con chrII)
  outData <- rbind(myDataFiltChrIIStatusZero, myData[which(myData$status == 1), ], overlappedWithControlsSegments)
} else {
  # save tables: rbind di le LOH (senza chrII), i segmenti in eterozigosi (senza chrII), le LOH filtrate perché presenti nei controlli (senza chrII)
  outData <- rbind(myDataFiltChrIIStatusZero, myData[which(myData$status == 1 & myData$chr != markerChrom), ], overlappedWithControlsSegments)
}
# sort by chrom and "first" position
outData <- outData[order(match(outData$chr, allChr), outData$first), ]
# sort by parent
outData <- outData[order(match(outData$ref, unique(outData$ref))), ]
# sort by sample
outData <- outData[order(match(outData$sample, unique(outData$sample))), ]
write.table(x = outData, file = fileTable, col.names = T, row.names = F, quote = F, append = F, sep = "\t", na = ".")

# data in RData for TestRearranged
save(outData, file = fileForTestRearranged)

# plot number of interstitial VS chromosome length
counterRef <- 1
for (indR in c(ref1, ref2)) {
  chromData <- c()
  counterChr <- 1
  for (indC in allChr) {
    interstitialChrData <- myDataFiltChrIIStatusZero[which(myDataFiltChrIIStatusZero$ref == indR
                                                           & myDataFiltChrIIStatusZero$chr == indC
                                                           & myDataFiltChrIIStatusZero$TI == "int"), ]
    sampCount <- c()
    for (indS in mySamples) {
      sampCount <- c(sampCount, length(which(interstitialChrData$sample == indS)))
    }
    chromData <- rbind(chromData, data.frame(mean(sampCount), sd(sampCount), sum(sampCount)))
    counterChr <- counterChr + 1
  }
  plotIntPerChrData <- data.frame(chromData, allChrLen[[counterRef]][, 2:3], allChr)
  colnames(plotIntPerChrData) <- c("mean", "sd", "sum", "chr_parent", "len", "chrom")
  # sort chromosome labels by length
  plotIntPerChrData$chrom <- factor(plotIntPerChrData$chrom, levels = plotIntPerChrData$chrom[order(plotIntPerChrData$len)])
  # filter out chromosomes with no interstitial event
  plotIntPerChrData <- plotIntPerChrData[which(plotIntPerChrData$sum != 0), ]
  # switch to Mb
  plotIntPerChrData$len <- plotIntPerChrData$len / 1E6
  
  # plotting, labels: indR, nSamp
  # mean number of interstitial events, per sample and length
  
  plotTolo <- ggplot(plotIntPerChrData, aes(x = len, y = (mean / len), color = chrom)) +
    labs(color = "Chromosomes\nby length") +
    ylab("MDIE [1/Mb]") +
    xlab("length [Mb]") +
    # scale_x_continuous(limits = c(0, 1.5)) +
    geom_pointrange(aes(ymin = (mean - sd) / len, ymax = (mean + sd) / len)) +
    ggtitle(paste0("Mean Density of Interstitial Events - ", indR),
            subtitle = paste("per Sample per Mb by chromosome\n",
                             "N samples = ", nSamp,
                             sep = ""))
  
  fileOutPlot <- file.path(dirOutPlot, paste("MDIE", indR, "pdf", sep = "."))
  pdf(file = fileOutPlot)
  print(plotTolo)
  dev.off()
  
  counterRef <- counterRef + 1
}


