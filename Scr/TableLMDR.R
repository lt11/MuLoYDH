# header ------------------------------------------------------------------

# OCCHIO
# non ho controllato che nel dataframe nucmerData
# le posizioni sul ref1 corrispondano a quelle sul ref2
# ma quando faccio l'unione non dovrebbero venire fuori casini
# a meno che e.g. non ci una cosa tipo:
# ref1 vs ref2: chrI 123 133 chrI
# ref2 vs ref1: chrI 123 143 chrI
# ma non dovrebbe accadere
# se anche accadesse, ce ne sbattiamo tanto qui usiamo
# le posizioni su ref1 e ref2 separatamente

# calculate low marker density regions
# ovvero le regioni con marker distanti più di 300 bp
# usando l'unione dei dati ref1 vs ref2 
# e ref2 vs ref1

# la correzione degli alleli nelle inversioni non è
# necessaria perché usiamo solo la posizione e il
# cromosoma

rm(list = ls())
options(stringsAsFactors = F)

SpecDec <- function(x, k) as.numeric(format(round(x, k), nsmall = k))

# settings ----------------------------------------------------------------

# arguments
argsVal <- commandArgs(trailingOnly = T)
readLen <- as.numeric(argsVal[1])
baseDir <- argsVal[2]
# commentami
# baseDir <- "/Users/lorenzo/Prog/Yeasts/DistanceMarkerSNV"

distThreshold <- readLen * 2

outDir <- file.path(baseDir, "MUMmer", "LowDensityRegions")
dir.create(path = outDir, showWarnings = F, recursive = T)

mumDir <- file.path(baseDir, "MUMmer")
inputDir <- list.files(mumDir, pattern = "_", full.names = T)

allChr <- c("chrI", "chrII", "chrIII", "chrIV", "chrV", "chrVI", "chrVII", "chrVIII", 
            "chrIX", "chrX", "chrXI", "chrXII", "chrXIII", "chrXIV", "chrXV", "chrXVI")
colNames <- c("posR1", "allR1", "allR2", "posR2", "buff", "dist", "lenR1", "lenR2", "frm", "tags", "chromR1", "chromR2")

# runner file
inputFilePath <- list.files(path = file.path(baseDir, "Scr"), pattern = "^MuLoYDH.*sh$", full.names = T)
inputFile <- readLines(con = inputFilePath)

# lettura argmenti che non "entrano" nella riga di comando
# (per evitare di farla lunga un km)
ref1Label <- gsub("^\"|\"$", "", 
                  unlist(strsplit(grep(pattern = "Ref1Label=", x = inputFile, value = T)[1], split = "=", fixed = T))[2])
ref2Label <- gsub("^\"|\"$", "", 
                  unlist(strsplit(grep(pattern = "Ref2Label=", x = inputFile, value = T)[1], split = "=", fixed = T))[2])
bothRef <- c(ref1Label, ref2Label)

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

# init out summary files by chromosome, core regions
outSummaryFileR1 <- file.path(outDir, paste0(bothRef[1], ".Summary.LMDR.Core.Union.All.txt"))
headerString <- c("chrom", "LMDR [bp]", "length [bp]", "LMDR fraction", "\n")
cat(x = headerString, file = outSummaryFileR1, append = F, sep = "\t")
outSummaryFileR2 <- file.path(outDir, paste0(bothRef[2], ".Summary.LMDR.Core.Union.All.txt"))
headerString <- c("chrom", "LMDR [bp]", "length [bp]", "LMDR fraction", "\n")
cat(x = headerString, file = outSummaryFileR2, append = F, sep = "\t")
# init out summary files by chromosome, whole chromosome
outWholeChrSummaryFileR1 <- file.path(outDir, paste0(bothRef[1], ".Summary.LMDR.WholeChr.Union.All.txt"))
cat(x = headerString, file = outWholeChrSummaryFileR1, append = F, sep = "\t")
outWholeChrSummaryFileR2 <- file.path(outDir, paste0(bothRef[2], ".Summary.LMDR.WholeChr.Union.All.txt"))
cat(x = headerString, file = outWholeChrSummaryFileR2, append = F, sep = "\t")

# input data
dirRes12 <- file.path(baseDir, "MUMmer", paste(ref1Label, ref2Label, sep = "_"))
file12 <- list.files(dirRes12, pattern = "prt\\.all$", full.names = T)

dirRes21 <- file.path(baseDir, "MUMmer", paste(ref2Label, ref1Label, sep = "_"))
file21 <- list.files(dirRes21, pattern = "prt\\.all$", full.names = T)

data12 <- read.table(file12, header = F, skip = 4)
colnames(data12) <- colNames

data21 <- read.table(file21, header = F, skip = 4)
colnames(data21) <- c("posR2", "allR2", "allR1", "posR1", "buff", "dist", "lenR2", "lenR1", 
                      "frm", "tags", "chromR2", "chromR1")

# collassa indels che sono riportate per esteso, tipo:
# ref1 -> ref2
# 8515	A	.	1581	1	1581	237661	215496	1	1	chrI	chrI
# 8516	T	.	1581	1	1581	237661	215496	1	1	chrI	chrI
# 8517	T	.	1581	1	1581	237661	215496	1	1	chrI	chrI
# 8518	A	.	1581	1	1581	237661	215496	1	1	chrI	chrI
# ref2 -> ref1
# 1581	.	A	8515	0	1581	215496	237661	1	1	chrI	chrI
# 1581	.	T	8516	0	1581	215496	237661	1	1	chrI	chrI
# 1581	.	T	8517	0	1581	215496	237661	1	1	chrI	chrI
# 1581	.	A	8518	0	1581	215496	237661	1	1	chrI	chrI

# collassa le indel nucmer ref1 vs ref2 -----------------------------------

# si parte con:
# posR1 allR1 allR2 posR2 buff  dist  lenR1  lenR2 frm tags chromR1 chromR2
# 773     .     T  7708    2   773 215496 237661   1    1    chrI    chrI
# 779     C     .  7713    4   779 215496 237661   1    1    chrI    chrI
# 860     .     T  7795    8   860 215496 237661   1    1    chrI    chrI
# 868     A     .  7802    8   868 215496 237661   1    1    chrI    chrI
# 1486     A     .  8419    9  1486 215496 237661   1    1    chrI    chrI
# 1520     .     T  8454   10  1520 215496 237661   1    1    chrI    chrI
# 1540     C     .  8473    3  1540 215496 237661   1    1    chrI    chrI
# 1581     .     A  8515    0  1581 215496 237661   1    1    chrI    chrI
# 1581     .     T  8516    0  1581 215496 237661   1    1    chrI    chrI
# 1581     .     T  8517    0  1581 215496 237661   1    1    chrI    chrI
# 1581     .     A  8518    0  1581 215496 237661   1    1    chrI    chrI
# 1581     .     C  8519    0  1581 215496 237661   1    1    chrI    chrI
# 1581     .     C  8520    0  1581 215496 237661   1    1    chrI    chrI
# 1581     .     G  8521    0  1581 215496 237661   1    1    chrI    chrI
# 1581     .     T  8522    0  1581 215496 237661   1    1    chrI    chrI
# 1581     .     C  8523    0  1581 215496 237661   1    1    chrI    chrI
# 1608     G     .  8549    5  1608 215496 237661   1    1    chrI    chrI
# 1613     T     .  8553    5  1613 215496 237661   1    1    chrI    chrI
# 1725     .     A  8666    5  1725 215496 237661   1    1    chrI    chrI
# 1730     .     T  8672    5  1730 215496 237661   1    1    chrI    chrI
# 1746     .     G  8689    4  1746 215496 237661   1    1    chrI    chrI
# 1785     .     A  8729    0  1785 215496 237661   1    1    chrI    chrI
# 1785     .     T  8730    0  1785 215496 237661   1    1    chrI    chrI
# 1785     .     T  8731    0  1785 215496 237661   1    1    chrI    chrI
# 2223     T     .  9168    4  2223 215496 237661   1    1    chrI    chrI
# 2235     .     G  9181    0  2235 215496 237661   1    1    chrI    chrI
# 2235     .     G  9182    0  2235 215496 237661   1    1    chrI    chrI
# 2235     .     G  9183    0  2235 215496 237661   1    1    chrI    chrI
# 2239     G     .  9186    1  2239 215496 237661   1    1    chrI    chrI
# 2240     T     .  9186    1  2240 215496 237661   1    1    chrI    chrI
# 2241     T     .  9186    1  2241 215496 237661   1    1    chrI    chrI
# 2353     .     T  9299    2  2353 215496 237661   1    1    chrI    chrI
# 2656     A     .  9601    4  2656 215496 237661   1    1    chrI    chrI
# 2667     A     .  9611    1  2667 215496 237661   1    1    chrI    chrI
# 2668     G     .  9611    1  2668 215496 237661   1    1    chrI    chrI
# 3887     A     . 10829    1  3887 215496 237661   1    1    chrI    chrI

# con il primo passaggio (indBon12Ref1) si arriva a:
# posR1 allR1 allR2 posR2 buff  dist  lenR1  lenR2 frm tags chromR1 chromR2
# 773     .     T  7708    2   773 215496 237661   1    1    chrI    chrI
# 779     C     .  7713    4   779 215496 237661   1    1    chrI    chrI
# 860     .     T  7795    8   860 215496 237661   1    1    chrI    chrI
# 868     A     .  7802    8   868 215496 237661   1    1    chrI    chrI
# 1486     A     .  8419    9  1486 215496 237661   1    1    chrI    chrI
# 1520     .     T  8454   10  1520 215496 237661   1    1    chrI    chrI
# 1540     C     .  8473    3  1540 215496 237661   1    1    chrI    chrI
# 1581     .     A  8515    0  1581 215496 237661   1    1    chrI    chrI
# 1608     G     .  8549    5  1608 215496 237661   1    1    chrI    chrI
# 1613     T     .  8553    5  1613 215496 237661   1    1    chrI    chrI
# 1725     .     A  8666    5  1725 215496 237661   1    1    chrI    chrI
# 1730     .     T  8672    5  1730 215496 237661   1    1    chrI    chrI
# 1746     .     G  8689    4  1746 215496 237661   1    1    chrI    chrI
# 1785     .     A  8729    0  1785 215496 237661   1    1    chrI    chrI
# 2223     T     .  9168    4  2223 215496 237661   1    1    chrI    chrI
# 2235     .     G  9181    0  2235 215496 237661   1    1    chrI    chrI
# 2239     G     .  9186    1  2239 215496 237661   1    1    chrI    chrI
# 2240     T     .  9186    1  2240 215496 237661   1    1    chrI    chrI
# 2241     T     .  9186    1  2241 215496 237661   1    1    chrI    chrI
# 2353     .     T  9299    2  2353 215496 237661   1    1    chrI    chrI
# 2656     A     .  9601    4  2656 215496 237661   1    1    chrI    chrI
# 2667     A     .  9611    1  2667 215496 237661   1    1    chrI    chrI
# 2668     G     .  9611    1  2668 215496 237661   1    1    chrI    chrI
# 3887     A     . 10829    1  3887 215496 237661   1    1    chrI    chrI
# 3888     C     . 10829    1  3888 215496 237661   1    1    chrI    chrI
# 3889     A     . 10829    1  3889 215496 237661   1    1    chrI    chrI

# e con il secondo (indBon12Ref2):
# posR1 allR1 allR2 posR2 buff dist  lenR1  lenR2 frm tags chromR1 chromR2
# 773     .     T  7708    2  773 215496 237661   1    1    chrI    chrI
# 779     C     .  7713    4  779 215496 237661   1    1    chrI    chrI
# 860     .     T  7795    8  860 215496 237661   1    1    chrI    chrI
# 868     A     .  7802    8  868 215496 237661   1    1    chrI    chrI
# 1486     A     .  8419    9 1486 215496 237661   1    1    chrI    chrI
# 1520     .     T  8454   10 1520 215496 237661   1    1    chrI    chrI
# 1540     C     .  8473    3 1540 215496 237661   1    1    chrI    chrI
# 1581     .     A  8515    0 1581 215496 237661   1    1    chrI    chrI
# 1608     G     .  8549    5 1608 215496 237661   1    1    chrI    chrI
# 1613     T     .  8553    5 1613 215496 237661   1    1    chrI    chrI
# 1725     .     A  8666    5 1725 215496 237661   1    1    chrI    chrI
# 1730     .     T  8672    5 1730 215496 237661   1    1    chrI    chrI
# 1746     .     G  8689    4 1746 215496 237661   1    1    chrI    chrI
# 1785     .     A  8729    0 1785 215496 237661   1    1    chrI    chrI
# 2223     T     .  9168    4 2223 215496 237661   1    1    chrI    chrI
# 2235     .     G  9181    0 2235 215496 237661   1    1    chrI    chrI
# 2239     G     .  9186    1 2239 215496 237661   1    1    chrI    chrI
# 2353     .     T  9299    2 2353 215496 237661   1    1    chrI    chrI
# 2656     A     .  9601    4 2656 215496 237661   1    1    chrI    chrI
# 2667     A     .  9611    1 2667 215496 237661   1    1    chrI    chrI
# 3887     A     . 10829    1 3887 215496 237661   1    1    chrI    chrI
# 3892     C     . 10831    3 3892 215496 237661   1    1    chrI    chrI

# estrai i SNM: sono tantini, non perdiamo tempo a fare i conti con loro
data12SNM <- data12[which(data12$allR1 != "." & data12$allR2 != "."), ]
# estrai indel
data12indel <- data12[which(data12$allR1 == "." | data12$allR2 == "."), ]

# posizioni genomiche indel ref1
str12PosChrRef1 <- paste(data12indel$posR1, data12indel$chromR1, sep = "_")
# posizioni genomiche start delezioni ref1 e tutte le posizioni insert
str12PosChrRef1Unq <- unique(str12PosChrRef1)
# indici start delezioni ref1 e tutte le posizioni insert
indBon12Ref1 <- match(str12PosChrRef1Unq, str12PosChrRef1)

# posizioni genomiche indel ref2
str12PosChrRef2 <- paste(data12indel$posR2, data12indel$chromR2, sep = "_")
# posizioni genomiche start delezioni ref2 e tutte le posizioni insert
str12PosChrRef2Unq <- unique(str12PosChrRef2)
# indici start delezioni ref2 e tutte le posizioni insert
indBon12Ref2 <- match(str12PosChrRef2Unq, str12PosChrRef2)

# con intersect restano solo gli start delle del (su ref1 e su ref2)
indelCollapsed12 <- data12indel[intersect(indBon12Ref1, indBon12Ref2), ]

# collassa le indel nucmer ref2 vs ref1 -----------------------------------

# estrai i SNM: sono tantini, non perdiamo tempo a fare i conti con loro
data21SNM <- data21[which(data21$allR1 != "." & data21$allR2 != "."), ]
# estrai indel
data21indel <- data21[which(data21$allR1 == "." | data21$allR2 == "."), ]

# posizioni genomiche indel ref1
str21PosChrRef1 <- paste(data21indel$posR1, data21indel$chromR1, sep = "_")
# posizioni genomiche start delezioni ref1 e tutte le posizioni insert
str21PosChrRef1Unq <- unique(str21PosChrRef1)
# indici start delezioni ref1 e tutte le posizioni insert
indBon21Ref1 <- match(str21PosChrRef1Unq, str21PosChrRef1)

# posizioni genomiche indel ref2
str21PosChrRef2 <- paste(data21indel$posR2, data21indel$chromR2, sep = "_")
# posizioni genomiche start delezioni ref2 e tutte le posizioni insert
str21PosChrRef2Unq <- unique(str21PosChrRef2)
# indici start delezioni ref2 e tutte le posizioni insert
indBon21Ref2 <- match(str21PosChrRef2Unq, str21PosChrRef2)

# con intersect restano solo gli start delle del (su ref1 e su ref2)
indelCollapsed21 <- data21indel[intersect(indBon21Ref1, indBon21Ref2), ]

# i df con le indel collassate

data12IndColl <- rbind(data12SNM, indelCollapsed12)
data21IndColl <- rbind(data21SNM, indelCollapsed21)

strRef1Dir12 <- paste(data12IndColl$posR1, data12IndColl$chromR1, sep = "_")
strRef1Dir21 <- paste(data21IndColl$posR1, data21IndColl$chromR1, sep = "_")
ref1Union <- union(strRef1Dir12, strRef1Dir21)

# length(strRef1Dir12)
# length(strRef1Dir21)
# length(ref1Union)

strRef2Dir12 <- paste(data12IndColl$posR2, data12IndColl$chromR2, sep = "_")
strRef2Dir21 <- paste(data21IndColl$posR2, data21IndColl$chromR2, sep = "_")
ref2Union <- union(strRef1Dir12, strRef1Dir21)

# length(strRef2Dir12)
# length(strRef2Dir21)
# length(ref2Union)

# df dell'unione dei dati di nucmer
posR1 <- as.numeric(sapply(strsplit(ref1Union, split = "_"), "[[", 1))
chromR1 <- sapply(strsplit(ref1Union, split = "_"), "[[", 2)
posR2 <- as.numeric(sapply(strsplit(ref2Union, split = "_"), "[[", 1))
chromR2 <- sapply(strsplit(ref2Union, split = "_"), "[[", 2)
nucmerData <- data.frame(chromR1, posR1, chromR2, posR2)

# remove MT markers
indNoMT <- which(nucmerData$chromR1 != "chrMT" & nucmerData$chromR2 != "chrMT")
nucmerDataFilt <- nucmerData[indNoMT, ]

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
  
  # whole chromosome --------------------------------------------------------
  
  # concatena 1, posizioni, length chromosome in posRefR1
  posRefR1SortedExtreme <- c(1, posRefR1Sorted, refR1ChromLen[counterChr, 2])
  # distance
  diffPosRefR1SortedExtreme <- diff(posRefR1SortedExtreme)
  # lunghezza regioni con marker troppo lontani
  tooFarReg <- diffPosRefR1SortedExtreme[which(diffPosRefR1SortedExtreme > distThreshold)]
  # correggo le lunghezza considerando le read che dagli estremi si allungano dentro la regione
  # così sovrastimo la correzione di readLen per le regioni (sub)telomeriche
  tooFarReg <- tooFarReg - 2*(readLen - 1)
  
  if (length(tooFarReg) != 0) {
    sumLowDenReg <- sum(tooFarReg)
    wChrLenR1 <- refR1ChromLen[counterChr, 2]
    outSummaryR1 <- data.frame(indC, sumLowDenReg, wChrLenR1, 
                               SpecDec(x = sumLowDenReg / wChrLenR1, 5))   
    write.table(x = outSummaryR1, file = outWholeChrSummaryFileR1, 
                append = T, row.names = F, col.names = F, quote = F, sep = "\t")
  } else {
    sumLowDenReg <- 0
    outSummaryR1 <- data.frame(indC, sumLowDenReg, wChrLenR1, 
                               SpecDec(x = sumLowDenReg / wChrLenR1, 5))
    write.table(x = outSummaryR1, file = outWholeChrSummaryFileR1, 
                append = T, row.names = F, col.names = F, quote = F, sep = "\t")
  }
  
  # core chromosome --------------------------------------------------------
  
  # concatena posizioni core
  posRefR1SortedExtreme <- posRefR1Sorted[which(posRefR1Sorted > lenSubTelR1[counterChr, 1] | 
                                                  posRefR1Sorted < lenSubTelR1[counterChr, 2])]
  # distance
  diffPosRefR1SortedExtreme <- diff(posRefR1SortedExtreme)
  # lunghezza regioni con marker troppo lontani
  tooFarReg <- diffPosRefR1SortedExtreme[which(diffPosRefR1SortedExtreme > distThreshold)]
  # correggo le lunghezza considerando le read che dagli estremi si allungano dentro la regione
  tooFarReg <- tooFarReg - 2*(readLen - 1)
  
  if (length(tooFarReg) != 0) {
    sumLowDenReg <- sum(tooFarReg)
    coreLenR1 <- refR1ChromLen[counterChr, 2] - 
      (lenSubTelR1[counterChr, 1] + refR1ChromLen[counterChr, 2] - lenSubTelR1[counterChr, 2])
    outSummaryR1 <- data.frame(indC, sumLowDenReg, coreLenR1, 
                               SpecDec(x = sumLowDenReg / coreLenR1, 5))   
    write.table(x = outSummaryR1, file = outSummaryFileR1, 
                append = T, row.names = F, col.names = F, quote = F, sep = "\t")
  } else {
    sumLowDenReg <- 0
    coreLenR1 <- refR1ChromLen[counterChr, 2] - 
      (lenSubTelR1[counterChr, 1] + refR1ChromLen[counterChr, 2] - lenSubTelR1[counterChr, 2])
    outSummaryR1 <- data.frame(indC, sumLowDenReg, coreLenR1, 
                               SpecDec(x = sumLowDenReg / coreLenR1, 5))
    write.table(x = outSummaryR1, file = outSummaryFileR1, 
                append = T, row.names = F, col.names = F, quote = F, sep = "\t")
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
  
  # whole chromosome --------------------------------------------------------
  
  # concatena 1, posizioni, length chromosome in posRefR2
  posRefR2SortedExtreme <- c(1, posRefR2Sorted, refR2ChromLen[counterChr, 2])
  # distance
  diffPosRefR2SortedExtreme <- diff(posRefR2SortedExtreme)
  # lunghezza regioni con marker troppo lontani
  tooFarReg <- diffPosRefR2SortedExtreme[which(diffPosRefR2SortedExtreme > distThreshold)]
  # correggo le lunghezza considerando le read che dagli estremi si allungano dentro la regione
  #  così sovrastimo la correzione di readLen per le regioni (sub)telomeriche
  tooFarReg <- tooFarReg - 2*(readLen - 1)
  
  if (length(tooFarReg) != 0) {
    sumLowDenReg <- sum(tooFarReg)
    wChrLenR2 <- refR2ChromLen[counterChr, 2]
    outSummaryR2 <- data.frame(indC, sumLowDenReg, wChrLenR2, 
                               SpecDec(x = sumLowDenReg / wChrLenR2, 5))   
    write.table(x = outSummaryR2, file = outWholeChrSummaryFileR2, 
                append = T, row.names = F, col.names = F, quote = F, sep = "\t")
  } else {
    sumLowDenReg <- 0
    outSummaryR2 <- data.frame(indC, sumLowDenReg, wChrLenR2, 
                               SpecDec(x = sumLowDenReg / wChrLenR2, 5))
    write.table(x = outSummaryR2, file = outWholeChrSummaryFileR2, 
                append = T, row.names = F, col.names = F, quote = F, sep = "\t")
  }
  
  # core chromosome --------------------------------------------------------
  
  # concatena posizioni core
  posRefR2SortedExtreme <- posRefR2Sorted[which(posRefR2Sorted > lenSubTelR2[counterChr, 1] | 
                                                  posRefR2Sorted < lenSubTelR2[counterChr, 2])]
  # distance
  diffPosRefR2SortedExtreme <- diff(posRefR2SortedExtreme)
  # lunghezza regioni con marker troppo lontani
  tooFarReg <- diffPosRefR2SortedExtreme[which(diffPosRefR2SortedExtreme > distThreshold)]
  # correggo le lunghezza considerando le read che dagli estremi si allungano dentro la regione
  tooFarReg <- tooFarReg - 2*(readLen - 1)
  
  if (length(tooFarReg) != 0) {
    sumLowDenReg <- sum(tooFarReg)
    coreLenR2 <- refR2ChromLen[counterChr, 2] - 
      (lenSubTelR2[counterChr, 1] + refR2ChromLen[counterChr, 2] - lenSubTelR2[counterChr, 2])
    outSummaryR2 <- data.frame(indC, sumLowDenReg, coreLenR2, 
                               SpecDec(x = sumLowDenReg / coreLenR2, 5))   
    write.table(x = outSummaryR2, file = outSummaryFileR2, 
                append = T, row.names = F, col.names = F, quote = F, sep = "\t")
  } else {
    sumLowDenReg <- 0
    coreLenR2 <- refR2ChromLen[counterChr, 2] - 
      (lenSubTelR2[counterChr, 1] + refR2ChromLen[counterChr, 2] - lenSubTelR2[counterChr, 2])
    outSummaryR2 <- data.frame(indC, sumLowDenReg, coreLenR2, 
                               SpecDec(x = sumLowDenReg / coreLenR2, 5))
    write.table(x = outSummaryR2, file = outSummaryFileR2, 
                append = T, row.names = F, col.names = F, quote = F, sep = "\t")
  }
  
  # increment chromosome
  counterChr <- counterChr + 1
} # chromosome loop


