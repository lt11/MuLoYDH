# header ------------------------------------------------------------------

# annotate summary table with interstitial/terminal labels
# make plots for interstitial and terminal events

rm(list = ls())
options(stringsAsFactors = F)
library(ggplot2)

# settings ----------------------------------------------------------------

# arguments
argsVal <- commandArgs(trailingOnly = T)
ref1Label <- argsVal[1]
ref2Label <- argsVal[2]
ref1 <- argsVal[3]
ref2 <- argsVal[4]
baseDir <- argsVal[5]

# paths
dirData <- baseDir
dirName <- file.path("LOH", "AllSegments")
annoDir <- file.path("LOH", "Annotations")
summaryFile <- "Annotations.SNV.txt"
outFileSummary <- "Annotations.SNV.IntTer.txt"
outDirSummary <- file.path(dirData, annoDir)
allEventsFolder <- file.path(dirData, dirName)
allEventsTable <- file.path(allEventsFolder, "AllEvents.csv")

binW <- 25

outDirHis <- file.path(dirData, dirName, paste("Bin", binW, "_HistoInterTerminal", sep = ""))
dir.create(path = outDirHis, showWarnings = F, recursive = T)

# init data
allChr <- c("chrI", "chrII", "chrIII", "chrIV", "chrV", "chrVI", "chrVII", "chrVIII", 
            "chrIX", "chrX", "chrXI", "chrXII", "chrXIII", "chrXIV", "chrXV", "chrXVI")
outTable <- "AllEvents.ZeroState.TI.RData"
strainLab <- c(ref1, ref2)
strainLabel <- c(ref1Label, ref2Label)
indTerminal <- c()

# load all events table
allEvents <- read.table(file = allEventsTable, sep = ",", header = T, na.strings = "antibodi")
# filter out non-0 state segments
indZeroState <- which(allEvents$status == 0)
allEventsZeroState <- allEvents[indZeroState, ]
rm(allEvents)
rownames(allEventsZeroState) <- seq(nrow(allEventsZeroState))

# by strain
for (indStrain in c(1:2)) {
  indEventsStrain <- which(allEventsZeroState$strain == strainLab[indStrain])
  allEventsZeroStateStrain <- allEventsZeroState[indEventsStrain, ]
  # terminal events
  indTermStart <- which(allEventsZeroStateStrain$start == 1)
  # load strain chromosome lengths
  chromLen <- read.table(file = file.path(baseDir, "CNV", "GCdata", strainLabel[indStrain], "LenChr.txt"), header = F, sep = "\t")[, 3]
  appCharEnd <- paste(allEventsZeroStateStrain$chr, allEventsZeroStateStrain$end, sep = "")
  intersectTermEnd <- intersect(paste(allChr, chromLen, sep = ""), appCharEnd)
  
  booTermEnd <- is.element(appCharEnd, intersectTermEnd)
  indTermEnd <- seq(length(booTermEnd))[booTermEnd]
  
  allEventsZeroStateStrainTerminal <- allEventsZeroStateStrain[sort(c(indTermStart, indTermEnd)), ]
  
  # terminal event lines in allEventsZeroState
  indTerminal <- c(indTerminal, as.numeric(rownames(allEventsZeroStateStrainTerminal)))
}

# add tag terminal
tagTerminalIterstitial <- rep("int", length = nrow(allEventsZeroState))
tagTerminalIterstitial[indTerminal] <- "ter"
allEventsZeroState <- data.frame(allEventsZeroState, tagTerminalIterstitial)
colnames(allEventsZeroState)[which(colnames(allEventsZeroState) == "tagTerminalIterstitial")] <- "TI"

# convertion to kb
allEventsZeroState$distES <- allEventsZeroState$distES / 1000
allEventsZeroState$distLF <- allEventsZeroState$distLF / 1000

save(allEventsZeroState, file = file.path(allEventsFolder, outTable))

# annotate TI to summary table and remove original summary table "Annotations.SNV.txt"
annoFile <- file.path(dirData, annoDir, summaryFile)
dfAnno <- read.table(file = annoFile, header = T, sep = "\t", na.strings = "antani")
unlink(annoFile)

m1 <- paste(dfAnno$sample, dfAnno$ref, dfAnno$chr, dfAnno$start, sep = "")
TI <- rep("-", length(m1))
terEvZeroAnn <- allEventsZeroState[which(allEventsZeroState$TI == "ter"), ]
m2 <- paste(terEvZeroAnn$hybrid, terEvZeroAnn$strain, terEvZeroAnn$chr, terEvZeroAnn$start, sep = "")
TI[is.element(m1, m2)] <- "ter"

intEvZeroAnn <- allEventsZeroState[which(allEventsZeroState$TI == "int"), ]
m2 <- paste(intEvZeroAnn$hybrid, intEvZeroAnn$strain, intEvZeroAnn$chr, intEvZeroAnn$start, sep = "")
TI[is.element(m1, m2)] <- "int"

allAnno <- data.frame(dfAnno, TI)
write.table(x = allAnno, file = file.path(outDirSummary, outFileSummary), sep = "\t", 
            quote = F, row.names = F)

# plotting ----------------------------------------------------------------

# RefA plots ----------------------------------------------------------------

# terminal events plots ------------------------------------------------------
plotEvents <- allEventsZeroState[which(allEventsZeroState$strain == ref1 & 
                                         allEventsZeroState$status == 0 & 
                                         allEventsZeroState$TI == "ter"), ]

nEv <- nrow(plotEvents)

# plot histogram with kernel density: terminal LF
pLFs <- ggplot(plotEvents, aes(x = distLF)) + 
  ggtitle(paste0("First-Last Marker Distance, ", ref1 , " strain, Terminal"), 
          subtitle = paste("N = ", nEv, "; ", 
                           "Bin = ", binW, " kb", 
                           sep = "")) + 
  geom_histogram(aes(y = ..density..), # histogram with density instead of count on y-axis
                 binwidth = binW, 
                 colour = "black", fill = "white") + 
  geom_density(alpha = .25, fill = "blue") + # overlay with transparent kernel density plot
  xlab("Length [kb]") + 
  ylab("Density") + 
  theme(
    # title style
    plot.title = element_text(size = 16, hjust = 1), 
    plot.subtitle = element_text(size = 12, face = "bold", hjust = 1), 
    # axis style
    axis.title = element_text(size = 16, face = "bold"), 
    axis.text = element_text(size = 16))
pPath <- file.path(outDirHis, paste0(ref1, ".LF.terminal.pdf"))
pdf(file = pPath)
print(pLFs)
dev.off()

# plot histogram: terminal LF with count of events per bin (no density)
pLFs <- ggplot(plotEvents, aes(x = distLF)) + 
  ggtitle(paste0("First-Last Marker Distance, ", ref1 , " strain, Terminal"), 
          subtitle = paste("N = ", nEv, "; ", 
                           "Bin = ", binW, " kb", 
                           sep = "")) + 
  geom_histogram(binwidth = binW, 
                 colour = "black", fill = "white") + 
  geom_density(alpha = .25, fill = "blue") + # overlay with transparent kernel density plot
  xlab("Length [kb]") + 
  ylab("Count") + 
  theme(
    # title style
    plot.title = element_text(size = 16, hjust = 1), 
    plot.subtitle = element_text(size = 12, face = "bold", hjust = 1), 
    # axis style
    axis.title = element_text(size = 16, face = "bold"), 
    axis.text = element_text(size = 16))
pPath <- file.path(outDirHis, paste0(ref1, ".LF.terminal.count.pdf"))
pdf(file = pPath)
print(pLFs)
dev.off()

# small events threshold [kb] & parameters
smallThres <- 25
binWlarge <- 25
binWsmall <- .5

indSmall <- which(plotEvents$distLF < smallThres)
plotSmallEvents <- plotEvents[indSmall, ]
plotLargeEvents <- plotEvents[-indSmall, ]

# small events plots ------------------------------------------------------
nEv <- nrow(plotSmallEvents)

# plot histogram with kernel density: small LF
pLFs <- ggplot(plotSmallEvents, aes(x = distLF)) + 
  ggtitle(paste0("First-Last Marker Distance, ", ref1 , " strain, Terminal"), 
          subtitle = paste("Length < ", smallThres, " kb; ", 
                           "N = ", nEv, "; ", 
                           "Bin = ", binWsmall, " kb", 
                           sep = "")) + 
  geom_histogram(aes(y = ..density..), # histogram with density instead of count on y-axis
                 binwidth = binWsmall, 
                 colour = "black", fill = "white") + 
  geom_density(alpha = .25, fill = "blue") + # overlay with transparent kernel density plot
  xlab("Length [kb]") + 
  ylab("Density") + 
  theme(
    # title style
    plot.title = element_text(size = 16, hjust = 1), 
    plot.subtitle = element_text(size = 12, face = "bold", hjust = 1), 
    # axis style
    axis.title = element_text(size = 16, face = "bold"), 
    axis.text = element_text(size = 16))
pPath <- file.path(outDirHis, paste0(ref1, ".LF.small.terminal.pdf"))
pdf(file = pPath)
print(pLFs)
dev.off()

# large events plots ------------------------------------------------------
nEv <- nrow(plotLargeEvents)

# plot histogram with kernel density: large LF
pFLl <- ggplot(plotLargeEvents, aes(x = distLF)) + 
  ggtitle(paste0("First-Last Marker Distance, ", ref1 , " strain, Terminal"), 
          # unicode character ≥ \u2265
          subtitle = paste("Length > ", smallThres, " kb; ", 
                           "N = ", nEv, "; ", 
                           "Bin = ", binWlarge, " kb", 
                           sep = "")) + 
  geom_histogram(aes(y = ..density..), # histogram with density instead of count on y-axis
                 binwidth = binWlarge, 
                 colour = "black", fill = "white") + 
  geom_density(alpha = .25, fill = "blue") + # overlay with transparent kernel density plot
  xlab("Length [kb]") + 
  ylab("Density") + 
  theme(
    # title style
    plot.title = element_text(size = 16, hjust = 1), 
    plot.subtitle = element_text(size = 12, face = "bold", hjust = 1), 
    # axis style
    axis.title = element_text(size = 16, face = "bold"), 
    axis.text = element_text(size = 16))
pPath <- file.path(outDirHis, paste0(ref1, ".LF.large.terminal.pdf"))
pdf(file = pPath)
print(pFLl)
dev.off()


# interstitial events plots ------------------------------------------------------
plotEvents <- allEventsZeroState[which(allEventsZeroState$strain == ref1 & 
                                         allEventsZeroState$status == 0 & 
                                         allEventsZeroState$TI == "int"), ]

nEv <- nrow(plotEvents)

# plot histogram with kernel density: interstitial LF
pLFs <- ggplot(plotEvents, aes(x = distLF)) + 
  ggtitle(paste0("First-Last Marker Distance, ", ref1 , " strain, Interstitial"), 
          subtitle = paste("N = ", nEv, "; ", 
                           "Bin = ", binW, " kb", 
                           sep = "")) + 
  geom_histogram(aes(y = ..density..), # histogram with density instead of count on y-axis
                 binwidth = binW, 
                 colour = "black", fill = "white") + 
  geom_density(alpha = .25, fill = "blue") + # overlay with transparent kernel density plot
  xlab("Length [kb]") + 
  ylab("Density") + 
  theme(
    # title style
    plot.title = element_text(size = 16, hjust = 1), 
    plot.subtitle = element_text(size = 12, face = "bold", hjust = 1), 
    # axis style
    axis.title = element_text(size = 16, face = "bold"), 
    axis.text = element_text(size = 16))
pPath <- file.path(outDirHis, paste0(ref1, ".LF.interstitial.pdf"))
pdf(file = pPath)
print(pLFs)
dev.off()

# plot histogram: interstitial LF with count of events per bin (no density)
pLFs <- ggplot(plotEvents, aes(x = distLF)) + 
  ggtitle(paste0("First-Last Marker Distance, ", ref1 , " strain, Interstitial"), 
          subtitle = paste("N = ", nEv, "; ", 
                           "Bin = ", binW, " kb", 
                           sep = "")) + 
  geom_histogram(binwidth = binW, 
                 colour = "black", fill = "white") + 
  geom_density(alpha = .25, fill = "blue") + # overlay with transparent kernel density plot
  xlab("Length [kb]") + 
  ylab("Count") + 
  theme(
    # title style
    plot.title = element_text(size = 16, hjust = 1), 
    plot.subtitle = element_text(size = 12, face = "bold", hjust = 1), 
    # axis style
    axis.title = element_text(size = 16, face = "bold"), 
    axis.text = element_text(size = 16))
pPath <- file.path(outDirHis, paste0(ref1, ".LF.interstitial.count.pdf"))
pdf(file = pPath)
print(pLFs)
dev.off()

# small events threshold [kb] & parameters
smallThres <- 25
binWlarge <- 25
binWsmall <- .5

indSmall <- which(plotEvents$distLF < smallThres)
plotSmallEvents <- plotEvents[indSmall, ]
plotLargeEvents <- plotEvents[-indSmall, ]

# small events plots ------------------------------------------------------
nEv <- nrow(plotSmallEvents)

# plot histogram with kernel density: small LF
pLFs <- ggplot(plotSmallEvents, aes(x = distLF)) + 
  ggtitle(paste0("First-Last Marker Distance, ", ref1 , " strain, Interstitial"), 
          subtitle = paste("Length < ", smallThres, " kb; ", 
                           "N = ", nEv, "; ", 
                           "Bin = ", binWsmall, " kb", 
                           sep = "")) + 
  geom_histogram(aes(y = ..density..), # histogram with density instead of count on y-axis
                 binwidth = binWsmall, 
                 colour = "black", fill = "white") + 
  geom_density(alpha = .25, fill = "blue") + # overlay with transparent kernel density plot
  xlab("Length [kb]") + 
  ylab("Density") + 
  theme(
    # title style
    plot.title = element_text(size = 16, hjust = 1), 
    plot.subtitle = element_text(size = 12, face = "bold", hjust = 1), 
    # axis style
    axis.title = element_text(size = 16, face = "bold"), 
    axis.text = element_text(size = 16))
pPath <- file.path(outDirHis, paste0(ref1, ".LF.small.interstitial.pdf"))
pdf(file = pPath)
print(pLFs)
dev.off()

# large events plots ------------------------------------------------------
nEv <- nrow(plotLargeEvents)

# plot histogram with kernel density: large LF
pFLl <- ggplot(plotLargeEvents, aes(x = distLF)) + 
  ggtitle(paste0("First-Last Marker Distance, ", ref1 , " strain, Interstitial"), 
          # unicode character ≥ \u2265
          subtitle = paste("Length > ", smallThres, " kb; ", 
                           "N = ", nEv, "; ", 
                           "Bin = ", binWlarge, " kb", 
                           sep = "")) + 
  geom_histogram(aes(y = ..density..), # histogram with density instead of count on y-axis
                 binwidth = binWlarge, 
                 colour = "black", fill = "white") + 
  geom_density(alpha = .25, fill = "blue") + # overlay with transparent kernel density plot
  xlab("Length [kb]") + 
  ylab("Density") + 
  theme(
    # title style
    plot.title = element_text(size = 16, hjust = 1), 
    plot.subtitle = element_text(size = 12, face = "bold", hjust = 1), 
    # axis style
    axis.title = element_text(size = 16, face = "bold"), 
    axis.text = element_text(size = 16))
pPath <- file.path(outDirHis, paste0(ref1, ".LF.large.interstitial.pdf"))
pdf(file = pPath)
print(pFLl)
dev.off()

# RefB plots ----------------------------------------------------------------

# terminal events plots ------------------------------------------------------
plotEvents <- allEventsZeroState[which(allEventsZeroState$strain == ref2 & 
                                         allEventsZeroState$status == 0 & 
                                         allEventsZeroState$TI == "ter"), ]

nEv <- nrow(plotEvents)

# plot histogram with kernel density: terminal LF
pLFs <- ggplot(plotEvents, aes(x = distLF)) + 
  ggtitle(paste0("First-Last Marker Distance, ", ref2 , " strain, Terminal"), 
          subtitle = paste("N = ", nEv, "; ", 
                           "Bin = ", binW, " kb", 
                           sep = "")) + 
  geom_histogram(aes(y = ..density..), # histogram with density instead of count on y-axis
                 binwidth = binW, 
                 colour = "black", fill = "white") + 
  geom_density(alpha = .25, fill = "red") + # overlay with transparent kernel density plot
  xlab("Length [kb]") + 
  ylab("Density") + 
  theme(
    # title style
    plot.title = element_text(size = 16, hjust = 1), 
    plot.subtitle = element_text(size = 12, face = "bold", hjust = 1), 
    # axis style
    axis.title = element_text(size = 16, face = "bold"), 
    axis.text = element_text(size = 16))
pPath <- file.path(outDirHis, paste0(ref2, ".LF.terminal.pdf"))
pdf(file = pPath)
print(pLFs)
dev.off()

# plot histogram: terminal LF with count of events per bin (no density)
pLFs <- ggplot(plotEvents, aes(x = distLF)) + 
  ggtitle(paste0("First-Last Marker Distance, ", ref2 , " strain, Terminal"), 
          subtitle = paste("N = ", nEv, "; ", 
                           "Bin = ", binW, " kb", 
                           sep = "")) + 
  geom_histogram(binwidth = binW, 
                 colour = "black", fill = "white") + 
  geom_density(alpha = .25, fill = "red") + # overlay with transparent kernel density plot
  xlab("Length [kb]") + 
  ylab("Count") + 
  theme(
    # title style
    plot.title = element_text(size = 16, hjust = 1), 
    plot.subtitle = element_text(size = 12, face = "bold", hjust = 1), 
    # axis style
    axis.title = element_text(size = 16, face = "bold"), 
    axis.text = element_text(size = 16))
pPath <- file.path(outDirHis, paste0(ref2, ".LF.terminal.count.pdf"))
pdf(file = pPath)
print(pLFs)
dev.off()

# small events threshold [kb] & parameters
smallThres <- 25
binWlarge <- 25
binWsmall <- .5

indSmall <- which(plotEvents$distLF < smallThres)
plotSmallEvents <- plotEvents[indSmall, ]
plotLargeEvents <- plotEvents[-indSmall, ]

# small events plots ------------------------------------------------------
nEv <- nrow(plotSmallEvents)

# plot histogram with kernel density: small LF
pLFs <- ggplot(plotSmallEvents, aes(x = distLF)) + 
  ggtitle(paste0("First-Last Marker Distance, ", ref2 , " strain, Terminal"), 
          subtitle = paste("Length < ", smallThres, " kb; ", 
                           "N = ", nEv, "; ", 
                           "Bin = ", binWsmall, " kb", 
                           sep = "")) + 
  geom_histogram(aes(y = ..density..), # histogram with density instead of count on y-axis
                 binwidth = binWsmall, 
                 colour = "black", fill = "white") + 
  geom_density(alpha = .25, fill = "red") + # overlay with transparent kernel density plot
  xlab("Length [kb]") + 
  ylab("Density") + 
  theme(
    # title style
    plot.title = element_text(size = 16, hjust = 1), 
    plot.subtitle = element_text(size = 12, face = "bold", hjust = 1), 
    # axis style
    axis.title = element_text(size = 16, face = "bold"), 
    axis.text = element_text(size = 16))
pPath <- file.path(outDirHis, paste0(ref2, ".LF.small.terminal.pdf"))
pdf(file = pPath)
print(pLFs)
dev.off()

# large events plots ------------------------------------------------------
nEv <- nrow(plotLargeEvents)

# plot histogram with kernel density: large LF
pFLl <- ggplot(plotLargeEvents, aes(x = distLF)) + 
  ggtitle(paste0("First-Last Marker Distance, ", ref2 , " strain, Terminal"), 
          # unicode character ≥ \u2265
          subtitle = paste("Length > ", smallThres, " kb; ", 
                           "N = ", nEv, "; ", 
                           "Bin = ", binWlarge, " kb", 
                           sep = "")) + 
  geom_histogram(aes(y = ..density..), # histogram with density instead of count on y-axis
                 binwidth = binWlarge, 
                 colour = "black", fill = "white") + 
  geom_density(alpha = .25, fill = "red") + # overlay with transparent kernel density plot
  xlab("Length [kb]") + 
  ylab("Density") + 
  theme(
    # title style
    plot.title = element_text(size = 16, hjust = 1), 
    plot.subtitle = element_text(size = 12, face = "bold", hjust = 1), 
    # axis style
    axis.title = element_text(size = 16, face = "bold"), 
    axis.text = element_text(size = 16))
pPath <- file.path(outDirHis, paste0(ref2, ".LF.large.terminal.pdf"))
pdf(file = pPath)
print(pFLl)
dev.off()


# interstitial events plots ------------------------------------------------------
plotEvents <- allEventsZeroState[which(allEventsZeroState$strain == ref2 & 
                                         allEventsZeroState$status == 0 & 
                                         allEventsZeroState$TI == "int"), ]

nEv <- nrow(plotEvents)

# plot histogram with kernel density: interstitial LF
pLFs <- ggplot(plotEvents, aes(x = distLF)) + 
  ggtitle(paste0("First-Last Marker Distance, ", ref2 , " strain, Interstitial"), 
          subtitle = paste("N = ", nEv, "; ", 
                           "Bin = ", binW, " kb", 
                           sep = "")) + 
  geom_histogram(aes(y = ..density..), # histogram with density instead of count on y-axis
                 binwidth = binW, 
                 colour = "black", fill = "white") + 
  geom_density(alpha = .25, fill = "red") + # overlay with transparent kernel density plot
  xlab("Length [kb]") + 
  ylab("Density") + 
  theme(
    # title style
    plot.title = element_text(size = 16, hjust = 1), 
    plot.subtitle = element_text(size = 12, face = "bold", hjust = 1), 
    # axis style
    axis.title = element_text(size = 16, face = "bold"), 
    axis.text = element_text(size = 16))
pPath <- file.path(outDirHis, paste0(ref2, ".LF.interstitial.pdf"))
pdf(file = pPath)
print(pLFs)
dev.off()

# plot histogram: interstitial LF with count of events per bin (no density)
pLFs <- ggplot(plotEvents, aes(x = distLF)) + 
  ggtitle(paste0("First-Last Marker Distance, ", ref2 , " strain, Interstitial"), 
          subtitle = paste("N = ", nEv, "; ", 
                           "Bin = ", binW, " kb", 
                           sep = "")) + 
  geom_histogram(binwidth = binW, 
                 colour = "black", fill = "white") + 
  geom_density(alpha = .25, fill = "red") + # overlay with transparent kernel density plot
  xlab("Length [kb]") + 
  ylab("Count") + 
  theme(
    # title style
    plot.title = element_text(size = 16, hjust = 1), 
    plot.subtitle = element_text(size = 12, face = "bold", hjust = 1), 
    # axis style
    axis.title = element_text(size = 16, face = "bold"), 
    axis.text = element_text(size = 16))
pPath <- file.path(outDirHis, paste0(ref2, ".LF.interstitial.count.pdf"))
pdf(file = pPath)
print(pLFs)
dev.off()

# small events threshold [kb] & parameters
smallThres <- 25
binWlarge <- 25
binWsmall <- .5

indSmall <- which(plotEvents$distLF < smallThres)
plotSmallEvents <- plotEvents[indSmall, ]
plotLargeEvents <- plotEvents[-indSmall, ]

# small events plots ------------------------------------------------------
nEv <- nrow(plotSmallEvents)

# plot histogram with kernel density: small LF
pLFs <- ggplot(plotSmallEvents, aes(x = distLF)) + 
  ggtitle(paste0("First-Last Marker Distance, ", ref2 , " strain, Interstitial"), 
          subtitle = paste("Length < ", smallThres, " kb; ", 
                           "N = ", nEv, "; ", 
                           "Bin = ", binWsmall, " kb", 
                           sep = "")) + 
  geom_histogram(aes(y = ..density..), # histogram with density instead of count on y-axis
                 binwidth = binWsmall, 
                 colour = "black", fill = "white") + 
  geom_density(alpha = .25, fill = "red") + # overlay with transparent kernel density plot
  xlab("Length [kb]") + 
  ylab("Density") + 
  theme(
    # title style
    plot.title = element_text(size = 16, hjust = 1), 
    plot.subtitle = element_text(size = 12, face = "bold", hjust = 1), 
    # axis style
    axis.title = element_text(size = 16, face = "bold"), 
    axis.text = element_text(size = 16))
pPath <- file.path(outDirHis, paste0(ref2, ".LF.small.interstitial.pdf"))
pdf(file = pPath)
print(pLFs)
dev.off()

# large events plots ------------------------------------------------------
nEv <- nrow(plotLargeEvents)

# plot histogram with kernel density: large LF
pFLl <- ggplot(plotLargeEvents, aes(x = distLF)) + 
  ggtitle(paste0("First-Last Marker Distance, ", ref2 , " strain, Interstitial"), 
          # unicode character ≥ \u2265
          subtitle = paste("Length > ", smallThres, " kb; ", 
                           "N = ", nEv, "; ", 
                           "Bin = ", binWlarge, " kb", 
                           sep = "")) + 
  geom_histogram(aes(y = ..density..), # histogram with density instead of count on y-axis
                 binwidth = binWlarge, 
                 colour = "black", fill = "white") + 
  geom_density(alpha = .25, fill = "red") + # overlay with transparent kernel density plot
  xlab("Length [kb]") + 
  ylab("Density") + 
  theme(
    # title style
    plot.title = element_text(size = 16, hjust = 1), 
    plot.subtitle = element_text(size = 12, face = "bold", hjust = 1), 
    # axis style
    axis.title = element_text(size = 16, face = "bold"), 
    axis.text = element_text(size = 16))
pPath <- file.path(outDirHis, paste0(ref2, ".LF.large.interstitial.pdf"))
pdf(file = pPath)
print(pFLl)
dev.off()


