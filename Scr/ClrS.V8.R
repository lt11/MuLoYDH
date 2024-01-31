# header ------------------------------------------------------------------

# use intersection of standard nucmer
# all chromosomes at once
# plot density of markers calculated by nucmer (after (sub)telomeric positions filtering)
# filters multiallelic positions here: remove non-matching GTs (i 2 non vengono contemplati)
# filters for markers: (sub)telomeric, quality, genotype, allele, deletions
# save vcf data (before sorting for LOH segment calls!!!)
# calculate segments of consecutive genotypes

rm(list = ls())
options(stringsAsFactors = F)
library(ggplot2)
library(scales)

# function(s) -------------------------------------------------------------

SpecDec <- function(x, k) as.numeric(format(round(x, k), nsmall = k))

ptmGlob <- proc.time()
ptmInit <- proc.time()

# settings ----------------------------------------------------------------

# arguments
argsVal <- commandArgs(trailingOnly = T)
ref1Label <- argsVal[1]
ref2Label <- argsVal[2]
ref1 <- argsVal[3]
ref2 <- argsVal[4]
baseDir <- argsVal[5]

# input folders
polDir <- file.path(baseDir, "ParsedVariants", "Results", "Markers")
cnvDir <- file.path(baseDir, "CNV")
dirMUMmer <- file.path(baseDir, "MUMmer")
# output folders
resDir <- file.path(baseDir, "LOH")

# set bin width
bWidth <- 1.0

# marker density output folder
outDirMUM <- file.path(resDir, "Markers")
dir.create(path = outDirMUM, showWarnings = F, recursive = T)

allChr <- c("chrI", "chrII", "chrIII", "chrIV", "chrV", "chrVI", "chrVII", "chrVIII", "chrIX", 
            "chrX", "chrXI", "chrXII", "chrXIII", "chrXIV", "chrXV", "chrXVI")

# label tables marker density
outLabelA <- "Density.Table.FltSubTel"
# label plots
outLabelB <- "Distance.FltSubTel"

# centromeri, (sub)telomeri, chromosomes length
chrLenRefA <- read.table(file = file.path(baseDir, "CNV", "GCdata", ref1Label, "LenChr.txt"), header = F, sep = "\t")[, 3]
chrLenRefB <- read.table(file = file.path(baseDir, "CNV", "GCdata", ref2Label, "LenChr.txt"), header = F, sep = "\t")[, 3]

centrSRefA <- as.numeric(read.table(file = file.path(baseDir, "Ref", "Ann", 
                                                   paste0(ref1Label, ".centromere.txt")), header = F)[, 1])
centrERefA <- as.numeric(read.table(file = file.path(baseDir, "Ref", "Ann", 
                                                   paste0(ref1Label, ".centromere.txt")), header = F)[, 2])
centrSRefB <- as.numeric(read.table(file = file.path(baseDir, "Ref", "Ann", 
                                                   paste0(ref2Label, ".centromere.txt")), header = F)[, 1])
centrERefB <- as.numeric(read.table(file = file.path(baseDir, "Ref", "Ann", 
                                                   paste0(ref2Label, ".centromere.txt")), header = F)[, 2])

subTelLRefA <- as.numeric(read.table(file = file.path(baseDir, "Ref", "Ann", paste0(ref1Label, ".subtel.txt")))[, 1])
subTelRRefA <- as.numeric(read.table(file = file.path(baseDir, "Ref", "Ann", paste0(ref1Label, ".subtel.txt")))[, 2])
subTelLRefB <- as.numeric(read.table(file = file.path(baseDir, "Ref", "Ann", paste0(ref2Label, ".subtel.txt")))[, 1])
subTelRRefB <- as.numeric(read.table(file = file.path(baseDir, "Ref", "Ann", paste0(ref2Label, ".subtel.txt")))[, 2])

myMarkers <- list.files(path = dirMUMmer, full.names = T, recursive = T, pattern = "intersect\\.snps$")

# shift
cPosRefA <- (centrSRefA + centrERefA) / 2
cPosRefB <- (centrSRefB + centrERefB) / 2
shiftRefBRefA <- cPosRefB - cPosRefA
# assembly di riferimento per plot label
labAssemblyRiferimento <- ref2

# retrieve MUMmer marker positions --------------------------------------------------

stdMUM <- read.table(file = myMarkers, skip = 4, header = F)
dfGG <- stdMUM[, c(11, 1, 12, 4, 2, 3, 9, 10)]
rm(stdMUM)
colnames(dfGG) <- c(paste0("chr", ref1), 
                    paste0("pos", ref1), 
                    paste0("chr", ref2), 
                    paste0("pos", ref2), 
                    paste0("ref", ref1), 
                    paste0("ref", ref2), 
                    paste0("strand", ref1), 
                    paste0("strand", ref2))

dfGG[, 2] <- as.numeric(dfGG[, 2])
dfGG[, 4] <- as.numeric(dfGG[, 4])

# save dfGG for de novo variants filtering
save(dfGG, file = file.path(outDirMUM, paste0("Markers.", ref1, "-", ref2, ".RData")))

# remove MT markers: da qui "non tornano" i row.names(dfGG)
outMT <- unique(c(which(dfGG[, 1] == "chrMT"), which(dfGG[, 3] == "chrMT")))
if (length(outMT) != 0) {
  dfGG <- dfGG[-outMT, ]
}

# plot ref1 MUMmer marker positions ---------------------------------------------------

dfSNP <- dfGG[, 1:2]
colnames(dfSNP) <- c("chr", "pos")

# make factor to set facet panels order
dfSNP$chr <- factor(dfSNP$chr, levels = allChr)

# add (sub)telomers and centromer to dataframe
nVar <- nrow(dfSNP)
chrCentrSRefA <- numeric(length = nVar)
chrCentrERefA <- numeric(length = nVar)
chrSubTelLRefA <- numeric(length = nVar)
chrSubTelRRefA <- numeric(length = nVar)
lenRefA <- numeric(length = nVar)

counterChr <- 1
for (indR in allChr) {
  indBon <- which(dfSNP$chr == indR)
  chrCentrSRefA[indBon] <- centrSRefA[counterChr]
  chrCentrERefA[indBon] <- centrERefA[counterChr]
  chrSubTelLRefA[indBon] <- subTelLRefA[counterChr]
  chrSubTelRRefA[indBon] <- subTelRRefA[counterChr]
  lenRefA[indBon] <- chrLenRefA[counterChr]
  counterChr <- counterChr + 1
}
dfSNP <- data.frame(dfSNP, lenRefA, chrCentrSRefA, chrCentrERefA, chrSubTelLRefA, chrSubTelRRefA)

# switching data to kb
dfSNP[, 2:7] <- dfSNP[, 2:7] / 1000
dfSNP <- dfSNP[order(match(dfSNP$chr, allChr)), ]

# plot
ref1StrXlab <- paste0("Position [kb]\n(",  ref1, ")")
pDenSNP <- ggplot(dfSNP) + 
  ggtitle("Marker Positions", subtitle = paste("Bin = ", bWidth, " [kb]", sep = "")) + 
  # centromere
  geom_rect(data = NULL, aes(xmin = chrCentrSRefA, xmax = chrCentrERefA, ymin = 1, ymax = Inf), 
            fill = "red") + 
  # left (sub)telomer
  geom_rect(data = NULL, aes(xmin = 1, xmax = chrSubTelLRefA, ymin = 1, ymax = Inf), 
            fill = "orange") + 
  # right (sub)telomer
  geom_rect(data = NULL, aes(xmin = chrSubTelRRefA, xmax = lenRefA, ymin = 1, ymax = Inf), 
            fill = "orange") + 
  geom_histogram(aes(x = pos), binwidth = bWidth) + 
  xlab(ref1StrXlab) + 
  ylab("Count") + 
  facet_wrap(facets = ~ chr, ncol = 2, scales = "free_x") + 
  theme(plot.margin = unit(c(1, 2, 1, 1), "cm"), 
        # title style
        plot.title = element_text(size = 16, face = "bold", hjust = .5), 
        plot.subtitle = element_text(size = 12, face = "bold", hjust = .5), 
        # axis style
        axis.title = element_text(size = 14, face = "bold"), 
        axis.text = element_text(size = 10), 
        # subtitle format
        strip.background = element_rect(fill = "white", colour = "white"), 
        strip.text.x = element_text(size = 10, face = "bold"), 
        # y labels
        axis.text.y = element_text(size = 6)) + 
  # scale_y_continuous(breaks = c(50, 100)) + 
  scale_x_continuous(labels = comma)

pPath <- file.path(outDirMUM, paste0("Density.", ref1, ".pdf"))
pdf(file = pPath)
print(pDenSNP)
dev.off()

# zoom plot
pDenSNP <- ggplot(dfSNP) + 
  ggtitle("Marker Positions", subtitle = paste("Bin = ", bWidth, " [kb]", sep = "")) + 
  # left (sub)telomer
  geom_rect(data = NULL, aes(xmin = 1, xmax = chrSubTelLRefA, ymin = 1, ymax = Inf), 
            fill = "orange") + 
  # right (sub)telomer
  geom_rect(data = NULL, aes(xmin = chrSubTelRRefA, xmax = lenRefA, ymin = 1, ymax = Inf), 
            fill = "orange") + 
  geom_histogram(aes(x = pos), binwidth = bWidth) + 
  # centromere
  geom_rect(data = NULL, aes(xmin = chrCentrSRefA, xmax = chrCentrERefA, ymin = 1, ymax = Inf), 
            fill = "red") + 
  xlab(ref1StrXlab) + 
  ylab("Count") + 
  facet_wrap(facets = ~ chr, ncol = 2, scales = "free_x") + 
  theme(plot.margin = unit(c(1, 2, 1, 1), "cm"), 
        # title style
        plot.title = element_text(size = 16, face = "bold", hjust = .5), 
        plot.subtitle = element_text(size = 12, face = "bold", hjust = .5), 
        # axis style
        axis.title = element_text(size = 14, face = "bold"), 
        axis.text = element_text(size = 10), 
        # subtitle format
        strip.background = element_rect(fill = "white", colour = "white"), 
        strip.text.x = element_text(size = 10, face = "bold"), 
        # y labels
        axis.text.y = element_text(size = 6)) + 
  # scale_y_continuous(breaks = c(50, 100)) + 
  coord_cartesian(ylim = c(0, 15)) + 
  scale_x_continuous(labels = comma)

pPath <- file.path(outDirMUM, paste0("Density.zm.", ref1, ".pdf"))
pdf(file = pPath)
print(pDenSNP)
dev.off()

# plot ref2 MUMmer marker positions ---------------------------------------------------

dfSNP <- dfGG[, 3:4]
colnames(dfSNP) <- c("chr", "pos")

# make factor to set facet panels order
dfSNP$chr <- factor(dfSNP$chr, levels = allChr)

# add (sub)telomers and centromer to dataframe
nVar <- nrow(dfSNP)
chrCentrSRefB <- numeric(length = nVar)
chrCentrERefB <- numeric(length = nVar)
chrSubTelLRefB <- numeric(length = nVar)
chrSubTelRRefB <- numeric(length = nVar)
lenRefB <- numeric(length = nVar)

counterChr <- 1
for (indR in allChr) {
  indBon <- which(dfSNP$chr == indR)
  chrCentrSRefB[indBon] <- centrSRefB[counterChr]
  chrCentrERefB[indBon] <- centrERefB[counterChr]
  chrSubTelLRefB[indBon] <- subTelLRefB[counterChr]
  chrSubTelRRefB[indBon] <- subTelRRefB[counterChr]
  lenRefB[indBon] <- chrLenRefB[counterChr]
  counterChr <- counterChr + 1
}
dfSNP <- data.frame(dfSNP, lenRefB, chrCentrSRefB, chrCentrERefB, chrSubTelLRefB, chrSubTelRRefB)

# switching data to kb
dfSNP[, 2:7] <- dfSNP[, 2:7] / 1000
dfSNP <- dfSNP[order(match(dfSNP$chr, allChr)), ]

# plot
ref2StrXlab <- paste0("Position [kb]\n(",  ref2, ")")
pDenSNP <- ggplot(dfSNP) + 
  ggtitle("Marker Positions", subtitle = paste("Bin = ", bWidth, " [kb]", sep = "")) + 
  # left (sub)telomer
  geom_rect(data = NULL, aes(xmin = 1, xmax = chrSubTelLRefB, ymin = 1, ymax = Inf), 
            fill = "orange") + 
  # right (sub)telomer
  geom_rect(data = NULL, aes(xmin = chrSubTelRRefB, xmax = lenRefB, ymin = 1, ymax = Inf), 
            fill = "orange") + 
  geom_histogram(aes(x = pos), binwidth = bWidth) + 
  # centromere
  geom_rect(data = NULL, aes(xmin = chrCentrSRefB, xmax = chrCentrERefB, ymin = 1, ymax = Inf), 
            fill = "red") + 
  xlab(ref2StrXlab) + 
  ylab("Count") + 
  facet_wrap(facets = ~ chr, ncol = 2, scales = "free_x") + 
  theme(plot.margin = unit(c(1, 2, 1, 1), "cm"), 
        # title style
        plot.title = element_text(size = 16, face = "bold", hjust = .5), 
        plot.subtitle = element_text(size = 12, face = "bold", hjust = .5), 
        # axis style
        axis.title = element_text(size = 14, face = "bold"), 
        axis.text = element_text(size = 10), 
        # subtitle format
        strip.background = element_rect(fill = "white", colour = "white"), 
        strip.text.x = element_text(size = 10, face = "bold"), 
        # y labels
        axis.text.y = element_text(size = 6)) + 
  # scale_y_continuous(breaks = c(50, 100)) + 
  scale_x_continuous(labels = comma)

pPath <- file.path(outDirMUM, paste0("Density.", ref2, ".pdf"))
pdf(file = pPath)
print(pDenSNP)
dev.off()

# zoom plot
pDenSNP <- ggplot(dfSNP) + 
  ggtitle("Marker Positions", subtitle = paste("Bin = ", bWidth, " [kb]", sep = "")) + 
  # centromere
  geom_rect(data = NULL, aes(xmin = chrCentrSRefB, xmax = chrCentrERefB, ymin = 1, ymax = Inf), 
            fill = "red") + 
  # left (sub)telomer
  geom_rect(data = NULL, aes(xmin = 1, xmax = chrSubTelLRefB, ymin = 1, ymax = Inf), 
            fill = "orange") + 
  # right (sub)telomer
  geom_rect(data = NULL, aes(xmin = chrSubTelRRefB, xmax = lenRefB, ymin = 1, ymax = Inf), 
            fill = "orange") + 
  geom_histogram(aes(x = pos), binwidth = bWidth) + 
  xlab(ref2StrXlab) + 
  ylab("Count") + 
  facet_wrap(facets = ~ chr, ncol = 2, scales = "free_x") + 
  theme(plot.margin = unit(c(1, 2, 1, 1), "cm"), 
        # title style
        plot.title = element_text(size = 16, face = "bold", hjust = .5), 
        plot.subtitle = element_text(size = 12, face = "bold", hjust = .5), 
        # axis style
        axis.title = element_text(size = 14, face = "bold"), 
        axis.text = element_text(size = 10), 
        # subtitle format
        strip.background = element_rect(fill = "white", colour = "white"), 
        strip.text.x = element_text(size = 10, face = "bold"), 
        # y labels
        axis.text.y = element_text(size = 6)) + 
  # scale_y_continuous(breaks = c(50, 100)) + 
  coord_cartesian(ylim = c(0, 15)) + 
  scale_x_continuous(labels = comma)

pPath <- file.path(outDirMUM, paste0("Density.zm.", ref2, ".pdf"))
pdf(file = pPath)
print(pDenSNP)
dev.off()

cat("Init time", "\n")
proc.time() - ptmInit

# remove (sub)telomeric markers -------------------------------------------

# numero del chr del marker ref1
indChrRef1 <- match(dfGG[, 1], allChr)
booRefA1 <- dfGG[, 2] > subTelLRefA[indChrRef1]
booRefA2 <- dfGG[, 2] < subTelRRefA[indChrRef1]
# boolean delle posiizoni dei marker nel core dei chr ref1
booRefA <- booRefA1 & booRefA2
# numero del chr del marker ref2
indChrRef2 <- match(dfGG[, 3], allChr)
booRefB1 <- dfGG[, 4] > subTelLRefB[indChrRef2]
booRefB2 <- dfGG[, 4] < subTelRRefB[indChrRef2]
# boolean delle posiizoni dei marker nel core dei chr ref2
booRefB <- booRefB1 & booRefB2
# boolean of non-(sub)telomeric markers
booRef1Ref2 <- booRefA & booRefB
# markers not within (sub)telomeric positions
dfGG <- dfGG[booRef1Ref2, ]

# plot filtered data (pasted from Markers.Dist che fa lo stesso sui dati non filtrati)
dfSNP <- dfGG

# chr and pos c(1, 2) for RefA (RefAsm1): sort by chr
dfSNPRefAsm1 <- dfSNP[order(match(dfSNP[, 1], allChr)), ][, c(1,2)]
# chr and pos c(3, 4) for RefB (RefAsm2): sort by chr
dfSNPRefAsm2 <- dfSNP[order(match(dfSNP[, 3], allChr)), ][, c(3,4)]

# header tables marker density
fileDenRefA <- file.path(outDirMUM, paste(outLabelA, ref1, "txt", sep = "."))
fileDenRefB <- file.path(outDirMUM, paste(outLabelA, ref2, "txt", sep = "."))

headerDensity <- data.frame("chr", "length", "# markers", "density")
write.table(x = headerDensity, file = fileDenRefA, col.names = F, row.names = F, quote = F, append = F, sep = "\t")
write.table(x = headerDensity, file = fileDenRefB, col.names = F, row.names = F, quote = F, append = F, sep = "\t")

accChrRefA <- c()
accChrRefB <- c()
accDiffRefA <- c()
accDiffRefB <- c()
counterChr <- 1
for (indC in allChr) {
  # RefAsm1 pos sorted
  posGG <- sort(dfSNPRefAsm1[which(dfSNPRefAsm1[, 1] == indC), 2])
  lenPos <- length(posGG)
  diffPosRefA <- diff(posGG)
  accChrRefA <- c(accChrRefA, rep(indC, lenPos - 1))
  accDiffRefA <- c(accDiffRefA, diffPosRefA)
  denRefA <- lenPos / chrLenRefA[counterChr]
  write.table(x = data.frame(indC, chrLenRefA[counterChr], lenPos, denRefA),  
              file = fileDenRefA, 
              col.names = F, row.names = F, quote = F, append = T, sep = "\t")
  # RefAsm1 all markers
  bWidth <- 100
  nEv <- length(diffPosRefA)
  diffPosRefA <- data.frame(diffPosRefA)
  colnames(diffPosRefA) <- c("dist")
  plot1RefAsm1 <- ggplot(diffPosRefA, aes(x = dist)) + 
    ggtitle(paste0("Marker distance, ", ref1 , " strain"), 
            subtitle = paste(indC, "; ", 
                             "N = ", nEv, "; ", 
                             "Bin = ", bWidth, " bp", 
                             sep = "")) + 
    geom_histogram(aes(y = ..density..), 
                   binwidth = bWidth, 
                   colour = "blue", fill = "white") + 
    geom_density(alpha = .2, fill = "blue") + 
    xlab("Distance [bp]") + 
    ylab("Density") + 
    theme(
      plot.margin = unit(c(1, 2, 1, 1), "cm"), 
      # title style
      plot.title = element_text(size = 16, hjust = 1), 
      plot.subtitle = element_text(size = 12, face = "bold", hjust = 1), 
      # axis style
      axis.title = element_text(size = 16, face = "bold"), 
      axis.text = element_text(size = 16))
  plo1PathRefA <- file.path(outDirMUM, paste(outLabelB, ref1, indC, "pdf", sep = "."))
  pdf(file = plo1PathRefA)
  print(plot1RefAsm1, width = 10, height = 26)
  dev.off()
  
  # < 500
  bWidth <- 10
  diffPosRefA500 <- diffPosRefA[which(diffPosRefA$dist < 500), ]
  nEv <- length(diffPosRefA500)
  diffPosRefA500 <- data.frame(diffPosRefA500)
  colnames(diffPosRefA500) <- c("dist")
  plot1RefAsm1 <- ggplot(diffPosRefA500, aes(x = dist)) + 
    ggtitle(paste0("Marker distance ( < 500 bp), ", ref1 , " strain"), 
            subtitle = paste(indC, "; ", 
                             "N = ", nEv, "; ", 
                             "Bin = ", bWidth, " bp", 
                             sep = "")) + 
    geom_histogram(aes(y = ..density..), 
                   binwidth = bWidth, 
                   colour = "blue", fill = "white") + 
    geom_density(alpha = .2, fill = "blue") + 
    xlab("Distance [bp]") + 
    ylab("Density") + 
    theme(
      plot.margin = unit(c(1, 2, 1, 1), "cm"), 
      # title style
      plot.title = element_text(size = 16, hjust = 1), 
      plot.subtitle = element_text(size = 12, face = "bold", hjust = 1), 
      # axis style
      axis.title = element_text(size = 16, face = "bold"), 
      axis.text = element_text(size = 16))
  plo1PathRefA <- file.path(outDirMUM, paste(outLabelB, ref1, indC, "500.pdf", sep = "."))
  pdf(file = plo1PathRefA)
  print(plot1RefAsm1, width = 10, height = 26)
  dev.off()
  
  # RefAsm2 pos sorted
  posGG <- sort(dfSNPRefAsm2[which(dfSNPRefAsm2[, 1] == indC), 2])
  lenPos <- length(posGG)
  diffPosRefB <- diff(posGG)
  accChrRefB <- c(accChrRefB, rep(indC, lenPos - 1))
  accDiffRefB <- c(accDiffRefB, diffPosRefB)
  denRefB <- lenPos / chrLenRefB[counterChr]
  write.table(x = data.frame(indC, chrLenRefB[counterChr], lenPos, denRefB), 
              file = fileDenRefB, 
              col.names = F, row.names = F, quote = F, append = T, sep = "\t")
  # RefAsm2 all markers
  bWidth <- 100
  nEv <- length(diffPosRefB)
  diffPosRefB <- data.frame(diffPosRefB)
  colnames(diffPosRefB) <- c("dist")
  plot1RefAsm2 <- ggplot(diffPosRefB, aes(x = dist)) + 
    ggtitle(paste0("Marker distance, ", ref2 , " strain"), 
            subtitle = paste(indC, "; ", 
                             "N = ", nEv, "; ", 
                             "Bin = ", bWidth, " bp", 
                             sep = "")) + 
    geom_histogram(aes(y = ..density..), 
                   binwidth = bWidth, 
                   colour = "blue", fill = "white") + 
    geom_density(alpha = .2, fill = "blue") + 
    xlab("Distance [bp]") + 
    ylab("Density") + 
    theme(
      plot.margin = unit(c(1, 2, 1, 1), "cm"), 
      # title style
      plot.title = element_text(size = 16, hjust = 1), 
      plot.subtitle = element_text(size = 12, face = "bold", hjust = 1), 
      # axis style
      axis.title = element_text(size = 16, face = "bold"), 
      axis.text = element_text(size = 16))
  plo1PathRefB <- file.path(outDirMUM, paste(outLabelB, ref2, indC, "pdf", sep = "."))
  pdf(file = plo1PathRefB)
  print(plot1RefAsm2, width = 10, height = 26)
  dev.off()
  
  # < 500
  bWidth <- 10
  diffPosRefB500 <- diffPosRefB[which(diffPosRefB$dist < 500), ]
  nEv <- length(diffPosRefB500)
  diffPosRefB500 <- data.frame(diffPosRefB500)
  colnames(diffPosRefB500) <- c("dist")
  plot1RefAsm2 <- ggplot(diffPosRefB500, aes(x = dist)) + 
    ggtitle(paste0("Marker distance ( < 500 bp), ", ref2 , " strain"), 
            subtitle = paste(indC, "; ", 
                             "N = ", nEv, "; ", 
                             "Bin = ", bWidth, " bp", 
                             sep = "")) + 
    geom_histogram(aes(y = ..density..), 
                   binwidth = bWidth, 
                   colour = "blue", fill = "white") + 
    geom_density(alpha = .2, fill = "blue") + 
    xlab("Distance [bp]") + 
    ylab("Density") + 
    theme(
      plot.margin = unit(c(1, 2, 1, 1), "cm"), 
      # title style
      plot.title = element_text(size = 16, hjust = 1), 
      plot.subtitle = element_text(size = 12, face = "bold", hjust = 1), 
      # axis style
      axis.title = element_text(size = 16, face = "bold"), 
      axis.text = element_text(size = 16))
  plo1PathRefB <- file.path(outDirMUM, paste(outLabelB, ref2, indC, "500.pdf", sep = "."))
  pdf(file = plo1PathRefB)
  print(plot1RefAsm2, width = 10, height = 26)
  dev.off()
  
  counterChr <- counterChr + 1
}

# RefAsm1 dataframe for ggplot 
dfGGDist4RefA <- data.frame(accChrRefA, accDiffRefA)
colnames(dfGGDist4RefA) <- c("chr", "dist")
# RefAsm2 dataframe for ggplot 
dfGGDist4RefB <- data.frame(accChrRefB, accDiffRefB)
colnames(dfGGDist4RefB) <- c("chr", "dist")

# make factor to set facet panels order
dfGGDist4RefA$chr <- factor(dfGGDist4RefA$chr, levels = allChr)
dfGGDist4RefB$chr <- factor(dfGGDist4RefB$chr, levels = allChr)
dfGGDist4RefA <- dfGGDist4RefA[order(match(dfGGDist4RefA$chr, allChr)), ]
dfGGDist4RefB <- dfGGDist4RefB[order(match(dfGGDist4RefB$chr, allChr)), ]

# RefAsm1 markers distant < 0.5 kb
dfGGDist4RefA500 <- dfGGDist4RefA[which(dfGGDist4RefA$dist < 500), ]
bWidth <- 10
nEv <- nrow(dfGGDist4RefA500)
plot3RefAsm1 <- ggplot(dfGGDist4RefA500, aes(x = dist)) + 
  ggtitle(paste0("Marker distance (< 500 bp), ", ref1 , " strain"), 
          subtitle = paste("N = ", nEv, "; ", 
                           "Bin = ", bWidth, " bp", 
                           sep = "")) + 
  geom_histogram(aes(y = ..density..), 
                 binwidth = bWidth, 
                 colour = "blue", fill = "white") + 
  geom_density(alpha = .2, fill = "blue") + 
  facet_wrap(facets = ~ chr, ncol = 2, scales = "free_x") + 
  xlab("Distance [bp]") + 
  ylab("Density") + 
  theme(
    plot.margin = unit(c(1, 2, 1, 1), "cm"), 
    # title style
    plot.title = element_text(size = 16, hjust = 1), 
    plot.subtitle = element_text(size = 12, face = "bold", hjust = 1), 
    # axis style
    axis.title = element_text(size = 16, face = "bold"), 
    axis.text = element_text(size = 16))
plo3PathRefA <- file.path(outDirMUM, paste(outLabelB, ref1, "500.pdf", sep = "."))
pdf(file = plo3PathRefA, width = 10, height = 26)
print(plot3RefAsm1)
dev.off()

# RefAsm2 markers distant < 0.5 kb
dfGGDist4RefB500 <- dfGGDist4RefB[which(dfGGDist4RefB$dist < 500), ]
bWidth <- 10
nEv <- nrow(dfGGDist4RefB500)
plot3RefAsm2 <- ggplot(dfGGDist4RefB500, aes(x = dist)) + 
  ggtitle(paste0("Marker distance (< 500 bp), ", ref2 , " strain"), 
          subtitle = paste("N = ", nEv, "; ", 
                           "Bin = ", bWidth, " bp", 
                           sep = "")) + 
  geom_histogram(aes(y = ..density..), 
                 binwidth = bWidth, 
                 colour = "red", fill = "white") + 
  geom_density(alpha = .2, fill = "red") + 
  facet_wrap(facets = ~ chr, ncol = 2, scales = "free_x") + 
  xlab("Distance [bp]") + 
  ylab("Density") + 
  theme(
    plot.margin = unit(c(1, 2, 1, 1), "cm"), 
    # title style
    plot.title = element_text(size = 16, hjust = 1), 
    plot.subtitle = element_text(size = 12, face = "bold", hjust = 1), 
    # axis style
    axis.title = element_text(size = 16, face = "bold"), 
    axis.text = element_text(size = 16))
plo3PathRefB <- file.path(outDirMUM, paste(outLabelB, ref2, "500.pdf", sep = "."))
pdf(file = plo3PathRefB, width = 10, height = 26)
print(plot3RefAsm2)
dev.off()

# sample loop ------------------------------------------------------

# per tutti gli ibridi
lFiles <- basename(list.files(path = polDir, recursive = T, pattern = "vcf\\.RData$"))
sampIDs <- unique(sapply(strsplit(lFiles, split = ".", fixed = T), "[[", 1))

for (ind1 in sampIDs) {
  # make path to output folder
  sampOutDir <- file.path(resDir, ind1)
  dir.create(sampOutDir, showWarnings = F, recursive = T)
  
  # quality filters ---------------------------------------------------------
  
  # load RefA vcf
  pathSampRefA <- file.path(polDir, ind1, grep(ref1Label, 
                                             grep(pattern = paste(ind1, ref1Label, sep = "."), 
                                                  lFiles, value = T), value = T))
  load(pathSampRefA) # load vcfData list
  vcfFiltRefA <- list(vcfData$meta, as.data.frame(vcfData$fix), as.data.frame(vcfData$format))
  names(vcfFiltRefA) <- c("meta", "fix", "format")
  
  # QUAL calculation
  vcfFiltRefA$fix[, 6] <- as.numeric(vcfFiltRefA$fix[, 6])
  mRefA <- round(mean(vcfFiltRefA$fix[, 6]))
  sRefA <- round(sd(vcfFiltRefA$fix[, 6]))
  cat(paste(ind1, " ", ref1, " mean(QUAL): ", mRefA, "\n", 
            ind1, " ", ref1, " sd(QUAL): ", sRefA, "\n", sep = ""), 
      file = file.path(sampOutDir, paste(ind1, ref1, "QUALdata.txt", sep = ".")))
  indQualRefA <- which(vcfFiltRefA$fix[, 6] > c(mRefA - sRefA))
  # save low quality vcf
  nCoreVarRefA <- length(vcfFiltRefA$fix$POS)
  indLowQualRefA <- setdiff(c(1:nCoreVarRefA), indQualRefA)
  lowQualRefAvcf <- list(vcfFiltRefA$meta, vcfFiltRefA$fix[indLowQualRefA, ], vcfFiltRefA$format[indLowQualRefA, ])
  names(lowQualRefAvcf) <- c("meta", "fix", "format")
  pathOutLowQualRefA <- file.path(sampOutDir, paste(ind1, "LQ", ref1, "RData", sep = "."))
  cat(c("Fraction of", ref1, "low quality variants (/# core positions)", paste(ind1, ":", sep = ""), 
        paste(length(indLowQualRefA), nCoreVarRefA, sep = "/"), "\n"))
  save(list = c("lowQualRefAvcf"), file = pathOutLowQualRefA)
  # QUAL filtering
  vcfFiltRefA$fix <- vcfFiltRefA$fix[indQualRefA, ]
  vcfFiltRefA$format <- vcfFiltRefA$format[indQualRefA, ]
  
  # load RefB vcf
  pathSampRefB <- file.path(polDir, ind1, grep(ref2Label, 
                                             grep(pattern = paste(ind1, ref2Label, sep = "."), 
                                                  lFiles, value = T), value = T))
  load(pathSampRefB) # load vcfData list
  vcfFiltRefB <- list(vcfData$meta, as.data.frame(vcfData$fix), as.data.frame(vcfData$format))
  names(vcfFiltRefB) <- c("meta", "fix", "format")
  
  # QUAL calculation
  vcfFiltRefB$fix[, 6] <- as.numeric(vcfFiltRefB$fix[, 6])
  mRefB <- round(mean(vcfFiltRefB$fix[, 6]))
  sRefB <- round(sd(vcfFiltRefB$fix[, 6]))
  cat(paste(ind1, " ", ref2, " mean(QUAL): ", mRefB, "\n", 
            ind1, " ", ref2, " sd(QUAL): ", sRefB, "\n", sep = ""), 
      file = file.path(sampOutDir, paste(ind1, ref2, "QUALdata.txt", sep = ".")))
  indQualRefB <- which(vcfFiltRefB$fix[, 6] > c(mRefB - sRefB))
  # save low quality vcf
  nCoreVarRefB <- length(vcfFiltRefB$fix$POS)
  indLowQualRefB <- setdiff(c(1:nCoreVarRefB), indQualRefB)
  lowQualRefBvcf <- list(vcfFiltRefB$meta, vcfFiltRefB$fix[indLowQualRefB, ], vcfFiltRefB$format[indLowQualRefB, ])
  names(lowQualRefBvcf) <- c("meta", "fix", "format")
  pathOutLowQualRefB <- file.path(sampOutDir, paste(ind1, "LQ", ref2, "RData", sep = "."))
  cat(c("Fraction of", ref2, "low quality variants (/# core positions)", paste(ind1, ":", sep = ""), 
        paste(length(indLowQualRefB), nCoreVarRefB, sep = "/"), "\n"))
  save(list = c("lowQualRefBvcf"), file = pathOutLowQualRefB)
  # QUAL filtering
  vcfFiltRefB$fix <- vcfFiltRefB$fix[indQualRefB, ]
  vcfFiltRefB$format <- vcfFiltRefB$format[indQualRefB, ]
  
  # hic fuit chr loop
  # chrGG <- dfGG[, c(2, 4)]
  # bad guy
  chrGG <- data.frame(paste(dfGG[, 1], dfGG[, 2], sep = ""), paste(dfGG[, 3], dfGG[, 4], sep = ""), 
                      dfGG[, 5], dfGG[, 6], dfGG[, 7], dfGG[, 8])
  names(chrGG) <- c("posRef1", "posRef2", "alleleRef1", "alleleRef2", "strandRef1", "strandRef2")
  
  dfVcfRefAfix <- vcfFiltRefA$fix
  dfVcfRefBfix <- vcfFiltRefB$fix
  chrVcfRefAformat <- vcfFiltRefA$format[, 2]
  chrVcfRefBformat <- vcfFiltRefB$format[, 2]
  
  # determine common chr variants
  # booleans 
  # T se la variante è nel vcf
  # F se non c'è
  pstChrPosRefA <- paste(dfVcfRefAfix$CHR, dfVcfRefAfix$POS, sep = "")
  pstChrPosRefB <- paste(dfVcfRefBfix$CHR, dfVcfRefBfix$POS, sep = "")
  ieRefA <- is.element(chrGG$posRef1, pstChrPosRefA)
  ieRefB <- is.element(chrGG$posRef2, pstChrPosRefB)
  goldBool <- ieRefA & ieRefB
  # va come una scheggia
  # a <- c("29068", "29125", "29129", "29153")
  # b <- c("29086", "29125", "29129")
  # pos <- c("29064", "29068", "29086", "29125", "29129", "29153")
  # ie1 <- is.element(pos, a)
  # ie2 <- is.element(pos, b)
  # ie1 & ie2
  # mappings present in both vcfs (gold)
  # Good Guys (filtered, i.e. core mapping) Gold
  chrGGG <- chrGG[goldBool, ]
  rm(chrGG)
  # pesca le posizioni GGG nei vcf di RefA e RefB
  # indici del secondo vettore da prendere: dfGGG vs vcf
  # ricordati che nei non-colineari, o in caso di inversioni
  # questi vettori sono diversi (ad esempio se il marker chrI1
  # del ref1 corrisponde al marker chrX1234, indPosRefA[1] sarà 1 ma
  # indPosRefB sarà N >> 1)
  indPosRefA <- match(chrGGG$posRef1, pstChrPosRefA)
  indPosRefB <- match(chrGGG$posRef2, pstChrPosRefB)
  # in questi df le varianti sono ordinate secondo la mappa dei marker chrGGG
  dfRefAfixGGG <- dfVcfRefAfix[indPosRefA, ]
  dfRefBfixGGG <- dfVcfRefBfix[indPosRefB, ]
  rm(dfVcfRefAfix, dfVcfRefBfix)
  
  chrRefAformatGGG <- chrVcfRefAformat[indPosRefA]
  chrRefBformatGGG <- chrVcfRefBformat[indPosRefB]
  rm(chrVcfRefAformat, chrVcfRefBformat)
  # spacchetta il campo sample
  # convert genotype to state ID integer
  # TODO: pre-allocate data structure
  dfRefAformatGGG <- data.frame(sapply(strsplit(chrRefAformatGGG, split = ":"), "[[", 1), 
                              sapply(strsplit(chrRefAformatGGG, split = ":"), "[[", 2), 
                              sapply(strsplit(chrRefAformatGGG, split = ":"), "[[", 3), 
                              sapply(strsplit(chrRefAformatGGG, split = ":"), "[[", 4), 
                              sapply(strsplit(chrRefAformatGGG, split = ":"), "[[", 5), 
                              sapply(strsplit(chrRefAformatGGG, split = ":"), "[[", 6), 
                              sapply(strsplit(chrRefAformatGGG, split = ":"), "[[", 7))
  rm(chrRefAformatGGG)
  dfRefBformatGGG <- data.frame(sapply(strsplit(chrRefBformatGGG, split = ":"), "[[", 1), 
                              sapply(strsplit(chrRefBformatGGG, split = ":"), "[[", 2), 
                              sapply(strsplit(chrRefBformatGGG, split = ":"), "[[", 3), 
                              sapply(strsplit(chrRefBformatGGG, split = ":"), "[[", 4), 
                              sapply(strsplit(chrRefBformatGGG, split = ":"), "[[", 5), 
                              sapply(strsplit(chrRefBformatGGG, split = ":"), "[[", 6), 
                              sapply(strsplit(chrRefBformatGGG, split = ":"), "[[", 7))
  rm(chrRefBformatGGG)
  colnames(dfRefAformatGGG) <- c("GT", "PL", "DP", "SP", "ADF", "ADR", "AD")
  colnames(dfRefBformatGGG) <- c("GT", "PL", "DP", "SP", "ADF", "ADR", "AD")
  
  # GT filters ---------------------------------------------------------
  
  # remove non-matching GTs
  nM1 <- which(dfRefAformatGGG$GT == "1/1" & dfRefBformatGGG$GT != "0/0")
  nM2 <- which(dfRefAformatGGG$GT == "0/0" & dfRefBformatGGG$GT != "1/1")
  nM3 <- which(dfRefAformatGGG$GT == "0/1" & dfRefBformatGGG$GT != "0/1")
  nM4 <- which(dfRefAformatGGG$GT != "0/1" & dfRefBformatGGG$GT == "0/1")
  nM5 <- which(dfRefBformatGGG$GT == "0/0" & dfRefAformatGGG$GT != "1/1")
  nM6 <- which(dfRefBformatGGG$GT == "1/1" & dfRefAformatGGG$GT != "0/0")
  
  nMatchGT <- unique(c(nM1, nM2, nM3, nM4, nM5, nM6))
  
  if (length(nMatchGT) > 0) {
    cat(c("Number of non-matching GT", paste(ind1, ":", sep = ""), length(nMatchGT), "\n"))
    
    unMtRefAfix <- dfRefAfixGGG[nMatchGT, ]
    dfRefAfixGGG <- dfRefAfixGGG[-nMatchGT, ]
    unMtRefBfix <- dfRefBfixGGG[nMatchGT, ]
    dfRefBfixGGG <- dfRefBfixGGG[-nMatchGT, ]
    
    unMtRefAformat <- dfRefAformatGGG[nMatchGT, ]
    dfRefAformatGGG <- dfRefAformatGGG[-nMatchGT, ]
    unMtRefBformat <- dfRefBformatGGG[nMatchGT, ]
    dfRefBformatGGG <- dfRefBformatGGG[-nMatchGT, ]
    
    unMtRefAvcf <- list(vcfFiltRefA$meta, unMtRefAfix, unMtRefAformat)
    pathOutUMRefA <- file.path(sampOutDir, paste(ind1, "UnMatchedGT", ref1, "RData", sep = "."))
    save(list = c("unMtRefAvcf"), file = pathOutUMRefA)
    
    unMtRefBvcf <- list(vcfFiltRefB$meta, unMtRefBfix, unMtRefBformat)
    pathOutUMRefB <- file.path(sampOutDir, paste(ind1, "UnMatchedGT", ref2, "RData", sep = "."))
    save(list = c("unMtRefBvcf"), file = pathOutUMRefB)
    
    unMtGGG <- chrGGG[nMatchGT, ]
    chrGGG <- chrGGG[-nMatchGT, ]
    pathOutUM <- file.path(sampOutDir, paste(ind1, "UnMatchedGT.RData", sep = "."))
    save(unMtGGG, file = pathOutUM)
  }
  
  # filtering non-matching alleles: la versione V4 non funziona per strain non colineari
  # perché nei ref riarrangiati le varianti nei vcf non sono
  # ordinate come nella mappa di nucmer
  # ovvero alla riga N della mappa compaiono due varianti che non necessariamente
  # si trovano alla riga N dei due vcf
  
  # inoltre il check degli alleli non funziona
  # nelle regioni invertite
  # perché BWA mappa tutte le reverse complement
  # quindi l'allele ALT del vcf1 è il complementare
  # dell'allele REF del vcf2
  
  # trick
  # V1 <- c("A","T","T","A","G")
  # V2 <- c("A","C","T","A","G")
  # lP <- length(V1)
  # pV1 <- paste(V1, rep(1:lP), sep = "")
  # pV2 <- paste(V2, rep(1:lP), sep = "")
  # !is.element(pV1, pV2)
  
  # non 0/0 positions in RefA
  indRefAvAltAll <- which(!is.na(dfRefAfixGGG$ALT))
  lPRefA <- length(indRefAvAltAll)
  # CHR, POS and ALT in RefA
  cPosAltRefA <- dfRefAfixGGG[indRefAvAltAll, c(1, 2, 5)]
  # CHR, POS and REF in RefB
  cPosRefRefB <- dfRefBfixGGG[indRefAvAltAll, c(1, 2, 4)]
  # it's a character trap!
  alleleBoolRefA <- !is.element(paste0(cPosAltRefA$ALT, rep(1:lPRefA)), 
                              paste0(cPosRefRefB$REF, rep(1:lPRefA)))
  posNoMatchRefA <- paste(cPosAltRefA$CHROM[alleleBoolRefA], cPosAltRefA$POS[alleleBoolRefA], sep = "")
  # mSuS1: indici delle varianti non 0/0 in RefA il cui allele ALT non
  # matcha con quello REF nella relativa posizione di RefB
  mSuS1 <- match(posNoMatchRefA, chrGGG$posRef1)
  # check for inversions in ref1
  complementAltRef1 <- chartr("ATCG", "TAGC", dfRefAfixGGG$ALT[mSuS1])
  refAlleleRef2 <- dfRefBfixGGG$REF[mSuS1]
  # T elements: la base complementare a quella ALT (del parent 1) matcha la base REF (del parent 2)
  booComp <- is.element(paste0(complementAltRef1, rep(1:length(complementAltRef1))), 
                        paste0(refAlleleRef2, rep(1:length(refAlleleRef2))))
  # T elements: i marker sono in una inversione (su uno qualsiasi dei reference)
  booInv <- (chrGGG[mSuS1, 5] == -1 | chrGGG[mSuS1, 6] == -1)
  # tieni i marker nelle inversioni con l'allele giusto 
  mSuS1 <- mSuS1[!(booComp & booInv)]
  
  # non 0/0 positions in RefB
  indRefBvAltAll <- which(!is.na(dfRefBfixGGG$ALT))
  lPRefB <- length(indRefBvAltAll)
  # CHR, POS and ALT in RefB
  cPosAltRefB <- dfRefBfixGGG[indRefBvAltAll, c(1, 2, 5)]
  # CHR, POS and REF in RefA
  cPosRefRefA <- dfRefAfixGGG[indRefBvAltAll, c(1, 2, 4)]
  # it's a character trap!
  alleleBoolRefB <- !is.element(paste0(cPosAltRefB$ALT, rep(1:lPRefB)), 
                              paste0(cPosRefRefA$REF, rep(1:lPRefB)))
  posNoMatchRefB <- paste(cPosAltRefB$CHROM[alleleBoolRefB], cPosAltRefB$POS[alleleBoolRefB], sep = "")
  # mSuS2: indici delle varianti non 0/0 in RefB il cui allele ALT non
  # matcha con quello REF della relativa posizione di RefA
  mSuS2 <- match(posNoMatchRefB, chrGGG$posRef2)
  # check for inversions in ref2
  complementAltRef2 <- chartr("ATCG", "TAGC", dfRefBfixGGG$ALT[mSuS2])
  refAlleleRef1 <- dfRefAfixGGG$REF[mSuS2]
  # T elements: la base complementare a quella ALT (del parent 2) matcha la base REF (del parent 1)
  booComp <- is.element(paste0(complementAltRef2, rep(1:length(complementAltRef2))), 
                        paste0(refAlleleRef1, rep(1:length(refAlleleRef1))))
  # T elements: i marker sono in una inversione (su uno qualsiasi dei reference)
  booInv <- (chrGGG[mSuS2, 5] == -1 | chrGGG[mSuS2, 6] == -1)
  # tieni i marker nelle inversioni con l'allele giusto 
  mSuS2 <- mSuS2[!(booComp & booInv)]
  
  # merge non-matching allele position 
  # and remove lines from chrGGG map, vcf fix and vcf format data
  nMatchAL <- unique(c(mSuS1, mSuS2))
  
  # alleles filters ---------------------------------------------------------
  
  if (length(nMatchAL) > 0) {
    cat(c("Number of non-matching alleles", paste(ind1, ":", sep = ""), length(nMatchAL), "\n"))
    
    unMtRefAfix <- dfRefAfixGGG[nMatchAL, ]
    dfRefAfixGGG <- dfRefAfixGGG[-nMatchAL, ]
    unMtRefBfix <- dfRefBfixGGG[nMatchAL, ]
    dfRefBfixGGG <- dfRefBfixGGG[-nMatchAL, ]
    
    unMtRefAformat <- dfRefAformatGGG[nMatchAL, ]
    dfRefAformatGGG <- dfRefAformatGGG[-nMatchAL, ]
    unMtRefBformat <- dfRefBformatGGG[nMatchAL, ]
    dfRefBformatGGG <- dfRefBformatGGG[-nMatchAL, ]
    
    unMtRefAvcf <- list(vcfFiltRefA$meta, unMtRefAfix, unMtRefAformat)
    pathOutUMRefA <- file.path(sampOutDir, paste(ind1, "UnMatchedAL", ref1, "RData", sep = "."))
    save(list = c("unMtRefAvcf"), file = pathOutUMRefA)
    
    unMtRefBvcf <- list(vcfFiltRefB$meta, unMtRefBfix, unMtRefBformat)
    pathOutUMRefB <- file.path(sampOutDir, paste(ind1, "UnMatchedAL", ref2, "RData", sep = "."))
    save(list = c("unMtRefBvcf"), file = pathOutUMRefB)
    
    unMtGGG <- chrGGG[nMatchAL, ]
    chrGGG <- chrGGG[-nMatchAL, ]
    pathOutUM <- file.path(sampOutDir, paste(ind1, "UnMatchedAL.RData", sep = "."))
    save(unMtGGG, file = pathOutUM)
  }
  
  # CNV filters ---------------------------------------------------------
  
  # load data
  inFileRefA <- list.files(path = file.path(cnvDir, "Results", ref1Label, ind1), 
                         pattern = "_CNVs.p.value.txt$", full.names = T)
  dfCopyNumberRefA <- read.table(file = inFileRefA, header = T, sep = "\t")
  inFileRefB <- list.files(path = file.path(cnvDir, "Results", ref2Label, ind1), 
                         pattern = "_CNVs.p.value.txt$", full.names = T)
  dfCopyNumberRefB<- read.table(file = inFileRefB, header = T, sep = "\t")
  # format chr
  dfCopyNumberRefA$chr <- gsub(pattern = "_.*", replacement = "", x = dfCopyNumberRefA$chr)
  dfCopyNumberRefB$chr <- gsub(pattern = "_.*", replacement = "", x = dfCopyNumberRefB$chr)
  # keep only robust loss events
  dfCopyNumberRefA <- dfCopyNumberRefA[which(dfCopyNumberRefA$status == "loss" & 
                                           dfCopyNumberRefA$WilcoxonRankSumTestPvalue < 0.01 & 
                                           dfCopyNumberRefA$KolmogorovSmirnovPvalue < 0.01), ]
  dfCopyNumberRefB <- dfCopyNumberRefB[which(dfCopyNumberRefB$status == "loss" & 
                                           dfCopyNumberRefB$WilcoxonRankSumTestPvalue < 0.01 & 
                                           dfCopyNumberRefB$KolmogorovSmirnovPvalue < 0.01), ]
  # if loss events were detected against both references
  if (nrow(dfCopyNumberRefA) != 0 & nrow(dfCopyNumberRefB) != 0) {
    # change chromosome encoding in CNV calls
    dfCopyNumberRefA$chr <- paste("chr", dfCopyNumberRefA$chr, sep = "")
    dfCopyNumberRefB$chr <- paste("chr", dfCopyNumberRefB$chr, sep = "")
    # prepare markers table
    dfGGG4CNV <- data.frame(gsub(pattern = "[[:digit:]].*$", "", chrGGG$posRef1), 
                            as.numeric(gsub(pattern = "^.*[[:alpha:]]", "", chrGGG$posRef1)), 
                            gsub(pattern = "[[:digit:]].*$", "", chrGGG$posRef2), 
                            as.numeric(gsub(pattern = "^.*[[:alpha:]]", "", chrGGG$posRef2)))
    colnames(dfGGG4CNV) <- c("chrRef1", "posRef1", "chrRef2", "posRef2")
    # check which markers fall within loss events
    stateCopyRefA <- character(length = nrow(dfGGG4CNV))
    for (indCNV in 1:nrow(dfCopyNumberRefA)) {
      indLossRefA <- which(dfGGG4CNV$chrRef1 == dfCopyNumberRefA$chr[indCNV] & 
                           dfGGG4CNV$posRef1 >= dfCopyNumberRefA$start[indCNV] & 
                           dfGGG4CNV$posRef1 <= dfCopyNumberRefA$end[indCNV])
      stateCopyRefA[indLossRefA] <- "loss"
    }
    stateCopyRefB <- character(length = nrow(dfGGG4CNV))
    for (indCNV in 1:nrow(dfCopyNumberRefB)) {
      indLossRefB <- which(dfGGG4CNV$chrRef2 == dfCopyNumberRefB$chr[indCNV] & 
                           dfGGG4CNV$posRef2 >= dfCopyNumberRefB$start[indCNV] & 
                           dfGGG4CNV$posRef2 <= dfCopyNumberRefB$end[indCNV])
      stateCopyRefB[indLossRefB] <- "loss"
    }
    # index of markers in loss events in one or both references
    indLoss <- which(stateCopyRefA == "loss" | stateCopyRefB == "loss")
    
    if (length(indLoss) > 0) {
      
      pChrPosVcfRefA <- paste(dfRefAfixGGG$CHROM, dfRefAfixGGG$POS, sep = "")
      booRefA <- is.element(pChrPosVcfRefA, chrGGG$posRef1[indLoss])
      inLossRefAfixGGG <- dfRefAfixGGG[booRefA, ]
      dfRefAfixGGG <- dfRefAfixGGG[!booRefA, ]
      inLossRefAformatGGG <- dfRefAformatGGG[booRefA, ]
      dfRefAformatGGG <- dfRefAformatGGG[!booRefA, ]
      
      inLossRefAvcf <- list(vcfFiltRefA$meta, inLossRefAfixGGG, inLossRefAformatGGG)
      pathOutLossRefA <- file.path(sampOutDir, paste(ind1, "Loss", ref1, "RData", sep = "."))
      names(inLossRefAvcf) <- c("meta", "fix", "format")
      save(list = c("inLossRefAvcf"), file = pathOutLossRefA)
      
      pChrPosVcfRefB <- paste(dfRefBfixGGG$CHROM, dfRefBfixGGG$POS, sep = "")
      booRefB <- is.element(pChrPosVcfRefB, chrGGG$posRef2[indLoss])
      inLossRefBfixGGG <- dfRefBfixGGG[booRefB, ]
      dfRefBfixGGG <- dfRefBfixGGG[!booRefB, ]
      inLossRefBformatGGG <- dfRefBformatGGG[booRefB, ]
      dfRefBformatGGG <- dfRefBformatGGG[!booRefB, ]
      
      inLossRefBvcf <- list(vcfFiltRefB$meta, inLossRefBfixGGG, inLossRefBformatGGG)
      pathOutLossRefB <- file.path(sampOutDir, paste(ind1, "Loss", ref2, "RData", sep = "."))
      names(inLossRefBvcf) <- c("meta", "fix", "format")
      save(list = c("inLossRefBvcf"), file = pathOutLossRefB)
      
      chrGGG <- chrGGG[-indLoss, ]
    }
  }
  
  # saving data ---------------------------------------------------------
  
  # save chrGGG
  
  # save filtered vcf used for segment detection
  lstRefAvcfGGG <- list(vcfFiltRefA$meta, dfRefAfixGGG, dfRefAformatGGG)
  names(lstRefAvcfGGG) <- c("meta", "fix", "format")
  pathOutRefAvcfGGG <- file.path(sampOutDir, paste(ind1, ref1 ,"VcfFilt.RData", sep = "."))
  save(list = c("lstRefAvcfGGG"), file = pathOutRefAvcfGGG)
  
  lstRefBvcfGGG <- list(vcfFiltRefB$meta, dfRefBfixGGG, dfRefBformatGGG)
  names(lstRefBvcfGGG) <- c("meta", "fix", "format")
  pathOutRefBvcfGGG <- file.path(sampOutDir, paste(ind1, ref2, "VcfFilt.RData", sep = "."))
  save(list = c("lstRefBvcfGGG"), file = pathOutRefBvcfGGG)
  
  # segment detection ---------------------------------------------------------
  
  # sort fix & format dataframes by chromosome & position
  sortRefAfix <- dfRefAfixGGG[order(as.numeric(dfRefAfixGGG$POS)), ]
  sortRefAformat <- dfRefAformatGGG[order(as.numeric(dfRefAfixGGG$POS)), ]
  dfRefAfixGGG <- sortRefAfix[order(match(sortRefAfix$CHROM, allChr)), ]
  dfRefAformatGGG <- sortRefAformat[order(match(sortRefAfix$CHROM, allChr)), ]
  
  sortRefBfix <- dfRefBfixGGG[order(as.numeric(dfRefBfixGGG$POS)), ]
  sortRefBformat <- dfRefBformatGGG[order(as.numeric(dfRefBfixGGG$POS)), ]
  dfRefBfixGGG <- sortRefBfix[order(match(sortRefBfix$CHROM, allChr)), ]
  dfRefBformatGGG <- sortRefBformat[order(match(sortRefBfix$CHROM, allChr)), ]
  
  lGen <- length(dfRefAformatGGG$GT)
  # apply
  # initialized with 0
  intGenRefA <- integer(length = lGen)
  intGenRefB <- integer(length = lGen)
  
  for (indG in 1:lGen) {
    intGenRefA[indG] <- switch(dfRefAformatGGG$GT[indG], 
                             "0/0" = 0, "0/1" = 1, "1/1" = 2)
    intGenRefB[indG] <- switch(dfRefBformatGGG$GT[indG],
                             "0/0" = 0, "0/1" = 1, "1/1" = 2)
  }
  
  for (indC in allChr) {
    # extract chromosome length
    chrStrRefA <- unlist(strsplit(grep(paste(indC, "_", sep = ""), lstRefAvcfGGG$meta, value = T), split = ","))[2]
    chromLenRefA <- as.numeric(gsub(pattern = "[[:alpha:][:punct:]]", x = chrStrRefA, replacement = ""))
    chrStrRefB <- unlist(strsplit(grep(paste(indC, "_", sep = ""), lstRefBvcfGGG$meta, value = T), split = ","))[2]
    chromLenRefB <- as.numeric(gsub(pattern = "[[:alpha:][:punct:]]", x = chrStrRefB, replacement = ""))
    
    indVarChrRefA <- which(dfRefAfixGGG$CHROM == indC)
    intGenChrRefA <- intGenRefA[indVarChrRefA]
    chrRefAfixGGG <- dfRefAfixGGG[indVarChrRefA, ]
    
    indVarChrRefB <- which(dfRefBfixGGG$CHROM == indC)
    intGenChrRefB <- intGenRefB[indVarChrRefB]
    chrRefBfixGGG <- dfRefBfixGGG[indVarChrRefB, ]
    
    resRefA <- rle(intGenChrRefA)
    resRefB <- rle(intGenChrRefB)
    
    nMarkerRefA <- sum(resRefA$lengths)
    nMarkerRefB <- sum(resRefB$lengths)
    if (nMarkerRefA == 0) {
      cat(paste("No marker", ref1, "found in chromosome", indC, "\n"))
      next()
    }
    if (nMarkerRefB == 0) {
      cat(paste("No marker", ref2, "found in chromosome", indC, "\n"))
      next()
    }
    # Het fraction RefA
    nHEt <- sum(resRefA$lengths[which(resRefA$values == 1)])
    ratHH <- nHEt / nMarkerRefA
    cat(paste("Het fraction", ind1, indC, paste0(ref1, ":")), SpecDec(ratHH, 5), "\n")
    # mrkRatRefA <- nMarkerRefA / chrLenRefA[counterChr] / vctDenRefA[counterChr]
    # cat("Marker ratio", ind1, indC, paste0(ref1, ":"), mrkRatRefA, nMarkerRefA, length(intGenChrRefA), "\n")
    # Het fraction RefB
    nHEt <- sum(resRefB$lengths[which(resRefB$values == 1)])
    ratHH <- nHEt / nMarkerRefB
    cat(paste("Het fraction", ind1, indC, paste0(ref2, ":")), SpecDec(ratHH, 5), "\n")
    # mrkRatRefB <- nMarkerRefB / chrLenRefB[counterChr] / vctDenRefB[counterChr]
    # cat("Marker ratio", ind1, indC, paste0(ref2, ":"), mrkRatRefB, nMarkerRefB, length(intGenChrRefB), "\n")
    # events map
    lEventsRefA <- length(resRefA$length)
    lEventsRefB <- length(resRefB$length)
    
    # make events RefA
    evRefA <- data.frame(
      chr = rep(indC, lEventsRefA), 
      start = numeric(length = lEventsRefA), 
      first = numeric(length = lEventsRefA), 
      last = numeric(length = lEventsRefA), 
      end = numeric(length = lEventsRefA), 
      status = resRefA$values, 
      len = as.numeric(resRefA$lengths), # number of SNPs
      denES = numeric(length = lEventsRefA), 
      evrES = numeric(length = lEventsRefA), 
      distES = numeric(length = lEventsRefA), 
      denLF = numeric(length = lEventsRefA), 
      evrLF = numeric(length = lEventsRefA), 
      distLF = numeric(length = lEventsRefA))
    evRefA$start[1] <- 1
    evRefA$first[1] <- as.numeric(chrRefAfixGGG$POS[1])
    evRefA$end[lEventsRefA] <- chromLenRefA
    evRefA$last[lEventsRefA] <- as.numeric(chrRefAfixGGG$POS[nrow(chrRefAfixGGG)])
    
    if (lEventsRefA > 1) {
      # indE: index of events
      for (indE in 2:lEventsRefA) {
        # ultimo snp dello stato precedente
        evRefA$last[indE - 1] <- as.numeric(chrRefAfixGGG$POS[sum(resRefA$lengths[1:c(indE - 1)])])
        # primo snp dello stato corrente
        evRefA$first[indE] <- as.numeric(chrRefAfixGGG$POS[sum(resRefA$lengths[1:c(indE - 1)]) + 1])
        # end of former state
        evRefA$end[indE - 1] <- floor((evRefA$first[indE] + evRefA$last[indE - 1]) / 2)
        # start of current state
        evRefA$start[indE] <- floor((evRefA$first[indE] + evRefA$last[indE - 1]) / 2) + 1
      }
    }
    # N snps / bp
    evRefA$distES <- c(evRefA$end - evRefA$start)
    evRefA$denES <- evRefA$len / evRefA$distES
    # +1 to take into account events with 1 SNP
    evRefA$distLF <- c(evRefA$last - evRefA$first + 1)
    evRefA$denLF <- evRefA$len / evRefA$distLF
    # bp / N snps
    evRefA$evrES <- evRefA$distES / evRefA$len
    evRefA$evrLF <- evRefA$distLF / evRefA$len
    # format events table
    evRefA$denES <- as.numeric(format(evRefA$denES, digits = 3))
    evRefA$evrES <- round(evRefA$evrES)
    evRefA$denLF <- as.numeric(format(evRefA$denLF, digits = 3))
    evRefA$evrLF <- round(evRefA$evrLF)
    
    # make events RefB
    evRefB <- data.frame(
      chr = rep(indC, lEventsRefB), 
      start = numeric(length = lEventsRefB), 
      first = numeric(length = lEventsRefB), 
      last = numeric(length = lEventsRefB), 
      end = numeric(length = lEventsRefB), 
      status = resRefB$values, 
      len = as.numeric(resRefB$lengths), # number of SNPs
      denES = numeric(length = lEventsRefB), 
      evrES = numeric(length = lEventsRefB), 
      distES = numeric(length = lEventsRefB), 
      denLF = numeric(length = lEventsRefB), 
      evrLF = numeric(length = lEventsRefB), 
      distLF = numeric(length = lEventsRefB))
    evRefB$start[1] <- 1
    evRefB$first[1] <- as.numeric(chrRefBfixGGG$POS[1])
    evRefB$end[lEventsRefB] <- chromLenRefB
    evRefB$last[lEventsRefB] <- as.numeric(chrRefBfixGGG$POS[nrow(chrRefBfixGGG)])
    
    if (lEventsRefB > 1) {
      # indE: index of events
      for (indE in 2:lEventsRefB) {
        # ultimo snp dello stato precedente
        evRefB$last[indE - 1] <- as.numeric(chrRefBfixGGG$POS[sum(resRefB$lengths[1:c(indE - 1)])])
        # primo snp dello stato corrente
        evRefB$first[indE] <- as.numeric(chrRefBfixGGG$POS[sum(resRefB$lengths[1:c(indE -1)]) + 1])
        # end of former state
        evRefB$end[indE - 1] <- floor((evRefB$first[indE] + evRefB$last[indE - 1]) / 2)
        # start of current state
        evRefB$start[indE] <- floor((evRefB$first[indE] + evRefB$last[indE - 1]) / 2) + 1
      }
    }
    # N snps / bp
    evRefB$distES <- c(evRefB$end - evRefB$start)
    evRefB$denES <- evRefB$len / evRefB$distES
    # +1 to take into account events with 1 SNP
    evRefB$distLF <- c(evRefB$last - evRefB$first + 1)
    evRefB$denLF <- evRefB$len / evRefB$distLF
    # bp / N snps
    evRefB$evrES <- evRefB$distES / evRefB$len
    evRefB$evrLF <- evRefB$distLF / evRefB$len
    # format events table
    evRefB$denES <- as.numeric(format(evRefB$denES, digits = 3))
    evRefB$evrES <- round(evRefB$evrES)
    evRefB$denLF <- as.numeric(format(evRefB$denLF, digits = 3))
    evRefB$evrLF <- round(evRefB$evrLF)
    
    # save data for AllSeg.R
    pathOutRefA <- file.path(sampOutDir, paste(ind1, indC, ref1, "Seg.RData", sep = "."))
    varOutName1 <- paste0("ev", ref1)
    assign(varOutName1, evRefA)
    save(list = varOutName1, file = pathOutRefA)
    pathOutRefB <- file.path(sampOutDir, paste(ind1, indC, ref2, "Seg.RData", sep = "."))
    varOutName2 <- paste0("ev", ref2)
    assign(varOutName2, evRefB)
    save(list = varOutName2, file = pathOutRefB)
    
    # plotting segments
    indChr <- which(allChr == indC)
    yRefA <- .25
    yRefB <- .75
    
    # invert!!!
    nrRefA <- nrow(evRefA)
    color <- character(length = nrRefA)
    for (indSw in 1:nrRefA) {
      color[indSw] <- switch(as.character(evRefA$status[indSw]), 
                             "0"="blue", 
                             # 0/1
                             "1"="darkgrey", 
                             # 1/1
                             "2"="red")
    }
    color <- factor(color)
    yVal <- rep(yRefA, nrRefA)
    evRefA <- data.frame(evRefA, color, yVal)
    
    # invert!!!
    nrRefB <- nrow(evRefB)
    color <- character(length = nrRefB)
    for (indSw in 1:nrRefB) {
      color[indSw] <- switch(as.character(evRefB$status[indSw]), 
                             "0"="red", 
                             # 0/1
                             "1"="darkgrey", 
                             # 1/1
                             "2"="blue")
    }
    color <- factor(color)
    yVal <- rep(yRefB, nrRefB)
    evRefB <- data.frame(evRefB, color, yVal)
    # all events
    strain <- c(rep(ref1, nrow(evRefA)), rep(ref2, nrow(evRefB)))
    allEv <- rbind(evRefA, evRefB)
    
    # single sample single chromosome plot
    
    # size markers
    szCentr <- 8
    szEvents <- 10
    szChrom <- 2
    
    # shifting RefA data
    indRefA <- which(allEv$yVal == .25)
    allEv[indRefA, c(2:5)] <- allEv[indRefA, c(2:5)] + shiftRefBRefA[indChr]
    
    # plot events from first SNP to last SNP
    pGG <- ggplot(allEv) + 
      # title
      ggtitle(ind1, subtitle = indC) + 
      coord_cartesian(ylim = c(.0, 1.)) + 
      scale_y_continuous(breaks = NULL, labels = NULL) + 
      annotate("text", x = 8, y = 0.95, label = ref2, color = "red", size = theme_get()$text[["size"]] / 1.5) + 
      annotate("text", x = 8, y = 0.45, label = ref1, color = "blue", size = theme_get()$text[["size"]] / 1.5) + 
      # whole chromosome
      geom_segment(aes(x = 1, y = yRefB, xend = chromLenRefB, yend = yRefB), 
                   colour = "black", size =  szChrom) + 
      geom_segment(aes(x = shiftRefBRefA[indChr], y = yRefA, xend = (chromLenRefA + shiftRefBRefA[indChr]), yend = yRefA), 
                   colour = "black", size =  szChrom) + 
      # (sub)telomers
      geom_segment(aes(x = 1, y = yRefB, xend = subTelLRefB[indChr], yend = yRefB), 
                   colour = "orange", size =  szChrom) + 
      geom_segment(aes(x = subTelRRefB[indChr], y = yRefB, xend = chromLenRefB, yend = yRefB), 
                   colour = "orange", size =  szChrom) + 
      geom_segment(aes(x = shiftRefBRefA[indChr], y = yRefA, xend = (subTelLRefA[indChr] + shiftRefBRefA[indChr]), yend = yRefA), 
                   colour = "orange", size =  szChrom) + 
      geom_segment(aes(x = (subTelRRefA[indChr] + shiftRefBRefA[indChr]), y = yRefA, 
                       xend = (chromLenRefA + shiftRefBRefA[indChr]), yend = yRefA), 
                   colour = "orange", size =  szChrom) + 
      # events
      geom_segment(aes(x = first, y = allEv$yVal, xend = last, yend = allEv$yVal), 
                   colour = allEv$color, size = szEvents) + 
      # centromere
      geom_point(shape = 1, aes(x = (centrSRefB[indChr] + centrERefB[indChr]) / 2, y = yRefB), size = szCentr) + 
      geom_point(shape = 1, aes(x = (centrSRefA[indChr] + centrERefA[indChr]) / 2 + shiftRefBRefA[indChr], y = yRefA), 
                 size = szCentr) + 
      # plain background
      theme(plot.margin = unit(c(1, 2, 1, 1), "cm"), 
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            panel.border = element_blank(), 
            panel.background = element_blank(), 
            # title style
            plot.title = element_text(size = 16, face = "bold", hjust = 1), 
            plot.subtitle = element_text(size = 12, face = "bold", hjust = 1), 
            # axis style
            axis.title = element_text(size = 16, face = "bold"), 
            axis.text = element_text(size = 16)) + 
      xlab(paste0("(", labAssemblyRiferimento, ") Position [bp]")) + 
      ylab(NULL) + 
      # integer x-axis values
      scale_x_continuous(labels = comma)
    pPath <- file.path(sampOutDir, paste(ind1, indC, "SegFL.pdf", sep = "."))
    pdf(file = pPath, width = 14, height = 4)
    print(pGG)
    dev.off()
    
    # plot events from start position to end position
    pGG <- ggplot(allEv) + 
      # title
      ggtitle(ind1, subtitle = indC) + 
      coord_cartesian(ylim = c(.0, 1.)) + 
      scale_y_continuous(breaks = NULL, labels = NULL) + 
      annotate("text", x = 8, y = 0.95, label = ref2, color = "red", size = theme_get()$text[["size"]] / 1.5) + 
      annotate("text", x = 8, y = 0.45, label = ref1, color = "blue", size = theme_get()$text[["size"]] / 1.5) + 
      # events
      geom_segment(aes(x = start, y = allEv$yVal, xend = end, yend = allEv$yVal), 
                   colour = allEv$color, size = szEvents) + 
      # whole chromosome
      geom_segment(aes(x = 1, y = yRefB, xend = chromLenRefB, yend = yRefB), 
                   colour = "black", size =  szChrom) + 
      geom_segment(aes(x = shiftRefBRefA[indChr], y = yRefA, xend = (chromLenRefA + shiftRefBRefA[indChr]), yend = yRefA), 
                   colour = "black", size =  szChrom) + 
      # (sub)telomers
      geom_segment(aes(x = 1, y = yRefB, xend = subTelLRefB[indChr], yend = yRefB), 
                   colour = "orange", size =  szChrom) + 
      geom_segment(aes(x = subTelRRefB[indChr], y = yRefB, xend = chromLenRefB, yend = yRefB), 
                   colour = "orange", size =  szChrom) + 
      geom_segment(aes(x = shiftRefBRefA[indChr], y = yRefA, xend = (subTelLRefA[indChr] + shiftRefBRefA[indChr]), yend = yRefA), 
                   colour = "orange", size =  szChrom) + 
      geom_segment(aes(x = (subTelRRefA[indChr] + shiftRefBRefA[indChr]), y = yRefA, 
                       xend = (chromLenRefA + shiftRefBRefA[indChr]), yend = yRefA), 
                   colour = "orange", size =  szChrom) + 
      # centromere
      geom_point(shape = 1, aes(x = (centrSRefB[indChr] + centrERefB[indChr]) / 2, y = yRefB), size = szCentr) + 
      geom_point(shape = 1, aes(x = (centrSRefA[indChr] + centrERefA[indChr]) / 2 + shiftRefBRefA[indChr], y = yRefA), 
                 size = szCentr) + 
      # plain background
      theme(plot.margin = unit(c(1, 2, 1, 1), "cm"), 
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            panel.border = element_blank(), 
            panel.background = element_blank(), 
            # title style
            plot.title = element_text(size = 16, face = "bold", hjust = 1), 
            plot.subtitle = element_text(size = 12, face = "bold", hjust = 1), 
            # axis style
            axis.title = element_text(size = 16, face = "bold"), 
            axis.text = element_text(size = 16)) + 
      xlab(paste0("(", labAssemblyRiferimento, ") Position [bp]")) + 
      ylab(NULL) + 
      # integer x-axis values
      scale_x_continuous(labels = comma)
    pPath <- file.path(sampOutDir, paste(ind1, indC, "SegSE.pdf", sep = "."))
    pdf(file = pPath, width = 14, height = 4)
    print(pGG)
    dev.off()
  } # chromosome loop
} # sample loop

cat("Global time", "\n")
proc.time() - ptmGlob


