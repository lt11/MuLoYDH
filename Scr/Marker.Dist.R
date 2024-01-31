# header ------------------------------------------------------------------

# riordina sulla posizione prima di diffPosRefA <- diff(posGG)
# perché se ci sono mapping inter-chr non torna una cippa!

# Density distribution of consecutive marker distances

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

# label tables marker density
outLabelA <- "Density.Table"
# label plots
outLabelB <- "Distance"

# load the dataframe of core markers, the good guys
load(file.path(baseDir, "LOH", "Markers", paste0("Markers.", ref1, "-", ref2, ".RData")))

# length of chromosomes
chrLenRefA <- read.table(file = file.path(baseDir, "CNV", "GCdata", ref1Label, "LenChr.txt"), header = F, sep = "\t")[, 3]
chrLenRefB <- read.table(file = file.path(baseDir, "CNV", "GCdata", ref2Label, "LenChr.txt"), header = F, sep = "\t")[, 3]

outDir <- file.path(baseDir, "LOH", "Markers")

allChr <- c("chrI", "chrII", "chrIII", "chrIV", "chrV", "chrVI", "chrVII", "chrVIII", 
            "chrIX", "chrX", "chrXI", "chrXII", "chrXIII", "chrXIV", "chrXV", "chrXVI")

# remove MT markers
outMT <- unique(which(dfGG$chrRefA == "chrMT"), which(dfGG$chrRefB == "chrMT"))

if (length(outMT) != 0) {
  dfSNP <- dfGG[-outMT, ]
} else {
  dfSNP <- dfGG
}

# chr and pos c(1, 2) for RefA (RefAsm1): sort by chr
dfSNPRefAsm1 <- dfSNP[order(match(dfSNP[, 1], allChr)), ][, c(1,2)]
# chr and pos c(3, 4) for RefB (RefAsm2): sort by chr
dfSNPRefAsm2 <- dfSNP[order(match(dfSNP[, 3], allChr)), ][, c(3,4)]

# header tables marker density
fileDenRefA <- file.path(outDir, paste(outLabelA, ref1, "txt", sep = "."))
fileDenRefB <- file.path(outDir, paste(outLabelA, ref2, "txt", sep = "."))

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
      plot.margin = unit(c(1, 1, 1, 1), "cm"), 
      # title style
      plot.title = element_text(size = 16, hjust = 1), 
      plot.subtitle = element_text(size = 12, face = "bold", hjust = 1), 
      # axis style
      axis.title = element_text(size = 16, face = "bold"), 
      axis.text = element_text(size = 16))
  plo1PathRefA <- file.path(outDir, paste(outLabelB, ref1, indC, "pdf", sep = "."))
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
      plot.margin = unit(c(1, 1, 1, 1), "cm"), 
      # title style
      plot.title = element_text(size = 16, hjust = 1), 
      plot.subtitle = element_text(size = 12, face = "bold", hjust = 1), 
      # axis style
      axis.title = element_text(size = 16, face = "bold"), 
      axis.text = element_text(size = 16))
  plo1PathRefA <- file.path(outDir, paste(outLabelB, ref1, indC, "500.pdf", sep = "."))
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
      plot.margin = unit(c(1, 1, 1, 1), "cm"), 
      # title style
      plot.title = element_text(size = 16, hjust = 1), 
      plot.subtitle = element_text(size = 12, face = "bold", hjust = 1), 
      # axis style
      axis.title = element_text(size = 16, face = "bold"), 
      axis.text = element_text(size = 16))
  plo1PathRefB <- file.path(outDir, paste(outLabelB, ref2, indC, "pdf", sep = "."))
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
      plot.margin = unit(c(1, 1, 1, 1), "cm"), 
      # title style
      plot.title = element_text(size = 16, hjust = 1), 
      plot.subtitle = element_text(size = 12, face = "bold", hjust = 1), 
      # axis style
      axis.title = element_text(size = 16, face = "bold"), 
      axis.text = element_text(size = 16))
  plo1PathRefB <- file.path(outDir, paste(outLabelB, ref2, indC, "500.pdf", sep = "."))
  pdf(file = plo1PathRefB)
  print(plot1RefAsm2, width = 10, height = 26)
  dev.off()
  
  counterChr <- counterChr + 1
}

# RefAsm1 dataframe for ggplot 
dfGGDistRefA <- data.frame(accChrRefA, accDiffRefA)
colnames(dfGGDistRefA) <- c("chr", "dist")
# RefAsm2 dataframe for ggplot 
dfGGDistRefB <- data.frame(accChrRefB, accDiffRefB)
colnames(dfGGDistRefB) <- c("chr", "dist")

# make factor to set facet panels order
dfGGDistRefA$chr <- factor(dfGGDistRefA$chr, levels = allChr)
dfGGDistRefB$chr <- factor(dfGGDistRefB$chr, levels = allChr)
dfGGDistRefA <- dfGGDistRefA[order(match(dfGGDistRefA$chr, allChr)), ]
dfGGDistRefB <- dfGGDistRefB[order(match(dfGGDistRefB$chr, allChr)), ]

# RefAsm1 markers distant < 0.5 kb
dfGGDistRefA500 <- dfGGDistRefA[which(dfGGDistRefA$dist < 500), ]
bWidth <- 10
nEv <- nrow(dfGGDistRefA500)
plot3RefAsm1 <- ggplot(dfGGDistRefA500, aes(x = dist)) + 
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
    plot.margin = unit(c(1, 1, 1, 1), "cm"), 
    # title style
    plot.title = element_text(size = 16, hjust = 1), 
    plot.subtitle = element_text(size = 12, face = "bold", hjust = 1), 
    # axis style
    axis.title = element_text(size = 16, face = "bold"), 
    axis.text = element_text(size = 16))
plo3PathRefA <- file.path(outDir, paste(outLabelB, ref1, "500.pdf", sep = "."))
pdf(file = plo3PathRefA, width = 10, height = 26)
print(plot3RefAsm1)
dev.off()

# RefAsm2 markers distant < 0.5 kb
dfGGDistRefB500 <- dfGGDistRefB[which(dfGGDistRefB$dist < 500), ]
bWidth <- 10
nEv <- nrow(dfGGDistRefB500)
plot3RefAsm2 <- ggplot(dfGGDistRefB500, aes(x = dist)) + 
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
    plot.margin = unit(c(1, 1, 1, 1), "cm"), 
    # title style
    plot.title = element_text(size = 16, hjust = 1), 
    plot.subtitle = element_text(size = 12, face = "bold", hjust = 1), 
    # axis style
    axis.title = element_text(size = 16, face = "bold"), 
    axis.text = element_text(size = 16))
plo3PathRefB <- file.path(outDir, paste(outLabelB, ref2, "500.pdf", sep = "."))
pdf(file = plo3PathRefB, width = 10, height = 26)
print(plot3RefAsm2)
dev.off()


