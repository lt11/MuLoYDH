# header ------------------------------------------------------------------

# plotting allele frequencies

# nota bene che AF1 è una 
# "Max-likelihood estimate of the first ALT allele frequency"
# quindi non va bene: 
# > unique(allFreq)
# [1] 0.5 1.0 0.0
# usiamo DP4: Number of high-quality ref-forward , ref-reverse, alt-forward and alt-reverse bases

rm(list = ls())
options(stringsAsFactors = F)
library(ggplot2)

# function(s)
SpecDec <- function(x, k) as.numeric(format(round(x, k), nsmall = k))

# settings ----------------------------------------------------------------

# arguments
argsVal <- commandArgs(trailingOnly = T)
ref1Label <- argsVal[1]
ref2Label <- argsVal[2]
ref1 <- argsVal[3]
ref2 <- argsVal[4]
baseDir <- argsVal[5]

# output folder
dirOutPlot <- file.path(baseDir, "LOH", "Summary", "Plot", "BAF")
dir.create(path = dirOutPlot, showWarnings = F, recursive = T)

# load markers positions
markerFile <- file.path(baseDir, "LOH", "Markers", paste0("Markers.", ref1, "-", ref2, ".RData"))
load(markerFile)

# leggo gli ID dei campioni dalla cartella Markers: ID.ref1.vcf.gz
markerCallsDir <- file.path(baseDir, "Markers")
allIDs <- sapply(strsplit(x = list.files(path = markerCallsDir, pattern = "\\.vcf\\.gz$"), split = "\\."), "[[", 1)
sampleIDs <- unique(grep(pattern = "^PAR", x = allIDs, invert = T, value = T))

for (indS in sampleIDs) {
  # load vcf data
  fileRef1 <- file.path(baseDir, "LOH", indS, paste(indS, ref1, "VcfFilt.RData", sep = "."))
  load(fileRef1)
  fileRef2 <- file.path(baseDir, "LOH", indS, paste(indS, ref2, "VcfFilt.RData", sep = "."))
  load(fileRef2)
  
  # parse depth data ref1
  unlstInfo <- unlist(strsplit(lstRefAvcfGGG$fix$INFO, split = ";", fixed = T))
  strDepth4 <- grep(pattern = "DP4", x = unlstInfo, value = T)
  commaDepth4 <- sapply(strsplit(strDepth4, split = "="), "[[", 2)
  nVar <- nrow(lstRefAvcfGGG$fix)
  numDepth4 <- matrix(unlist(strsplit(commaDepth4, split = ",")), nrow = nVar, byrow = T)
  numDepth4 <- apply(numDepth4, 2, as.numeric)
  # colnames(numDepth4) <- c("dpRefFwd", "dpRefRev", "dpAltFwd", "dpAltRev")
  allDepth <- rowSums(numDepth4)
  betaAlleleFractionRef1 <- (numDepth4[, 3] + numDepth4[, 4]) / allDepth
  
  # parse depth data ref2
  unlstInfo <- unlist(strsplit(lstRefBvcfGGG$fix$INFO, split = ";", fixed = T))
  strDepth4 <- grep(pattern = "DP4", x = unlstInfo, value = T)
  commaDepth4 <- sapply(strsplit(strDepth4, split = "="), "[[", 2)
  nVar <- nrow(lstRefBvcfGGG$fix)
  numDepth4 <- matrix(unlist(strsplit(commaDepth4, split = ",")), nrow = nVar, byrow = T)
  numDepth4 <- apply(numDepth4, 2, as.numeric)
  # colnames(numDepth4) <- c("dpRefFwd", "dpRefRev", "dpAltFwd", "dpAltRev")
  allDepth <- rowSums(numDepth4)
  betaAlleleFractionRef2 <- (numDepth4[, 3] + numDepth4[, 4]) / allDepth
  
  # make dataframe for ggplot
  dfBAF <- data.frame(betaAlleleFractionRef1, betaAlleleFractionRef2, 
                      lstRefAvcfGGG$fix$CHROM, lstRefAvcfGGG$fix$POS, lstRefAvcfGGG$format$GT, 
                      lstRefBvcfGGG$fix$CHROM, lstRefBvcfGGG$fix$POS, lstRefBvcfGGG$format$GT)
  colnames(dfBAF) <- c( "BAF1", "BAF2", "chromRef1", "posRef1", "GT1", "chromRef2", "posRef2", "GT2")
  
  # colori (basta usare un reference, tanto sono già controllati da ClrS)
  # ma se proprio vuoi controllare
  # length(indRed) + length(indBlue) + length(indGrey)
  indGrey <- which(dfBAF$GT1 == "0/1")
  indBlue <- which(dfBAF$GT1 == "0/0")
  indRed <- which(dfBAF$GT1 == "1/1")
  
  vctColor <- character(length = nVar)
  vctColor[indGrey] <- "grey70"
  vctColor[indBlue] <- "blue2"
  vctColor[indRed] <- "red3"
  vctColor <- factor(vctColor)
  
  dfBAF <- data.frame(dfBAF, vctColor)
  
  plotTolo <- ggplot(dfBAF, aes(x = BAF2, y = BAF1, shape = "a")) +
    # labs(color = "Chromosomes\nby length") +
    xlab(paste0("BAF ", ref2)) + 
    ylab(paste0("BAF ", ref1)) + 
    scale_shape_discrete(solid = F) + 
    geom_point(color = vctColor, show.legend = F, alpha = .4) + 
    theme_minimal() + 
    geom_smooth(method = lm, color = "black") + 
    annotate(x = .75, y = .75, 
             label = paste("R = ", SpecDec(cor(dfBAF$BAF1, dfBAF$BAF2), 3)), 
             geom = "text", size = 8, color = "black") + 
    theme(axis.title = element_text(size = 16, face = "bold"), 
          axis.text = element_text(size = 16), 
          axis.title.x = element_text(colour = "red"), 
          axis.title.y = element_text(colour = "blue"))
  
  fileOutPlot <- file.path(dirOutPlot, paste0(indS, ".BAF.pdf"))
  pdf(file = fileOutPlot)
  print(plotTolo)
  dev.off()
}


