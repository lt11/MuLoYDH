# header ------------------------------------------------------------------

# annotation of genes involved in all segments

rm(list = ls())
options(stringsAsFactors = F)
library(IRanges)

# settings ----------------------------------------------------------------

# arguments
argsVal <- commandArgs(trailingOnly = T)
ref1Label <- argsVal[1]
ref2Label <- argsVal[2]
ref1 <- argsVal[3]
ref2 <- argsVal[4]
baseDir <- argsVal[5]

# paths
pathRefA <- file.path(baseDir, "Ref", "Ann", paste0(ref1Label, ".all_feature.gff"))
pathRefB <- file.path(baseDir, "Ref", "Ann", paste0(ref2Label, ".all_feature.gff"))
segmeDir <- "AllSegments"
allEvFile <- "AllEvents.csv"
outFile <- "Annotations.Genes.txt"
headerString <- c("chr\tstart\tfirst\tlast\tend\tstatus\tlen\tdenES\tevrES\tdistES\tdenLF\tevrLF\tdistLF\tstrain\thybrid\tgeneF\tgeneL\tgenesO")
outFile <- "Annotations.Genes.txt"
eventsDir <- file.path(baseDir, "LOH")
outDirGenes <- file.path(eventsDir, "Annotations")
write.table(x = headerString, file = file.path(outDirGenes, outFile), 
            append = F, sep = "", row.names = F, col.names = F, quote = F)

# load all events
allEvents <- read.table(file = file.path(eventsDir, segmeDir, allEvFile), 
                        header = T, na.strings = "antani", sep = ",")
# load references annotations
lstAnnoGene <- list()
ref1Anno <- read.table(file = pathRefA, header = F)
ref2Anno <- read.table(file = pathRefB, header = F)
refAnnoS <- list(ref1Anno, ref2Anno)
for (indA in c(1:2)) {
  dfAnno <- refAnnoS[[indA]]
  dfAnnoGene <- dfAnno[which(dfAnno$V3 == "gene"), ]
  V9 <- gsub(pattern = "Name=", replacement = "", 
             sapply(strsplit(dfAnnoGene$V9, split = ";"), "[[", 2))
  dfAnnoGeneName <- data.frame(dfAnnoGene[, c(1:8)], V9)
  # list with reference gene annotation dataframes
  lstAnnoGene[[indA]] <- dfAnnoGeneName
}

myRef <- c(ref1, ref2)

for (indR in c(1:2)) {
  eventsRef <- allEvents[which(allEvents$strain == myRef[indR]), ]
  for (indE in c(1:nrow(eventsRef))) {
    # broken gene first marker
    indBroken <- which(lstAnnoGene[[indR]][, 1] == eventsRef[indE, 1] & 
                         lstAnnoGene[[indR]][, 4] < eventsRef[indE, 3] & 
                         lstAnnoGene[[indR]][, 5] > eventsRef[indE, 3])
    brokenGeneF <- "-"
    if (length(indBroken) == 1) {
      brokenGeneF <- lstAnnoGene[[indR]]$V9[indBroken]
    }
    # broken gene last marker
    indBroken <- which(lstAnnoGene[[indR]][, 1] == eventsRef[indE, 1] & 
                         lstAnnoGene[[indR]][, 4] < eventsRef[indE, 4] & 
                         lstAnnoGene[[indR]][, 5] > eventsRef[indE, 4])
    brokenGeneE <- "-"
    if (length(indBroken) == 1) {
      brokenGeneE <- lstAnnoGene[[indR]]$V9[indBroken]
    }
    # overlapping genes
    eventQuery <- IRanges(eventsRef[indE, 3], eventsRef[indE, 4])
    lstAnnoGeneChr <- lstAnnoGene[[indR]][which(lstAnnoGene[[indR]][, 1] == eventsRef[indE, 1]), ]
    geneSubects <- IRanges(lstAnnoGeneChr[, 4], lstAnnoGeneChr[, 5])
    overlapQS <- findOverlaps(eventQuery, geneSubects)
    vctGenes <- lstAnnoGeneChr[as.matrix(overlapQS)[, 2], 9]
    genesOverlap <- "-"
    if (length(vctGenes) != 0) {
      genesOverlap <- paste(vctGenes, collapse = ",")
    }
    evAnnoGene <- data.frame(eventsRef[indE, 1:15], brokenGeneF, brokenGeneE, genesOverlap)
    write.table(x = evAnnoGene, file = file.path(outDirGenes, outFile), 
                sep = "\t", row.names = F, col.names = F, quote = F, append = T)
  }
}



