# header ------------------------------------------------------------------

# mancano le annotazioni delle features upstream al primo evento

# annotation of genes involved in all segments
# from allEvents file
# V2 annota anche start e end delle feature
# (sia quelle broken che quelle overlap)
# whole-chromosome segments are NOT annotated for start and end coordinates of the features
# tanto non servono a una cippa

# scrive nelle cartelle generate da StatLOH.V1.R

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
inputDir <- "AllSegments"
allEvFile <- "AllEvents.csv"
eventsDir <- file.path(baseDir, "LOH")
outDir <- file.path(eventsDir, "Summary", "Tables")
dir.create(outDir, recursive = T, showWarnings = F)
outFile <- file.path(outDir, "LOH.features.txt")
fieldsNameGff <- c("chrom", "ref", "feature", "start", "end", "score", "strand", "frame", "attribute")

allChr <- c("chrI", "chrII", "chrIII", "chrIV", "chrV", "chrVI", "chrVII", "chrVIII", 
            "chrIX", "chrX", "chrXI", "chrXII", "chrXIII", "chrXIV", "chrXV", "chrXVI")

# load all events
allEvents <- read.table(file = file.path(eventsDir, inputDir, allEvFile), 
                        header = T, na.strings = "antani", sep = ",")
# format allEvents table
allEventsFormat <- data.frame(allEvents$hybrid,
                              allEvents$strain,
                              allEvents[, 1:5],
                              allEvents[, 6:13],
                              allEvents[, 16:17])

colnames(allEventsFormat) <- c("samples", "ref", colnames(allEvents[, 1:5]), 
                               colnames(allEvents[, 6:13]), colnames( allEvents[, 16:17]))
# sort table
# sort by chrom and "first" position
allEventsFormat <- allEventsFormat[order(match(allEventsFormat$chr, allChr), allEventsFormat$first), ]
# sort by parent
allEventsFormat <- allEventsFormat[order(match(allEventsFormat$ref, unique(allEventsFormat$ref))), ]
# sort by sample
allEventsFormat <- allEventsFormat[order(match(allEventsFormat$sample, unique(allEventsFormat$sample))), ]
# rebuild table
allEvents <- allEventsFormat

# init output file with header
cat(colnames(allEvents), file = outFile, append = F, sep = "\t")
cat("", "upstream feature(s)", "upstream f. start", "upstream f. end", 
    "downstream feature(s)", "downstream f. start", "downstream f. end", 
    "overlap feature(s)", "overlap f. start", "overlap f. end", 
    file = outFile, append = T, sep = "\t")
cat("\n", file = outFile, append = T)

# load references annotations
lstAnnoGene <- list()
ref1Anno <- read.table(file = pathRefA, header = F)
colnames(ref1Anno) <- fieldsNameGff
ref2Anno <- read.table(file = pathRefB, header = F)
colnames(ref2Anno) <- fieldsNameGff
refAnnoS <- list(ref1Anno, ref2Anno)
for (indA in c(1:2)) {
  dfAnno <- refAnnoS[[indA]]
  # togli le colonne ridondanti delle annotazioni "gene"
  dfAnnoGene <- dfAnno[-which(dfAnno$feature == "CDS" | dfAnno$feature == "exon" | dfAnno$feature == "mRRefB"), ]
  bonAttribute <- gsub(pattern = "Name=", replacement = "", 
                       sapply(strsplit(dfAnnoGene$attribute, split = ";"), "[[", 2))
  dfAnnoGeneName <- data.frame(dfAnnoGene[, c(1:8)], bonAttribute)
  # list with reference gene annotation dataframes
  lstAnnoGene[[indA]] <- dfAnnoGeneName
}

myRef <- c(ref1, ref2)

for (indR in c(1:2)) {
  # gli eventi contro un reference
  eventsRef <- allEvents[which(allEvents$ref == myRef[indR]), ]
  for (indS in unique(eventsRef$sample)) {
    # gli eventi di un campione
    eventsRefSample <- eventsRef[which(eventsRef$sample == indS), ]
    for (indC in unique(eventsRefSample$chr)) {
      # le features di un cromosoma
      lstAnnoGeneChr <- lstAnnoGene[[indR]][which(lstAnnoGene[[indR]]$chrom == indC), ]
      featuresQueries <- IRanges(lstAnnoGeneChr$start, lstAnnoGeneChr$end)
      # gli eventi di un cromosoma
      eventsRefSampleChrom <- eventsRefSample[which(eventsRefSample$chr == indC), ]
      numEvents <- nrow(eventsRefSampleChrom)
      if (numEvents == 1 & eventsRefSampleChrom$status[1] == 1) {
        brokenGeneF <- "-"
        brokenGeneFstart <- "-"
        brokenGeneFend <- "-"
        brokenGeneL <- "-"
        brokenGeneLstart <- "-"
        brokenGeneLend <- "-"
        genesOverlap <- paste(lstAnnoGeneChr$bonAttribute, collapse = ",")
        genesOverlapStart <- "all"
        genesOverlapEnd <- "all"
        evAnnoGene <- data.frame(eventsRefSampleChrom[1, ], brokenGeneF, brokenGeneFstart, brokenGeneFend, 
                                 brokenGeneL, brokenGeneLstart, brokenGeneLend, 
                                 genesOverlap, genesOverlapStart, genesOverlapEnd)
        write.table(x = evAnnoGene, file = outFile, 
                    sep = "\t", row.names = F, col.names = F, quote = F, append = T)
      } else if (numEvents == 1 & eventsRefSampleChrom$status[1] == 2) {
        brokenGeneF <- "-"
        brokenGeneFstart <- "-"
        brokenGeneFend <- "-"
        brokenGeneL <- "-"
        brokenGeneLstart <- "-"
        brokenGeneLend <- "-"
        genesOverlap <- paste(lstAnnoGeneChr$bonAttribute, collapse = ",")
        genesOverlapStart <- "all"
        genesOverlapEnd <- "all"
        evAnnoGene <- data.frame(eventsRefSampleChrom[1, ], brokenGeneF, brokenGeneFstart, brokenGeneFend, 
                                 brokenGeneL, brokenGeneLstart, brokenGeneLend, 
                                 genesOverlap, genesOverlapStart, genesOverlapEnd)
        write.table(x = evAnnoGene, file = outFile, 
                    sep = "\t", row.names = F, col.names = F, quote = F, append = T)
      } else if (numEvents == 1 & eventsRefSampleChrom$status[1] == 0) {
        brokenGeneF <- "-"
        brokenGeneFstart <- "-"
        brokenGeneFend <- "-"
        brokenGeneL <- "-"
        brokenGeneLstart <- "-"
        brokenGeneLend <- "-"
        genesOverlap <- paste(lstAnnoGeneChr$bonAttribute, collapse = ",")
        genesOverlapStart <- "all"
        genesOverlapEnd <- "all"
        evAnnoGene <- data.frame(eventsRefSampleChrom[1, ], brokenGeneF, brokenGeneFstart, brokenGeneFend, 
                                 brokenGeneL, brokenGeneLstart, brokenGeneLend, 
                                 genesOverlap, genesOverlapStart, genesOverlapEnd)
        write.table(x = evAnnoGene, file = outFile, 
                    sep = "\t", row.names = F, col.names = F, quote = F, append = T)
      } else {
        for (indE in c(1:numEvents)) {
          # inizializza i valori in caso non ci siano features in overlap
          # upstream o downstream con la LOH
          # o se il segmento è in eterozigosi
          brokenGeneF <- "-"
          brokenGeneFstart <- "-"
          brokenGeneFend <- "-"
          brokenGeneL <- "-"
          brokenGeneLstart <- "-"
          brokenGeneLend <- "-"
          flagStatus <- eventsRefSampleChrom$status[indE]
          if (flagStatus == 0 | flagStatus == 2) {
            # overlap features-segment
            eventSubject <- IRanges(eventsRefSampleChrom$first[indE], eventsRefSampleChrom$last[indE])
            overlapQS <- findOverlaps(featuresQueries, eventSubject, type = "within")
            vctGenes <- lstAnnoGeneChr[as.matrix(overlapQS)[, 1], 9]
            vctStart <- lstAnnoGeneChr[as.matrix(overlapQS)[, 1], 4]
            vctEnd <- lstAnnoGeneChr[as.matrix(overlapQS)[, 1], 5]
            genesOverlap <- "-"
            genesOverlapStart <- "-"
            genesOverlapEnd <- "-"
            if (length(vctGenes) != 0) {
              genesOverlap <- paste(vctGenes, collapse = ",")
              genesOverlapStart <- paste(vctStart, collapse = ",")
              genesOverlapEnd <- paste(vctEnd, collapse = ",")
            }
            # initial terminal LOH: allo start del cromosoma
            if (indE == 1) {
              # overlap downstream-features
              eventSubject <- IRanges(eventsRefSampleChrom$last[indE], eventsRefSampleChrom$first[indE + 1])
              overlapQS <- findOverlaps(featuresQueries, eventSubject, type = "any")
              vctGenes <- lstAnnoGeneChr[as.matrix(overlapQS)[, 1], 9]
              vctStart <- lstAnnoGeneChr[as.matrix(overlapQS)[, 1], 4]
              vctEnd <- lstAnnoGeneChr[as.matrix(overlapQS)[, 1], 5]
              if (length(vctGenes) != 0) {
                brokenGeneL <- paste(vctGenes, collapse = ",")
                brokenGeneLstart <- paste(vctStart, collapse = ",")
                brokenGeneLend <- paste(vctEnd, collapse = ",")
              }
              # final terminal LOH: alla fine del cromosoma
            } else if (indE == numEvents) {
              # overlap upstream-features
              eventSubject <- IRanges(eventsRefSampleChrom$last[indE - 1], eventsRefSampleChrom$first[indE])
              overlapQS <- findOverlaps(featuresQueries, eventSubject, type = "any")
              vctGenes <- lstAnnoGeneChr[as.matrix(overlapQS)[, 1], 9]
              vctStart <- lstAnnoGeneChr[as.matrix(overlapQS)[, 1], 4]
              vctEnd <- lstAnnoGeneChr[as.matrix(overlapQS)[, 1], 5]
              if (length(vctGenes) != 0) {
                brokenGeneF <- paste(vctGenes, collapse = ",")
                brokenGeneFstart <- paste(vctStart, collapse = ",")
                brokenGeneFend <- paste(vctEnd, collapse = ",")
              }
              # interstitial LOH  
            } else {
              # overlap upstream-features
              eventSubject <- IRanges(eventsRefSampleChrom$last[indE - 1], eventsRefSampleChrom$first[indE])
              overlapQS <- findOverlaps(featuresQueries, eventSubject, type = "any")
              vctGenes <- lstAnnoGeneChr[as.matrix(overlapQS)[, 1], 9]
              vctStart <- lstAnnoGeneChr[as.matrix(overlapQS)[, 1], 4]
              vctEnd <- lstAnnoGeneChr[as.matrix(overlapQS)[, 1], 5]
              if (length(vctGenes) != 0) {
                brokenGeneF <- paste(vctGenes, collapse = ",")
                brokenGeneFstart <- paste(vctStart, collapse = ",")
                brokenGeneFend <- paste(vctEnd, collapse = ",")
              }
              # overlap downstream-features
              eventSubject <- IRanges(eventsRefSampleChrom$last[indE], eventsRefSampleChrom$first[indE + 1])
              overlapQS <- findOverlaps(featuresQueries, eventSubject, type = "any")
              vctGenes <- lstAnnoGeneChr[as.matrix(overlapQS)[, 1], 9]
              vctStart <- lstAnnoGeneChr[as.matrix(overlapQS)[, 1], 4]
              vctEnd <- lstAnnoGeneChr[as.matrix(overlapQS)[, 1], 5]
              if (length(vctGenes) != 0) {
                brokenGeneL <- paste(vctGenes, collapse = ",")
                brokenGeneLstart <- paste(vctStart, collapse = ",")
                brokenGeneLend <- paste(vctEnd, collapse = ",")
              }
            }
            evAnnoGene <- data.frame(eventsRefSampleChrom[indE, ], brokenGeneF, brokenGeneFstart, brokenGeneFend, 
                                     brokenGeneL, brokenGeneLstart, brokenGeneLend, 
                                     genesOverlap, genesOverlapStart, genesOverlapEnd)
            write.table(x = evAnnoGene, file = outFile, 
                        sep = "\t", row.names = F, col.names = F, quote = F, append = T)
          } else {
            # eterozygous segments
            # overlap features-segment
            eventSubject <- IRanges(eventsRefSampleChrom$first[indE], eventsRefSampleChrom$last[indE])
            overlapQS <- findOverlaps(featuresQueries, eventSubject, type = "within")
            vctGenes <- lstAnnoGeneChr[as.matrix(overlapQS)[, 1], 9]
            vctStart <- lstAnnoGeneChr[as.matrix(overlapQS)[, 1], 4]
            vctEnd <- lstAnnoGeneChr[as.matrix(overlapQS)[, 1], 5]
            genesOverlap <- "-"
            genesOverlapStart <- "-"
            genesOverlapEnd <- "-"
            if (length(vctGenes) != 0) {
              genesOverlap <- paste(vctGenes, collapse = ",")
              genesOverlapStart <- paste(vctStart, collapse = ",")
              genesOverlapEnd <- paste(vctEnd, collapse = ",")
            }
            evAnnoGene <- data.frame(eventsRefSampleChrom[indE, ], brokenGeneF, brokenGeneFstart, brokenGeneFend, 
                                     brokenGeneL, brokenGeneLstart, brokenGeneLend, 
                                     genesOverlap, genesOverlapStart, genesOverlapEnd)
            write.table(x = evAnnoGene, file = outFile, 
                        sep = "\t", row.names = F, col.names = F, quote = F, append = T)
          }
        }
      }
    }
  }
}


