# header ------------------------------------------------------------------

# this makes the bed files of the homo-/heterozygous regions
# to call small variants accordingly

rm(list = ls())
options(stringsAsFactors = F)

# settings ----------------------------------------------------------------

# arguments
argsVal <- commandArgs(trailingOnly = T)
ref1Label <- argsVal[1]
ref2Label <- argsVal[2]
ref1 <- argsVal[3]
ref2 <- argsVal[4]
baseDir <- argsVal[5]

allChr <- c("chrI", "chrII", "chrIII", "chrIV", "chrV", "chrVI", "chrVII", "chrVIII", 
            "chrIX", "chrX", "chrXI", "chrXII", "chrXIII", "chrXIV", "chrXV", "chrXVI")

# output folder
outDir <- file.path(baseDir, "VariantCalls", "HomoHeteroReg")
dir.create(outDir, showWarnings = F, recursive = T)

# load input segments
inFileEvents <- file.path(baseDir, "LOH", "AllSegments", "AllEvents.csv")
tabEvents <- read.table(inFileEvents, header = T, sep = ",", na.strings = "antani")
allSamples <- unique(tabEvents$hybrid)

# filter status 2
tabEvents <- tabEvents[which(tabEvents$status != 2), ]

# change reference label
tabEvents$strain <- gsub(pattern = ref1, replacement = ref1Label, tabEvents$strain)
tabEvents$strain <- gsub(pattern = ref2, replacement = ref2Label, tabEvents$strain)

for (indS in allSamples) {
  # sample table sorting
  tabSamp <- tabEvents[which(tabEvents$hybrid == indS), ]
  tabSampSorted <- tabSamp[order(match(tabSamp$strain, c(ref1Label, ref2Label)), match(tabSamp$chr, allChr), 
                                 tabSamp$start), ]
  coordOut <- paste(tabSampSorted$chr, tabSampSorted$strain, sep = "_")
  tabOut <- data.frame(coordOut, tabSampSorted$start - 1, tabSampSorted$end, tabSampSorted$status)
  indHetero <- which(tabOut[, 4] == 1)
  indHomo <- which(tabOut[, 4] == 0)
  tabOutHetero <- tabOut[indHetero, c(1:3)]
  tabOutHomo <- tabOut[indHomo, c(1:3)]
  write.table(tabOutHetero, file = file.path(outDir, paste0("Hetero.", indS, ".bed")), sep = "\t", 
              row.names = F, col.names = F, quote = F)
  write.table(tabOutHomo, file = file.path(outDir, paste0("Homo.", indS, ".bed")), sep = "\t", 
              row.names = F, col.names = F, quote = F)
  # se non ci sono regioni in LOH (cosa molto possibile nel ParHyb)
  # faccio una riga fake per assicurarmi che venga generato un vcf
  if (nrow(tabOutHomo) == 0) {
    tabOutHomo <- data.frame(paste0(allChr[1], "_", ref1Label), c(0), c(1))
    write.table(tabOutHomo, file = file.path(outDir, paste0("Homo.", indS, ".bed")), sep = "\t", 
                row.names = F, col.names = F, quote = F)
  }
}


