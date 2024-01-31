# header ------------------------------------------------------------------

# parser for markers
# this is also the guy who removes MT variants

rm(list = ls())
library(vcfR)
library(ape)

# settings ----------------------------------------------------------------

# arguments
argsVal <- commandArgs(trailingOnly = T)
workDir <- argsVal[1]

# dirs & paths
outDir <- file.path(workDir, "ParsedVariants")
dir.create(outDir, showWarnings = F)
vcfMarkers <- file.path(workDir, "Markers")
allVarDirs <- c(vcfMarkers)

# run ---------------------------------------------------------------------

for (vcfDir in allVarDirs) {
  listFiles <- list.files(vcfDir, pattern = "\\.vcf\\.gz$")
  outDirRes <- file.path(outDir, "Results", basename(vcfDir))
  
  for (ind2 in listFiles) {
    myVcf <- read.vcfR(file.path(vcfDir, ind2), 
                       limit = 1e+07, nrows = -1, skip = 0, cols = NULL,
                       convertNA = T, verbose = T)
    # retrieve reference path
    refPath <- unlist(strsplit(grep("reference", myVcf@meta, value = T), 
                               split = "##reference=file://", fixed = T))[2]
    # convert path to original reference
    # l'unica differenza tra quelli in Ref e quelli in Mod è il nome dei cromosomi
    refPath <- gsub(".chrref", "", x = gsub("Mod/", "", x = refPath, fixed = T), fixed = T)
    refName <- unlist(strsplit(basename(refPath), split = "\\."))[[1]]
    # filter out MT data
    indOut <- grep("chrMT_", myVcf@fix[, 1])
    if (length(indOut) > 0) {
      myVcf@fix <- myVcf@fix[-indOut, ]
      myVcf@gt <- myVcf@gt[-indOut, ]
    }
    # remove _STRAIN from marker vcf
    myVcf@fix[, 1] <- gsub(pattern = "_.*$", "", myVcf@fix[, 1])
    loopChr <- unique(myVcf@fix[, 1])
    # add sample folder
    splitIn <- unlist(strsplit(ind2, split = "\\."))
    sampleID <- splitIn[1]
    outDirResSamp <- file.path(outDirRes, sampleID)
    dir.create(outDirResSamp, recursive = T, showWarnings = F)
    # save RData
    outRData <- file.path(outDirResSamp, paste(sampleID, refName, "vcf", "RData", sep = "."))
    vcfData <- list(myVcf@meta, myVcf@fix, myVcf@gt)
    names(vcfData) <- c("meta", "fix", "format")
    save(list = c("vcfData"), file = outRData)
  }
}



