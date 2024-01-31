# header ------------------------------------------------------------------

rm(list = ls())
options(stringsAsFactors = F)

# settings ----------------------------------------------------------------

argsVal <- commandArgs(trailingOnly = T)
baseDir <- argsVal[1]

dirDepthDdp <- file.path(baseDir, "DepthDdp")
dirDepthRaw <- file.path(baseDir, "DepthRaw")

allDir <- c(dirDepthDdp, dirDepthRaw)

for (indD in allDir) { 
  filesDepth <- list.files(indD, pattern = "\\.mean\\.txt$", full.names = T)
  expName <- c()
  mappingName <- c()
  finalData <- c()
  
  for (indF in filesDepth) {
    # increase row with exp data
    depthData <- read.table(indF, sep = "=", header = FALSE)[, 2]
    finalData <- rbind(finalData, depthData)
    # mapping name row 
    fileUnlist <- unlist(strsplit(basename(indF), "\\."))
    expName <- c(expName, fileUnlist[1])
    mappingName <- c(mappingName, fileUnlist[2])
  } # files loop
  dfOut <- data.frame(expName, mappingName, finalData, row.names = NULL)
  colnames(dfOut) <- c("Experiment", "Mapping", "Covered Reference Sites [bp]", 
                       "% Reference Covered", "Reference Length [bp]", "Average Depth", "Depth SD")
  outName <- paste("Table", basename(indD),"txt", sep = ".")
  pathOut <- file.path(baseDir, basename(indD), outName)
  write.table(dfOut, pathOut, sep ="\t", row.names = F, col.names = T, quote = F)
} # folders loop



