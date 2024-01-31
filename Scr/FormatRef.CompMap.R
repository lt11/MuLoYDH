# header ------------------------------------------------------------------

# prepare Mod folder with references
# with modified contig names
# and concatenated for competitive mapping

rm(list = ls())
options(warn = 1)
options(stringsAsFactors = F)

# settings ----------------------------------------------------------------

# arguments
argsVal <- commandArgs(trailingOnly = T)
refName1 <- argsVal[1]
refName2 <- argsVal[2]
refDir <- argsVal[3]

# dirs & paths
outDir <- file.path(refDir, "Mod")
unlink(outDir, recursive = T)
dir.create(outDir, showWarnings = F)
refS <- c(refName1, refName2)

# run ---------------------------------------------------------------------

for (refName in refS) {
  inRef <- file.path(refDir, paste(refName, ".genome.fa", sep = ""))
  linesRef <- readLines(inRef)
  linesRef[grep(">", linesRef)] <- paste(linesRef[grep(">", linesRef)], refName, sep = "_")
  write.table(linesRef, file = file.path(refDir, "Mod", paste(refName, ".genome.chrref.fa", sep = "")), 
              quote = F, col.names = F, row.names = F)
  write.table(linesRef, file = file.path(refDir, "Mod", paste(paste(refName1, refName2, sep = "+"), 
                                                       ".genome.chrref.fa", sep = "")), 
              quote = F, col.names = F, row.names = F, append = T)
}



