# header ------------------------------------------------------------------

# legge i vcf nelle cartelle HybSubtracted di Ploidy1, Ploidy2 e Unphased
# fa il merge
# salva il vcf mergiato in VariantCalls/Merged/HybSubtracted
# fa il parsing di questi e salva in file.path(baseDir, "ParsedVariants", "Results", "HybSubtracted")

rm(list = ls())
options(stringsAsFactors = F)
library(vcfR)
library(ape)

# settings ----------------------------------------------------------------

# arguments
argsVal <- commandArgs(trailingOnly = T)
baseDir <- argsVal[1]
# commentami
# baseDir <- "/Users/lorenzo/Prog/Yeasts/DevMuLo4"

allChr <- c("chrI", "chrII", "chrIII", "chrIV", "chrV", "chrVI", "chrVII", "chrVIII", 
            "chrIX", "chrX", "chrXI", "chrXII", "chrXIII", "chrXIV", "chrXV", "chrXVI")

vcfClasses <- c("character", "numeric", "character", "character", "character", 
                "numeric", "character", "character", "character", "character")

# merger ------------------------------------------------------------------

# output folder
outDir <- file.path(baseDir, "VariantCalls", "Merged", "HybSubtracted")
dir.create(outDir, showWarnings = F, recursive = T)
# trova gli ID dei campioni
ploidy1Dir <- file.path(baseDir, "VariantCalls", "Ploidy1", "HybSubtracted")
ploidy1Files <- list.files(ploidy1Dir)
ploidy2Dir <- file.path(baseDir, "VariantCalls", "Ploidy2", "HybSubtracted")
ploidy2Files <- list.files(ploidy1Dir)
unphasedDir <- file.path(baseDir, "VariantCalls", "Unphased", "HybSubtracted")
unphasedFiles <- list.files(unphasedDir)

sampleIDs <- c(sapply(strsplit(ploidy1Files, split = "\\."), "[[", 1), 
               sapply(strsplit(ploidy2Files, split = "\\."), "[[", 1), 
               sapply(strsplit(unphasedFiles, split = "\\."), "[[", 1))
sampleIDs <- unique(sampleIDs)

stringTag <- "##INFO=<ID=PH,Number=1,Type=Character,Description=\"Boolean specifiing whether the variant is phased (T) or not phased (F)\">"

mapID <- grep(pattern = "+", 
              unlist(strsplit( ploidy1Files[1], split = "\\.")), 
              fixed = T, value = T)
outSuffix <- paste0(".", mapID, ".vcf")

ref1 <- unlist(strsplit(mapID, split = "\\+"))[1]
ref2 <- unlist(strsplit(mapID, split = "\\+"))[2]

allChrSort <- c(paste(allChr, ref1, sep = "_"), 
                paste(allChr, ref2, sep = "_"))

for (indS in sampleIDs) {
  outFile <- file.path(outDir, paste0(indS, outSuffix))
  ploidy1File <- list.files(path = ploidy1Dir, 
                            pattern = paste0(indS, "\\.", ref1, "\\+", ref2, "\\.vcf\\.gz"), 
                            full.names = T)
  ploidy2File <- list.files(path = ploidy2Dir, 
                            pattern = paste0(indS, "\\.", ref1, "\\+", ref2, "\\.vcf\\.gz"), 
                            full.names = T)
  unphasedRef1File <- list.files(path = unphasedDir, 
                                 pattern = paste0(indS, "\\.", ref1, "\\.vcf\\.gz"), 
                                 full.names = T)
  unphasedRef2File <- list.files(path = unphasedDir, 
                                 pattern = paste0(indS, "\\.", ref2, "\\.vcf\\.gz"), 
                                 full.names = T)
  # aggiungo all'header la descrizione del TAG PH
  allVcf <- c(ploidy1File, ploidy2File, unphasedRef1File, unphasedRef2File)
  lsVcf <- list()
  for (indH in 1:length(allVcf)) {
    lsVcf[[indH]] <- grep(pattern = "^#", readLines(con = allVcf[[indH]]), value = T)
  }
  # prendo l'header più lungo
  longHeader <- which.max(lapply(lsVcf, function(x) length(x)))
  # cat header
  inputHeader <- lsVcf[[longHeader]]
  outputHeader <- inputHeader[1:length(inputHeader) - 1]
  outputHeader <- c(outputHeader, stringTag, inputHeader[length(inputHeader)])
  cat(outputHeader, sep = "\n", file = outFile)
  
  wholeVcf <- c()
  if (class(try(read.table(ploidy1File), silent = T)) != 'try-error') {
    ploidy1Vcf <- read.table(ploidy1File, colClasses = vcfClasses)
    # aggiungo il TAG PH
    ploidy1Vcf[, 8] <- paste0(ploidy1Vcf[, 8], ";PH=T")
    wholeVcf <- rbind(wholeVcf, ploidy1Vcf)
  }
  if (class(try(read.table(ploidy2File), silent = T)) != 'try-error') {
    ploidy2Vcf <- read.table(ploidy2File, colClasses = vcfClasses)
    # aggiungo il TAG PH
    ploidy2Vcf[, 8] <- paste0(ploidy2Vcf[, 8], ";PH=F")
    wholeVcf <- rbind(wholeVcf, ploidy2Vcf)
    wholeVcf <- wholeVcf[order(match(wholeVcf$V1, allChrSort), wholeVcf$V2), ]
  }
  if (class(try(read.table(unphasedRef1File), silent = T)) != 'try-error') {
    unphasedRef1Vcf <- read.table(unphasedRef1File, colClasses = vcfClasses)
    # aggiungo il TAG PH
    unphasedRef1Vcf[, 8] <- paste0(unphasedRef1Vcf[, 8], ";PH=F")
    wholeVcf <- rbind(wholeVcf, unphasedRef1Vcf)
    wholeVcf <- wholeVcf[order(match(wholeVcf$V1, allChrSort), wholeVcf$V2), ]
  }
  if (class(try(read.table(unphasedRef2File), silent = T)) != 'try-error') {
    unphasedRef2Vcf <- read.table(unphasedRef2File, colClasses = vcfClasses)
    # aggiungo il TAG PH
    unphasedRef2Vcf[, 8] <- paste0(unphasedRef2Vcf[, 8], ";PH=F")
    wholeVcf <- rbind(wholeVcf, unphasedRef2Vcf)
    wholeVcf <- wholeVcf[order(match(wholeVcf$V1, allChrSort), wholeVcf$V2), ]
  }
  
  write.table(wholeVcf, file = outFile, append = T, row.names = F, col.names = F, quote = F, sep = "\t")
}

# parser ------------------------------------------------------------------

# dirs & paths
outDir <- file.path(baseDir, "ParsedVariants")
dir.create(outDir, showWarnings = F)
vcfMerged <- file.path(baseDir, "VariantCalls", "Merged", "HybSubtracted")
allVarDirs <- c(vcfMerged)

# run ---------------------------------------------------------------------

for (vcfDir in allVarDirs) {
  listFiles <- list.files(vcfDir, pattern = "\\.vcf$")
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


