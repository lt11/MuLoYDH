# header ------------------------------------------------------------------

# WARN: crossref!

# saves "quality" & "region" vcf
# using START-END marker coordinates
# annotates LOHs & non-LOHs with SNV & indels (quality-filtered, sdFactor*sdPol) from competitive mappings
# variants within telomers are filtered
# since V5 it does not consider CnC & CbC positions for filtering

rm(list = ls())
options(stringsAsFactors = F)
# disable scientific notation
options(scipen=999)

# function(s) -------------------------------------------------------------

SpecDec <- function(x, k) as.numeric(format(round(x, k), nsmall = k))

ptmGlob <- proc.time()

# settings ----------------------------------------------------------------

# arguments
argsVal <- commandArgs(trailingOnly = T)
ref1Label <- argsVal[1]
ref2Label <- argsVal[2]
ref1 <- argsVal[3]
ref2 <- argsVal[4]
baseDir <- argsVal[5]
# commentami
# ref1Label <- "SK1"
# ref2Label <- "S288C"
# ref1 <- "SK1"
# ref2 <- "BY"
# baseDir <- "/Users/lorenzo/Prog/Yeasts/DevMuLo4"

# init
sdFactor <- 1

# paths
varCmpDir <- file.path(baseDir, "ParsedVariants", "Results", "HybSubtracted")
outFile <- "Annotations.SNV.txt"
eventsDir <- file.path(baseDir, "LOH")
mummerDir <- file.path(baseDir, "MUMmer")
outDirSmallVarAnn <- file.path(baseDir, "LOH", "Annotations")
dir.create(path = outDirSmallVarAnn, showWarnings = F, recursive = T)
outDirQual <- file.path(baseDir, "VariantCalls", "Merged", "Quality")
dir.create(path = outDirQual, showWarnings = F, recursive = T)
outDirReg <- file.path(baseDir, "VariantCalls", "Merged", "Regions")
dir.create(path = outDirReg, showWarnings = F, recursive = T)


headerSummary <- data.frame("sample", "ref", "chr", "start", "first", "last", "end", "status", 
                            "len", "denES", "evrES", "distES", "denLF", "evrLF", "distLF", "SNV", 
                            "denSNV", "InDel", "denInDels", "hetSNV", "hetInDel")
# init outFile
write.table(x = headerSummary, file = file.path(outDirSmallVarAnn, outFile), 
            col.names = F, row.names = F, quote = F, append = F, sep = "\t")

allChr <- c("chrI", "chrII", "chrIII", "chrIV", "chrV", "chrVI", "chrVII", "chrVIII", "chrIX", 
            "chrX", "chrXI", "chrXII", "chrXIII", "chrXIV", "chrXV", "chrXVI")

subTelLRefA <- as.numeric(read.table(file = file.path(baseDir, "Ref", "Ann", paste0(ref1Label, ".subtel.txt")))[, 1])
subTelRRefA <- as.numeric(read.table(file = file.path(baseDir, "Ref", "Ann", paste0(ref1Label, ".subtel.txt")))[, 2])
subTelLRefB <- as.numeric(read.table(file = file.path(baseDir, "Ref", "Ann", paste0(ref2Label, ".subtel.txt")))[, 1])
subTelRRefB <- as.numeric(read.table(file = file.path(baseDir, "Ref", "Ann", paste0(ref2Label, ".subtel.txt")))[, 2])

corePosFile <- paste("Markers", paste(ref1, ref2, sep = "-"),"RData", sep = ".")
load(file.path(eventsDir, "Markers", corePosFile))
sampIDs <- list.files(varCmpDir)

formatNames <- c("GT", "PL", "DP", "SP", "ADF", "ADR", "AD")
fixNames <- c("CHR", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")

# sample loop -------------------------------------------------------------
for (indS in sampIDs) {
  sampFileCMP <- list.files(file.path(varCmpDir, indS), full.names = T)
  load(sampFileCMP) # load CMP vcfData list
  vcfDf <- list(vcfData$meta, 
                data.frame(vcfData$fix[, 1], as.numeric(vcfData$fix[, 2]), 
                           vcfData$fix[, 3:5], as.numeric(vcfData$fix[, 6]), vcfData$fix[, 7:8]), 
                data.frame(sapply(strsplit(vcfData$format[, 2], split = ":"), "[[", 1), 
                           sapply(strsplit(vcfData$format[, 2], split = ":"), "[[", 2), 
                           sapply(strsplit(vcfData$format[, 2], split = ":"), "[[", 3), 
                           sapply(strsplit(vcfData$format[, 2], split = ":"), "[[", 4), 
                           sapply(strsplit(vcfData$format[, 2], split = ":"), "[[", 5), 
                           sapply(strsplit(vcfData$format[, 2], split = ":"), "[[", 6), 
                           sapply(strsplit(vcfData$format[, 2], split = ":"), "[[", 7)))
  
  names(vcfDf) <- c("meta", "fix", "format")
  colnames(vcfDf$fix) <- fixNames
  colnames(vcfDf$format) <- formatNames
  
  # RefA filters --------------------------------------------------------------
  
  # filter for reference
  vcfRefA <- vcfDf
  indRef <- grep(pattern = ref1Label, x = vcfRefA$fix$CHR)
  vcfRefA$fix <- vcfDf$fix[indRef, ]
  vcfRefA$fix$CHR <- gsub(pattern = "_.*$", x = vcfRefA$fix$CHR, replacement = "")
  vcfRefA$format <- vcfRefA$format[indRef, ]
  
  # carica le qualità dei marker RefA
  qPolRefA <- read.table(file = file.path(eventsDir, indS, paste(indS, ref1, "QUALdata.txt", sep = ".")), 
                       na.strings = " ")
  qualPolRefA <- qPolRefA$V4[1]
  sdPolRefA <- qPolRefA$V4[2]
  qThresh <- qualPolRefA - sdFactor*sdPolRefA
  
  # quality filters 
  indBonSnv <- which(vcfRefA$fix$QUAL > qThresh)
  vcfRefAQual <- list(vcfRefA$meta, vcfRefA$fix[indBonSnv, ], vcfRefA$format[indBonSnv, ])
  names(vcfRefAQual) <- c("meta", "fix", "format")
  rm(list = "vcfRefA")
  
  # save vcf quality filtered
  vcfFileOutRefA <- file.path(outDirQual, paste(indS, ref1, "vcf", sep = "."))
  write.table(x = vcfRefAQual$meta, file = vcfFileOutRefA, col.names = F, row.names = F, quote = F)
  # source bam file
  chrSamCoommand <- grep(pattern = "samtoolsCommand", x = vcfRefAQual$meta, value = T)
  ulsSamCommand <- unlist(strsplit(chrSamCoommand, split = " "))
  sourceBam <- ulsSamCommand[length(ulsSamCommand)]
  # header string
  headerString <- c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", sourceBam, "\n")
  cat(x = headerString, file = vcfFileOutRefA, append = T, sep = "\t")
  # merge dataframes
  nRowVcf <- nrow(vcfRefAQual$format)
  write.table(x = data.frame(vcfRefAQual$fix, rep(paste(formatNames, collapse = ":"), nRowVcf), 
                             apply(vcfRefAQual$format[, c(1:ncol(vcfRefAQual$format))], 1, paste, collapse = ":" )), 
              file = vcfFileOutRefA, col.names = F, row.names = F, quote = F, 
              append = T, sep = "\t", na = ".")
  
  # write header of regions filtered vcf file
  vcfFileOutRegRefA <- file.path(outDirReg, paste(indS, ref1, "vcf", sep = "."))
  write.table(x = vcfRefAQual$meta, file = vcfFileOutRegRefA, col.names = F, row.names = F, quote = F)
  # header string
  cat(x = headerString, file = vcfFileOutRegRefA, append = T, sep = "\t")
  
  # separate small indels from SNV
  indSid <- grep(pattern = "INDEL", x = (vcfRefAQual$fix$INFO))
  vcfSidQual <- list(vcfRefAQual$meta, vcfRefAQual$fix[indSid, ], vcfRefAQual$format[indSid, ])
  names(vcfSidQual) <- c("meta", "fix", "format")
  indSnv <- grep(pattern = "INDEL", x = (vcfRefAQual$fix$INFO), invert = T)
  vcfSnvQual <- list(vcfRefAQual$meta, vcfRefAQual$fix[indSnv, ], vcfRefAQual$format[indSnv, ])
  names(vcfSnvQual) <- c("meta", "fix", "format")
  
  # filter out markers
  # crossref: prendi CHR_POS_ALT di ref1 nelle varianti e nella mappa dei marker
  stringVarRefA <- paste(vcfSnvQual$fix$CHR, vcfSnvQual$fix$POS, vcfSnvQual$fix$ALT, sep = "")
  stringPolRefA <- paste(dfGG$chrRefA, dfGG$posRefA, dfGG$refRefB, sep = "")
  intVarPol <- intersect(stringVarRefA, stringPolRefA)
  booBonRefA <- !is.element(stringVarRefA, intVarPol)
  vcfSnvQualFilt <- list(vcfSnvQual$meta, 
                         vcfSnvQual$fix[booBonRefA, ], 
                         vcfSnvQual$format[booBonRefA, ])
  names(vcfSnvQualFilt) <- c("meta", "fix", "format")
  rm(vcfSnvQual)
  # dummy line to keep variable name concordant with SNVs
  vcfSidQualFilt <- vcfSidQual
  rm(vcfSidQual)
  
  # # calculate indel / SNV ratio
  # nSid <- length(indSid)
  # ratIndelSnv <- nSid / (length(vcfRefAQual$fix$INFO) - nSid)
  # ratIndelSnv <- SpecDec(ratIndelSnv, 3)
  # cat(indS, "(RefA) indel / SNV ratio: ", ratIndelSnv, "\n")
  
  # RefA data (ref1Label) -----------------------------------------------------
  counterChr <- 1
  for (indC in allChr) {
    # accumulatore varianti by cromosoma
    accuVarRegChrRefA <- c()
    
    samplePath <- file.path(eventsDir, indS)
    fileSeg <- list.files(samplePath, 
                          pattern = paste(indS, indC, ref1, "Seg.RData", sep = "."), 
                          full.names = T)
    if (length(fileSeg) == 0) {
      counterChr <- counterChr + 1
      next()
    }
    # variabile segmenti: evRefA
    load(fileSeg)
    assign("evRefA", get(paste0("ev", ref1)))
    
    # se non ci sono indel o snv si pianta
    # NB se aggiungi un next, RICORDATI di incrementare counterChr
    indChrSnv <- which(vcfSnvQualFilt$fix$CHR == indC)
    indChrSid <- which(vcfSidQualFilt$fix$CHR == indC)
    vcfSnvQualChr <- list(vcfSnvQualFilt$fix[indChrSnv, ], vcfSnvQualFilt$format[indChrSnv, ])
    names(vcfSnvQualChr) <- c("fix", "format")
    vcfSidQualChr <- list(vcfSidQualFilt$fix[indChrSid, ], vcfSidQualFilt$format[indChrSid, ])
    names(vcfSidQualChr) <- c("fix", "format")
    
    # filter (sub)telomeric SNVs
    indTelVar <- which(vcfSnvQualChr$fix$POS < subTelLRefA[counterChr] | 
                         vcfSnvQualChr$fix$POS > subTelRRefA[counterChr])
    if (length(indTelVar) != 0) {
      vcfSnvQualChr$fix <- vcfSnvQualChr$fix[-indTelVar, ]
      vcfSnvQualChr$format <- vcfSnvQualChr$format[-indTelVar, ]
    }
    # filter (sub)telomeric indels
    indTelSid <- which(vcfSidQualChr$fix$POS < subTelLRefA[counterChr] | 
                         vcfSidQualChr$fix$POS > subTelRRefA[counterChr])
    if (length(indTelSid) != 0) {
      vcfSidQualChr$fix <- vcfSidQualChr$fix[-indTelSid, ]
      vcfSidQualChr$format <- vcfSidQualChr$format[-indTelSid, ]
    }
    
    # LOH events --------------------------------------------------------------
    
    indLOH <- which(evRefA$status == 0)
    evLOH <- evRefA[indLOH, ]
    numLOH <- length(indLOH)
    if (numLOH != 0) {
      nR <- nrow(evLOH)
      indSnvBon <- replicate(nR, list())
      indSidBon <- replicate(nR, list())
      lohIdentifierSnv <- replicate(nR, list())
      lohIdentifierSid <- lohIdentifierSnv
      # quali SNV & indel cadono nei segmenti tra FIRST e LAST marker
      for (indSeg in 1:nR) {
        indSnvBon[[indSeg]] <- which(vcfSnvQualChr$fix$POS > evLOH$start[indSeg] & 
                                       vcfSnvQualChr$fix$POS < evLOH$end[indSeg])
        indSidBon[[indSeg]] <- which(vcfSidQualChr$fix$POS > evLOH$start[indSeg] & 
                                       vcfSidQualChr$fix$POS < evLOH$end[indSeg])
        lohIdentifierSnv[[indSeg]] <- rep(indSeg, length(indSnvBon[[indSeg]]))
        lohIdentifierSid[[indSeg]] <- rep(indSeg, length(indSidBon[[indSeg]]))
      }
      lossHetIdenSnv <- unlist(lohIdentifierSnv)
      lossHetIdenSid <- unlist(lohIdentifierSid)
      snVar <- unlist(lapply(indSnvBon, length))
      idVar <- unlist(lapply(indSidBon, length))
      # annota tabella LOH con densità varianti
      evLOH <- data.frame(evLOH, snVar, SpecDec(c(snVar / evLOH$distLF), 7))
      colnames(evLOH)[14:15] <- c("SNV", "denSNV")
      evLOH <- data.frame(evLOH, idVar, SpecDec(c(idVar / evLOH$distLF), 7))
      colnames(evLOH)[16:17] <- c("InDel", "denInDel")
      
      # # numero regioni LOH
      # numLOH <- nrow(evLOH)
      # cat(c("STAT SNV LOH RefA", indS, indC, "# events:", numLOH, "\n"), sep = " ")
      # # Total lenght (LF)
      # sumLen <- sum(evLOH$distLF)
      # cat(c("STAT SNV LOH RefA", indS, indC, "Tot. length FL [bp]:", 
      #       sumLen, "\n"), sep = " ")
      # # Total lenght (ES)
      # sumLenES <- sum(evLOH$distES)
      # cat(c("STAT SNV LOH RefA", indS, indC, "Tot. length SE [bp]:", 
      #       sumLenES, "\n"), sep = " ")
      # # quanti marker
      # sumPol <- sum(evLOH$len)
      # cat(c("STAT SNV LOH RefA", indS, indC, "# marker:", 
      #       sumPol, "\n"), sep = " ")
      
      # annota il vcf delle SNV con l'ID del segmento
      indSnvChr <- unlist(indSnvBon)
      # make sample output folder
      sampOutDir <- file.path(outDirSmallVarAnn, indS)
      dir.create(path = sampOutDir, showWarnings = F, recursive = T)
      # annota quante SNV het nel file degli eventi
      annHetSNV <- rep(0, numLOH)
      
      if (length(indSnvChr) != 0) {
        dfFormat <- data.frame(vcfSnvQualChr$format[indSnvChr, ], lossHetIdenSnv)
        colnames(dfFormat)[8] <- c("lohID")
        chrLoh <- list(vcfSidQualFilt$meta, 
                       vcfSnvQualChr$fix[indSnvChr, ], 
                       dfFormat)
        names(chrLoh) <- c("meta", "fix", "format")
        # save vcf RData SNV
        varName <- paste("vcfSNV", indS, indC, ref1, sep = "")
        assign(varName, chrLoh)
        save(list = varName, file = file.path(sampOutDir, 
                                              paste(indS, indC, ref1, "Vcf.SNV.LOH.RData", sep = ".")))
        # accumulate variants for regions folder
        regVar <- data.frame(vcfSnvQualChr$fix[indSnvChr, ], rep(paste(formatNames, collapse = ":"), nrow(dfFormat)), 
                             apply(dfFormat[, c(1:ncol(dfFormat)-1)], 1, paste, collapse = ":" ))
        accuVarRegChrRefA <- rbind(accuVarRegChrRefA, regVar)
        
        # numero LOH con SNV
        withSNV <- which(evLOH$SNV != 0)
        lenEvSNV <- length(withSNV)
        cat(c("STAT SNV LOH", ref1, indS, indC, "# SNV events:", lenEvSNV, "\n"), sep = " ")
        # status of events with SNV
        statusSNV <- evLOH$status[withSNV]
        cat(c("STAT SNV LOH", ref1, indS, indC, "status of events with SNV:", 
              paste(statusSNV, collapse = ","), "\n"), sep = " ")
        # quante SNV hom + SNV het
        sumSNV <- sum(evLOH$SNV)
        cat(c("STAT SNV LOH", ref1, indS, indC, "# SNV var:", sumSNV, "\n"), sep = " ")
        # quante SNV het
        numHetSNV <- length(which(chrLoh$format$GT == "0/1"))
        cat(c("STAT SNV LOH", ref1, indS, indC, "# SNV het:", numHetSNV, "\n"), sep = " ")
        
        # annota quante SNV het nel file degli eventi
        wSNV <- unique(chrLoh$format$lohID)
        for (indU in wSNV) {
          annHetSNV[indU] <- length(which(chrLoh$format$lohID == indU & 
                                            chrLoh$format$GT == "0/1"))
        }
        evLOH <- data.frame(evLOH, annHetSNV)
        colnames(evLOH)[which(colnames(evLOH) == "annHetSNV")] <- "hetSNV"
        rm(chrLoh)
        rm(list = varName)
      } else {
        evLOH <- data.frame(evLOH, annHetSNV)
        colnames(evLOH)[which(colnames(evLOH) == "annHetSNV")] <- "hetSNV"
      }
      
      # annota il vcf delle Sid con l'ID del segmento
      indSidChr <- unlist(indSidBon)
      # annota quante Sid het nel file degli eventi
      annHetSid <- rep(0, numLOH)
      
      if (length(indSidChr) != 0) {
        dfFormat <- data.frame(vcfSidQualChr$format[indSidChr, ], lossHetIdenSid)
        colnames(dfFormat)[8] <- c("lohID")
        chrLoh <- list(vcfSidQualFilt$meta, 
                       vcfSidQualChr$fix[indSidChr, ], 
                       dfFormat)
        names(chrLoh) <- c("meta", "fix", "format")
        # save vcf InDel
        varName <- paste("vcfInDel", indS, indC, ref1, sep = "")
        assign(varName, chrLoh)
        save(list = varName, file = file.path(sampOutDir, 
                                              paste(indS, indC, ref1, "Vcf.InDel.LOH.RData", sep = ".")))
        # accumulate variants for regions folder
        regVar <- data.frame(vcfSidQualChr$fix[indSidChr, ], rep(paste(formatNames, collapse = ":"), nrow(dfFormat)), 
                             apply(dfFormat[, c(1:ncol(dfFormat)-1)], 1, paste, collapse = ":" ))
        accuVarRegChrRefA <- rbind(accuVarRegChrRefA, regVar)
        
        # numero LOH con InDel
        withSid <- which(evLOH$InDel != 0)
        lenEvSid <- length(withSid)
        cat(c("STAT InDel LOH", ref1, indS, indC, "# InDel events:", lenEvSid, "\n"), sep = " ")
        # status of events with InDel
        statusSid <- evLOH$status[withSid]
        cat(c("STAT InDel LOH", ref1, indS, indC, "status of events with InDel:", 
              paste(statusSid, collapse = ","), "\n"), sep = " ")
        # quante InDel hom + InDel het
        sumSid <- sum(evLOH$InDel)
        cat(c("STAT InDel LOH", ref1, indS, indC, "# InDel var:", sumSid, "\n"), sep = " ")
        # quante InDel het
        numHetSid <- length(which(chrLoh$format$GT == "0/1"))
        cat(c("STAT InDel LOH", ref1, indS, indC, "# InDel het:", numHetSid, "\n"), sep = " ")
        
        # annota quante InDel het nel file degli eventi
        wSid <- unique(chrLoh$format$lohID)
        for (indU in wSid) {
          annHetSid[indU] <- length(which(chrLoh$format$lohID == indU & 
                                            chrLoh$format$GT == "0/1"))
        }
        evLOH <- data.frame(evLOH, annHetSid)
        colnames(evLOH)[which(colnames(evLOH) == "annHetSid")] <- "hetInDels"
        rm(chrLoh)
        rm(list = varName)
      } else {
        evLOH <- data.frame(evLOH, annHetSid)
        colnames(evLOH)[which(colnames(evLOH) == "annHetSid")] <- "hetInDels"
      }
      
      # save annotated data: evLOH
      varName <- paste(indS, indC, ref1, "evLOH", sep = "")
      # assign annotated LOH table to e.g. A315R52chrVIevLOH
      assign(varName, evLOH)
      save(list = varName, file = file.path(sampOutDir, paste(indS, indC, ref1, "LOH.RData", sep = ".")))
      # appending data to summary table by events
      write.table(x = data.frame(rep(indS, nR), rep(ref1, nR), get(varName)), 
                  file = file.path(outDirSmallVarAnn, outFile), 
                  col.names = F, row.names = F, quote = F, append = T, sep = "\t")
      
      rm(evLOH)
      rm(list = varName)
    }
    
    # HET events --------------------------------------------------------------
    
    indHET <- which(evRefA$status == 1)
    evHET <- evRefA[indHET, ]
    numHET <- length(indHET)
    if (numHET != 0) {
      nR <- nrow(evHET)
      indSnvBon <- replicate(nR, list())
      indSidBon <- replicate(nR, list())
      lohIdentifierSnv <- replicate(nR, list())
      lohIdentifierSid <- lohIdentifierSnv
      # quali SNV & indel cadono nei segmenti tra FIRST e LAST marker
      for (indSeg in 1:nR) {
        indSnvBon[[indSeg]] <- which(vcfSnvQualChr$fix$POS > evHET$start[indSeg] & 
                                       vcfSnvQualChr$fix$POS < evHET$end[indSeg])
        indSidBon[[indSeg]] <- which(vcfSidQualChr$fix$POS > evHET$start[indSeg] & 
                                       vcfSidQualChr$fix$POS < evHET$end[indSeg])
        lohIdentifierSnv[[indSeg]] <- rep(indSeg, length(indSnvBon[[indSeg]]))
        lohIdentifierSid[[indSeg]] <- rep(indSeg, length(indSidBon[[indSeg]]))
      }
      lossHetIdenSnv <- unlist(lohIdentifierSnv)
      lossHetIdenSid <- unlist(lohIdentifierSid)
      snVar <- unlist(lapply(indSnvBon, length))
      idVar <- unlist(lapply(indSidBon, length))
      # annota tabella HET con densità varianti
      evHET <- data.frame(evHET, snVar, SpecDec(c(snVar / evHET$distLF), 7))
      colnames(evHET)[14:15] <- c("SNV", "denSNV")
      evHET <- data.frame(evHET, idVar, SpecDec(c(idVar / evHET$distLF), 7))
      colnames(evHET)[16:17] <- c("InDel", "denInDel")
      
      # # numero regioni HET
      # numHET <- nrow(evHET)
      # cat(c("STAT SNV HET RefA", indS, indC, "# events:", numHET, "\n"), sep = " ")
      # # Total lenght (LF)
      # sumLen <- sum(evHET$distLF)
      # cat(c("STAT SNV HET RefA", indS, indC, "Tot. length FL [bp]:", 
      #       sumLen, "\n"), sep = " ")
      # # Total lenght (ES)
      # sumLenES <- sum(evHET$distES)
      # cat(c("STAT SNV HET RefA", indS, indC, "Tot. length SE [bp]:", 
      #       sumLenES, "\n"), sep = " ")
      # # quanti marker
      # sumPol <- sum(evHET$len)
      # cat(c("STAT SNV HET RefA", indS, indC, "# marker:", 
      #       sumPol, "\n"), sep = " ")
      
      # annota il vcf delle SNV con l'ID del segmento
      indSnvChr <- unlist(indSnvBon)
      # make sample output folder
      sampOutDir <- file.path(outDirSmallVarAnn, indS)
      dir.create(path = sampOutDir, showWarnings = F, recursive = T)
      # annota quante SNV het nel file degli eventi
      annHetSNV <- rep(0, numHET)
      
      if (length(indSnvChr) != 0) {
        dfFormat <- data.frame(vcfSnvQualChr$format[indSnvChr, ], lossHetIdenSnv)
        colnames(dfFormat)[8] <- c("lohID")
        chrLoh <- list(vcfSidQualFilt$meta, 
                       vcfSnvQualChr$fix[indSnvChr, ], 
                       dfFormat)
        names(chrLoh) <- c("meta", "fix", "format")
        # save vcf RData SNV
        varName <- paste("vcfSNV", indS, indC, ref1, sep = "")
        assign(varName, chrLoh)
        save(list = varName, file = file.path(sampOutDir, 
                                              paste(indS, indC, ref1, "Vcf.SNV.HET.RData", sep = ".")))
        # accumulate variants for regions folder
        regVar <- data.frame(vcfSnvQualChr$fix[indSnvChr, ], rep(paste(formatNames, collapse = ":"), nrow(dfFormat)), 
                             apply(dfFormat[, c(1:ncol(dfFormat)-1)], 1, paste, collapse = ":" ))
        accuVarRegChrRefA <- rbind(accuVarRegChrRefA, regVar)
        
        # numero HET con SNV
        withSNV <- which(evHET$SNV != 0)
        lenEvSNV <- length(withSNV)
        cat(c("STAT SNV HET", ref1, indS, indC, "# SNV events:", lenEvSNV, "\n"), sep = " ")
        # status of events with SNV
        statusSNV <- evHET$status[withSNV]
        cat(c("STAT SNV HET", ref1, indS, indC, "status of events with SNV:", 
              paste(statusSNV, collapse = ","), "\n"), sep = " ")
        # quante SNV hom + SNV het
        sumSNV <- sum(evHET$SNV)
        cat(c("STAT SNV HET", ref1, indS, indC, "# SNV var:", sumSNV, "\n"), sep = " ")
        # quante SNV het
        numHetSNV <- length(which(chrLoh$format$GT == "0/1"))
        cat(c("STAT SNV HET", ref1, indS, indC, "# SNV het:", numHetSNV, "\n"), sep = " ")
        
        # annota quante SNV het nel file degli eventi
        wSNV <- unique(chrLoh$format$lohID)
        for (indU in wSNV) {
          annHetSNV[indU] <- length(which(chrLoh$format$lohID == indU & 
                                            chrLoh$format$GT == "0/1"))
        }
        evHET <- data.frame(evHET, annHetSNV)
        colnames(evHET)[which(colnames(evHET) == "annHetSNV")] <- "hetSNV"
        rm(chrLoh)
        rm(list = varName)
      } else {
        evHET <- data.frame(evHET, annHetSNV)
        colnames(evHET)[which(colnames(evHET) == "annHetSNV")] <- "hetSNV"
      }
      
      # annota il vcf delle Sid con l'ID del segmento
      indSidChr <- unlist(indSidBon)
      # annota quante Sid het nel file degli eventi
      annHetSid <- rep(0, numHET)
      
      if (length(indSidChr) != 0) {
        dfFormat <- data.frame(vcfSidQualChr$format[indSidChr, ], lossHetIdenSid)
        colnames(dfFormat)[8] <- c("lohID")
        chrLoh <- list(vcfSidQualFilt$meta, 
                       vcfSidQualChr$fix[indSidChr, ], 
                       dfFormat)
        names(chrLoh) <- c("meta", "fix", "format")
        # save vcf InDel
        varName <- paste("vcfInDel", indS, indC, ref1, sep = "")
        assign(varName, chrLoh)
        save(list = varName, file = file.path(sampOutDir, 
                                              paste(indS, indC, ref1, "Vcf.InDel.HET.RData", sep = ".")))
        # accumulate variants for regions folder
        regVar <- data.frame(vcfSidQualChr$fix[indSidChr, ], rep(paste(formatNames, collapse = ":"), nrow(dfFormat)), 
                             apply(dfFormat[, c(1:ncol(dfFormat)-1)], 1, paste, collapse = ":" ))
        accuVarRegChrRefA <- rbind(accuVarRegChrRefA, regVar)
        
        # numero HET con InDel
        withSid <- which(evHET$InDel != 0)
        lenEvSid <- length(withSid)
        cat(c("STAT InDel HET", ref1, indS, indC, "# InDel events:", lenEvSid, "\n"), sep = " ")
        # status of events with InDel
        statusSid <- evHET$status[withSid]
        cat(c("STAT InDel HET", ref1, indS, indC, "status of events with InDel:", 
              paste(statusSid, collapse = ","), "\n"), sep = " ")
        # quante InDel hom + InDel het
        sumSid <- sum(evHET$InDel)
        cat(c("STAT InDel HET", ref1, indS, indC, "# InDel var:", sumSid, "\n"), sep = " ")
        # quante InDel het
        numHetSid <- length(which(chrLoh$format$GT == "0/1"))
        cat(c("STAT InDel HET", ref1, indS, indC, "# InDel het:", numHetSid, "\n"), sep = " ")
        
        # annota quante InDel het nel file degli eventi
        wSid <- unique(chrLoh$format$lohID)
        for (indU in wSid) {
          annHetSid[indU] <- length(which(chrLoh$format$lohID == indU & 
                                            chrLoh$format$GT == "0/1"))
        }
        evHET <- data.frame(evHET, annHetSid)
        colnames(evHET)[which(colnames(evHET) == "annHetSid")] <- "hetInDels"
        rm(chrLoh)
        rm(list = varName)
      } else {
        evHET <- data.frame(evHET, annHetSid)
        colnames(evHET)[which(colnames(evHET) == "annHetSid")] <- "hetInDels"
      }
      
      # save annotated data: evHET
      varName <- paste(indS, indC, ref1, "evHET", sep = "")
      # assign annotated HET table to e.g. A315R52chrVIevHET
      assign(varName, evHET)
      save(list = varName, file = file.path(sampOutDir, paste(indS, indC, ref1, "HET.RData", sep = ".")))
      # appending data to summary table by events
      write.table(x = data.frame(rep(indS, nR), rep(ref1, nR), get(varName)), 
                  file = file.path(outDirSmallVarAnn, outFile), 
                  col.names = F, row.names = F, quote = F, append = T, sep = "\t")
      
      rm(evHET)
      rm(list = varName)
    }
    if (length(accuVarRegChrRefA$POS) != 0) {
      # sort acummulatore varianti
      srtAccuVarRegChrRefA <- accuVarRegChrRefA[order(accuVarRegChrRefA$POS), ]
      # scrivi l'accumulatore sorted
      write.table(x = srtAccuVarRegChrRefA, file = vcfFileOutRegRefA, col.names = F, row.names = F, quote = F, 
                  append = T, sep = "\t", na = ".")
    }
    
    counterChr <- counterChr + 1
  } # chromosome loop RefA
  
  # RefB filters --------------------------------------------------------------
  
  # filter for reference
  vcfRefB <- vcfDf
  indRef <- grep(pattern = ref2Label, x = vcfRefB$fix$CHR)
  vcfRefB$fix <- vcfDf$fix[indRef, ]
  vcfRefB$fix$CHR <- gsub(pattern = "_.*$", x = vcfRefB$fix$CHR, replacement = "")
  vcfRefB$format <- vcfRefB$format[indRef, ]
  
  # carica le qualità dei marker RefB
  qPolRefB <- read.table(file = file.path(eventsDir, indS, paste(indS, ref2, "QUALdata.txt", sep = ".")),
                       na.strings = " ")
  qualPolRefB <- qPolRefB$V4[1]
  sdPolRefB <- qPolRefB$V4[2]
  qThresh <- qualPolRefB - sdFactor*sdPolRefB
  
  # quality filters
  indBonSnv <- which(vcfRefB$fix$QUAL > qThresh)
  vcfRefBQual <- list(vcfRefB$meta, vcfRefB$fix[indBonSnv, ], vcfRefB$format[indBonSnv, ])
  names(vcfRefBQual) <- c("meta", "fix", "format")
  rm(list = "vcfRefB")
  
  # save vcf quality filtered
  vcfFileOutRefB <- file.path(outDirQual, paste(indS, ref2, "vcf", sep = "."))
  write.table(x = vcfRefBQual$meta, file = vcfFileOutRefB, col.names = F, row.names = F, quote = F)
  # header string
  cat(x = headerString, file = vcfFileOutRefB, append = T, sep = "\t")
  # merge dataframes
  nRowVcf <- nrow(vcfRefBQual$format)
  write.table(x = data.frame(vcfRefBQual$fix, rep(paste(formatNames, collapse = ":"), nRowVcf), 
                             apply(vcfRefBQual$format[, c(1:ncol(vcfRefBQual$format))], 1, paste, collapse = ":" )),
              file = vcfFileOutRefB, col.names = F, row.names = F, quote = F,
              append = T, sep = "\t", na = ".")
  
  # write header of regions filtered vcf file
  vcfFileOutRegRefB <- file.path(outDirReg, paste(indS, ref2, "vcf", sep = "."))
  write.table(x = vcfRefBQual$meta, file = vcfFileOutRegRefB, col.names = F, row.names = F, quote = F)
  # header string
  cat(x = headerString, file = vcfFileOutRegRefB, append = T, sep = "\t")
  
  # separate small indels from SNV
  indSid <- grep(pattern = "INDEL", x = (vcfRefBQual$fix$INFO))
  vcfSidQual <- list(vcfRefBQual$meta, vcfRefBQual$fix[indSid, ], vcfRefBQual$format[indSid, ])
  names(vcfSidQual) <- c("meta", "fix", "format")
  indSnv <- grep(pattern = "INDEL", x = (vcfRefBQual$fix$INFO), invert = T)
  vcfSnvQual <- list(vcfRefBQual$meta, vcfRefBQual$fix[indSnv, ], vcfRefBQual$format[indSnv, ])
  names(vcfSnvQual) <- c("meta", "fix", "format")
  
  # filter out markers
  # crossref: prendi CHR_POS_ALT di ref2 nelle varianti e nella mappa dei marker
  stringVarRefB <- paste(vcfSnvQual$fix$CHR, vcfSnvQual$fix$POS, vcfSnvQual$fix$ALT, sep = "")
  stringPolRefB <- paste(dfGG$chrRefB, dfGG$posRefB, dfGG$refXX, sep = "")
  intVarPol <- intersect(stringVarRefB, stringPolRefB)
  booBonRefB <- !is.element(stringVarRefB, intVarPol)
  vcfSnvQualFilt <- list(vcfSnvQual$meta,
                         vcfSnvQual$fix[booBonRefB, ],
                         vcfSnvQual$format[booBonRefB, ])
  names(vcfSnvQualFilt) <- c("meta", "fix", "format")
  rm(vcfSnvQual)
  # dummy line to keep variable name concordant with SNVs
  vcfSidQualFilt <- vcfSidQual
  rm(vcfSidQual)
  
  # # calculate indel / SNV ratio
  # nSid <- length(indSid)
  # ratIndelSnv <- nSid / (length(vcfRefBQual$fix$INFO) - nSid)
  # ratIndelSnv <- SpecDec(ratIndelSnv, 3)
  # cat(indS, "(RefB) indel / SNV ratio: ", ratIndelSnv, "\n")
  
  # RefB data (ref2Label) -----------------------------------------------------
  counterChr <- 1
  for (indC in allChr) {
    accuVarRegChrRefB <- c()
    
    samplePath <- file.path(eventsDir, indS)
    fileSeg <- list.files(samplePath,
                          pattern = paste(indS, indC, ref2, "Seg.RData", sep = "."),
                          full.names = T)
    if (length(fileSeg) == 0) {
      counterChr <- counterChr + 1
      next()
    }
    
    # variabile segmenti: evRefB
    load(fileSeg)
    assign("evRefB", get(paste0("ev", ref2)))
    
    # se non ci sono indel o snv si pianta
    # NB se aggiungi un next, RICORDATI di incrementare counterChr
    indChrSnv <- which(vcfSnvQualFilt$fix$CHR == indC)
    indChrSid <- which(vcfSidQualFilt$fix$CHR == indC)
    vcfSnvQualChr <- list(vcfSnvQualFilt$fix[indChrSnv, ], vcfSnvQualFilt$format[indChrSnv, ])
    names(vcfSnvQualChr) <- c("fix", "format")
    vcfSidQualChr <- list(vcfSidQualFilt$fix[indChrSid, ], vcfSidQualFilt$format[indChrSid, ])
    names(vcfSidQualChr) <- c("fix", "format")
    
    # filter (sub)telomeric SNVs
    indTelVar <- which(vcfSnvQualChr$fix$POS < subTelLRefB[counterChr] |
                         vcfSnvQualChr$fix$POS > subTelRRefB[counterChr])
    if (length(indTelVar) != 0) {
      vcfSnvQualChr$fix <- vcfSnvQualChr$fix[-indTelVar, ]
      vcfSnvQualChr$format <- vcfSnvQualChr$format[-indTelVar, ]
    }
    # filter (sub)telomeric indels
    indTelSid <- which(vcfSidQualChr$fix$POS < subTelLRefB[counterChr] |
                         vcfSidQualChr$fix$POS > subTelRRefB[counterChr])
    if (length(indTelSid) != 0) {
      vcfSidQualChr$fix <- vcfSidQualChr$fix[-indTelSid, ]
      vcfSidQualChr$format <- vcfSidQualChr$format[-indTelSid, ]
    }
    
    # LOH events --------------------------------------------------------------
    
    indLOH <- which(evRefB$status == 0)
    evLOH <- evRefB[indLOH, ]
    numLOH <- length(indLOH)
    if (numLOH != 0) {
      nR <- nrow(evLOH)
      indSnvBon <- replicate(nR, list())
      indSidBon <- replicate(nR, list())
      lohIdentifierSnv <- replicate(nR, list())
      lohIdentifierSid <- lohIdentifierSnv
      # quali SNV & indel cadono nei segmenti tra FIRST e LAST marker
      for (indSeg in 1:nR) {
        indSnvBon[[indSeg]] <- which(vcfSnvQualChr$fix$POS > evLOH$start[indSeg] &
                                       vcfSnvQualChr$fix$POS < evLOH$end[indSeg])
        indSidBon[[indSeg]] <- which(vcfSidQualChr$fix$POS > evLOH$start[indSeg] &
                                       vcfSidQualChr$fix$POS < evLOH$end[indSeg])
        lohIdentifierSnv[[indSeg]] <- rep(indSeg, length(indSnvBon[[indSeg]]))
        lohIdentifierSid[[indSeg]] <- rep(indSeg, length(indSidBon[[indSeg]]))
      }
      lossHetIdenSnv <- unlist(lohIdentifierSnv)
      lossHetIdenSid <- unlist(lohIdentifierSid)
      snVar <- unlist(lapply(indSnvBon, length))
      idVar <- unlist(lapply(indSidBon, length))
      # annota tabella LOH con densità varianti
      evLOH <- data.frame(evLOH, snVar, SpecDec(c(snVar / evLOH$distLF), 7))
      colnames(evLOH)[14:15] <- c("SNV", "denSNV")
      evLOH <- data.frame(evLOH, idVar, SpecDec(c(idVar / evLOH$distLF), 7))
      colnames(evLOH)[16:17] <- c("InDel", "denInDel")
      
      # # numero regioni LOH
      # numLOH <- nrow(evLOH)
      # cat(c("STAT SNV LOH RefB", indS, indC, "# events:", numLOH, "\n"), sep = " ")
      # # Total lenght (LF)
      # sumLen <- sum(evLOH$distLF)
      # cat(c("STAT SNV LOH RefB", indS, indC, "Tot. length FL [bp]:",
      #       sumLen, "\n"), sep = " ")
      # # Total lenght (ES)
      # sumLenES <- sum(evLOH$distES)
      # cat(c("STAT SNV LOH RefB", indS, indC, "Tot. length SE [bp]:",
      #       sumLenES, "\n"), sep = " ")
      # # quanti marker
      # sumPol <- sum(evLOH$len)
      # cat(c("STAT SNV LOH RefB", indS, indC, "# marker:",
      #       sumPol, "\n"), sep = " ")
      
      # annota il vcf delle SNV con l'ID del segmento
      indSnvChr <- unlist(indSnvBon)
      # make sample output folder
      sampOutDir <- file.path(outDirSmallVarAnn, indS)
      dir.create(path = sampOutDir, showWarnings = F, recursive = T)
      # annota quante SNV het nel file degli eventi
      annHetSNV <- rep(0, numLOH)
      
      if (length(indSnvChr) != 0) {
        dfFormat <- data.frame(vcfSnvQualChr$format[indSnvChr, ], lossHetIdenSnv)
        colnames(dfFormat)[8] <- c("lohID")
        chrLoh <- list(vcfSidQualFilt$meta,
                       vcfSnvQualChr$fix[indSnvChr, ],
                       dfFormat)
        names(chrLoh) <- c("meta", "fix", "format")
        # save vcf RData SNV
        varName <- paste("vcfSNV", indS, indC, ref2, sep = "")
        assign(varName, chrLoh)
        save(list = varName, file = file.path(sampOutDir,
                                              paste(indS, indC, ref2, "Vcf.SNV.LOH.RData", sep = ".")))
        # accumulate variants for regions folder
        regVar <- data.frame(vcfSnvQualChr$fix[indSnvChr, ], rep(paste(formatNames, collapse = ":"), nrow(dfFormat)), 
                             apply(dfFormat[, c(1:ncol(dfFormat)-1)], 1, paste, collapse = ":" ))
        accuVarRegChrRefB <- rbind(accuVarRegChrRefB, regVar)
        
        # numero LOH con SNV
        withSNV <- which(evLOH$SNV != 0)
        lenEvSNV <- length(withSNV)
        cat(c("STAT SNV LOH", ref2, indS, indC, "# SNV events:", lenEvSNV, "\n"), sep = " ")
        # status of events with SNV
        statusSNV <- evLOH$status[withSNV]
        cat(c("STAT SNV LOH", ref2, indS, indC, "status of events with SNV:",
              paste(statusSNV, collapse = ","), "\n"), sep = " ")
        # quante SNV hom + SNV het
        sumSNV <- sum(evLOH$SNV)
        cat(c("STAT SNV LOH", ref2, indS, indC, "# SNV var:", sumSNV, "\n"), sep = " ")
        # quante SNV het
        numHetSNV <- length(which(chrLoh$format$GT == "0/1"))
        cat(c("STAT SNV LOH", ref2, indS, indC, "# SNV het:", numHetSNV, "\n"), sep = " ")
        
        # annota quante SNV het nel file degli eventi
        wSNV <- unique(chrLoh$format$lohID)
        for (indU in wSNV) {
          annHetSNV[indU] <- length(which(chrLoh$format$lohID == indU &
                                            chrLoh$format$GT == "0/1"))
        }
        evLOH <- data.frame(evLOH, annHetSNV)
        colnames(evLOH)[which(colnames(evLOH) == "annHetSNV")] <- "hetSNV"
        rm(chrLoh)
        rm(list = varName)
      } else {
        evLOH <- data.frame(evLOH, annHetSNV)
        colnames(evLOH)[which(colnames(evLOH) == "annHetSNV")] <- "hetSNV"
      }
      
      # annota il vcf delle Sid con l'ID del segmento
      indSidChr <- unlist(indSidBon)
      # annota quante Sid het nel file degli eventi
      annHetSid <- rep(0, numLOH)
      
      if (length(indSidChr) != 0) {
        dfFormat <- data.frame(vcfSidQualChr$format[indSidChr, ], lossHetIdenSid)
        colnames(dfFormat)[8] <- c("lohID")
        chrLoh <- list(vcfSidQualFilt$meta,
                       vcfSidQualChr$fix[indSidChr, ],
                       dfFormat)
        names(chrLoh) <- c("meta", "fix", "format")
        # save vcf InDel
        varName <- paste("vcfInDel", indS, indC, ref2, sep = "")
        assign(varName, chrLoh)
        save(list = varName, file = file.path(sampOutDir,
                                              paste(indS, indC, ref2, "Vcf.InDel.LOH.RData", sep = ".")))
        # accumulate variants for regions folder
        regVar <- data.frame(vcfSidQualChr$fix[indSidChr, ], rep(paste(formatNames, collapse = ":"), nrow(dfFormat)), 
                             apply(dfFormat[, c(1:ncol(dfFormat)-1)], 1, paste, collapse = ":" ))
        accuVarRegChrRefB <- rbind(accuVarRegChrRefB, regVar)
        
        # numero LOH con InDel
        withSid <- which(evLOH$InDel != 0)
        lenEvSid <- length(withSid)
        cat(c("STAT InDel LOH", ref2, indS, indC, "# InDel events:", lenEvSid, "\n"), sep = " ")
        # status of events with InDel
        statusSid <- evLOH$status[withSid]
        cat(c("STAT InDel LOH", ref2, indS, indC, "status of events with InDel:",
              paste(statusSid, collapse = ","), "\n"), sep = " ")
        # quante InDel hom + InDel het
        sumSid <- sum(evLOH$InDel)
        cat(c("STAT InDel LOH", ref2, indS, indC, "# InDel var:", sumSid, "\n"), sep = " ")
        # quante InDel het
        numHetSid <- length(which(chrLoh$format$GT == "0/1"))
        cat(c("STAT InDel LOH", ref2, indS, indC, "# InDel het:", numHetSid, "\n"), sep = " ")
        
        # annota quante InDel het nel file degli eventi
        wSid <- unique(chrLoh$format$lohID)
        for (indU in wSid) {
          annHetSid[indU] <- length(which(chrLoh$format$lohID == indU &
                                            chrLoh$format$GT == "0/1"))
        }
        evLOH <- data.frame(evLOH, annHetSid)
        colnames(evLOH)[which(colnames(evLOH) == "annHetSid")] <- "hetInDels"
        rm(chrLoh)
        rm(list = varName)
      } else {
        evLOH <- data.frame(evLOH, annHetSid)
        colnames(evLOH)[which(colnames(evLOH) == "annHetSid")] <- "hetInDels"
      }
      
      # save annotated data: evLOH
      varName <- paste(indS, indC, ref2, "evLOH", sep = "")
      # assign annotated LOH table to e.g. A315R52chrVIevLOH
      assign(varName, evLOH)
      save(list = varName, file = file.path(sampOutDir, paste(indS, indC, ref2, "LOH.RData", sep = ".")))
      # appending data to summary table by events
      write.table(x = data.frame(rep(indS, nR), rep(ref2, nR), get(varName)),
                  file = file.path(outDirSmallVarAnn, outFile),
                  col.names = F, row.names = F, quote = F, append = T, sep = "\t")
      
      rm(evLOH)
      rm(list = varName)
    }
    
    # HET events --------------------------------------------------------------
    
    indHET <- which(evRefB$status == 1)
    evHET <- evRefB[indHET, ]
    numHET <- length(indHET)
    if (numHET != 0) {
      nR <- nrow(evHET)
      indSnvBon <- replicate(nR, list())
      indSidBon <- replicate(nR, list())
      lohIdentifierSnv <- replicate(nR, list())
      lohIdentifierSid <- lohIdentifierSnv
      # quali SNV & indel cadono nei segmenti tra FIRST e LAST marker
      for (indSeg in 1:nR) {
        indSnvBon[[indSeg]] <- which(vcfSnvQualChr$fix$POS > evHET$start[indSeg] &
                                       vcfSnvQualChr$fix$POS < evHET$end[indSeg])
        indSidBon[[indSeg]] <- which(vcfSidQualChr$fix$POS > evHET$start[indSeg] &
                                       vcfSidQualChr$fix$POS < evHET$end[indSeg])
        lohIdentifierSnv[[indSeg]] <- rep(indSeg, length(indSnvBon[[indSeg]]))
        lohIdentifierSid[[indSeg]] <- rep(indSeg, length(indSidBon[[indSeg]]))
      }
      lossHetIdenSnv <- unlist(lohIdentifierSnv)
      lossHetIdenSid <- unlist(lohIdentifierSid)
      snVar <- unlist(lapply(indSnvBon, length))
      idVar <- unlist(lapply(indSidBon, length))
      # annota tabella HET con densità varianti
      evHET <- data.frame(evHET, snVar, SpecDec(c(snVar / evHET$distLF), 7))
      colnames(evHET)[14:15] <- c("SNV", "denSNV")
      evHET <- data.frame(evHET, idVar, SpecDec(c(idVar / evHET$distLF), 7))
      colnames(evHET)[16:17] <- c("InDel", "denInDel")
      
      # # numero regioni HET
      # numHET <- nrow(evHET)
      # cat(c("STAT SNV HET RefB", indS, indC, "# events:", numHET, "\n"), sep = " ")
      # # Total lenght (LF)
      # sumLen <- sum(evHET$distLF)
      # cat(c("STAT SNV HET RefB", indS, indC, "Tot. length FL [bp]:",
      #       sumLen, "\n"), sep = " ")
      # # Total lenght (ES)
      # sumLenES <- sum(evHET$distES)
      # cat(c("STAT SNV HET RefB", indS, indC, "Tot. length SE [bp]:",
      #       sumLenES, "\n"), sep = " ")
      # # quanti marker
      # sumPol <- sum(evHET$len)
      # cat(c("STAT SNV HET RefB", indS, indC, "# marker:",
      #       sumPol, "\n"), sep = " ")
      
      # annota il vcf delle SNV con l'ID del segmento
      indSnvChr <- unlist(indSnvBon)
      # make sample output folder
      sampOutDir <- file.path(outDirSmallVarAnn, indS)
      dir.create(path = sampOutDir, showWarnings = F, recursive = T)
      # annota quante SNV het nel file degli eventi
      annHetSNV <- rep(0, numHET)
      
      if (length(indSnvChr) != 0) {
        dfFormat <- data.frame(vcfSnvQualChr$format[indSnvChr, ], lossHetIdenSnv)
        colnames(dfFormat)[8] <- c("lohID")
        chrLoh <- list(vcfSidQualFilt$meta,
                       vcfSnvQualChr$fix[indSnvChr, ],
                       dfFormat)
        names(chrLoh) <- c("meta", "fix", "format")
        # save vcf RData SNV
        varName <- paste("vcfSNV", indS, indC, ref2, sep = "")
        assign(varName, chrLoh)
        save(list = varName, file = file.path(sampOutDir,
                                              paste(indS, indC, ref2, "Vcf.SNV.HET.RData", sep = ".")))
        # accumulate variants for regions folder
        regVar <- data.frame(vcfSnvQualChr$fix[indSnvChr, ], rep(paste(formatNames, collapse = ":"), nrow(dfFormat)), 
                             apply(dfFormat[, c(1:ncol(dfFormat)-1)], 1, paste, collapse = ":" ))
        accuVarRegChrRefB <- rbind(accuVarRegChrRefB, regVar)
        
        # numero HET con SNV
        withSNV <- which(evHET$SNV != 0)
        lenEvSNV <- length(withSNV)
        cat(c("STAT SNV HET", ref2, indS, indC, "# SNV events:", lenEvSNV, "\n"), sep = " ")
        # status of events with SNV
        statusSNV <- evHET$status[withSNV]
        cat(c("STAT SNV HET", ref2, indS, indC, "status of events with SNV:",
              paste(statusSNV, collapse = ","), "\n"), sep = " ")
        # quante SNV hom + SNV het
        sumSNV <- sum(evHET$SNV)
        cat(c("STAT SNV HET", ref2, indS, indC, "# SNV var:", sumSNV, "\n"), sep = " ")
        # quante SNV het
        numHetSNV <- length(which(chrLoh$format$GT == "0/1"))
        cat(c("STAT SNV HET", ref2, indS, indC, "# SNV het:", numHetSNV, "\n"), sep = " ")
        
        # annota quante SNV het nel file degli eventi
        wSNV <- unique(chrLoh$format$lohID)
        for (indU in wSNV) {
          annHetSNV[indU] <- length(which(chrLoh$format$lohID == indU &
                                            chrLoh$format$GT == "0/1"))
        }
        evHET <- data.frame(evHET, annHetSNV)
        colnames(evHET)[which(colnames(evHET) == "annHetSNV")] <- "hetSNV"
        rm(chrLoh)
        rm(list = varName)
      } else {
        evHET <- data.frame(evHET, annHetSNV)
        colnames(evHET)[which(colnames(evHET) == "annHetSNV")] <- "hetSNV"
      }
      
      # annota il vcf delle Sid con l'ID del segmento
      indSidChr <- unlist(indSidBon)
      # annota quante Sid het nel file degli eventi
      annHetSid <- rep(0, numHET)
      
      if (length(indSidChr) != 0) {
        dfFormat <- data.frame(vcfSidQualChr$format[indSidChr, ], lossHetIdenSid)
        colnames(dfFormat)[8] <- c("lohID")
        chrLoh <- list(vcfSidQualFilt$meta,
                       vcfSidQualChr$fix[indSidChr, ],
                       dfFormat)
        names(chrLoh) <- c("meta", "fix", "format")
        # save vcf InDel
        varName <- paste("vcfInDel", indS, indC, ref2, sep = "")
        assign(varName, chrLoh)
        save(list = varName, file = file.path(sampOutDir,
                                              paste(indS, indC, ref2, "Vcf.InDel.HET.RData", sep = ".")))
        # accumulate variants for regions folder
        regVar <- data.frame(vcfSidQualChr$fix[indSidChr, ], rep(paste(formatNames, collapse = ":"), nrow(dfFormat)), 
                             apply(dfFormat[, c(1:ncol(dfFormat)-1)], 1, paste, collapse = ":" ))
        accuVarRegChrRefB <- rbind(accuVarRegChrRefB, regVar)
        
        # numero HET con InDel
        withSid <- which(evHET$InDel != 0)
        lenEvSid <- length(withSid)
        cat(c("STAT InDel HET", ref2, indS, indC, "# InDel events:", lenEvSid, "\n"), sep = " ")
        # status of events with InDel
        statusSid <- evHET$status[withSid]
        cat(c("STAT InDel HET", ref2, indS, indC, "status of events with InDel:",
              paste(statusSid, collapse = ","), "\n"), sep = " ")
        # quante InDel hom + InDel het
        sumSid <- sum(evHET$InDel)
        cat(c("STAT InDel HET", ref2, indS, indC, "# InDel var:", sumSid, "\n"), sep = " ")
        # quante InDel het
        numHetSid <- length(which(chrLoh$format$GT == "0/1"))
        cat(c("STAT InDel HET", ref2, indS, indC, "# InDel het:", numHetSid, "\n"), sep = " ")
        
        # annota quante InDel het nel file degli eventi
        wSid <- unique(chrLoh$format$lohID)
        for (indU in wSid) {
          annHetSid[indU] <- length(which(chrLoh$format$lohID == indU &
                                            chrLoh$format$GT == "0/1"))
        }
        evHET <- data.frame(evHET, annHetSid)
        colnames(evHET)[which(colnames(evHET) == "annHetSid")] <- "hetInDels"
        rm(chrLoh)
        rm(list = varName)
      } else {
        evHET <- data.frame(evHET, annHetSid)
        colnames(evHET)[which(colnames(evHET) == "annHetSid")] <- "hetInDels"
      }
      
      # save annotated data: evHET
      varName <- paste(indS, indC, ref2, "evHET", sep = "")
      # assign annotated HET table to e.g. A315R52chrVIevHET
      assign(varName, evHET)
      save(list = varName, file = file.path(sampOutDir, paste(indS, indC, ref2, "HET.RData", sep = ".")))
      # appending data to summary table by events
      write.table(x = data.frame(rep(indS, nR), rep(ref2, nR), get(varName)),
                  file = file.path(outDirSmallVarAnn, outFile),
                  col.names = F, row.names = F, quote = F, append = T, sep = "\t")
      
      rm(evHET)
      rm(list = varName)
    }
    if (length(accuVarRegChrRefB$POS) != 0) {
      # sort acummulatore varianti
      srtAccuVarRegChrRefB <- accuVarRegChrRefB[order(accuVarRegChrRefB$POS), ]
      # scrivi l'accumulatore sorted
      write.table(x = srtAccuVarRegChrRefB, file = vcfFileOutRegRefB, col.names = F, row.names = F, quote = F,
                  append = T, sep = "\t", na = ".")
    }
    
    counterChr <- counterChr + 1
  } # chromosome loop RefB
} # sample loop

cat("Global time", "\n")
proc.time() - ptmGlob


