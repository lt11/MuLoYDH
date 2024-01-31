# header ------------------------------------------------------------------

# elimina i loci CBS432 che compaiono 
# nell'allineamento contro N17

# se finisce in pipeline va messo un if 
# che capisca se esiste già 
# il file CBS432_YPS128.intersect.snps.bck

rm(list = ls())
options(warn = 1)
options(stringsAsFactors = F)

argsVal <- commandArgs(trailingOnly = T)
ref1Label <- argsVal[1]
ref2Label <- argsVal[2]
refAux <- argsVal[3]
baseDir <- argsVal[4]

# ref1Label <- "CBS432"
# ref2Label <- "YPS128"
# refAux <- "N17"
# baseDir <- "/Users/lorenzo/Prog/Yeasts/SummaryRTG1/N17cbs432NA"

# il file da correggere
fileMark <- file.path(baseDir, "MUMmer", paste(ref1Label, ref2Label, sep = "_"), 
                      paste(paste(ref1Label, ref2Label, sep = "_"), "intersect", "snps", sep = "."))
# qui salvo una copia dell'originale
fileBackup <- file.path(baseDir, "MUMmer", paste(ref1Label, ref2Label, sep = "_"), 
                        paste(paste(ref1Label, ref2Label, sep = "_"), "intersect", "snps", "bck", sep = "."))

tableMark <- read.table(file = fileMark, skip = 4, header = F)
# le righe dell'header
headerOutFile <- readLines(con = fileMark, n = 4)

headerData <- c("pos1", "allele1", "allele2", "pos2", "buff", "dist", "lenR", "lenQ", "frm", "strand", "chrom1", "chrom2")
colnames(tableMark) <- headerData
# correttore 1
filePar <- file.path(baseDir, "MUMmer", "Correction", paste(ref1Label, refAux, sep = "_"), 
                     paste(paste(ref1Label, refAux, sep = "_"), "prt", "snps", sep = "."))
tablePar <- read.table(file = filePar, skip = 4, header = F)
colnames(tablePar) <- headerData
# correttore 2
fileRev <- file.path(baseDir, "MUMmer", "Correction", paste(refAux, ref1Label, sep = "_"), 
                     paste(paste(refAux, ref1Label, sep = "_"), "prt", "snps", sep = "."))
tableRev <- read.table(file = fileRev, skip = 4, header = F)
colnames(tableRev) <- headerData

# stringhe con pos_chrom, 9683_chrI, di
# tabella marker (CBS432 vs NA),
# tabella per correzioni CBS432 vs N17)
# tabella per correzioni (N17 vs CBS432)
stringMark <- paste(tableMark$pos1, tableMark$chrom1, sep = "_")
stringPar <- paste(tablePar$pos1, tablePar$chrom1, sep = "_")
stringRev <- paste(tableRev$pos2, tableRev$chrom2, sep = "_")

booTrue <- is.element(stringMark, stringPar) | is.element(stringMark, stringRev)

filtMark <- tableMark[!booTrue, ]

file.rename(fileMark, fileBackup)
# header del file
cat(headerOutFile, file = fileMark, sep = "\n")
# tabella filtrate
write.table(filtMark, file = fileMark, append = T, quote = F, col.names = F, row.names = F, sep = "\t")


