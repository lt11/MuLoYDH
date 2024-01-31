# header ------------------------------------------------------------------

# filtre le varianti condivise da N campioni
# secondo la soglia shareThreshold
# il calcolo delle varianti buone è fatto su un vcf
# multisample, poi le varianti buone vengono
# ricercate nei vcf singoli per recuperare 
# le QUAL

rm(list = ls())
options(stringsAsFactors = F)

# function(s) -------------------------------------------------------------

# calcola il nunmero di campioni nel vcf che hanno una variante
# i campioni dove la variante è assente riportano ".:.:.:.:.:.:."
# mentre se la variante c'è "1:114,13:4:0:0,0:1,3:1,3" abbiamo un
# sacco di bei numerini
PresVar <- function(x) length(grep("[[:digit:]]+", x))

# settings ----------------------------------------------------------------

argsVal <- commandArgs(trailingOnly = T)
subDir <- argsVal[1]
baseDir <- argsVal[2]
shareThreshold <- 0.95
# commentami
# subDir <- "Ploidy1"
# baseDir <- "/Users/lorenzo/Prog/Yeasts/DevMuLo4"

inputDir <- file.path(baseDir, "VariantCalls", subDir, "ParSubtracted")
outDir <- file.path(baseDir, "VariantCalls", subDir, "HybSubtracted")

mergedVcf <- read.table(file = file.path(inputDir, "Merged.txt.gz"), sep = "\t", comment.char = "{")
# chissà perché R deve scassare il cazzo
# si rifiuta di leggere l'header
headerStr <- mergedVcf[1, ]
headerStr[1] <- "CHROM"
mergedVcf <- mergedVcf[-1, ]
colnames(mergedVcf) <- headerStr

nCol <- ncol(mergedVcf)
nSamples <- nCol - 9

# numero di campioni con la variante
nSampVar <- apply(mergedVcf[, c(10:nCol)], 1, PresVar)
# frazione di campioni con una variante
fractSsamples <- nSampVar / nSamples
filtVcf <- mergedVcf[which(fractSsamples <= shareThreshold), ]
# filtra i vcf singoli così recuperi le QUAL
# NB la colonna 10 è quella del primo sample
for (indS in 10:nCol) {
  # nome del vcf single sample
  vcfName <- paste0(paste(unlist(strsplit(as.character(headerStr[indS]), split = "\\."))[1:2], collapse = "."), ".vcf")
  vcfPath <- file.path(inputDir, paste0(vcfName, ".gz"))
  headerFileOut <- readLines(con = vcfPath)
  headerFileOut <- headerFileOut[grep("^#", headerFileOut)]
  # file output
  outVcf <- file.path(outDir, vcfName)
  # write header
  cat(headerFileOut, file = outVcf, sep = "\n")

  # prendi le varianti del campione
  # presenti nel vcf multisample
  indBonFiltVcf <- grep(pattern = "[[:digit:]]+", x = filtVcf[, indS])
  # stringa varianti buone da confrontare con la stringa da vcf single sample
  strBonFiltVcf <- paste(filtVcf$CHROM[indBonFiltVcf], filtVcf$POS[indBonFiltVcf], sep = "_")
  
  # read vcf data
  singSampleVcf <- read.table(vcfPath)
  strSingSampVcf <- paste(singSampleVcf$V1, singSampleVcf$V2, sep = "_")
  # filter single sample vcf
  singSampleVcfFilt <- singSampleVcf[is.element(strSingSampVcf, strBonFiltVcf), ]
  
  # append al file dove ho stampato l'header
  write.table(singSampleVcfFilt, file = outVcf, append = T, col.names = F, row.names = F, quote = F, sep = "\t")
}



