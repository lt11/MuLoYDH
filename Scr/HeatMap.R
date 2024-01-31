# header ------------------------------------------------------------------

# l'overlap viene fatto usando le coordinate start-end
# (anche se di solito usiamo first-last)
# perché sennò gli eventi con un solo marker
# (ma anche quelli con 2 o pochi marker vicini)
# non si vedono

# inoltre così si hanno più sfumature e il plot è
# più bellino esteticamente (e questo è un plot tendente 
# all'illustrativo)

# plot heat map of LOH events by chromosome
# using start (of LOH) and end-of-chromosome coordinate
# tutte le brutte cose che accadono in questo script
# (tipo plottare fino a end-of-chromosome)
# servono a evitare
# che si vedano dei gap nel plot, ahimè è un limite di ggplot

# detto meglio: il loop for (indE in 2:nLim) vede le regioni con zero eventi alla fine di un cromosoma
# ma non all'inizio perché conta le LOH in overlap (ovviamente) usando eventi a status ZERO

# nei cromosomi che finiscono con LOH a count > 0 si crea un evento
# a count ZERO che però non si vede nel plot
# (perché START = END)

# nei cromosomi che finiscono con LOH a count = 0 si deve
# modificare l'end dell'evento
# (che nel loop for (indE in 2:nLim) resta uguale a ZERO
# come da inizializzazione)

rm(list = ls())
options(stringsAsFactors = F)
library(ggplot2)
library(scales)

SpecDec <- function(x, k) as.numeric(format(round(x, k), nsmall = k))

# settings ----------------------------------------------------------------

# arguments
argsVal <- commandArgs(trailingOnly = T)
refBonLabel <- argsVal[1]
refBon <- argsVal[2]
baseDir <- argsVal[3]
# commentami
# refBonLabel <- "UWOPS034614"
# refBon <- "MA"
# baseDir <- "/Users/lorenzo/Prog/Yeasts/SummaryRTG/MARefB_A347R6"

# init
# input segments folder
inputSegmentsDir <- file.path(baseDir, "LOH", "AllSegments")

allChr <- c("chrI", "chrII", "chrIII", "chrIV", "chrV", "chrVI", "chrVII", "chrVIII", 
            "chrIX", "chrX", "chrXI", "chrXII", "chrXIII", "chrXIV", "chrXV", "chrXVI")

# centromeri, (sub)telomeri, chromosomes length
centrSRefA <- as.numeric(read.table(file = file.path(baseDir, "Ref", "Ann", 
                                                   paste0(refBonLabel, ".centromere.txt")), header = F)[, 1])
centrERefA <- as.numeric(read.table(file = file.path(baseDir, "Ref", "Ann", 
                                                   paste0(refBonLabel, ".centromere.txt")), header = F)[, 2])

subTelLRefA <- as.numeric(read.table(file = file.path(baseDir, "Ref", "Ann", paste0(refBonLabel, ".subtel.txt")))[, 1])
subTelRRefA <- as.numeric(read.table(file = file.path(baseDir, "Ref", "Ann", paste0(refBonLabel, ".subtel.txt")))[, 2])

chromLenRefA <- read.table(file = file.path(baseDir, "CNV", "GCdata", refBonLabel, "LenChr.txt"), header = F, sep = "\t")[, 3]

# switch to kb
centrSRefA <- centrSRefA / 1000
centrERefA <- centrERefA / 1000
subTelLRefA <- subTelLRefA / 1000
subTelRRefA <- subTelRRefA / 1000
chromLenRefA <- chromLenRefA / 1000

# plotting ----------------------------------------------------------------

olapAllChr <- c()
counterChr <- 1

for (indC in allChr) {
  # load allEv dataframe
  load(file.path(inputSegmentsDir, paste("Events", indC, "RData", sep = ".")))
  # remove ParHyb data
  allEv <- allEv[grep(pattern = "^ParHyb", allEv$hybrid, invert = T), ]
  # select strain to plot
  allEv <- allEv[which(allEv$strain == refBon), ]
  # remove heterozygous segments
  evRTG <- allEv[which(allEv$status == 0), ]
  # fake olapDF: un evento fake per coprire eventuali regioni
  # start senza LOH nel left arm
  # (che altrimenti verrebbero fuori nel plot
  # come dei buchi perché poi, vedi "eventi veri",
  # contiamo solo gli eventi LOH)
  olapDf <- data.frame(indC,
                       1, 
                       # * 1000 perché l'end degli eventi
                       # viene riscalato dopo, sigh
                       chromLenRefA[counterChr] * 1000, 
                       0, 
                       .0, 
                       subTelLRefA[counterChr], 
                       subTelRRefA[counterChr], 
                       (centrSRefA[counterChr] + centrERefA[counterChr]) / 2, 
                       chromLenRefA[counterChr], 
                       counterChr)
  colnames(olapDf) <- c("rep.indC..nOlap.", "sOut", "eOut", "countEv", "percEv", 
                        "subTelL", "subTelR", "centromere", "chromLen", "yPlotCoord")
  olapAllChr <- rbind.data.frame(olapAllChr, olapDf)
  
  # qui invece ci sono gli eventi veri
  if (nrow(evRTG) > 0) {
    # number of hybrids
    nHyb <- length(levels(evRTG$hybrid))
    # numero dei limiti di eventi in omozigosi
    nEv <- length(evRTG$start)
    nLim <- 2 * nEv
    # data ordering
    limEv <- c(evRTG$start, evRTG$end)
    # start label: +1, end label -1
    limLb <- c(rep(1, nEv), rep(-1, nEv))
    # sort label according to limits ordering
    limLb <- limLb[order(limEv)]
    # sort(limEv) = limEv[order(limEv)]
    limEv <- sort(limEv)
    
    countEv <- integer(length = nLim)
    sOut <- integer(length = nLim)
    eOut <- integer(length = nLim)
    countEv[1] <- 1
    sOut[1] <- limEv[1]
    for (indE in 2:nLim) {
      countEv[indE] <- countEv[indE - 1] + limLb[indE]
      sOut[indE] <- limEv[indE]
      eOut[indE - 1] <- sOut[indE] - 1
    }
    olapDf <- data.frame(sOut, eOut, countEv)
    # l'end dell'ultimo evento è ZERO, quindi va cambiato
    olapDf$eOut[nrow(olapDf)] <- chromLenRefA[counterChr] * 1E3
    # filter lines due to overlapping limits
    indOverlapping <- which(olapDf$eOut - olapDf$sOut < 0)
    if (length(indOverlapping) > 0) {
      olapDf <- olapDf[-indOverlapping, ]
    }

    percEv <- (olapDf$countEv / nHyb) * 100
    percEv <- SpecDec(x = percEv, k = 2)
    nOlap <- nrow(olapDf)
    cS <- centrSRefA[counterChr]
    cE <- centrERefA[counterChr]
    cPos <- (cS + cE) / 2
    centromere <- rep(cPos, nOlap)
    subTelL <- rep(subTelLRefA[counterChr], nOlap)
    subTelR <- rep(subTelRRefA[counterChr], nOlap)
    chromLen <- rep(chromLenRefA[counterChr], nOlap)
    # e qui si vola: gli end vanno fino in fondo per evitare i gap nei plot
    # l'ultimo va già fino in fondo (è stato sviluppato così)
    olapDf$eOut[1:(length(olapDf$eOut) - 1)] <- chromLenRefA[counterChr] * 1000
    yPlotCoord <- rep(counterChr, nOlap)
    olapDf <- data.frame(rep(indC, nOlap) , olapDf, 
                         percEv, subTelL, subTelR, centromere, chromLen, yPlotCoord)
    rm(percEv)
    olapAllChr <- rbind.data.frame(olapAllChr, olapDf)
  }
  counterChr <- counterChr + 1
}
# set names of column 1
colnames(olapAllChr)[1] <- "chr"
colnames(olapAllChr)[2] <- "start"
colnames(olapAllChr)[3] <- "end"

# converti in kb gli eventi
olapAllChr$start <- olapAllChr$start / 1000
olapAllChr$end <- olapAllChr$end / 1000
# lo start (cromosomi, telomero)
allStart <- 1 / 1000

# save plot data

# esempio: se plotti da start a chromLen
# risolvi il problema delle linee di discontinuità
# tra gli "eventi" e torna tutto perché il segmento
# a fine cromosoma è l'ultimo plottato (gli altri sono 
# sotto e non si vedono)

# plotting heatmap
# size markers
szCentr <- 2
szEvents <- 4
szChrom <- .5
szTel <- 2

# commentami, artefatti per vedere se i colori rendono
# olapAllChr$percEv[which(olapAllChr$chr == "chrVIII")[1] - 1] <- 50.0
# olapAllChr$percEv[which(olapAllChr$chr == "chrVIII")[1] - 2] <- 20.0

pHeatRainbow <- ggplot(olapAllChr) + 
  coord_cartesian(ylim = c(.5, 16.5)) + 
  scale_y_continuous(breaks = olapAllChr$yPlotCoord, labels = olapAllChr$chr) + 
  # heat maps
  geom_segment(aes(x = olapAllChr$start, y = olapAllChr$yPlotCoord, 
                   xend = olapAllChr$end, yend = olapAllChr$yPlotCoord, color = olapAllChr$percEv), 
               size = szEvents, lineend = "butt") + 
  scale_color_gradient2(low = "blue", mid = "yellow", 
                        high = "red", midpoint = 50, 
                        name = "% of samples\nbearing LOH") + 
  # (sub)telomers
  geom_segment(aes(x = allStart, y = olapAllChr$yPlotCoord,
                   xend = olapAllChr$subTelL, yend = olapAllChr$yPlotCoord),
               colour = "black", size =  szTel) + 
  geom_segment(aes(x = olapAllChr$subTelR, y = olapAllChr$yPlotCoord,
                   xend = olapAllChr$chromLen, yend = olapAllChr$yPlotCoord),
               colour = "black", size =  szTel) + 
  # centromere
  geom_point(shape = 1, aes(x = olapAllChr$centromere, y = olapAllChr$yPlotCoord), 
             size = szCentr) + 
  # plain background
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        panel.background = element_blank(), 
        # title style
        # plot.title = element_text(size = 16, face="bold", hjust = 1), 
        # plot.subtitle = element_text(size = 12, face="bold", hjust = 1), 
        
        # axis style
        axis.title = element_text(size = 16, face = "bold"), 
        axis.text = element_text(size = 16)) + 
  xlab(paste0("Position [kb]\n(", refBon, ")")) + 
  ylab("Chromosome") + 
  # integer x-axis values
  scale_x_continuous(labels = comma)

# save pdf
hmPath <- file.path(inputSegmentsDir, paste("HeatMap", refBon, "pdf", sep = "."))
pdf(file = hmPath)
print(pHeatRainbow)
dev.off()
# save jpg
hmPath <- file.path(inputSegmentsDir, paste("HeatMap", refBon, "jpg", sep = "."))
jpeg(file = hmPath)
print(pHeatRainbow)
dev.off()



