# header ------------------------------------------------------------------

# makes segments plots for all samples by chromosome
# and saves all samples data by chromosome (for heat maps) 
# [salva anche gli eventi status 2 ma li filtro dopo negli Heat.R]
# creates folder "AllSegments" for output pdfs
# and SE, FL subfolders for plots
# if plotSeg is TRUE
# makes all samples plots
# questa riga non serve a nulla

# very peculiar stuff:
# length(chromStartRefB) -> 1
# length(chromStartRefA) -> 16

rm(list = ls())
options(stringsAsFactors = F)
library(ggplot2)
library(scales)

# settings ----------------------------------------------------------------

# arguments
argsVal <- commandArgs(trailingOnly = T)
ref1Label <- argsVal[1]
ref2Label <- argsVal[2]
ref1 <- argsVal[3]
ref2 <- argsVal[4]
baseDir <- argsVal[5]

# init
# to plot or not to plot? this is the question
plotSeg <- F
# input folder
dirData <- file.path(baseDir, "LOH")
# main output folder
dirName <- "AllSegments"

outDirHis <- file.path(dirData, dirName, "HistoLength")
outDirData <- file.path(dirData, dirName)
if (plotSeg) {
  outDirFL <- file.path(dirData, dirName, "FL")
  outDirSE <- file.path(dirData, dirName, "SE")
  dir.create(path = outDirFL, showWarnings = F, recursive = T)
  dir.create(path = outDirSE, showWarnings = F, recursive = T)
}
dir.create(path = outDirHis, showWarnings = F, recursive = T)

summaryCsv <- file.path(outDirData, paste("AllEvents", "csv", sep = "."))
if (file.exists(summaryCsv)) file.remove(summaryCsv)

allChr <- c("chrI", "chrII", "chrIII", "chrIV", "chrV", "chrVI", "chrVII", "chrVIII", 
            "chrIX", "chrX", "chrXI", "chrXII", "chrXIII", "chrXIV", "chrXV", "chrXVI")

# centromeri, (sub)telomeri, chromosomes length
chromLenRefA <- read.table(file = file.path(baseDir, "CNV", "GCdata", ref1Label, "LenChr.txt"), header = F, sep = "\t")[, 3]
chromLenRefB <- read.table(file = file.path(baseDir, "CNV", "GCdata", ref2Label, "LenChr.txt"), header = F, sep = "\t")[, 3]

centrSRefA <- as.numeric(read.table(file = file.path(baseDir, "Ref", "Ann", 
                                                   paste0(ref1Label, ".centromere.txt")), header = F)[, 1])
centrERefA <- as.numeric(read.table(file = file.path(baseDir, "Ref", "Ann", 
                                                   paste0(ref1Label, ".centromere.txt")), header = F)[, 2])
centrSRefB <- as.numeric(read.table(file = file.path(baseDir, "Ref", "Ann", 
                                                   paste0(ref2Label, ".centromere.txt")), header = F)[, 1])
centrERefB <- as.numeric(read.table(file = file.path(baseDir, "Ref", "Ann", 
                                                   paste0(ref2Label, ".centromere.txt")), header = F)[, 2])

subTelLRefA <- as.numeric(read.table(file = file.path(baseDir, "Ref", "Ann", paste0(ref1Label, ".subtel.txt")))[, 1])
subTelRRefA <- as.numeric(read.table(file = file.path(baseDir, "Ref", "Ann", paste0(ref1Label, ".subtel.txt")))[, 2])
subTelLRefB <- as.numeric(read.table(file = file.path(baseDir, "Ref", "Ann", paste0(ref2Label, ".subtel.txt")))[, 1])
subTelRRefB <- as.numeric(read.table(file = file.path(baseDir, "Ref", "Ann", paste0(ref2Label, ".subtel.txt")))[, 2])

# switch to kb
centrSRefA <- centrSRefA / 1000
centrERefA <- centrERefA / 1000
centrSRefB <- centrSRefB / 1000
centrERefB <- centrERefB / 1000
subTelLRefB <- subTelLRefB / 1000
subTelLRefA <- subTelLRefA / 1000
subTelRRefB <- subTelRRefB / 1000
subTelRRefA <- subTelRRefA / 1000
chromLenRefB <- chromLenRefB / 1000
chromLenRefA <- chromLenRefA / 1000

# add chromosome start (for shifting coordinates)
chromStartRefA <- 0.001
chromStartRefB <- 0.001

# centromere position and shift
cPosRefA <- (centrSRefA + centrERefA) / 2
cPosRefB <- (centrSRefB + centrERefB) / 2
shiftRefBRefA <- cPosRefB - cPosRefA

# shifting RefA (sub)telomers, chromosome start and end
subTelLRefA <- subTelLRefA + shiftRefBRefA
subTelRRefA <- subTelRRefA + shiftRefBRefA
chromLenRefA <- chromLenRefA + shiftRefBRefA
chromStartRefA <- chromStartRefA + shiftRefBRefA

# shift centromere
cPosRefA <- cPosRefA + shiftRefBRefA

ptmChr <- proc.time()
ptmGlob <- proc.time()

# all segments files path
allFiles <- list.files(pattern = ".Seg.RData$", path = dirData, recursive = T, full.names = T)

# init variable to accumulate events @ create chr*Events
allEvVar <- c()

# write header string
headerString <- c("chr,start,first,last,end,status,len,denES,evrES,distES,denLF,evrLF,distLF,strain,hybrid,color,yVal")
write.table(x = headerString, 
            file = summaryCsv, append = F, sep = "", row.names = F, col.names = F, quote = F)

# plotting ----------------------------------------------------------------

# carica i dati di tutti i campioni, un solo chr
chrCount <- 1
for (indC in allChr) {
  allEv <- c()
  myChrFiles <- grep(pattern = paste(".", indC, ".", sep = ""), x = allFiles, fixed = T, value = T)
  # if control sample has one whole-chromosome deletion e.g. chrV, all *.chrV.*.Seg.RData files do not exist
  if (length(myChrFiles) == 0) {
    chrCount <- chrCount + 1
    next()
  }
  
  # files of segments by reference, chromosome, sample in sample folder under Events
  # e.g. Events/A452R77/A452R77.chrI.RefB.Seg.RData
  for (indD in myChrFiles) {
    load(indD) 
    bName <- basename(indD)
    nameEl <- unlist(strsplit(x = bName, split = "\\."))
    sampID <- nameEl[1]
    refName <- nameEl[3]
    # reconstruct dataframe name
    evDF <- paste("ev", refName, sep = "")
    # get: call an R object using a character string
    nrDF <- nrow(get(evDF))
    appEv <- data.frame(get(evDF), rep(refName, nrDF), rep(sampID, nrDF))
    colnames(appEv)[14:15] <- c("strain", "hybrid")
    allEv <- rbind(allEv, appEv)
  }
  allEv$hybrid <- factor(x = allEv$hybrid, levels = unique(allEv$hybrid))
  # plot
  # add color and yVal
  nrEv <- nrow(allEv)
  color <- character(length = nrEv)
  color[which(allEv$status == 1)] <- "darkgrey"
  color[which(allEv$status == 0 & allEv$strain == ref1)] <- "blue"
  color[which(allEv$status == 2 & allEv$strain == ref1)] <- "red"
  color[which(allEv$status == 0 & allEv$strain == ref2)] <- "red"
  color[which(allEv$status == 2 & allEv$strain == ref2)] <- "blue"
  color <- factor(color)
  
  yRefA <- .25
  yRefB <- .75
  
  yVal <- numeric(length = nrEv)
  yVal[which(allEv$strain == ref1)] <- yRefA
  yVal[which(allEv$strain == ref2)] <- yRefB
  allEv <- data.frame(allEv, color, yVal)
  
  # saving data
  save(allEv, file = file.path(outDirData, paste("Events", indC, "RData", sep = ".")))
  write.table(x = allEv, file = summaryCsv, append = T, sep = ",", row.names = F, col.names = F, quote = F)

  # size markers
  szCentr <- 4
  szEvents <- 10
  szChromTel <- 2
  
  # switch data to kb
  allEv$start <- allEv$start / 1000
  allEv$end <- allEv$end / 1000
  allEv$first <- allEv$first / 1000
  allEv$last <- allEv$last / 1000
  
  # shifting RefA data
  indRefA <- which(allEv$strain == ref1)
  allEv[indRefA, c(2:5)] <- allEv[indRefA, c(2:5)] + shiftRefBRefA[chrCount]
  
  # plot first-last marker
  if (plotSeg) {
    pAllSeg <- ggplot(allEv) + 
      # title
      ggtitle(indC) + 
      coord_cartesian(ylim = c(.0, 1.)) + 
      scale_y_continuous(breaks = allEv$yVal, labels = allEv$strain) + 
      # events
      geom_segment(aes(x = allEv$first, y = allEv$yVal, xend = allEv$last, yend = allEv$yVal), 
                   colour = allEv$color, size = szEvents) + 
      # whole chromosome
      geom_segment(aes(x = chromStartRefB, y = yRefB, xend = chromLenRefB[chrCount], yend = yRefB), 
                   colour = "black", size =  szChromTel) + 
      geom_segment(aes(x = chromStartRefA[chrCount], y = yRefA, xend = chromLenRefA[chrCount], yend = yRefA), 
                   colour = "black", size =  szChromTel) + 
      # (sub)telomers
      geom_segment(aes(x = chromStartRefB, y = yRefB, xend = subTelLRefB[chrCount], yend = yRefB), 
                   colour = "orange", size =  szChromTel) + 
      geom_segment(aes(x = subTelRRefB[chrCount], y = yRefB, xend = chromLenRefB[chrCount], yend = yRefB), 
                   colour = "orange", size =  szChromTel) + 
      geom_segment(aes(x = chromStartRefA[chrCount], y = yRefA, xend = subTelLRefA[chrCount], yend = yRefA), 
                   colour = "orange", size =  szChromTel) + 
      geom_segment(aes(x = subTelRRefA[chrCount], y = yRefA, xend = chromLenRefA[chrCount], yend = yRefA), 
                   colour = "orange", size =  szChromTel) + 
      # centromere
      geom_point(shape = 1, aes(x = cPosRefB[chrCount], y = yRefB), size = szCentr) + 
      geom_point(shape = 1, aes(x = cPosRefA[chrCount], y = yRefA), size = szCentr) + 
      # plain background
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            panel.border = element_blank(), 
            panel.background = element_blank(), 
            # title style
            plot.title = element_text(size = 16, face = "bold", hjust = 1), 
            plot.subtitle = element_text(size = 12, face = "bold", hjust = 1), 
            # axis style
            axis.title = element_text(size = 16, face = "bold"),
            axis.text = element_text(size = 16), 
            # subtitle format
            strip.background = element_rect(fill = "white", colour = "white"), 
            strip.text.x = element_text(size = 12, face = "bold"), 
            plot.margin = unit(c(1, 1, 1, 1), "cm")) + 
      # position shifted to ref2
      xlab(paste0("Position [kb]\n(", ref2, ")")) + 
      ylab(NULL) + 
      facet_wrap(facets = ~ hybrid, ncol = 4, scales = "free_x") + 
      # integer x-axis values
      scale_x_continuous(labels = comma)
    
    pPath <- file.path(outDirFL, paste("AllSeg", indC, "pdf", sep = "."))
    pdf(file = pPath, width = 26, height = 10)
    print(pAllSeg)
    dev.off()
  }
  
  # plot start-end event
  if (plotSeg) {
    pAllSeg <- ggplot(allEv) + 
      # title
      ggtitle(indC) + 
      coord_cartesian(ylim = c(.0, 1.)) + 
      scale_y_continuous(breaks = allEv$yVal, labels = allEv$strain) + 
      # events
      geom_segment(aes(x = allEv$start, y = allEv$yVal, xend = allEv$end, yend = allEv$yVal), 
                   colour = allEv$color, size = szEvents) + 
      # whole chromosome
      geom_segment(aes(x = chromStartRefB, y = yRefB, xend = chromLenRefB[chrCount], yend = yRefB), 
                   colour = "black", size =  szChromTel) + 
      geom_segment(aes(x = chromStartRefA[chrCount], y = yRefA, xend = chromLenRefA[chrCount], yend = yRefA), 
                   colour = "black", size =  szChromTel) + 
      # (sub)telomers
      geom_segment(aes(x = chromStartRefB, y = yRefB, xend = subTelLRefB[chrCount], yend = yRefB), 
                   colour = "orange", size =  szChromTel) + 
      geom_segment(aes(x = subTelRRefB[chrCount], y = yRefB, xend = chromLenRefB[chrCount], yend = yRefB), 
                   colour = "orange", size =  szChromTel) + 
      geom_segment(aes(x = chromStartRefA[chrCount], y = yRefA, xend = subTelLRefA[chrCount], yend = yRefA), 
                   colour = "orange", size =  szChromTel) + 
      geom_segment(aes(x = subTelRRefA[chrCount], y = yRefA, xend = chromLenRefA[chrCount], yend = yRefA), 
                   colour = "orange", size =  szChromTel) + 
      # centromere
      geom_point(shape = 1, aes(x = cPosRefB[chrCount], y = yRefB), size = szCentr) + 
      geom_point(shape = 1, aes(x = cPosRefA[chrCount], y = yRefA), size = szCentr) + 
      # plain background
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            panel.border = element_blank(), 
            panel.background = element_blank(), 
            # title style
            plot.title = element_text(size = 16, face="bold", hjust = 1), 
            plot.subtitle = element_text(size = 12, face="bold", hjust = 1), 
            # axis style
            axis.title = element_text(size = 16, face = "bold"),
            axis.text = element_text(size = 16), 
            # subtitle format
            strip.background = element_rect(fill = "white", colour = "white"), 
            strip.text.x = element_text(size = 12, face = "bold"), 
            plot.margin = unit(c(1, 1, 1, 1), "cm")) + 
      # position shifted to ref2
      xlab(paste0("Position [kb]\n(", ref2, ")")) + 
      ylab(NULL) + 
      facet_wrap(facets = ~ hybrid, ncol = 4, scales = "free_x") + 
      # integer x-axis values
      scale_x_continuous(labels = comma)
    
    pPath <- file.path(outDirSE, paste("AllSeg", indC, "pdf", sep = "."))
    pdf(file = pPath, width = 26, height = 10)
    print(pAllSeg)
    dev.off()
  }
  
  chrCount <- chrCount + 1
  
  # create chr*Events
  evVar <- paste(indC, "Events", sep = "")
  assign(evVar, value = allEv)
  allEvVar <- c(allEvVar, evVar)
} # chromosome loop
proc.time() - ptmChr

# prepare data for events density plots -----------------------------------
concEvents <- paste("globEvents <- rbind.data.frame(", paste(allEvVar, collapse = ", "), ")", sep = "")
# eval -> execute and create globEvents variable
eval(parse(text = concEvents))

# convertion to kb
globEvents$distES <- globEvents$distES / 1000
globEvents$distLF <- globEvents$distLF / 1000

# ref1 plots ----------------------------------------------------------------

plotEvents <- globEvents[which(globEvents$strain == ref1 & globEvents$status == 0), ]

# small events threshold [kb] & parameters
smallThres <- 25
binWlarge <- 25
binWsmall <- .5

indSmall <- which(plotEvents$distLF < smallThres)
plotSmallEvents <- plotEvents[indSmall, ]
plotLargeEvents <- plotEvents[-indSmall, ]

# small events plots ------------------------------------------------------
nEv <- nrow(plotSmallEvents)

# plot histogram with kernel density: small LF
pLFs <- ggplot(plotSmallEvents, aes(x = distLF)) + 
  ggtitle(paste0("First-Last Marker Distance, ", ref1 , " strain"), 
          subtitle = paste("Length < ", smallThres, " kb; ", 
                           "N = ", nEv, "; ", 
                           "Bin = ", binWsmall, " kb", 
                           sep = "")) + 
  geom_histogram(aes(y = ..density..), # histogram with density instead of count on y-axis
                 binwidth = binWsmall, 
                 colour = "black", fill = "white") + 
  geom_density(alpha = .25, fill = "blue") + # overlay with transparent kernel density plot
  xlab("Length [kb]") + 
  ylab("Density") + 
  theme(
    # title style
    plot.title = element_text(size = 16, hjust = 1), 
    plot.subtitle = element_text(size = 12, face = "bold", hjust = 1), 
    # axis style
    axis.title = element_text(size = 16, face = "bold"), 
    axis.text = element_text(size = 16))
pPath <- file.path(outDirHis, paste0(ref1, ".LF.small.pdf"))
pdf(file = pPath)
print(pLFs)
dev.off()
pPath <- file.path(outDirHis, paste0(ref1, ".LF.small.jpg"))
jpeg(file = pPath)
print(pLFs)
dev.off()

# rifare filtri sulla distES del dataframe, ora sono solo su indSmall <- which(plotEvents$distLF < smallThres)
# plot histogram with kernel density: small ES
# pESs <- ggplot(plotSmallEvents, aes(x = distES)) + 
#   ggtitle("Start-End Position Distance, RefA strain", 
#           subtitle = paste("Length < ", smallThres, " kb; ", 
#                            "N = ", nEv, "; ", 
#                            "Bin = ", binWsmall, " kb", 
#                            sep = "")) + 
#   geom_histogram(aes(y = ..density..), # histogram with density instead of count on y-axis
#                  binwidth = binWsmall, 
#                  colour = "black", fill = "white") + 
#   geom_density(alpha = .25, fill = "blue") + # overlay with transparent kernel density plot
#   xlab("Length [kb]") + 
#   ylab("Density") + 
#   theme(
#     # title style
#     plot.title = element_text(size = 16, hjust = 1), 
#     plot.subtitle = element_text(size = 12, face = "bold", hjust = 1), 
#     # axis style
#     axis.title = element_text(size = 16, face = "bold"), 
#     axis.text = element_text(size = 16))
# pPath <- file.path(outDirHis, "RefA.ES.small.pdf")
# pdf(file = pPath)
# print(pESs)
# dev.off()

# large events plots ------------------------------------------------------
nEv <- nrow(plotLargeEvents)

# plot histogram with kernel density: large LF
pFLl <- ggplot(plotLargeEvents, aes(x = distLF)) + 
  ggtitle(paste0("First-Last Marker Distance, ", ref1 , " strain"), 
          # unicode character ≥ \u2265
          subtitle = paste("Length > ", smallThres, " kb; ", 
                           "N = ", nEv, "; ", 
                           "Bin = ", binWlarge, " kb", 
                           sep = "")) + 
  geom_histogram(aes(y = ..density..), # histogram with density instead of count on y-axis
                 binwidth = binWlarge, 
                 colour = "black", fill = "white") + 
  geom_density(alpha = .25, fill = "blue") + # overlay with transparent kernel density plot
  xlab("Length [kb]") + 
  ylab("Density") + 
  theme(
    # title style
    plot.title = element_text(size = 16, hjust = 1), 
    plot.subtitle = element_text(size = 12, face = "bold", hjust = 1), 
    # axis style
    axis.title = element_text(size = 16, face = "bold"), 
    axis.text = element_text(size = 16))
pPath <- file.path(outDirHis, paste0(ref1, ".LF.large.pdf"))
pdf(file = pPath)
print(pFLl)
dev.off()
pPath <- file.path(outDirHis, paste0(ref1, ".LF.large.jpg"))
jpeg(file = pPath)
print(pFLl)
dev.off()

# rifare filtri sulla distES del dataframe, ora sono solo su indSmall <- which(plotEvents$distLF < smallThres)
# plot histogram with kernel density: large ES
# pESl <- ggplot(plotLargeEvents, aes(x = distES)) + 
#   ggtitle("Start-End Position Distance, RefA strain", 
#           # unicode character ≥ \u2265
#           subtitle = paste("Length > ", smallThres, " kb; ", 
#                            "N = ", nEv, "; ", 
#                            "Bin = ", binWlarge, " kb", 
#                            sep = "")) + 
#   geom_histogram(aes(y = ..density..), # histogram with density instead of count on y-axis
#                  binwidth = binWlarge, 
#                  colour = "black", fill = "white") + 
#   geom_density(alpha = .25, fill = "blue") + # overlay with transparent kernel density plot
#   xlab("Length [kb]") + 
#   ylab("Density") + 
#   theme(
#     # title style
#     plot.title = element_text(size = 16, hjust = 1), 
#     plot.subtitle = element_text(size = 12, face = "bold", hjust = 1), 
#     # axis style
#     axis.title = element_text(size = 16, face = "bold"), 
#     axis.text = element_text(size = 16))
# pPath <- file.path(outDirHis, "RefA.ES.large.pdf")
# pdf(file = pPath)
# print(pESl)
# dev.off()

# ref2 plots ----------------------------------------------------------------

plotEvents <- globEvents[which(globEvents$strain == ref2 & globEvents$status == 0), ]

# small events threshold [kb] & parameters
smallThres <- 25
binWlarge <- 25
binWsmall <- .5

indSmall <- which(plotEvents$distLF < smallThres)
plotSmallEvents <- plotEvents[indSmall, ]
plotLargeEvents <- plotEvents[-indSmall, ]

# small events plots ------------------------------------------------------
nEv <- nrow(plotSmallEvents)

# plot histogram with kernel density: small LF
pLFs <- ggplot(plotSmallEvents, aes(x = distLF)) + 
  ggtitle(paste0("First-Last Marker Distance, ", ref2 , " strain"), 
          subtitle = paste("Length < ", smallThres, " kb; ", 
                           "N = ", nEv, "; ", 
                           "Bin = ", binWsmall, " kb", 
                           sep = "")) + 
  geom_histogram(aes(y = ..density..), # histogram with density instead of count on y-axis
                 binwidth = binWsmall, 
                 colour = "black", fill = "white") + 
  geom_density(alpha = .25, fill = "red") + # overlay with transparent kernel density plot
  xlab("Length [kb]") + 
  ylab("Density") + 
  theme(
    # title style
    plot.title = element_text(size = 16, hjust = 1), 
    plot.subtitle = element_text(size = 12, face = "bold", hjust = 1), 
    # axis style
    axis.title = element_text(size = 16, face = "bold"), 
    axis.text = element_text(size = 16))
pPath <- file.path(outDirHis, paste0(ref2, ".LF.small.pdf"))
pdf(file = pPath)
print(pLFs)
dev.off()
pPath <- file.path(outDirHis, paste0(ref2, ".LF.small.jpg"))
jpeg(file = pPath)
print(pLFs)
dev.off()

# rifare filtri sulla distES del dataframe, ora sono solo su indSmall <- which(plotEvents$distLF < smallThres)
# plot histogram with kernel density: small ES
# pESs <- ggplot(plotSmallEvents, aes(x = distES)) + 
#   ggtitle("Start-End Position Distance, RefB strain", 
#           subtitle = paste("Length < ", smallThres, " kb; ", 
#                            "N = ", nEv, "; ", 
#                            "Bin = ", binWsmall, " kb", 
#                            sep = "")) + 
#   geom_histogram(aes(y = ..density..), # histogram with density instead of count on y-axis
#                  binwidth = binWsmall, 
#                  colour = "black", fill = "white") + 
#   geom_density(alpha = .25, fill = "red") + # overlay with transparent kernel density plot
#   xlab("Length [kb]") + 
#   ylab("Density") + 
#   theme(
#     # title style
#     plot.title = element_text(size = 16, hjust = 1), 
#     plot.subtitle = element_text(size = 12, face = "bold", hjust = 1), 
#     # axis style
#     axis.title = element_text(size = 16, face = "bold"), 
#     axis.text = element_text(size = 16))
# pPath <- file.path(outDirHis, "RefB.ES.small.pdf")
# pdf(file = pPath)
# print(pESs)
# dev.off()

# large events plots ------------------------------------------------------
nEv <- nrow(plotLargeEvents)

# plot histogram with kernel density: large LF
pFLl <- ggplot(plotLargeEvents, aes(x = distLF)) + 
  ggtitle(paste0("First-Last Marker Distance, ", ref2 , " strain"), 
          # unicode character ≥ \u2265
          subtitle = paste("Length > ", smallThres, " kb; ", 
                           "N = ", nEv, "; ", 
                           "Bin = ", binWlarge, " kb", 
                           sep = "")) + 
  geom_histogram(aes(y = ..density..), # histogram with density instead of count on y-axis
                 binwidth = binWlarge, 
                 colour = "black", fill = "white") + 
  geom_density(alpha = .25, fill = "red") + # overlay with transparent kernel density plot
  xlab("Length [kb]") + 
  ylab("Density") + 
  theme(
    # title style
    plot.title = element_text(size = 16, hjust = 1), 
    plot.subtitle = element_text(size = 12, face = "bold", hjust = 1), 
    # axis style
    axis.title = element_text(size = 16, face = "bold"), 
    axis.text = element_text(size = 16))
pPath <- file.path(outDirHis, paste0(ref2, ".LF.large.pdf"))
pdf(file = pPath)
print(pFLl)
dev.off()
pPath <- file.path(outDirHis, paste0(ref2, ".LF.large.jpg"))
jpeg(file = pPath)
print(pFLl)
dev.off()

# plot histogram with kernel density: large ES
# pESl <- ggplot(plotLargeEvents, aes(x = distES)) + 
#   ggtitle("Start-End Position Distance, RefB strain", 
#           # unicode character ≥ \u2265
#           subtitle = paste("Length > ", smallThres, " kb; ", 
#                            "N = ", nEv, "; ", 
#                            "Bin = ", binWlarge, " kb", 
#                            sep = "")) + 
#   geom_histogram(aes(y = ..density..), # histogram with density instead of count on y-axis
#                  binwidth = binWlarge, 
#                  colour = "black", fill = "white") + 
#   geom_density(alpha = .25, fill = "red") + # overlay with transparent kernel density plot
#   xlab("Length [kb]") + 
#   ylab("Density") + 
#   theme(
#     # title style
#     plot.title = element_text(size = 16, hjust = 1), 
#     plot.subtitle = element_text(size = 12, face = "bold", hjust = 1), 
#     # axis style
#     axis.title = element_text(size = 16, face = "bold"), 
#     axis.text = element_text(size = 16))
# pPath <- file.path(outDirHis, "RefB.ES.large.pdf")
# pdf(file = pPath)
# print(pESl)
# dev.off()

proc.time() - ptmGlob


