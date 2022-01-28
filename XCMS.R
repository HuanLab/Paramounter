###############################################################
directory <- "F:/Jian_Guo/Parameter_optimize_20201105/mzXML"
###############################################################
library(xcms)
library(MSnbase)
library(dplyr)
setwd(directory)
rawfile <- list.files(pattern = ".mzXML")
data <- readMSData(rawfile, mode = "onDisk")
cwp <- CentWaveParam(ppm=10,
                    peakwidth=c(5, 24),
                    mzdiff = -0.01,
                    snthresh = 1.065,
                    integrate = 2,
                    prefilter = c(3, 134),
                    noise = 134)
data <- findChromPeaks(data, param = cwp)
data_filtered <- filterMsLevel(data, msLevel = 1L)
xset <- as(data_filtered, 'xcmsSet') 
#ALIGNMENT  
xset@peaks <- xset@peaks[order(xset@peaks[,11]),]
table <- as.data.frame(xset@peaks)
for(n in (1:length(rawfile))){
    sampleOutput <- table[table$sample == n, ]
    sampleOutput <- sampleOutput[order(sampleOutput[,1]),]
    colnames(sampleOutput)[9] <- "intMax"
    row.names(sampleOutput) <- 1:nrow(sampleOutput)
    # write.csv(sampleOutput, file = paste(n,"featuretable.csv",sep = "_"))
}
xset <- group(xset, bw = 5, minfrac = 0.5, mzwid = 0.012, minsamp = 1, max = 100)
xset <- retcor(xset, method = "obiwarp", profStep = 1)
xset <- group(xset, bw = 5, minfrac = 0.5, mzwid = 0.012, minsamp = 1, max = 100)
xset <- fillPeaks(xset)
XCMt <- data.frame(xset@groups)
xcmI <- groupval(xset, value = "maxo")
featureTable <- cbind(XCMt$mzmed, XCMt$rtmed, XCMt$rtmin, XCMt$rtmax, xcmI)
colnames(featureTable)[1:4] <- c("mz", "rt", "rtmin", "rtmax")
featureTable <- featureTable[order(featureTable[,1]),]
featureTable <- as.data.frame(featureTable)
#Output
write.csv(featureTable, file = "paramounter.csv", row.names = FALSE)
