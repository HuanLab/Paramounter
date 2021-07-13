###############################################################
#This is the script to perform parameters estimation
#Jian Guo, Tao Huan 2021-07-28
#Copyright @ University of British Columbia
###############################################################

library(xcms)
library(MSnbase)
library(dplyr)
library(ggplot2)
library(gridExtra)

# User input the directory and software to optimize parameters for (XCMS, MSDIAL, or MZMINE2)
directory <- "F:/Jian_Guo/Parameter_optimize_20201105/ThermoQExactiveMTBLS201_HILIC-"
################################################################################################
setwd(directory)
filename <- list.files(pattern = ".mzXML")
start_time <- Sys.time()
mzDiff <- c()
ppm <- c()
mzDiff2D <- as.data.frame(matrix(ncol = 3, nrow = 1))
colnames(mzDiff2D) <- c("mz", "rt", "mzdiff")
ppm2D <- as.data.frame(matrix(ncol = 3, nrow = 1))
colnames(ppm2D) <- c("mz", "rt", "ppm")

for (q in 1:(length(filename))){
  # Parameter setting
  ms1data <- readMSData(files = filename[q], mode = "onDisk", msLevel. = 1)
  mzRange <- c(min(unlist(mz(ms1data))), max(unlist(mz(ms1data))))
  ROI <- seq(mzRange[1], mzRange[2], 0.05)
  mzData <- mz(ms1data)
  intData <- intensity(ms1data)
  rtime <- rtime(ms1data)
  ppm2Ddist <- as.data.frame(matrix(ncol = 3, nrow = 1))
  colnames(ppm2Ddist) <- c("mz", "rt", "ppm")
  mzdiff2Ddist <- as.data.frame(matrix(ncol = 3, nrow = 1))
  colnames(mzdiff2Ddist) <- c("mz", "rt", "mzdiff")
  
  # ROI detection and universal parameter estimation 
  for(i in 1:(length(ROI) - 1)) {
    # Obtain data lists in each m/z bin 
    currmzRange <- c(ROI[i], ROI[i+1])
    tmpMZdata <- mzData
    tmpINTdata <- intData
    for(j in 1:length(mzData)){
      index <- which(tmpMZdata[[j]] >= currmzRange[1] & tmpMZdata[[j]] < currmzRange[2])
      tmpMZdata[[j]] <- tmpMZdata[[j]][index]
      tmpINTdata[[j]] <- tmpINTdata[[j]][index]
    }
    # Extract the intensity vectors from each m/z bin 
    eicINT <- c()
    eicRT <- c()
    for(k in 1:length(mzData)){
      if(length(tmpINTdata[[k]]) > 0){
        eicINT[k] <- mean(tmpINTdata[[k]])
      }else{
        eicINT[k] <- 0
      }
      eicRT[k] <- rtime[k]
    }
    if(sum(eicINT != 0) == 0) next()
    # Sort the intensity vectors from each m/z bin, estimate the noise cut off and average
    
    eicNon0 <- sort(eicINT[eicINT > 0])
    if(length(eicNon0) > 10){
      for(x in seq(10,length(eicNon0), 10)){
        sd <- sd(eicNon0[1:x])
        blk <- sum(eicNon0[1:x])/x
        thres <- blk + 3*sd
        if(x+1 <= length(eicNon0)){
          if(eicNon0[x+1] >= thres) break()
        }
      }
      cutOFF <- eicNon0[x]
    }else{
      cutOFF <- max(eicNon0)
    }

    aboveTHindex <- which(eicINT > cutOFF)
    if(length(aboveTHindex) == 0) next()
    candidateSegInd <- split(aboveTHindex, cumsum(c(1, diff(aboveTHindex) != 1)))
    peakInd <- c()
    for(x in 1:length(candidateSegInd)){
      peakInd[x] <- which(eicINT[candidateSegInd[[x]]] == max(eicINT[candidateSegInd[[x]]]))[1] + min(candidateSegInd[[x]]) - 1
    }
    refMZvec <- c()
    for(y in 1:length(peakInd)){
      highestINT <- which(tmpINTdata[[peakInd[y]]] == max(tmpINTdata[[peakInd[y]]]))[1]
      refMZvec[y] <- tmpMZdata[[peakInd[y]]][highestINT]
    }
    
    # Estimate the universal parameters (mass tolerance, peak height, and peak width) for each m/z bin
    ppmDiff <- c()
    for(z in 1:length(peakInd)){
      currPeakInd <- peakInd[z]
      currRefMz <- refMZvec[z]
      currSamePeakMass <- c()
      currSamePeakMass <- c(currSamePeakMass, currRefMz)
      leftInd <- currPeakInd-1
      rightInd <- currPeakInd+1
      if(leftInd > 0){
        while (length(tmpMZdata[[leftInd]]) > 0 & mean(tmpINTdata[[leftInd]]) >= cutOFF) {
          if (length(tmpMZdata[[leftInd]]) == 1){
            currSamePeakMass <- c(currSamePeakMass, tmpMZdata[[leftInd]])
          } else {
            abvector <- abs(tmpMZdata[[leftInd]] - currRefMz)
            NearInd <- which(abvector == min(abvector))[1]
            currSamePeakMass <- c(currSamePeakMass, tmpMZdata[[leftInd]][NearInd])
          }
          leftInd <- leftInd-1
          if(leftInd <= 0) break()
        }
      }
      if(rightInd <= length(tmpMZdata)){
        while (length(tmpMZdata[[rightInd]]) > 0 & mean(tmpINTdata[[rightInd]]) >= cutOFF) {
          if (length(tmpMZdata[[rightInd]]) == 1){
            currSamePeakMass <- c(currSamePeakMass, tmpMZdata[[rightInd]])
          } else {
            abvector <- abs(tmpMZdata[[rightInd]] - currRefMz)
            NearInd <- which(abvector == min(abvector))[1]
            currSamePeakMass <- c(currSamePeakMass, tmpMZdata[[rightInd]][NearInd])
          }
          rightInd <- rightInd+1
          if(rightInd > length(tmpMZdata)) break()
        }
      }
      
      if(length(currSamePeakMass) > 1){
        ppmDiff[z] <- (sd(currSamePeakMass))/currRefMz * 1e6
        ppm2Ddist <- rbind(ppm2Ddist, c(currRefMz, rtime[[peakInd[z]]], ppmDiff[z]))
      }
    }
  }

  ppm2D <- rbind(ppm2D, ppm2Ddist)
}

ppm2D <- ppm2D[complete.cases(ppm2D),]
ppm2D <- ppm2D[order(ppm2D[,3]),]
ppm2D <- ppm2D[1:round(nrow(ppm2D)*0.97),]
plot(ppm2D$mz, ppm2D$ppm, ylab = "ppm", xlab = "m/z", pch=1, cex.main=4, cex.lab=1.7, cex.axis=2)
message("Please find the cutoff line in the generated ppm distribution, and run Paramounter part 2 using the ppm cutoff")

