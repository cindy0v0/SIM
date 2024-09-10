library(gtools)
library(rlang)
library(msdata)
library(MSnbase)
library(xcms)


#mzfiles <- list.files(path = ".", pattern="CD-6KUCT.mzML")
#mzfiles <- list.files(path = ".", pattern= "CD-9OS5Y.mzML")
# object <- readMSData(mzfiles, msLevel = 1)


# SMOOTHING 
setwd('C:/Users/cindy/Downloads/lab 2021-2022/mzML-sim/data/0415_5_intensity1-2_linear')
mzfiles <- list.files(path = ".", pattern= "1.mzML")
file_1 <- readMSData(mzfiles, msLevel = 1, mode = "onDisk")
rawInt <- intensity(file_1)
rawMz <- mz(file_1) # rt = 5-20 min -> 300-1200 sec
rtime <- rtime(file_1)
peak <- numeric(length(rawMz))
for (i in 1:length(rawMz)){
  peak[i] <- 0
  temp <- c()
  for (j in 1:length(rawMz[[i]])){
    if (rawMz[[i]][[j]] >= 283.082 && rawMz[[i]][[j]] <= 283.086){ # 283.0403 - 283.6205
      temp <- append(temp, rawInt[[i]][[j]])
    }
  }
  if (length(temp) > 0) {
    peak[i] <- mean(temp)
  }
}
plot(rtime, peak, type = "l", xlim=c(350, 450))
#smoothedInt <- peak_smooth(rawInt$F1.S001, 2)
smoothed_peak <- peak_smooth(peak, 4)
tiff("rt_2.tiff", units="in", width=5, height=5, res=300)
plot(rtime, smoothed_peak, type = "l", xlim=c(350, 450), ylab="intensity", xlab="retention time (sec)")
dev.off()
# par(mar=c(1,1,1,1))
#mz = c(280, 280.18)
#rt = c(475, 570)
ggplot(aes(rtime,smoothed_peak)) + geom_point() + 
  geom_label(data = . %>% group_by(z) %>% filter(y == max(y)), aes(label = sprintf('%0.2f', y)), hjust = -0.5) +
  facet_wrap(~z)


# PEAK PICKING
cwp <- CentWaveParam(noise = 100, peakwidth = c(10, 60), prefilter = c(3, 100), 
                     snthresh = 6, integrate = 2, ppm = 20, mzdiff = 0.001)
xdata <- findChromPeaks(object.ondisk, param = cwp)
regions <- data.frame(chromPeaks(xdata))
# regions <- regions[order(regions$mzmin),]
write.csv(regions, "fold2_features.csv")

# CHROMATOGRAM VIEWING
chrs <- chromatogram(xdata, rt = c(300, 500), mz = c(283.082712515440, 283.08505722), ) # 283.0403 - 383.6205  +- 0.01
chrs
plot(chrs)
bis <- chromatogram(object.ondisk, aggregationFun = "max")
bpis <- chromatogram(object.ondisk, rt = c(480, 600), mz = c(249, 251), aggregationFun = "max")
bpis <- chromatogram(subset, aggregationFun = "max")
plot(bpis)


# ALIGNMENT
dir <- '20_gauss_(1-2)_pca_0414'
setwd(paste('C:/Users/cindy/Downloads/lab 2021-2022/mzML-sim/data/', dir, sep = ""))
mzfiles <- list.files(path = ".", pattern=".ML")
data0 <- readMSData(mzfiles, msLevel = 1, mode = "onDisk")
cwp <- CentWaveParam(noise = 100, peakwidth = c(10, 60), prefilter = c(3, 100), 
                     snthresh = 6, integrate = 2, ppm = 20, mzdiff = 0.001)
xdata <- findChromPeaks(data0, param = cwp)
# xdata_filtered <- filterMsLevel(xdata, msLevel = 1)
xset <- as(xdata, 'xcmsSet')

xset@peaks <- xset@peaks[order(xset@peaks[,11]),] # order by samples
table <- as.data.frame(xset@peaks)
for(n in (1:length(mzfiles))){
  sampleOutput <- table[table$sample == n, ]
  sampleOutput <- sampleOutput[order(sampleOutput[,1]),]
  colnames(sampleOutput)[9] <- "intMax"
  row.names(sampleOutput) <- 1:nrow(sampleOutput)
}
xset <- group(xset, bw = 10, minfrac = 0.5, mzwid = 0.006, minsamp = 1, max = 100)
xset <- retcor(xset, method = "obiwarp", profStep = 1)
xset <- group(xset, bw = 10, minfrac = 0.5, mzwid = 0.006, minsamp = 1, max = 100)
xset <- fillPeaks(xset)
XCMt <- data.frame(xset@groups)
xcmI <- groupval(xset, value = "into")
featureTable <- cbind(XCMt$mzmed, XCMt$rtmed, XCMt$rtmin, XCMt$rtmin, xcmI)
featureTable <- featureTable[, -c(1:4)]
# colnames(featureTable)[0] <- "features"
# featureTable <- featureTable[order(featureTable[,1]),]
featureTable <- as.data.frame(featureTable)
featureTable <- featureTable[ , mixedsort(names(featureTable))]
#featureTable <- featureTable[,c(ncol(featureTable),1:(ncol(featureTable)-1))]
colnames(featureTable)[1:dim(featureTable)[2]] <- 1:dim(featureTable)[2]-1
featureTable <- rbind(sample.int(2, dim(featureTable)[2], replace=TRUE), featureTable)
rownames(featureTable)[1] <- "groups"

write.csv(featureTable, file = paste("XCMS_ALIGNMENT_", "guass", ".csv", sep = ""))

# GET DATA FILES
mzfiles <- list.files(path = ".", pattern="CD-6KUCT_0.mzML")
data2 <- readMSData(mzfiles, msLevel = 1, mode = "onDisk")
xdata <- findChromPeaks(data2, param = cwp)
regions <- data.frame(chromPeaks(xdata))
regions <- regions[order(regions$mzmin),]
write.csv(regions, "xcms_result_2_CD6.csv")

mzfiles <- list.files(path = ".", pattern="*.mzML")
data.rand <- readMSData(mzfiles, msLevel = 1, mode = "onDisk")
xdata <- findChromPeaks(data.rand, param = cwp)
regions.rand <- chromPeaks(xdata)
write.csv(regions.rand, "xcms_result_randshift.csv")

# GET RAW DATA
mzfiles <- list.files(path = ".", pattern= "2fold_1.mzXML")
object.ondisk <- readMSData(mzfiles, msLevel = 1, mode = "onDisk")
intensity <- intensity(object.ondisk)
intensity <- plyr::ldply(intensity, rbind)
intensity.df <- data.frame(intensity)
drops <- c(".id")
intensity.df <- intensity.df[ , !(names(intensity.df) %in% drops)]
write.csv(intensity.df, "2fold_1_intensity.csv")

mz <- mz(object.ondisk)
mz <- plyr::ldply(mz, rbind)
mz.df <- data.frame(mz)
mz.df <- mz.df[ , !(names(mz.df) %in% drops)]
write.csv(mz.df, "2fold_1_mz.csv")

rt <- rtime(object.ondisk)
rt.df <- data.frame(rt)
rt.df <- rt.df[ , !(names(rt.df) %in% drops)]
write.csv(rt.df, "2fold_1_rt.csv")



# get raw data (version 1)

mzfiles <- list.files(path = ".", pattern="*.mzML")
ms_fl <- openMSfile(mzfiles, backend = "pwiz")
pks <- spectra(ms_fl)
intensity <- pks
for (v in 1:length(pks)){
  intensity[[v]] <- intensity[[v]][,2]
}
intensity <- plyr::ldply(intensity, rbind)
intensity.df <- data.frame(intensity)
write.csv(intensity.df, "2fold_1_intensity.csv")

mz <- pks
for (v in 1:length(pks)){
  mz[[v]] <- mz[[v]][,1]
}
mz <- plyr::ldply(mz, rbind)
mz.df <- data.frame(mz)
write.csv(mz.df, "2fold_1_mz.csv")

hdr <- header(ms_fl)
rt <- hdr$retentionTime
write.csv(rt, "2fold_1_rt.csv")

data.rand <- readMSData(mzfiles, msLevel = 1, mode = "onDisk")
xdata <- findChromPeaks(data.rand, param = cwp)
regions <- chromPeaks(xdata)
write.csv(regions, "2fold_1_features.csv")

# SMOOTHING FUNCTION
peak_smooth <- function(x,level=smooth){
  n <- level
  if(length(x) < 2*n){
    return(x)
  } else if(length(unique(x))==1){
    return(x)
  } else{
    y <- vector(length=length(x))
    for(i in 1:n){
      y[i] <- sum(c((n-i+2):(n+1),n:1)*x[1:(i+n)])/sum(c((n-i+2):(n+1),n:1))
    }
    for(i in (n+1):(length(y)-n)){
      y[i] <-  sum(c(1:(n+1),n:1)*x[(i-n):(i+n)])/sum(c(1:(n+1),n:1))
    }
    for(i in (length(y)-n+1):length(y)){
      y[i] <- sum(c(1:n,(n+1):(n+i-length(x)+1))*x[(i-n):length(x)])/sum(c(1:n,(n+1):(n+i-length(x)+1)))
    }
    return(y)
  }
}
