rm(list=ls(all=TRUE))
assign("last.warning", NULL, envir = baseenv())

library(raster)
library(rgdal)
library(biomod2)
library(parallel)
library(foreach)
library(doParallel)
library(speedglm)
library(gnm)
library(speedglm)
library(ROCR)
library(ranger)

setwd("C:\\Roberto\\Saeed Harati\\Analysis infestation data - cumulative\\New analysis Roberto")

###################################################
# Some functions are found here.
source("MPB functions V14.R")

# Initial parameters.

threshold.initial <- pnorm(-1)
year.calculation <- 2000
i.start <- 2002 # Starting year
one.year <- T
years <- ifelse(one.year,"one year only","two years")
time.step <- 3
w.halfsize <- switch(time.step,
                     "1"=50,
                     "2"=100,
                     "3"=150)
load(paste("Regression V7 - ",year.calculation," - time step ",time.step," - w_halfsize ",w.halfsize," - ",years,sep=""))

# Raster stack with predictor variables. Make lulc a separate raster.
r <- get.predictor.rasters()
lulc <- r[["lulc.mask"]]
r <- dropLayer(r,which("lulc.mask"==names(r)))

cat(paste("Prediction ",i.start+time.step,"...\n",sep=""))
a <- read.cumkill(i.start,time.step,threshold.initial)
a <- dropLayer(a,3)
if (one.year) a <- dropLayer(a,1)
proj4string(a) <- proj4string(r)

# Now we compute neighborhood with "focal" function.
cat(paste("Calculating neighborhood with focal function...\n"))
rr <- stack(r,focal.parallel(raster(a,1),w.halfsize,T))
if (!one.year) rr <- stack(rr,focal.parallel(raster(a,2),w.halfsize,T))

# Building a mask that must fulfill four conditions:
# 1.- only locations that are not infested in t=1 are valid.
# 2.- only locations that are close to already infested locations are valid.
# 3.- only locations with appropriate lulc classes can be infested.
# 4.- only locations where other variables do not have NA are valid.
cat("   Creating the mask...\n")
mask <- if (one.year) a==0 else raster(a,2)==0
mask <- mask & rr[[5]]>0
mask <- mask & lulc==1
mask <- mask & calc(rr,fun=function(x) all(!is.na(x)))

# Create the dataset for the fit.
cat(paste("   Selecting training and test datasets...\n"))
dat <- as.data.frame(rr)
dat <- dat[which(as.data.frame(mask)==1),]

for (j in 1:4) {
  file.name <- switch(j,
                      "1"="Binomial",
                      "2"="Binomial parabollic elevation",
                      "3"="Binomial all interactions",
                      "4"="Random forest")
  file.name <- paste(file.name," - year calculation ",year.calculation,sep="")
  cat(paste("Predictions ",file.name,"...\n",sep=""))
  reg <- switch(j,
               "1"=binomial,
               "2"=binomial.parabollic.elevation,
               "3"=binomial.interactions,
               "4"=rf)
  pr <- if(j==4) predict(reg$reg,data=dat) else predict(reg$reg,newdata=dat,type="response")
  b <- if (one.year) a else raster(a,2)
  b[Which(mask==1)] <- if (j==4) pr$predictions else pr
  file.name <- paste(".\\Maps\\",file.name,", start ",i.start," - ",
                     "final ",i.start+time.step," - ",
                     years,sep="")
  youden.cutoff <- reg$youden.cutoff
  kappa.cutoff <- reg$kappa.cutoff
  writeRaster(b,
              filename=paste(file.name," - probability map year ",".tif",sep=""),
              format="GTiff",overwrite=T)
  writeRaster((b>youden.cutoff)*1,
              filename=paste(file.name," - TSS cutoff - ",".tif",sep=""),
              format="GTiff",overwrite=T)
  writeRaster((b>kappa.cutoff)*1,
              filename=paste(file.name," - Kappa cutoff - ",".tif",sep=""),
              format="GTiff",overwrite=T)
}

quit(save="no")
