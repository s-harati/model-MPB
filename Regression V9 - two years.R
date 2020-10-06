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
library(reghelper)

setwd("C:\\Roberto\\Saeed Harati\\Analysis infestation data - cumulative\\New analysis Roberto")

###################################################
# Some functions are found here.
source("MPB functions V14.R")

# Initial parameters.
pi <- 3.141592653
threshold.initial <- pnorm(-1)
time.step <- 1 # this simulation
w.halfsize <- switch(time.step,
                     "1"=50,
                     "2"=100,
                     "3"=150)
###################################################

# Raster stack with predictor variables. Make lulc a separate raster.
r <- get.predictor.rasters()
lulc <- r[["lulc.mask"]]
r <- dropLayer(r,which("lulc.mask"==names(r)))
# 
# # To store the results.
# binomial <- list()
# binomial.interactions <- list()
# binomial.parabollic.elevation <- list()
# rf <- list()

# for (i in c(2000,2004,2008,2011)) {
  for (i in 2012) {
  cat(paste("Reading ",i-1,", ",i," and ",i+time.step," rasters...\n",sep=""))
  a <- read.cumkill(i,time.step,threshold.initial)
  proj4string(a) <- proj4string(r)

  # Now we compute neighborhood with "focal" function.
  cat(paste("   Calculating neighborhood with focal function...\n"))
  cat(paste("      Neighborhood at t=-1\n"))
  rr <- stack(r,focal.parallel(raster(a,1),w.halfsize,T))
  cat(paste("      Neighborhood at t=0\n"))
  rr <- stack(rr,focal.parallel(raster(a,2),w.halfsize,T))
  rr <- stack(rr,raster(a,3))
  names(rr)[dim(rr)[3]] <- "y"

  # Building a mask that must fulfill four conditions:
  # 1.- only locations that are not infested in t=1 are valid.
  # 2.- only locations that are close to already infested locations are valid.
  # 3.- only locations with appropriate lulc classes can be infested.
  # 4.- only locations where other variables do not have NA are valid.
  cat("   Creating the mask...\n")
  mask <- raster(a,2)==0
  mask <- mask & rr[["identity"]]>0
  mask <- mask & lulc==1
  mask <- mask & calc(rr,fun=function(x) all(!is.na(x)))

# Create the dataset for the fit.
  cat(paste("   Selecting training and test datasets...\n"))
  dat <- as.data.frame(rr)
  dat <- dat[which(as.data.frame(mask)==1),]
  
  data.train <- dat
  data.test <- dat

  cat(paste("   Binomial fit...\n"))
  binomial <- fit.model(data.train,data.test,type="binomial")
  cat(paste("   Binomial fit with parabollic elevation...\n"))
  binomial.parabollic.elevation <- fit.model(data.train,data.test,type="binomial.parabollic.elevation")
  cat(paste("   Binomial fit with interactions...\n"))
  binomial.interactions <- fit.model(data.train,data.test,type="binomial.interactions")
  cat(paste("   Random forest...\n"))
  rf <- fit.model(data.train,data.test,type="randomforest",n.threads=7)

 # Saving. # no need to do the saving of regression results. we use them presently.
  cat("   Saving...\n")
  save(time.step,w.halfsize,binomial,binomial.parabollic.elevation,binomial.interactions,rf,compress=T,
       file=paste("Regression V9 - ",i," - time step ",time.step," - w_halfsize ",w.halfsize," - two years",sep=""))
  }
