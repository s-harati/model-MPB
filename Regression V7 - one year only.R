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
w.halfsize <- 150
time.step <- 3 # this simulation
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
  for (i in 2000) {
  cat(paste("Reading ",i-1,", ",i," and ",i+time.step," rasters...\n",sep=""))
  a <- read.cumkill(i,time.step,threshold.initial)
  a <- dropLayer(a,1)
  proj4string(a) <- proj4string(r)

  # Now we compute neighborhood with "focal" function.
  cat(paste("   Calculating neighborhood with focal function...\n"))
  rr <- stack(r,focal.parallel(raster(a,1),w.halfsize,T))
  rr <- stack(rr,raster(a,2))
  names(rr)[dim(rr)[3]] <- "y"

  # Building a mask that must fulfill four conditions:
  # 1.- only locations that are not infested in t=1 are valid.
  # 2.- only locations that are close to already infested locations are valid.
  # 3.- only locations with appropriate lulc classes can be infested.
  # 4.- only locations where other variables do not have NA are valid.
  cat("   Creating the mask...\n")
  mask <- raster(a,1)==0
  mask <- mask & rr[[5]]>0
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
  save(time.step,binomial,binomial.parabollic.elevation,binomial.interactions,rf,compress=T,
       file=paste("Regression V7 - ",i," - time step ",time.step," - w_halfsize ",w.halfsize," - one year only",sep=""))
}

#############################

# browser()
# 
# library(reghelper)
# load("Regression V5")
# 
# j <- 0
# for (i in c(2001,2005,2009,2014)) {
#   j <- j+1
#   write.table(coef(beta(binomial[[j]]$reg)),
#               file=paste("Binomial beta coefficients, year ",i,".csv",sep=""),col.names=NA,row.names=T,sep=";",dec=".")
#   write.table(coef(beta(binomial.parabollic.elevation[[j]]$reg)),
#               file=paste("Binomial beta coefficients, parabollic elevation, year ",i,".csv",sep=""),col.names=NA,row.names=T,sep=";",dec=".")
#   write.table(coef(beta(binomial.interactions[[j]]$reg)),
#               file=paste("Binomial beta coefficients, all pairwise interactions, year ",i,".csv",sep=""),col.names=NA,row.names=T,sep=";",dec=".")
# }
# 
# im <- sapply(1:length(rf),function(i) importance(rf[[i]]$reg)/max(importance(rf[[i]]$reg))*100)
# colnames(im) <- c(2001,2005,2009,2014)
# write.table(im,file="Random forest importance V4.csv",col.names=NA,row.names=T,sep=";",dec=".")
# 
#   
# #############################
# 
# r <- get.predictor.rasters()
# lulc <- r[["lulc.mask"]]
# r <- dropLayer(r,which("lulc.mask"==names(r)))
# 
# 
# i <- 2000
# time.step <- 1
# load(paste("Regression V5 - ",i," - time step ",time.step,sep=""))
# 
# 
# 
# #########################
# library(pdp)
# load("Regression V5 - 2000 - time step 1 - w_halfsize 50")
# pdp <- partial(rf$reg,pred.var=c("elevation"))
