# Calculating the focal weight matrices.
weighting.window <- function(window.type="identity",w.halfsize=25) {
  window.type <- tolower(window.type)
  if (!any(window.type==c("identity","inverse","squared","linear"))) stop("Wrong weighting window type")
  w.identity <- matrix(1/(2*w.halfsize+1)^2,2*w.halfsize+1,2*w.halfsize+1)
  w <- outer((-w.halfsize:w.halfsize)^2,(-w.halfsize:w.halfsize)^2,"+")
  w.identity[w>w.halfsize^2] <- 0
  w.identity[w.halfsize+1,w.halfsize+1] <- 0
  w.identity <- w.identity/sum(w.identity)
  if (window.type=="identity") {
    return(w.identity)
  } else switch(window.type,
                inverse={
                  w.inverse <- matrix(0,2*w.halfsize+1,2*w.halfsize+1)
                  w.inverse[w.identity!=0] <- 1/sqrt(w[w.identity!=0])
                  return(w.inverse)
                },
                squared={
                  w.squared <- matrix(0,2*w.halfsize+1,2*w.halfsize+1)
                  w.squared[w.identity!=0] <- 1/w[w.identity!=0]
                  return(w.squared)
                },
                linear={
                  w.linear <- matrix(0,2*w.halfsize+1,2*w.halfsize+1)
                  w.linear[w.identity!=0] <- max(sqrt(w))-max(sqrt(w))/w.halfsize*sqrt(w[w.identity!=0])
                  return(w.linear)
                })
}

# It computes the Kappa index for two arrays.
kappa.index <- function(a,b) {
  nab <- sum(!is.na(b))
  a0 <- !a
  b0 <- !b
  p.a <- (sum(a0 & b0) + sum(a & b))/nab
  p.e <- sum(a0)/nab*sum(b0)/nab+sum(a)/nab*sum(b)/nab
  return((p.a-p.e)/(1-p.e))
}

# Model fit.
fit.model <- function(data.train,data.test,type="randomforest",n.threads=4,verbose=T) {
# Regressions.
  if (type=="randomforest") {
    reg <- ranger(dependent.variable.name="y",data=data.train,verbose=T,importance="permutation",
                                    write.forest=T,num.threads=n.threads,num.trees=500,classification=F)
    pred <- predict(reg,data=data.test,verbose=F)$predictions
  } else if (type=="binomial") {
    reg <- glm(y~.,family=binomial(link="logit"),data=data.train,x=F,y=F)
  } else if (type=="binomial.interactions") {
    reg <- glm(y~.^2,family=binomial(link="logit"),data=data.train,x=F,y=F)
  } else if (type=="binomial.parabollic.elevation") {
    reg <- glm(y~.+I(elevation^2),family=binomial(link="logit"),data=data.train,x=F,y=F)
  } else stop("Wrong fit option")

  if (type!="randomforest") {
    if (verbose) print(summary(reg))
    sum.model <- summary(reg)$coefficients
    sum.beta <- beta(reg)$coefficients
    reg <- stripGlmLR(reg)
    pred <- predict(reg,newdata=data.test,type="response")
  } else {
    if (verbose) print(reg)
    sum.model <- c()
    sum.beta <- c()
  }
  roc.pred <- prediction(pred,data.test$y)
  
  # Youden's J.
  youden.perf <- performance(roc.pred,"sens","spec")
  x.youden <- youden.perf@alpha.values[[1]]
  y.youden <- youden.perf@y.values[[1]]+youden.perf@x.values[[1]]-1

  # Kappa.
  x.kappa <- seq(0.001,.999,length=1000)
  y.kappa <- sapply(x.kappa,function(x) kappa.index(pred>=x,data.test$y))

  return(list(reg=reg,
              sum.model=sum.model,
              sum.beta=sum.beta,
              auc=performance(roc.pred,"auc")@y.values,
              y.youden=max(y.youden),
              youden.cutoff=x.youden[which.max(y.youden)],
              y.kappa=max(y.kappa),
              kappa.cutoff=x.kappa[which.max(y.kappa)]))
}


# Write binary image onto disk. Cutoff is chosen to maximize Kappa or Youden' J.
write.cutoff.image <- function(name.output,regres,image.to.match,j.index,type="youden",
                               extension=ext.image,projection=proj.image,dataset=dataset) {
  cat(paste("   Image ",name.output,"\n",sep=""))
  if (type=="youden") {
    cat(paste("      Max. Youden J = ",max(regres$youden$youden),"\n",sep=""))
    cutoff <- regres$youden$cutoff[which.max(regres$youden$youden)]
    cat(paste("      Cutoff Youden J = ",cutoff,"\n",sep=""))
  } else {
    cat(paste("      Max. Kappa = ",max(regres$kappa$kappa),"\n",sep=""))
    cutoff <- regres$kappa$cutoff[which.max(regres$kappa$kappa)]
    cat(paste("      Cutoff Kappa = ",cutoff,"\n",sep=""))
  }
  if (class(regres$reg)[1]=="speedglm" | class(regres$reg)[1]=="glm") {
    pred <- predict(regres$reg,newdata=dataset,type="response")
  } else {
    pred <- predict(regres$reg,data=dataset,type="response")$predictions
  }
  image.to.match[j.index] <- (pred>cutoff)*1
  image.to.match <- raster(image.to.match)
  projection(image.to.match) <- proj.image
  extent(image.to.match) <- ext.image
  # setMinMax(image.to.match) <- c(0,1)
  writeRaster(x=image.to.match,format="ascii",overwrite=T,filename=name.output)
}

focal.parallel <- function(a,w.halfsize,parallel.focal=T) {
  wlabel <- c("identity","linear","inverse","squared")
  if (parallel.focal) {
    cl <- makeCluster(4)
    registerDoParallel(cl)
      out <- foreach (i=1:4,.packages=c("raster"),.export="weighting.window") %dopar% 
        focal(a,weighting.window(wlabel[i],w.halfsize),na.rm=T,pad=T)
    stopCluster(cl)
  } else out <- lapply(1:4,function(i) as.matrix(focal(a,weighting.window(wlabel[i],w.halfsize),na.rm=T,pad=T)))
  r <- stack(out[[1]])
  for (i in 2:4) r <- stack(r,out[[i]])
  names(r) <- wlabel
  return(r)
}

get.predictor.rasters <- function(neighbors=8,w.halfsize=25) {
  # Aspect and elevation datasets.
  
  cat("Reading elevation, aspect and mask data...\n")
  elevation <- raster(readGDAL("..\\..\\Infestation data - increment\\dem\\w001001.adf",silent=T))
  proj.images <- proj4string(elevation)
  ruggedness <- terrain(elevation,opt="TPI",neighbors=neighbors)
  slope <- terrain(elevation,opt="slope",neighbors=neighbors)
  aspect <- raster(readGDAL("..\\..\\Infestation data - increment\\aspect\\w001001.adf",silent=T))
  aspect.sin <- sin(2*pi/360*aspect)
  aspect.cos <- cos(2*pi/360*aspect)
  lulc.mask <- raster(readGDAL("..\\..\\Infestation data - increment\\Masks (urban, crops, lakes...)\\mask_water_lc_xposedland\\mask_wr_lc_xp.asc",silent=T))

  cat("Clipping rasters and converting to matrices...\n")
  a <- raster(readGDAL("..\\..\\Infestation data - cumulative\\cumkill1999.asc",silent=T))
  proj4string(a) <- proj.images
  elevation <- crop(elevation,a)
  ruggedness <- crop(ruggedness,a)
  aspect.sin <- crop(aspect.sin,a)
  aspect.cos <- crop(aspect.cos,a)
  lulc.mask <- crop(lulc.mask,a)
  slope <- crop(slope,a)

  x <- stack(elevation,ruggedness,aspect.sin,aspect.cos,slope,lulc.mask)
  names(x) <- c("elevation","ruggedness","aspect.sin","aspect.cos","slope","lulc.mask")
  return(x)
}

read.cumkill <- function(i,time.step,th,silent=T) {
  f <- function(aa) {
    aa[is.na(aa)] <- 0
    return(aa)
  }
  x <- "..\\..\\Infestation data - cumulative\\cumkill"
  a <- f(raster(readGDAL(paste(x,i-1,".asc",sep=""),silent=silent))>th)
  a <- stack(a,f(raster(readGDAL(paste(x,i,".asc",sep=""),silent=silent))>th))
  a <- stack(a,f(raster(readGDAL(paste(x,i+time.step,".asc",sep=""),silent=silent))>th))
  names(a) <- c("i-1","i","i+time.step")
  return(a)
}


stripGlmLR = function(cm) {
  cm$y = c()
  cm$model = c()
  
  cm$residuals = c()
  cm$fitted.values = c()
  cm$effects = c()
  cm$qr$qr = c()  
  cm$linear.predictors = c()
  cm$weights = c()
  cm$prior.weights = c()
  cm$data = c()
  
  
  cm$family$variance = c()
  cm$family$dev.resids = c()
  cm$family$aic = c()
  cm$family$validmu = c()
  cm$family$simulate = c()
  attr(cm$terms,".Environment") = c()
  attr(cm$formula,".Environment") = c()
  
  cm
}


compute.projection <- function(a,reg.object,b=NULL,w.halfsize=w.halfsize) {

  # Raster stack with predictor variables. Make lulc a separate raster.
  cat(paste("Get predictors...\n"))
  r <- get.predictor.rasters()
  lulc <- r[["lulc.mask"]]
  r <- dropLayer(r,which("lulc.mask"==names(r)))

  
  # Now we compute neighborhood with "focal" function.
  cat(paste("Calculating neighborhood with focal function...\n"))
  proj4string(a) <- proj4string(r)
  r <- stack(r,focal.parallel(a,w.halfsize,T))
  if (!is.null(b)) r <- stack(r,focal.parallel(b,w.halfsize,T))

  # Building the mask.
  cat("Creating the mask...\n")
  mask <- if (is.null(b)) a==0 else b==0
  mask <- mask & r[[5]]>0
  mask <- mask & lulc==1
  mask <- mask & calc(r,fun=function(x) all(!is.na(x)))

  # Create the dataset for the fit.
  dat <- as.data.frame(r)
  dat <- dat[which(as.data.frame(mask)==1),]
  pr <- predict(reg.object,newdata=dat,type="response")
  if (!is.null(b)) a <- b
  a[mask==1] <- pr
  
  return(a)
}



