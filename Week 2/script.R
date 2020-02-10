# TMB group week 2
setwd("/Users/jgao/Documents/Work/Teach/TMB\ study\ group/Week 2")
#-------------------------------------------------------------------------------
# problem 1 Setting bounds on parameters
library(TMB)
compile("p3.cpp",flags="-Wno-Wint-in-bool-context")
dyn.load(dynlib("p3"))
data <- list()
data$X <- c(128,158,92,122)
param <- list()
param$alpha <- c(0,0,0)
obj <- MakeADFun(data, param, DLL="p3", silent=TRUE)
opt <- nlminb(obj$par, obj$fn, obj$gr)  
summary(sdreport(obj))
#-------------------------------------------------------------------------------
# problem 2 map argument
library(TMB)
compile("map.cpp",flags="-Wno-Wint-in-bool-context") 
dyn.load(dynlib("map"))
param <- list(logAlpha=rep(0,nlevels(InsectSprays$spray)))
obj <- MakeADFun(InsectSprays, param, DLL="map") 
opt <- nlminb(obj$par, obj$fn, obj$gr)

map1=list(logAlpha=factor(c(1,1,2,3,4,1)))
obj1 <- MakeADFun(InsectSprays, param, map=map1, silent=TRUE, DLL="map") 
opt1 <- nlminb(obj1$par, obj1$fn, obj1$gr) 
1-pchisq(2*(opt1$obj-opt$obj),2)

map2=list(logAlpha=factor(c(NA,NA,2,3,4,NA)))
param2<-param
param2$logAlpha[c(1,2,6)]<-log(15)
obj2 <- MakeADFun(InsectSprays, param2, map=map2, silent=TRUE, DLL="map")
opt2 <- nlminb(obj2$par, obj2$fn, obj2$gr) 
1-pchisq(2*(opt2$obj-opt1$obj),1)

#-------------------------------------------------------------------------------
# problem 3 Beverton-Holt Curve fitting
library(TMB)

compile("bevholt.cpp",flags="-Wno-Wint-in-bool-context")
dyn.load(dynlib("bevholt"))

dat <- read.table("bevholt.dat", header=TRUE)
data <- list(SSB=dat$ssb,logR=dat$logR) 
parameters <- list(logA=0, logB=0)

obj <- MakeADFun(data,parameters,DLL="bevholt")
obj$env$beSilent() # silences console output
obj$fn()
obj$gr()

### --------------------------------------------------
### Exercise answers
opt <- nlminb(obj$par,obj$fn,obj$gr)
## The MLEs:
logA <- opt$par[1]
logB <- opt$par[2]
plot(dat$ssb,dat$logR)
ssb <- seq(0, 300000, len=1000)
ssb <- data$SSB
## Predicted logR:
logR <- logA+log(ssb)-log(1+exp(logB)*ssb)
lines(ssb, logR, lwd=3)

## Reoptimize model from a grid of starting points
x1 <- seq(-5, 10, len=100)
x2 <- seq(-15,0, len=100)
inits <- expand.grid(logA=x1,logB=x2)
nlls <- evals <- rep(NA, len=nrow(inits))
for(i in 1:nrow(inits)){
  opt.temp <- nlminb(inits[i,],obj$fn,obj$gr)
  evals[i] <- opt.temp$iterations
  nlls[i] <- opt.temp$objective
}
## These should all be the same:
summary(nlls)
summary(evals)

## Make a contour of the negative log-likelihood surface
o <- 6 # number of SEs away from MLE
logA.se <- sqrt(diag(sdreport(obj)$cov.fixed))[1]
logB.se <- sqrt(diag(sdreport(obj)$cov.fixed))[2]
x1 <- seq(logA-o*logA.se, logA+o*logA.se, len=100)
x2 <- seq(logB-o*logB.se, logB+o*logB.se, len=100)
z <- sapply(x1, function(logA) sapply(x2, function(logB)
  obj$fn(c(logA,logB))))
contour(x1,x2,z, nlevels=100)
points(logA, logB, pch=16, col=2)
covar <- sdreport(obj)$cov.fixed
## Add the confidence ellipse and standard errors
require(ellipse)
CI.region <- ellipse(x=covar, centre=c(logA, logB))
lines(CI.region, col=2)
lines(x=c(logA-1.96*logA.se, logA+1.96*logA.se), y=c(logB,logB), col=4)
lines(y=c(logB-1.96*logB.se, logB+1.96*logB.se), x=c(logA,logA), col=4)

## Simulation testing the BH model
Nsim <- 2000
mles <- ses <- coverage <- matrix(NA, nrow=Nsim, ncol=2)
logR.expected <- logA+log(data$SSB)-log(1+exp(logB)*data$SSB)
## Note: sigma=0.4 is fixed in the template and matches here
for(i in 1:Nsim){
  error <- rnorm(n=length(logR.expected), mean=0, sd=0.4)
  data$logR <- logR.expected+error
  obj <- MakeADFun(data,parameters,DLL="bevholt")
  obj$env$beSilent()
  opt <- nlminb(obj$par, obj$fn, obj$gr)
  mles[i,] <- opt$par
  ses[i,] <- sqrt(diag(sdreport(obj)$cov.fixed))
  ## Coverage is whether the confidence interval covers the true
  ## parameter.
  coverage[i,1] <-
    mles[i,1]-1.96*ses[i,1] < logA &
    logA < mles[i,1]+1.96*ses[i,1]
  coverage[i,2] <-
    mles[i,2]-1.96*ses[i,2] < logB &
    logB < mles[i,2]+1.96*ses[i,2]
}

par(mfrow=c(1,2))
breaks <- 20
hist(mles[,1], main='logA', breaks=breaks)
abline(v=logA, col=2)
hist(mles[,2], main='logB', breaks=breaks)
abline(v=logB, col=2)

## Check coverage
100*apply(coverage, 2, mean)

#-------------------------------------------------------------------------------
# problem 4 common Kalman filter

library(TMB)
compile("kf.cpp")
dyn.load("kf.so")
data <- list(y = scan("rw.dat"))
plot(seq(1, length(data$y), 1), data$y, pch=1, type='o')
parameters <- list(logSdRw=0, logSdObs=0, lam0=0)
obj <- MakeADFun(data, parameters, DLL="kf")
obj$fn()
obj$gr()
opt <- nlminb(obj$par, obj$fn, obj$gr)
rep <- obj$report()
save(rep, file = "kf.RData")

#png(filename="kf_fit.png")
plot(seq(1, length(data$y), 1), data$y, pch=1, 'o', xlab = 'Time', ylab = '')
lines(rep$lamSmooth, col='red')
lines(1.96*sqrt(rep$lamSmoothVar)+rep$lamSmooth, col = 'red', lty=2)
lines(-1.96*sqrt(rep$lamSmoothVar)+rep$lamSmooth, col = 'red', lty=2)
#dev.off()
#img = png::readPNG("kf_fit.png")
#grid::grid.raster(img)

# TMB+laplace
library(TMB)
#compile("rw.cpp", flags="-w")
#compile("rw.cpp", flags="-Wall")
compile("rw.cpp", flags="-Wno-Wint-in-bool-context")
dyn.load(dynlib("rw"))
data <- list(y=scan("rw.dat"))
parameters <- list(
  logSdRw=0,
  logSdObs=0,
  lam0=0,
  lam=rep(0,length(data$y))
  )

obj <- MakeADFun(data,parameters,random="lam",DLL="rw")

obj$fn()
obj$gr()
opt<-nlminb(obj$par,obj$fn,obj$gr)

sdr<-sdreport(obj)
pl <- as.list(sdr,"Est")
plsd <- as.list(sdr,"Std")
save(pl,plsd,file="rw.RData")

lines(1.96*plsd$lam+pl$lam, col = 'blue', lty=2)
lines(-1.96**plsd$lam+pl$lam, col = 'blue', lty=2)







