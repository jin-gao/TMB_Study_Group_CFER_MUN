library(TMB)

## Exercise 1: using TMB derivative
compile("polynomial.cpp")
dyn.load(dynlib("polynomial"))
obj <- MakeADFun(data=list(), parameters=list(x=0), DLL='polynomial')
obj$fn(5)
obj$gr(5)

## Solution
x <- seq(-3,2, len=200)
y <- g <- rep(NA, 200)
for(i in 1:length(x)) {
  y[i] <- obj$fn(x[i])
  g[i] <- obj$gr(x[i])
}
png("polynomial.png", width=5.5, height=2.5, units='in', res=500)
par(mfrow=c(1,2), tck=-.02, mgp=c(1.5,.25,0), mar=c(3,3,1,.1))
plot(x,y, type='l', lwd=2, main='Height', ylab='fn')
plot(x,g, type='l', lwd=2, main='Derivative', ylab='gr')
abline(h=0)
dev.off()

## Exercise 2 solution: Poisson likelihood
k <- 4
loglike <- function(lambda) k*log(lambda)-lambda-log(factorial(k))
loglike(5.5)
dpois(x=4, lambda=5.5, log=TRUE)
## lgamma is better because it is more stable numerically
lgamma(k+1)
log(factorial(k))
ll <- seq(0.1, 15, len=1000)
plot(ll, -loglike(ll), type='l', lwd=2)


## Excercise 3 solution: Fitting linear model: y=a+b*x
x <- c(1.87, 1.96, 1.39, 2.24, 2.33, 2.24, 2.67, 2.47, 1.35, 2.00)
y <- c(2.47, 2.42, 2.2, 2.72, 2.65, 2.5, 2.85, 2.77, 2.28, 2.45)
plot(x,y)

## Use the built-in lm function
fit1 <- lm(y~x)
abline(fit1, lwd=2)
coef(fit1)

## Manual NLL function
nll <- function(pars){
  ## Extract parameters from vector
  intercept <- pars[1]
  slope <- pars[2]
  sd <- exp(pars[3])
  ## Predict y given parameters
  mu <- intercept + slope*x
  ## Calculate log likelihood by hand
  nll <- -sum(dnorm(y, mu, sd, log=T))
  return(nll)
}

## Third way is to use TMB to calculate likelihood
compile("linmod.cpp")
dyn.load(dynlib("linmod"))
obj <- MakeADFun(data=list(x=x,y=y),
                 parameters=list(intercept=1, slope=0.5, logsd=-2),
                 DLL='linmod')

## Check that all three ways matches other ways
pars <- c(1.2, 0.6, -2)
nll(pars)
obj$fn(pars)
## TMB also provides gradient using automatic differentiation
obj$gr(pars)
fit.R <- nlminb(start=pars, objective=nll)
fit.tmb <- nlminb(start=pars, objective=obj$fn, gradient=obj$gr)
coef(fit1)
(mle <- fit.R$par)
## Should be close to zero:
obj$gr(fit.tmb$par)

## Show the improvement in NLL from intitial to MLE
obj$fn(pars)
obj$fn(mle)
plot(x,y)
abline(a=pars[1], b=pars[2], col='blue', lty=2, lwd=2)
abline(a=mle[1], b=mle[2], col='red', lty=3, lwd=2)

## What about the variance estimate?
summary(lm(y~x))
exp(fit.tmb$par[3])
## The sigma parameters do not match because lm() uses REML. We can
## recreate that by integrating out the other fixed effects except the
## variance term.
obj <- MakeADFun(data=list(x=x,y=y), random=c("intercept", "slope"),
                 parameters=list(intercept=1, slope=0.5, logsd=-2),
                 DLL='linmod')
test <- nlminb(start=-2, obj$fn)
exp(test$par)
