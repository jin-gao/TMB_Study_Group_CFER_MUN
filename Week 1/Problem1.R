library(TMB)

## : using TMB derivative
compile("Week 1/polynomial.cpp")
dyn.load(dynlib("Week 1/polynomial"))
obj <- MakeADFun(data=list(), parameters=list(x=0), DLL='polynomial')
obj$fn(5)
obj$gr(5)

