# a scripte for viewing MCMC outputs and diagnostics

### THE LIBRARIES:
library(tidyverse)
library(dplyr)
library(rjags)
library(coda)
library(MCMCvis)
library(ecoforecastR)


### OUT object:
# model name
jagmod = j.pests
# iterations and thinning
niter = 100000
jthin = 1
# variable names object
jagvars <- c("beta0", "alpha0",
             "beta[1]", "alpha[1]",
             "R",
             "pa0")

# run from model in SM_JAGS_Model_script
jpout<-coda.samples(jagmod,
                    variable.names = jagvars,
                    n.iter = niter,
                    thin = jthin)

### Let's do some diagnostics:
# plot trace and density plots
plot(jpout)

# GBR 
gelman.diag(jpout)

# GBR plot
gelman.plot(jpout)

# discarding burn in
# set burn in based on GBR
burnin <- 15000
# remove burn in
jburn <- window(jpout, start = burnin)
# plot burn in output
plot(jburn)

# autocorrelation plots
acfplot(jburn)  # ask Mike about this one

# effective sample size
effectiveSize(jburn)

# summary
summary(jburn)

# correlation plots
# convert MCMC to matrix
out <- as.matrix(jburn)
# pairs plot
pairs(out)
# correlations summary
cor(out)


### Model selection step:
# DIC
DIC <- dic.samples(j.pests, n,iter = 10000)
sum <- sum(DIC$deviance, DIC$penalty)



### JAGS Output Figure Code ###
# For when it is time to make figures
# getting proper names for mcmc outputs
##' @param w mcmc object containing matrix outputs
##' @param pre prefix (variable name) for the matrix variable to be extracted
##' @param numeric boolean, whether to coerce class to numeric
parse.MatrixNames <- function(w, pre = "x", numeric = FALSE) {
  w <- sub(pre, "", w)
  w <- sub("[", "", w, fixed = TRUE)
  w <- sub("]", "", w, fixed = TRUE)
  w <- matrix(unlist(strsplit(w, ",")), nrow = length(w), byrow = TRUE)
  if (numeric) {
    class(w) <- "numeric"
  }
  colnames(w) <- c("row", "col")
  return(as.data.frame(w))
}