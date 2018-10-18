# ANLP-SB 
ANLP-SB is a Bayesian sparse multivariate regression model with variable selection using asymetric nonlocal priors (NLPs) for count data.  Details of the model specification and its implementation are described in the accompanying paper *A Bayesian Sparse Multivariate Regression Model with
Asymmetric Nonlocal Priors for Microbiome Data*.

## Installation
After installing devtools run:
```r
library(devtools)
install_github("kurtis-s/anlpsb")
```

## Example Code
The code below shows how to run the model on the simulated count dataset included in the package:
```r
rm(list=ls())

set.seed(389228)

library(anlpsb)

### Simulated data
ALL_dat <- readRDS(system.file("extdata", "SIM_dat.rds", package="anlpsb"))

Y.mat <- ALL_dat$Y
Z.mat <- as.matrix(ALL_dat$Z[, -c(1:2)]) # Z matrix without time index and days
rep.K.cumsum <- c(0, cumsum(ALL_dat$rep_K)[-ALL_dat$n])
colnames(Y.mat) <- paste("OTU", 1:ncol(Y.mat), sep="")

r.hat.vec <- rowSums(Y.mat)/sum(Y.mat)
r.tilde.hat.vec <- log(r.hat.vec)
alpha0.hat.vec <-  log(colMeans(apply(Y.mat, 2, function(col) col/r.hat.vec)))

# Settings ----------------------------------------------------------------
SETTINGS <- list()
### MCMC
SETTINGS$n.mcmc <- 10
SETTINGS$print.mcmc <- SETTINGS$n.mcmc/100
SETTINGS$n.thin <- 2
### Kernal
M.knots <- 70 # Number of knot points
SETTINGS$u.m.basis <- seq(min(ALL_dat$Z$days)-0.5, max(ALL_dat$Z$days)+0.5, length.out=M.knots)
SETTINGS$phi <- (2*ALL_dat$n/M.knots)
### s_j
SETTINGS$s.j.tilde.proposal.sd <- 0.1
### r_tk
SETTINGS$L.r.trunc <- 50 # Truncation level for variance inflation factors
SETTINGS$w.ell.r.proposal.sd <- 0.04
SETTINGS$r.tk.tilde.proposal.sd <- 0.075
### alpha0_j
SETTINGS$L.alpha0.trunc <- 50 # Truncation level for baseline OTU abundance
SETTINGS$w.ell.alpha0.proposal.sd <- 0.04
SETTINGS$alpha0.j.proposal.sd <- 0.5
### theta_mj
SETTINGS$theta.mj.proposal.sd <- 0.5
### beta_jp
SETTINGS$beta.jp.proposal.sd <- 0.1
SETTINGS$beta.jp.joint.proposal.sd <- 1
SETTINGS$beta.jp.joint.otu.proposal.sd <- 0.05
## iota_p
SETTINGS$iota.p.proposal.sd <- 0.1
## iota_p and sigma_p_2 joint
SETTINGS$iota.p.joint.proposal.sd <- 0.05
SETTINGS$sigma.p.2.joint.proposal.sd <- 0.05

# Prior -------------------------------------------------------------------
PRIOR <- list()
### s_j (h and var.sig)
PRIOR$a.h <- -10.0
PRIOR$b.h.2 <- 10.0^2
PRIOR$a.kappa.2 <- 1e-05
PRIOR$b.kappa.2 <- 1e-05
### r_tk
PRIOR$a.psi.r <- rep(1.0, SETTINGS$L.r.trunc)
PRIOR$a.w.r <- 0.5
PRIOR$b.w.r <- 0.5
PRIOR$u.r.2 <- 0.1
PRIOR$upsilon.r <- mean(r.tilde.hat.vec)
PRIOR$b.eta.r.2 <- 0.3
# Alpha0_j
PRIOR$a.psi.alpha0 <- rep(1.0, SETTINGS$L.alpha0.trunc)
PRIOR$a.w.alpha0 <- 0.5
PRIOR$b.w.alpha0 <- 0.5
PRIOR$u.alpha0.2 <- 0.1
PRIOR$upsilon.alpha0 <- mean(alpha0.hat.vec)
PRIOR$b.eta.alpha0.2 <- 1
# tau.j.2
PRIOR$a.tau <- 1.0
PRIOR$b.tau <- 1.0
PRIOR$a.pi.0 <- 1
PRIOR$b.pi.0 <- ALL_dat$P
PRIOR$a.pi.1 <- 5
PRIOR$b.pi.1 <- 5
# sigma.p.2
PRIOR$a.sigma <- 1.0
PRIOR$b.sigma <- 1.0
# iota
PRIOR$a.iota <- 2.5
PRIOR$b.iota <- 10

Y.mat <- ALL_dat$Y
time <- ALL_dat$Z$days


# Missing values setup ----------------------------------------------------
missing.vars <- colnames(ALL_dat$Miss_ind[, -c(1, 2)])[colSums(ALL_dat$Miss_ind[, -c(1, 2)]) > 0] # First two columns are the index and the day

n.missing.vars <- length(missing.vars) # Number of missing variables, NOT using dummy encoding e.g. "bloom" is only 1 variable, not 4
max.miss.cat.vec <- sapply(missing.vars, function(missing.var) { # Number of categories for each of the missing variables{
    sum(startsWith(colnames(Z.mat), missing.var))
})

first.p.missing.idx1.vec <- sapply(missing.vars, function(missing.var) { # Number of categories for each of the missing variables{
    min((1:ncol(Z.mat))[startsWith(colnames(Z.mat), missing.var)])
})

Miss.var.ind.mat <- ALL_dat$Miss_ind[, -c(1, 2)]
Miss.var.ind.mat <- ALL_dat$Miss_ind[, missing.vars, drop=FALSE] # Remove 'Index', 'days', and non-missing vars


# Cleanup -----------------------------------------------------------------
rm(list=ls()[!(ls() %in% c("PRIOR", "SETTINGS", "ALL_dat", "Z.mat", "Y.mat", "time", "Miss.var.ind.mat", "max.miss.cat.vec", "first.p.missing.idx1.vec"))])
mod <- anlpsbm(Y.mat, Z.mat, time, ALL_dat$rep_K, Miss.var.ind.mat, max.miss.cat.vec, first.p.missing.idx1.vec, PRIOR=PRIOR, SETTINGS=SETTINGS)
```
