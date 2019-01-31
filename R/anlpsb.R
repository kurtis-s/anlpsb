#' @useDynLib anlpsb
#' @importFrom Rcpp sourceCpp
#' @importFrom RcppTN rtn
NULL

#' Bayesian sparse multivariate regression model with variable selection using asymetric nonlocal priors (NLPs)
#'
#' ANLP-SB models are fit using the \code{\link{anlpsbm}} function
#'
#' @author Kurtis Shuler
#' @references \emph{A Bayesian Sparse Multivariate Regression Model with Asymmetric Nonlocal Priors for Microbiome Data}
#'
#' @docType package
#' @name anlpsb
#' @title ANLP-SB
NULL

#' Fits ANLP-SB model on count data.  Details of the model description are available in the companion paper listed in the references.
#'
#' @param Y.mat \eqn{N x J} matrix of counts
#' @param X.mat \eqn{n x P} design matrix
#' @param time \eqn{n}-dimensional vector of observation timepoints
#' @param rep.K \eqn{n}-dimensional vector of the number of replicates at each timepoint
#' @param Miss.var.ind.mat \eqn{n x Q} (0, 1) matrix where \eqn{Q} is the number of covariates with missing values.  A value of "1" in the matrix indicates the covariate value is missing at that timepoint, and a "0" indicates that covariate value was observed.
#' @param max.miss.cat.vec \eqn{Q}-dimensional vector where the \eqn{q}th element has the number of factor levels in X.mat for the missing covariates.
#' @param first.p.missing.idx1.vec \eqn{Q}-dimensional vector with the \eqn{q}th element indicating the column index where the missing covariate occurs in X.mat
#' @param PRIOR List of the model hyperparameters with the following elements:
#' \itemize{
#'     \item a.h; normal mean for h
#'     \item b.h.2; normal variance for h
#'     \item a.kappa.2; normal mean for kappa.2
#'     \item b.kappa.2; normal variance for kappa.2
#'     \item a.psi.r; dirichlet concentration parameter for psi.r
#'     \item a.w.r; beta shape parameter for w.r
#'     \item b.w.r; beta shape parameter for w.r
#'     \item u.r.w; normal variance for r.tk mixture
#'     \item upsilon.r; r.tk mean constraint
#'     \item b.eta.r.2; normal variance for eta.r
#'     \item a.psi.alpha0; Dirichlet concentration parameter for psi.alpha0
#'     \item a.w.alpha0; beta shape parameter for w.alpha0
#'     \item b.w.alpha0; beta shape parameter for w.alpha0
#'     \item u.alpha0.2; normal variance for alpha.0j mixture
#'     \item upsilon.alpha0; alpha.0j mean constraint
#'     \item b.eta.alpha0.2; normal variance for eta.alpha
#'     \item a.tau; inverse-gamma shape paramter for tau.2
#'     \item b.tau; inverse-gamma scale parameter for tau.2
#'     \item a.pi.0; beta shape parameter for pi.0
#'     \item b.pi.0; beta shape parameter for pi.0
#'     \item a.pi.1; beta shape paramter for pi.1
#'     \item b.pi.1; beta shape paramter for pi.1
#'     \item a.sigma; inverse-gamma shape parameter for sigma.2.p
#'     \item b.sigma; inverse-gamma scale parameter for sigma.2.p
#'     \item a.iota; gamma shape parameter for iota.p
#'     \item b.iota; gamma rate parameter for iota.p
#' }
#'     The companion paper in the references has more details on the prior specification.
#' @param SETTINGS List of model settings with the following elements:
#' \itemize{
#'     \item n.mcmc; number of MCMC iterations
#'     \item print.mcmc; print MCMC progress every print.mcmc iterations.  Leave NULL to supress progress output.
#'     \item n.thin; MCMC thinning interval
#'     \item u.m.basis; knot points for process convolutions
#'     \item phi; variance/range parameter for process convolutions kernel
#'     \item s.j.tilde.proposal.sd; SD for Metropolis proposal of s.j
#'     \item L.r.trunc; number of mixture components for r.tk prior
#'     \item w.ell.r.proposal.sd; SD for Metropolis proposal of w.ell.r
#'     \item r.tk.tilde.proposal.sd; SD for Metropolis proposal of r.tk
#'     \item L.alpha0.trunc; number of mixture components for alpha0j prior
#'     \item w.ell.alpha0.proposal.sd; SD for Metropolis proposal of w.ell.alpha
#'     \item alpha0.j.proposal.sd; SD for Metropolis proposal of alpha.0.j
#'     \item theta.mj.proposal.sd; SD for Metropolis proposal of theta.mj
#'     \item beta.jp.proposal.sd; SD for Metropolis proposal of beta.jp
#'     \item beta.jp.joint.proposal.sd; step-size for joint Metropolis proposal of beta.jp for all P
#'     \item beta.jp.joint.otu.proposal.sd; step-size for joint Metropolis proposal of beta.jp for all J
#'     \item iota.p.proposal.sd; SD for Metropolis proposal of iota.p
#'     \item iota.p.joint.proposal.sd; step-size for joint Metropolis proposal of iota.p
#'     \item sigma.p.2.joint.proposal.sd; step-size for joint Metropolis proposal of sigma.p.2
#' }
#'
#' @return List with elements:
#' \itemize{
#'     \item samps; MCMC samples
#'     \item prior; prior specification
#'     \item settings; model settings
#' }
#'
#' @title ANLP-SB Model
#' @references \emph{A Bayesian Sparse Multivariate Regression Model with Asymmetric Nonlocal Priors for Microbiome Data}
#'
#' @export
anlpsbm <- function(Y.mat, X.mat, time, rep.K, Miss.var.ind.mat=NULL, max.miss.cat.vec=NULL, first.p.missing.idx1.vec=NULL, PRIOR, SETTINGS) {
    J.dim <- ncol(Y.mat) # Number of OTUs
    P.dim <- ncol(X.mat) # Number of covariates
    n.dim <- length(time) # Number of time points
    Tot_N <- nrow(Y.mat)
    SETTINGS$M <- length(SETTINGS$u.m.basis)

    # Change prior for pi to mean/sample size parametrization
    PRIOR$kappa.pi1 <- PRIOR$a.pi.1 + PRIOR$b.pi.1
    PRIOR$xi.pi1 <- PRIOR$a.pi.1/(PRIOR$a.pi.1 + PRIOR$b.pi.1)
    PRIOR$kappa.pi0 <- PRIOR$a.pi.0 + PRIOR$b.pi.0
    PRIOR$xi.pi0 <- PRIOR$a.pi.0/(PRIOR$a.pi.0 + PRIOR$b.pi.0)

    # Old naming convention for kappa.2
    PRIOR$a.var.sig <- PRIOR$a.kappa.2
    PRIOR$b.var.sig <- PRIOR$b.kappa.2

    # Estimated params --------------------------------------------------------
    rep.K.cumsum <- c(0, cumsum(rep.K)[-n.dim])
    r.hat.vec <- rowSums(Y.mat)/sum(Y.mat)
    r.tilde.hat.vec <- log(r.hat.vec)
    alpha0.hat.vec <-  log(colMeans(apply(Y.mat, 2, function(col) col/r.hat.vec)))

    time.idx <- rep(1:n.dim, times=rep.K) # Time index for each of the N=158 obs
    Y.frame <- data.frame(cbind(factor(time.idx), Y.mat))
    colnames(Y.frame) <- c("TimeIdx", colnames(Y.mat))
    Mu.tjk.hat.frame <- data.frame(cbind(time.idx, apply(Y.mat, 2, function(otu.col) otu.col/r.hat.vec)))
    colnames(Mu.tjk.hat.frame) <- c("TimeIdx", paste("Y", 1:ncol(Y.mat), sep=""))
    Mu.hat.mat <- as.matrix(aggregate(. ~ TimeIdx, Mu.tjk.hat.frame, mean)[,-1])
    Mu.hat.mat <- Mu.hat.mat + 1 # Add one just so mu != 0

    # Initialization ----------------------------------------------------------
    K.mat <- t(dnorm(sapply(time, function(t) t-SETTINGS$u.m.basis), 0, SETTINGS$phi))
    Mean.func.design.mat <- cbind(K.mat, X.mat)
    Mean.func.hat.mat <- sapply(1:J.dim, function(j) coef(arm::bayesglm(log(Mu.hat.mat)[, j] ~ Mean.func.design.mat, prior.df=Inf)))
    alpha0.hat.vec2 <- Mean.func.hat.mat[1, ]
    Theta.hat.mat <- Mean.func.hat.mat[2:(SETTINGS$M+1), ]
    Beta.hat.mat <- t(Mean.func.hat.mat[(SETTINGS$M+2):nrow(Mean.func.hat.mat), ])
    sigma.2.hat.vec <- apply(Beta.hat.mat, 2, var)
    Beta.hat.mat.standardized <- scale(Beta.hat.mat, center=FALSE, scale=TRUE)
    Beta.hat.mat <- ifelse((-0.2 < Beta.hat.mat.standardized) & (Beta.hat.mat.standardized < 0.2), 0, Beta.hat.mat)
    Gamma.hat.mat <- matrix(ifelse(Beta.hat.mat > 0, 1, -1), nrow=J.dim, ncol=P.dim)
    Gamma.hat.mat[Beta.hat.mat==0] <- 0

    # MCMC Initialization --------------------------------------------------
    ## r.tk
    r.vec <- r.hat.vec
    r.tilde.vec <- log(r.hat.vec)
    c.r.vec <- .bincode(r.vec, quantile(r.vec, (0:5)/5), right=TRUE, include.lowest=TRUE)
    d.r.vec <- tabulate(c.r.vec, nbins=SETTINGS$L.r.trunc)
    eta.r.vec <- rnorm(SETTINGS$L.r.trunc, 0, 1)
    w.r.vec <- rbeta(SETTINGS$L.r.trunc, 0.5, 0.5)
    lambda.r.vec <- rbinom(Tot_N, size=1, prob=w.r.vec[c.r.vec])
    psi.r.vec <- fnSampPsi(d.r.vec, PRIOR$a.psi.r)
    ## alpha0.j
    alpha0.vec <- alpha0.hat.vec2
    c.alpha0.vec <- .bincode(alpha0.vec, quantile(alpha0.vec, (0:5)/5), right=TRUE, include.lowest=TRUE)
    d.alpha0.vec <- tabulate(c.alpha0.vec, nbins=SETTINGS$L.alpha0.trunc)
    eta.alpha0.vec <- rnorm(SETTINGS$L.alpha0.trunc, 0, 1)
    w.alpha0.vec <- rbeta(SETTINGS$L.alpha0.trunc, 0.5, 0.5)
    lambda.alpha0.vec <- rbinom(J.dim, size=1, prob=w.alpha0.vec[c.alpha0.vec])
    psi.alpha0.vec <- fnSampPsi(d.alpha0.vec, PRIOR$a.psi.alpha0)
    ## theta.mj
    Theta.mat <- Theta.hat.mat
    tau.2.vec <- fnSampTau2(Theta.mat, PRIOR, SETTINGS, J.dim)
    ## beta.jp
    Beta.mat <- Beta.hat.mat
    Gamma.mat <- Gamma.hat.mat
    xi.pi0 <- PRIOR$xi.pi0
    xi.pi1 <- PRIOR$xi.pi1
    kappa.pi0 <- PRIOR$kappa.pi0
    kappa.pi1 <- PRIOR$kappa.pi1
    pi0.vec <- fnSampPi0Vec(kappa.pi0, xi.pi0, Gamma.mat, J.dim, P.dim)
    pi1.vec <- fnSampPi1Vec(kappa.pi1, xi.pi1, Gamma.mat, P.dim)
    sigma.2.vec <- sigma.2.hat.vec
    tmp.for.iota <- abs(Beta.mat)
    tmp.for.iota[tmp.for.iota==0] <- NA
    tmp.for.iota <- apply(tmp.for.iota, 2, min, na.rm=TRUE)
    iota.vec <- tmp.for.iota/sqrt(sigma.2.vec) - 10e-10
    ## s.j
    h.scal <- -0.7
    var.sig.2 <- 1/6
    s.tilde.vec <- rnorm(J.dim, h.scal, sqrt(var.sig.2))
    s.vec <- exp(s.tilde.vec)
    ## Mu.mat
    Mu.mat <- fnGetMuMat(alpha0.vec, Theta.mat, Beta.mat, K.mat, X.mat, n.dim, J.dim)
    ## Cleanup
    rm(Mu.hat.mat)
    rm(r.hat.vec)
    rm(r.tilde.hat.vec)
    rm(alpha0.hat.vec)
    rm(Mean.func.design.mat)
    rm(Mean.func.hat.mat)
    rm(alpha0.hat.vec2)
    rm(Theta.hat.mat)
    rm(Beta.hat.mat)
    rm(Beta.hat.mat.standardized)
    rm(Gamma.hat.mat)
    rm(tmp.for.iota)
    rm(sigma.2.hat.vec)

    # MCMC --------------------------------------------------------------------
    SAMPS <- list()
    ## s.j
    SAMPS$s.vec <- matrix(numeric(0), nrow=SETTINGS$n.mcmc, ncol=J.dim)
    SAMPS$h.scal <- matrix(numeric(0), nrow=SETTINGS$n.mcmc, ncol=1)
    SAMPS$var.sig.2 <- matrix(numeric(0), nrow=SETTINGS$n.mcmc, ncol=1)
    ## r.tk
    SAMPS$psi.r.vec <- matrix(numeric(0), nrow=SETTINGS$n.mcmc, ncol=SETTINGS$L.r.trunc)
    SAMPS$w.r.vec <- matrix(numeric(0), nrow=SETTINGS$n.mcmc, ncol=SETTINGS$L.r.trunc)
    SAMPS$eta.r.vec <- matrix(numeric(0), nrow=SETTINGS$n.mcmc, ncol=SETTINGS$L.r.trunc)
    SAMPS$c.r.vec <- matrix(numeric(0), nrow=SETTINGS$n.mcmc, ncol=Tot_N)
    SAMPS$lambda.r.vec <- matrix(numeric(0), nrow=SETTINGS$n.mcmc, ncol=Tot_N)
    SAMPS$r.vec <- matrix(numeric(0), nrow=SETTINGS$n.mcmc, ncol=Tot_N)
    ## alpha0.j
    SAMPS$psi.alpha0.vec <- matrix(numeric(0), nrow=SETTINGS$n.mcmc, ncol=SETTINGS$L.alpha0.trunc)
    SAMPS$w.alpha0.vec <- matrix(numeric(0), nrow=SETTINGS$n.mcmc, ncol=SETTINGS$L.alpha0.trunc)
    SAMPS$eta.alpha0.vec <- matrix(numeric(0), nrow=SETTINGS$n.mcmc, ncol=SETTINGS$L.alpha0.trunc)
    SAMPS$c.alpha0.vec <- matrix(numeric(0), nrow=SETTINGS$n.mcmc, ncol=J.dim)
    SAMPS$lambda.alpha0.vec <- matrix(numeric(0), nrow=SETTINGS$n.mcmc, ncol=J.dim)
    SAMPS$alpha0.vec <- matrix(numeric(0), nrow=SETTINGS$n.mcmc, ncol=J.dim)
    ## theta.mj
    SAMPS$Theta.mat <- matrix(numeric(0), nrow=SETTINGS$n.mcmc, ncol=SETTINGS$M*J.dim)
    SAMPS$tau.2.vec <- matrix(numeric(0), nrow=SETTINGS$n.mcmc, ncol=J.dim)
    ## beta.jp
    SAMPS$Beta.mat <- matrix(numeric(0), nrow=SETTINGS$n.mcmc, ncol=J.dim*P.dim)
    SAMPS$Gamma.mat <- matrix(numeric(0), nrow=SETTINGS$n.mcmc, ncol=J.dim*P.dim)
    SAMPS$iota.vec <- matrix(numeric(0), nrow=SETTINGS$n.mcmc, ncol=P.dim)
    SAMPS$sigma.2.vec <- matrix(numeric(0), nrow=SETTINGS$n.mcmc, ncol=P.dim)
    ## pi0/pi1
    SAMPS$pi0.vec <- matrix(numeric(0), nrow=SETTINGS$n.mcmc, ncol=P.dim)
    SAMPS$pi1.vec <- matrix(numeric(0), nrow=SETTINGS$n.mcmc, ncol=P.dim)
    ## Mu.ij
    SAMPS$Mu.mat <- matrix(numeric(0), nrow=SETTINGS$n.mcmc, ncol=n.dim*J.dim)
    for(b in 1:SETTINGS$n.mcmc) {
        for(thin.iter in 1:SETTINGS$n.thin) {
            ## Impute missing covariates
            if(!is.null(Miss.var.ind.mat)) {
                fnImputeZ(X.mat, Miss.var.ind.mat, max.miss.cat.vec, first.p.missing.idx1.vec, Mu.mat, Beta.mat, s.vec, r.vec, rep.K, rep.K.cumsum, Y.mat)
            }

            ## s.j
            fnSampSVec(s.vec, s.tilde.vec, r.vec, Mu.mat, rep.K, Y.mat, h.scal, var.sig.2, SETTINGS$s.j.tilde.proposal.sd)
            h.scal <- fnSampH(s.tilde.vec, var.sig.2, PRIOR, J.dim)
            var.sig.2 <- fnSampVarSig2(s.tilde.vec, h.scal, PRIOR, J.dim)
            ## r.tk
            psi.r.vec <- fnSampPsi(d.r.vec, PRIOR$a.psi.r)
            fnSampWVec(w.r.vec, r.tilde.vec, lambda.r.vec,
                       c.r.vec, eta.r.vec, PRIOR$upsilon.r, PRIOR$a.w.r,
                       PRIOR$b.w.r, PRIOR$u.r.2, SETTINGS$w.ell.r.proposal.sd)
            fnSampEtaVec(eta.r.vec, r.tilde.vec, w.r.vec, lambda.r.vec, c.r.vec,
                         PRIOR$u.r.2, PRIOR$upsilon.r, PRIOR$b.eta.r.2)
            fnSampCVec(c.r.vec, d.r.vec, r.tilde.vec, psi.r.vec, w.r.vec, eta.r.vec,
                       PRIOR$upsilon.r, PRIOR$u.r.2)
            fnSampLambdaVec(lambda.r.vec, r.tilde.vec, c.r.vec, w.r.vec,
                            eta.r.vec, PRIOR$upsilon.r, PRIOR$u.r.2)
            fnSampRVec(r.vec, r.tilde.vec, Mu.mat, s.vec, lambda.r.vec, c.r.vec,
                       eta.r.vec, w.r.vec, rep.K, PRIOR$upsilon.r, PRIOR$u.r.2,
                       Y.mat, SETTINGS$r.tk.tilde.proposal.sd)

            ## alpha0.j
            psi.alpha0.vec <- fnSampPsi(d.alpha0.vec, PRIOR$a.psi.alpha0)
            fnSampWVec(w.alpha0.vec, alpha0.vec, lambda.alpha0.vec,
                       c.alpha0.vec, eta.alpha0.vec, PRIOR$upsilon.alpha0, PRIOR$a.w.alpha0,
                       PRIOR$b.w.alpha0, PRIOR$u.alpha0.2, SETTINGS$w.ell.alpha0.proposal.sd)
            fnSampEtaVec(eta.alpha0.vec, alpha0.vec, w.alpha0.vec, lambda.alpha0.vec, c.alpha0.vec,
                         PRIOR$u.alpha0.2, PRIOR$upsilon.alpha0, PRIOR$b.eta.alpha0.2)
            fnSampCVec(c.alpha0.vec, d.alpha0.vec, alpha0.vec, psi.alpha0.vec, w.alpha0.vec, eta.alpha0.vec,
                       PRIOR$upsilon.alpha0, PRIOR$u.alpha0.2)
            fnSampLambdaVec(lambda.alpha0.vec, alpha0.vec, c.alpha0.vec, w.alpha0.vec,
                            eta.alpha0.vec, PRIOR$upsilon.alpha0, PRIOR$u.alpha0.2)
            fnSampAlpha0Vec(alpha0.vec, r.vec, Mu.mat, s.vec, lambda.alpha0.vec,
                            c.alpha0.vec, eta.alpha0.vec, w.alpha0.vec,
                            PRIOR$upsilon.alpha0, PRIOR$u.alpha0.2, rep.K,
                            Y.mat, SETTINGS$alpha0.j.proposal.sd)

            ## theta.mj
            fnSampThetaMat(Theta.mat, tau.2.vec, s.vec, r.vec, Mu.mat, rep.K, K.mat, Y.mat, SETTINGS$theta.mj.proposal.sd)
            tau.2.vec <- fnSampTau2(Theta.mat, PRIOR, SETTINGS, J.dim)

            ## beta.jp
            fnSampBetaMat(Beta.mat, Gamma.mat, Mu.mat, s.vec, r.vec, X.mat, sigma.2.vec, rep.K, Y.mat, iota.vec, SETTINGS$beta.jp.proposal.sd)
            fnSampBetaMatGammaMatJoint(Beta.mat, Gamma.mat, Mu.mat, X.mat, sigma.2.vec, iota.vec, s.vec, r.vec, rep.K, pi0.vec, pi1.vec, Y.mat, SETTINGS$beta.jp.joint.proposal.sd)
            fnBetaGammaAddSwpDel(Beta.mat, Gamma.mat, sigma.2.vec, iota.vec, Mu.mat, X.mat, s.vec, r.vec, rep.K, pi0.vec, pi1.vec, Y.mat)

            ## iota.p
            fnSampIotaVec(iota.vec, Beta.mat, Gamma.mat, sigma.2.vec, PRIOR$a.iota, PRIOR$b.iota, SETTINGS$iota.p.proposal.sd)

            ## pi0/pi1
            pi0.vec <- fnSampPi0Vec(kappa.pi0, xi.pi0, Gamma.mat, J.dim, P.dim)
            pi1.vec <- fnSampPi1Vec(kappa.pi1, xi.pi1, Gamma.mat, P.dim)
            ## sigma.p.2
            fnSampSigma2Vec(sigma.2.vec, Beta.mat, Gamma.mat, iota.vec, rtigamma, PRIOR$a.sigma, PRIOR$b.sigma)

            ## beta.jp, sigma.p.2, and iota.p joint update
            fnSampIotaVecSigma2VecJoint(sigma.2.vec, iota.vec, Beta.mat, Gamma.mat, Mu.mat, X.mat, s.vec, r.vec, rep.K, pi0.vec, pi1.vec, PRIOR$a.sigma, PRIOR$b.sigma, PRIOR$a.iota, PRIOR$b.iota, Y.mat, SETTINGS$iota.p.joint.proposal.sd, SETTINGS$sigma.p.2.joint.proposal.sd)

            ## beta.jp joint update for the entire OTU
            fnSampBetaMatGammaMatJointWholeOtu(Beta.mat, Gamma.mat, Mu.mat, X.mat, sigma.2.vec, iota.vec, s.vec, r.vec, rep.K, pi0.vec, pi1.vec, Y.mat, SETTINGS$beta.jp.joint.otu.proposal.sd)
        }

        ### Save samples
        ## s.j
        SAMPS$s.vec[b, ] <- s.vec
        SAMPS$h.scal[b, ] <- h.scal
        SAMPS$var.sig.2[b, ] <- var.sig.2
        ## r.tk
        SAMPS$psi.r.vec[b, ] <- psi.r.vec
        SAMPS$w.r.vec[b, ] <- w.r.vec
        SAMPS$eta.r.vec[b, ] <- eta.r.vec
        SAMPS$c.r.vec[b, ] <- c.r.vec
        SAMPS$lambda.r.vec[b, ] <- lambda.r.vec
        SAMPS$r.vec[b, ] <- r.vec
        ## alpha0.j
        SAMPS$psi.alpha0.vec[b, ] <- psi.alpha0.vec
        SAMPS$w.alpha0.vec[b, ] <- w.alpha0.vec
        SAMPS$eta.alpha0.vec[b, ] <- eta.alpha0.vec
        SAMPS$c.alpha0.vec[b, ] <- c.alpha0.vec
        SAMPS$lambda.alpha0.vec[b, ] <- lambda.alpha0.vec
        SAMPS$alpha0.vec[b, ] <- alpha0.vec
        ## theta.mj
        SAMPS$Theta.mat[b, ] <- Theta.mat
        SAMPS$tau.2.vec[b, ] <- tau.2.vec
        ## beta.jp
        SAMPS$Beta.mat[b, ] <- Beta.mat
        SAMPS$Gamma.mat[b, ] <- Gamma.mat
        SAMPS$iota.vec[b, ] <- iota.vec
        SAMPS$sigma.2.vec[b, ] <- sigma.2.vec
        ## pi0/pi1
        SAMPS$pi0.vec[b, ] <- pi0.vec
        SAMPS$pi1.vec[b, ] <- pi1.vec
        ## mu.ij
        SAMPS$Mu.mat[b, ] <- Mu.mat

        if(!is.null(SETTINGS$print.mcmc) && ((b %% SETTINGS$print.mcmc) == 0)) {
            cat(format(Sys.time()), " b=", b, "\n", sep="")
        }
    }

    return(list(samps=SAMPS, prior=PRIOR, settings=SETTINGS))
}
