fnSampH <- function(s.tilde.vec, var.sig.2, PRIOR, J.dim) {
    s.bar <- mean(s.tilde.vec)
    precision <- 1/PRIOR$b.h.2 + J.dim/var.sig.2
    draw.mean <- ( (1/PRIOR$b.h.2)*PRIOR$a.h + (J.dim/var.sig.2)*s.bar )/precision
    draw.var <- precision^(-1)

    return(rnorm(1, draw.mean, sqrt(draw.var)))
}

fnSampVarSig2 <- function(s.tilde.vec, h.scal, PRIOR, J.dim) {
    alpha.param <- PRIOR$a.var.sig + J.dim/2
    beta.param <- PRIOR$b.var.sig + (1/2)*sum( (s.tilde.vec - h.scal)^2)

    return(1/rgamma(1, alpha.param, beta.param))
}

fnSampPsi <- function(d.vec, a.psi) {
    dirichlet.param <- a.psi + d.vec

    return(c(rdirichlet(1, dirichlet.param)))
}

fnSampTau2 <- function(Theta.mat, PRIOR, SETTINGS, J.dim) {
    alpha.params <- PRIOR$a.tau + SETTINGS$M/2 # Only length 1, but is recycled in the IG draw
    beta.params <- PRIOR$b.tau + (1/2)*colSums(Theta.mat^2) # Length J

    1/rgamma(J.dim, alpha.params, beta.params)
}

fnSampPi0Vec <- function(kappa.pi0, xi.pi0, Gamma.mat, J.dim, P.dim) {
    gamma.zero.sums <- colSums(Gamma.mat == 0)
    alpha.params <- kappa.pi0 * xi.pi0 + (J.dim - gamma.zero.sums) # Length P, '(J - gamma.zero.sums)' is the sum of gamma_jp over J where gamma_jp != 0
    beta.params <- kappa.pi0 * (1-xi.pi0) + gamma.zero.sums # Length P

    rbeta(P.dim, alpha.params, beta.params)
}

fnSampPi1Vec <- function(kappa.pi1, xi.pi1, Gamma.mat, P.dim) {
    alpha.params <- kappa.pi1 * xi.pi1 + colSums(Gamma.mat == 1)
    beta.params <- kappa.pi1 * (1- xi.pi1) + colSums(Gamma.mat == -1)

    rbeta(P.dim, alpha.params, beta.params)
}

fnXiLogLik <- function(xi.pi, kappa.pi, pi.vec, a.xi, b.xi) {
    if( (xi.pi < 0) || (xi.pi > 1) ) {
        return(-Inf)
    }

    prior_part_1 <- (b.xi*a.xi - 1) * log(xi.pi)
    prior_part_2 <- (b.xi*(1-a.xi) - 1) * log(1-xi.pi)
    likelihood_parts <- -lgamma(kappa.pi*xi.pi) - lgamma(kappa.pi*(1-xi.pi)) +
        (kappa.pi*xi.pi - 1) * log(pi.vec) + (kappa.pi*(1-xi.pi) - 1) * log(1-pi.vec)

    return(prior_part_1 + prior_part_2 + sum(likelihood_parts))
}

fnSampXi <- function(xi.pi, kappa.pi, pi.vec, a.xi, b.xi, xi.pi.proposal.sd) {
    xi.pi.prop <- xi.pi + rnorm(1, 0, xi.pi.proposal.sd)
    curr.log.lik <- fnXiLogLik(xi.pi, kappa.pi, pi.vec, a.xi, b.xi)
    prop.log.lik <- fnXiLogLik(xi.pi.prop, kappa.pi, pi.vec, a.xi, b.xi)

    if(acceptProposal(curr.log.lik, prop.log.lik)) {
        return(xi.pi.prop)
    }
    else {
        return(xi.pi)
    }
}

fnKappaLogLik <- function(kappa.pi, xi.pi, pi.vec, a.kappa, b.kappa) {
    if(kappa.pi < 0) {
        return(-Inf)
    }

    prior_part_1 <- (a.kappa - 1) * log(kappa.pi)
    prior_part_2 <- -b.kappa * kappa.pi
    likelihood_parts <- -lgamma(kappa.pi*xi.pi) - lgamma(kappa.pi*(1-xi.pi)) + lgamma(kappa.pi) +
        (kappa.pi*xi.pi - 1) * log(pi.vec) + (kappa.pi*(1-xi.pi) - 1) * log(1-pi.vec)

    return(prior_part_1 + prior_part_2 + sum(likelihood_parts))
}

fnSampKappa <- function(kappa.pi, xi.pi, pi.vec, a.kappa, b.kappa, kappa.pi.proposal.sd) {
    kappa.pi.prop <- kappa.pi + rnorm(1, 0, kappa.pi.proposal.sd)
    curr.log.lik <- fnKappaLogLik(kappa.pi, xi.pi, pi.vec, a.kappa, b.kappa)
    prop.log.lik <- fnKappaLogLik(kappa.pi.prop, xi.pi, pi.vec, a.kappa, b.kappa)

    if(acceptProposal(curr.log.lik, prop.log.lik)) {
        return(kappa.pi.prop)
    }
    else {
        return(kappa.pi)
    }
}

fnGetMuMat <- function(alpha0.vec, Theta.mat, Beta.mat, K.mat, X.mat, n.dim, J.dim) {
    Mu.mat <- matrix(0, nrow=n.dim, ncol=J.dim)
    for(i in 1:n.dim) {
        for(j in 1:J.dim) {
            Mu.mat[i, j] <- exp(alpha0.vec[j] + t(K.mat[i,]) %*% Theta.mat[, j] + t(X.mat[i, ]) %*% Beta.mat[j, ])
        }
    }

    return(Mu.mat)
}

# Random generators -------------------------------------------------------
rdirichlet <- function (n, alpha) # Taken from MCMCpack
{
    l <- length(alpha)
    x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
    sm <- x %*% rep(1, l)
    return(x/as.vector(sm))
}

#' @export
rtigamma <- function(n, alpha, beta, ubound) {
    beta/qgamma(runif(n) * pgamma(beta / ubound, alpha, 1, lower.tail = FALSE), alpha, 1, lower.tail = FALSE)
}
