library(Rcpp)
library(RcppEigen)
library(numDeriv)
library(expm)

Sys.setenv("PKG_CPPFLAGS" = "-march=native")
sourceCpp(file = "Z:/network/CoxCNet/SimSet.cpp", verbose = TRUE, rebuild = TRUE)

n = 20
p = 1
Nij = matrix(0, nrow = n, ncol = n)
Zij = matrix(rnorm(n*n), nrow = n, ncol = n)
M = matrix(0, nrow = (n+p-1), ncol = 3)

for (i in 1:8) {
  M[i, 1] = 1
}

for (i in 9:19) {
  M[i, 2] = 1
}

M[20, 3] = 1


beta = rep(0, n + p - 1)
betax = rep(0, 3)


lik <- function(beta, Nij, Zij, n, p) {

  Pij = matrix(0, nrow = n, ncol = n)

  x = c(0, beta[1:(n-1)])
  y = beta[n:(n+p-1)]
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      Pij[i, j] = x[i] + x[j]
      for (z in 1:p) {
        Pij[i, j] = Pij[i, j] + Zij[i, j] * y[z]
      }
      Pij[i, j] = exp(Pij[i, j])
      Pij[j, i] = Pij[i, j]
    }
  }

  P0 = sum(Pij) / 2
  Pi = colSums(Pij)

  N0 = sum(Nij)

  l = 0
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      l = l + Nij[i, j] * log(Pij[i, j])
    }
  }
  l = l - N0 * log(P0)

  return(-l)
}


likMx <- function(betax, M, Nij, Zij, n, p) {

  beta = M %*% betax
  Pij = matrix(0, nrow = n, ncol = n)

  x = c(0, beta[1:(n-1)])
  y = beta[n:(n+p-1)]
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      Pij[i, j] = x[i] + x[j]
      for (z in 1:p) {
        Pij[i, j] = Pij[i, j] + Zij[i, j] * y[z]
      }
      Pij[i, j] = exp(Pij[i, j])
      Pij[j, i] = Pij[i, j]
    }
  }

  P0 = sum(Pij) / 2
  Pi = colSums(Pij)

  N0 = sum(Nij)

  l = 0
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      l = l + Nij[i, j] * log(Pij[i, j])
    }
  }
  l = l - N0 * log(P0)

  return(-l)
}

# Test functions

lik(beta, Nij, Zij, n, p)
likMx(betax, M, Nij, Zij, n, p)

g1 = grad(lik, beta, Nij = Nij, Zij = Zij, n = n, p = p)
g2 = grad(likMx, betax, M = M, Nij = Nij, Zij = Zij, n = n, p = p)

h1 = hessian(lik, beta, Nij = Nij, Zij = Zij, n = n, p = p)
h2 = hessian(likMx, betax, M = M, Nij = Nij, Zij = Zij, n = n, p = p)


# Simulation --------------------------------------------------------------

betaAll = matrix(0, 100, 3)
betaAllz = matrix(0, 100, 3)

for (z in 1:100) {

  set.seed(z)
  Zij = matrix(rnorm(n*n), nrow = n, ncol = n)
  trail_sim = SimSetC(n, 0.5, array(Zij, c(n, n, 1)))
  nn = length(trail_sim[[1]])

  Nij = matrix(0, nrow = n, ncol = n)
  for (i in 1:nn) {
    p1 = trail_sim[[1]][i]
    q1 = trail_sim[[2]][i]
    Nij[p1, q1] = Nij[p1, q1] + 1
  }

  # # FULL
  # beta0 = rep(0, n + p - 1)
  # for (i in 1:10) {
  #
  #   g1 = grad(lik, beta0, Nij = Nij, Zij = Zij, n = n, p = p)
  #   h1 = hessian(lik, beta0, Nij = Nij, Zij = Zij, n = n, p = p)
  #
  #   beta1 = beta0 - solve(h1, g1)
  #   beta0 = beta1
  #
  #   cat(i, "\n")
  # }

  # GROUP
  beta0 = rep(0, 3)
  for (i in 1:10) {

    g2 = grad(likMx, beta0, M = M, Nij = Nij, Zij = Zij, n = n, p = p)
    h2 = hessian(likMx, beta0, M = M, Nij = Nij, Zij = Zij, n = n, p = p)

    beta1 = beta0 - solve(h2, g2)
    beta0 = beta1

  }

  betaAll[z, ] = beta0
  betaAllz[z, ] = sqrtm(-h2) %*% (beta0 -  c(0, 1, 0.2))
  cat(z, "\n")

}



t(M) %*% h1 %*% M

beta = c(rep(0, 8), rep(1, 11), 0.5)
betax = c(0, 1, 0.5)
betax = beta0
betax[3] = 0.5

likMx(betax, M, Nij, Zij, n, p)


# ADMM SCAD ---------------------------------------------------------------

st <- function(t, lam) {
  return(sign(t)*(abs(t) - lam)*((abs(t) - lam) > 0))
}

stscad <- function(muij, lam, v, gam, n) {
  if (abs(muij) <= lam + lam/v) {
    return(st(muij, lam/v))
  } else if (lam + lam/v < abs(muij) && abs(muij) < gam*lam) {
    return(st(muij, gam*lam/((gam-1)*v))/ (1 - 1/((gam-1)*v)))
  } else {
    return(muij)
  }
}


likP <- function(beta, Nij, Zij, n, p, vij, etaij, v) {

  Pij = matrix(0, nrow = n, ncol = n)

  x = c(0, beta[1:(n-1)])
  y = beta[n:(n+p-1)]
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      Pij[i, j] = x[i] + x[j]
      for (z in 1:p) {
        Pij[i, j] = Pij[i, j] + Zij[i, j] * y[z]
      }
      Pij[i, j] = exp(Pij[i, j])
      Pij[j, i] = Pij[i, j]
    }
  }

  P0 = sum(Pij) / 2
  Pi = colSums(Pij)

  N0 = sum(Nij)

  l = 0
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      l = l + Nij[i, j] * log(Pij[i, j])
    }
  }
  l = l - N0 * log(P0)

  l = -l
  for (i in 1:(n-2)) {
    for (j in (i+1):(n-1)) {
    l = l + vij[i, j] * (beta[i] - beta[j] - etaij[i, j]) + v/2 * (beta[i] - beta[j] - etaij[i, j])^2
    }
  }
  return(l)
}

likP(beta0, Nij, Zij, n, p, vij, etaij, v)
likP(beta1, Nij, Zij, n, p, vij, etaij, v)
likP(beta2, Nij, Zij, n, p, vij, etaij, v)

beta2 = beta0 + solve(h1, g1)/100
lik(beta2, Nij, Zij, n, p)

admm_scad <- function(beta, Nij, Zij, n, p, lam, v) {

  vij = matrix(0, nrow = n - 1, ncol = n - 1)
  etaij = matrix(0, nrow = n - 1, ncol = n - 1)
  dij = matrix(0, nrow = n - 1, ncol = n - 1)
  one = matrix(1, nrow = n - 1, ncol = n - 1)
  Iij = diag(1, n-1)

  beta0 = rep(0.1, n+p-1)
  beta = beta0[1:(n-1)]
  sum_beta = sum(beta)

  for (i1 in 1:10) {
    for (i2 in 1:10) {
      # g1 = grad(lik, beta0, Nij = Nij, Zij = Zij, n = n, p = p)
      # g12 = v * ((n-1)*beta - sum_beta - colSums(etaij) - rowSums(etaij)) + colSums(vij) + rowSums(vij)
      # g12 = c(g12, 0)
      # g1 = g1 + g12

      # h1 = hessian(lik, beta0, Nij = Nij, Zij = Zij, n = n, p = p)
      # h1[1:(n-1), 1:(n-1)] = h1[1:(n-1), 1:(n-1)] + (n-1)*v*Iij - v*one

      g13 = grad(likP, beta0, Nij = Nij, Zij = Zij, n = n, p = p, vij = vij, etaij = etaij, v = v)
      h13 = hessian(likP, beta0, Nij = Nij, Zij = Zij, n = n, p = p, vij = vij, etaij = etaij, v = v)

      beta1 = beta0 - solve(h13, g13)
      beta0 = beta1
    }

    beta = beta0[1:(n-1)]
    sum_beta = sum(beta)

    for (i in 1:(n-2)) {
      for (j in (i+1):(n-1)) {
        dij[i, j] = beta[i] - beta[j] + vij[i, j]/v
        etaij[i, j] = stscad(dij[i, j], lam, v, gam, n)
        vij[i, j] = vij[i, j] + v*(beta[i] - beta[j] - etaij[i, j])
      }
    }

    cat(round(beta, 2), "\n")
  }
}
