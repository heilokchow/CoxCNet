library(Rcpp)
library(RcppEigen)
library(numDeriv)

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

  return(l)
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

  return(l)
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
  cat(z, "\n")

}



t(M) %*% h1 %*% M

beta = c(rep(0, 8), rep(1, 11), 0.5)
betax = c(0, 1, 0.5)
betax = beta0
betax[3] = 0.5

likMx(betax, M, Nij, Zij, n, p)

