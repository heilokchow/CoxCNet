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

trail_sim = SimSetC(n, 0.5, array(Zij, c(n, n, 1)))
nn = length(trail_sim[[1]])

for (i in 1:nn) {
  p1 = trail_sim[[1]][i]
  q1 = trail_sim[[2]][i]
  Nij[p1, q1] = Nij[p1, q1] + 1
  Nij[q1, p1] = Nij[q1, p1] + 1
}



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

  N0 = sum(Nij) / 2

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

  N0 = sum(Nij) / 2

  l = 0
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      l = l + Nij[i, j] * log(Pij[i, j])
    }
  }
  l = l - N0 * log(P0)

  return(l)
}


lik(beta, Nij, Zij, n, p)
likMx(betax, M, Nij, Zij, n, p)

g1 = grad(lik, beta, Nij = Nij, Zij = Zij, n = n, p = p)
g2 = grad(likMx, betax, M = M, Nij = Nij, Zij = Zij, n = n, p = p)

h1 = hessian(lik, beta, Nij = Nij, Zij = Zij, n = n, p = p)
h2 = hessian(likMx, betax, M = M, Nij = Nij, Zij = Zij, n = n, p = p)

t(M) %*% h1 %*% M

beta = c(rep(0, 8), rep(1, 11), 0.5)
betax = c(0, 1, 0.5)
