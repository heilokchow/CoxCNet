library(Rcpp)
library(RcppEigen)
library(numDeriv)
library(expm)
library(latex2exp)

Sys.setenv("PKG_CPPFLAGS" = "-march=native")
sourceCpp(file = "Z:/network/CoxCNet/SimSet.cpp", verbose = TRUE, rebuild = TRUE)
sourceCpp(file = "Z:/network/CoxCNet/NT.cpp", verbose = TRUE, rebuild = TRUE)

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

test = GradCN(beta, Nij, Zij, n, p, vij, etaij, v)

# Simulation --------------------------------------------------------------

betaAll = matrix(0, 100, 3)
betaAllz = matrix(0, 100, 3)

for (z in 1:100) {

  set.seed(z)
  Zij = matrix(rnorm(n*n), nrow = n, ncol = n)
  Zij = (Zij + t(Zij)) / 2
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

for (i in 1:(n-1)) {
  for (j in (i+1):n) {
    etaij[i, j] = 1
  }
}

grad(likP, beta, Nij = Nij, Zij = Zij, n = n, p = p, vij = vij, etaij = etaij, v = v)
hessian(likP, beta, Nij = Nij, Zij = Zij, n = n, p = p, vij = vij, etaij = etaij, v = v)
test = GradCN(beta, Nij, Zij, n, p, vij, etaij, v)

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
  for (i in 1:(n-1)) {
    for (j in (i+1):(n)) {
    l = l + vij[i, j] * (x[i] - x[j] - etaij[i, j]) + v/2 * (x[i] - x[j] - etaij[i, j])^2
    }
  }
  return(l)
}

likP(beta0, Nij, Zij, n, p, vij, etaij, v)
likP(beta1, Nij, Zij, n, p, vij, etaij, v)
likP(beta2, Nij, Zij, n, p, vij, etaij, v)

beta2 = beta0 + solve(h1, g1)/100
lik(beta2, Nij, Zij, n, p)

admm_scad <- function(beta, Nij, Zij, n, p, lam, v, gam) {

  vij = matrix(0, nrow = n, ncol = n)
  etaij = matrix(0, nrow = n, ncol = n)
  dij = matrix(0, nrow = n, ncol = n)
  one = matrix(1, nrow = n, ncol = n)
  Iij = diag(1, n-1)

  beta0 = rep(0.1, n+p-1)
  beta = beta0[1:(n-1)]
  sum_beta = sum(beta)


  beta00 = beta0
  etaij0 = etaij
  diff1 = 3
  cont1 = 0

  while (cont1 < 100 && diff1 > 1e-3) {

    diff = 10
    cont = 1
    while (cont < 20 && diff > 1e-3) {
      # g1 = grad(lik, beta0, Nij = Nij, Zij = Zij, n = n, p = p)
      # g12 = v * ((n-1)*beta - sum_beta - colSums(etaij) - rowSums(etaij)) + colSums(vij) + rowSums(vij)
      # g12 = c(g12, 0)
      # g1 = g1 + g12

      # h1 = hessian(lik, beta0, Nij = Nij, Zij = Zij, n = n, p = p)
      # h1[1:(n-1), 1:(n-1)] = h1[1:(n-1), 1:(n-1)] + (n-1)*v*Iij - v*one

      # g13 = grad(likP, beta0, Nij = Nij, Zij = Zij, n = n, p = p, vij = vij, etaij = etaij, v = v)
      # h13 = hessian(likP, beta0, Nij = Nij, Zij = Zij, n = n, p = p, vij = vij, etaij = etaij, v = v)

      test = GradCN(beta0, Nij, Zij, n, p, vij, etaij, v, 1)

      g13 = test$g
      h13 = test$h

      beta1 = beta0 - solve(h13, g13)

      diff = sum(abs(beta1 - beta0))
      cont = cont + 1

      beta0 = beta1
    }

    beta = beta0[1:(n-1)]
    sum_beta = sum(beta)
    x = c(0, beta[1:(n-1)])

    for (i in 1:(n-1)) {
      for (j in (i+1):n) {
        dij[i, j] = x[i] - x[j] + vij[i, j]/v
        etaij[i, j] = stscad(dij[i, j], lam, v, gam, n)
        vij[i, j] = vij[i, j] + v*(x[i] - x[j] - etaij[i, j])
      }
    }

    diff1 = sum(abs(beta00 - beta0)) + sum(abs(etaij0 - etaij))
    cont1 = cont1 + 1

    beta00 = beta0
    etaij0 = etaij

    cat(round(beta, 2), "\n")
  }

  return(beta)
}


# Simulation 1 (Oracle) ---------------------------------------------------

sourceCpp(file = "Z:/network/CoxCNet/SimSet.cpp", verbose = TRUE, rebuild = TRUE)
sourceCpp(file = "Z:/network/CoxCNet/NT.cpp", verbose = TRUE, rebuild = TRUE)


# n = 45 -------------------------------------------------------------------

n = 45
p = 1
Nij = matrix(0, nrow = n, ncol = n)
Zij = matrix(rnorm(n*n), nrow = n, ncol = n)
M = matrix(0, nrow = (n+p-1), ncol = 4)

for (i in 1:14) {
  M[i, 1] = 1
}

for (i in 15:29) {
  M[i, 2] = 1
}

for (i in 30:44) {
  M[i, 3] = 1
}

M[45, 4] = 1

beta = rep(0, n + p - 1)
betax = rep(0, 4)

betaAll = matrix(0, 1000, 4)
betaAllz = matrix(0, 1000, 4)
betaAllse = matrix(0, 1000, 4)

vij = matrix(0, nrow = n, ncol = n)
etaij = matrix(0, nrow = n, ncol = n)
v = 2

for (z in 1:1000) {

  set.seed(z)
  Zij = matrix(rnorm(n*n), nrow = n, ncol = n)
  Zij = (Zij + t(Zij)) / 2
  trail_sim = SimSetC(n, 0.5, array(Zij, c(n, n, 1)))
  nn = length(trail_sim[[1]])

  Nij = matrix(0, nrow = n, ncol = n)
  for (i in 1:nn) {
    p1 = trail_sim[[1]][i]
    q1 = trail_sim[[2]][i]
    Nij[p1, q1] = Nij[p1, q1] + 1
  }

  beta0 = rep(0, 4)

  diff = 10
  cont = 1
  while (cont < 10 && diff > 1e-3) {

    beta01 = M %*% beta0

    test = GradCN(beta01, Nij, Zij, n, p, vij, etaij, v, 0)

    g2 = test$g
    h2 = test$h

    g22 = t(M) %*% g2
    h22 = t(M) %*% h2 %*% M

    beta1 = beta0 - solve(h22, g22)
    cont = cont + 1
    diff = sum(abs(beta1 - beta0))

    beta0 = beta1
  }

  betaAll[z, ] = beta0
  betaAllz[z, ] = sqrtm(h22) %*% (beta0 -  c(0, 0.5, 1, 0.2))
  betaAllse[z, ] = sqrt(diag(solve(h22)))
  cat(z, "\n")

}

apply(betaAll, 2, mean)
apply(betaAll, 2, sd)
apply(betaAllse, 2, mean)
1-apply(abs(betaAllz) > 1.96, 2, sum)/1000

hist(betaAllz[,4])


# n = 100 -------------------------------------------------------------------

n = 100
p = 1
Nij = matrix(0, nrow = n, ncol = n)
Zij = matrix(rnorm(n*n), nrow = n, ncol = n)
M = matrix(0, nrow = (n+p-1), ncol = 5)

for (i in 1:24) {
  M[i, 1] = 1
}

for (i in 25:49) {
  M[i, 2] = 1
}

for (i in 50:74) {
  M[i, 3] = 1
}

for (i in 75:99) {
  M[i, 4] = 1
}

M[100, 5] = 1

beta = rep(0, n + p - 1)
betax = rep(0, 5)

betaAll = matrix(0, 1000, 5)
betaAllz = matrix(0, 1000, 5)
betaAllse = matrix(0, 1000, 5)

vij = matrix(0, nrow = n, ncol = n)
etaij = matrix(0, nrow = n, ncol = n)
v = 2

for (z in 1:1000) {

  set.seed(z)
  Zij = matrix(rnorm(n*n), nrow = n, ncol = n)
  Zij = (Zij + t(Zij)) / 2
  trail_sim = SimSetC(n, 0.5, array(Zij, c(n, n, 1)))
  nn = length(trail_sim[[1]])

  Nij = matrix(0, nrow = n, ncol = n)
  for (i in 1:nn) {
    p1 = trail_sim[[1]][i]
    q1 = trail_sim[[2]][i]
    Nij[p1, q1] = Nij[p1, q1] + 1
  }

  beta0 = rep(0, 5)

  diff = 10
  cont = 1
  while (cont < 10 && diff > 1e-3) {

    beta01 = M %*% beta0

    test = GradCN(beta01, Nij, Zij, n, p, vij, etaij, v, 0)

    g2 = test$g
    h2 = test$h

    g22 = t(M) %*% g2
    h22 = t(M) %*% h2 %*% M

    beta1 = beta0 - solve(h22, g22)
    cont = cont + 1
    diff = sum(abs(beta1 - beta0))

    beta0 = beta1
  }

  betaAll[z, ] = beta0
  betaAllz[z, ] = sqrtm(h22) %*% (beta0 -  c(0, 0.3, 0.7, 1, 0.2))
  betaAllse[z, ] = sqrt(diag(solve(h22)))
  cat(z, "\n")

}

apply(betaAll, 2, mean)
apply(betaAll, 2, sd)
apply(betaAllse, 2, mean)
1-apply(abs(betaAllz) > 1.96, 2, sum)/1000


# n = 240 -------------------------------------------------------------------

n = 240
p = 1
Nij = matrix(0, nrow = n, ncol = n)
Zij = matrix(rnorm(n*n), nrow = n, ncol = n)
M = matrix(0, nrow = (n+p-1), ncol = 7)

for (i in 1:39) {
  M[i, 1] = 1
}

for (i in 40:79) {
  M[i, 2] = 1
}

for (i in 80:119) {
  M[i, 3] = 1
}

for (i in 120:159) {
  M[i, 4] = 1
}

for (i in 160:199) {
  M[i, 5] = 1
}

for (i in 200:239) {
  M[i, 6] = 1
}


M[240, 7] = 1

beta = rep(0, n + p - 1)
betax = rep(0, 7)

betaAll = matrix(0, 1000, 7)
betaAllz = matrix(0, 1000, 7)
betaAllse = matrix(0, 1000, 7)

vij = matrix(0, nrow = n, ncol = n)
etaij = matrix(0, nrow = n, ncol = n)
v = 2

for (z in 1:1000) {

  set.seed(z)
  Zij = matrix(rnorm(n*n), nrow = n, ncol = n)
  Zij = (Zij + t(Zij)) / 2
  trail_sim = SimSetC(n, 0.5, array(Zij, c(n, n, 1)))
  nn = length(trail_sim[[1]])

  Nij = matrix(0, nrow = n, ncol = n)
  for (i in 1:nn) {
    p1 = trail_sim[[1]][i]
    q1 = trail_sim[[2]][i]
    Nij[p1, q1] = Nij[p1, q1] + 1
  }

  beta0 = rep(0, 7)

  diff = 10
  cont = 1
  while (cont < 10 && diff > 1e-3) {

    beta01 = M %*% beta0

    test = GradCN(beta01, Nij, Zij, n, p, vij, etaij, v, 0)

    g2 = test$g
    h2 = test$h

    g22 = t(M) %*% g2
    h22 = t(M) %*% h2 %*% M

    beta1 = beta0 - solve(h22, g22)
    cont = cont + 1
    diff = sum(abs(beta1 - beta0))

    beta0 = beta1
  }

  betaAll[z, ] = beta0
  betaAllz[z, ] = sqrtm(h22) %*% (beta0 -  c(0, 0.2, 0.4, 0.6, 0.8, 1, 0.2))
  betaAllse[z, ] = sqrt(diag(solve(h22)))
  cat(z, "\n")

}

apply(betaAll, 2, mean)
apply(betaAll, 2, sd)
apply(betaAllse, 2, mean)
1-apply(abs(betaAllz) > 1.96, 2, sum)/1000


hist(betaAllz[,1], seq(-5, 5, 0.5), ylim = c(0, 0.45), freq = FALSE, main = "", xlab = TeX('Standardized $\\widehat{x}_1$'))
curve(dnorm(x), xlab = "", ylab = "", add = T, lwd = 2.0)

hist(betaAllz[,3], seq(-5, 5, 0.5), ylim = c(0, 0.45), freq = FALSE, main = "", xlab = TeX('Standardized $\\widehat{x}_3$'))
curve(dnorm(x), xlab = "", ylab = "", add = T, lwd = 2.0)

hist(betaAllz[,5], seq(-5, 5, 0.5), ylim = c(0, 0.45), freq = FALSE, main = "", xlab = TeX('Standardized $\\widehat{x}_5$'))
curve(dnorm(x), xlab = "", ylab = "", add = T, lwd = 2.0)

hist(betaAllz[,6], seq(-5, 5, 0.5), ylim = c(0, 0.45), freq = FALSE, main = "", xlab = TeX('Standardized $\\widehat{\\gamma}$(t)'))
curve(dnorm(x), xlab = "", ylab = "", add = T, lwd = 2.0)



# simulation SCAD ---------------------------------------------------------

n = 45
p = 1
Nij = matrix(0, nrow = n, ncol = n)
Zij = matrix(rnorm(n*n), nrow = n, ncol = n)
M = matrix(0, nrow = (n+p-1), ncol = 4)

for (i in 1:14) {
  M[i, 1] = 1
}

for (i in 15:29) {
  M[i, 2] = 1
}

for (i in 30:44) {
  M[i, 3] = 1
}

M[45, 4] = 1

beta = rep(0, n + p - 1)
betax = rep(0, 4)

betaAll = matrix(0, 1000, 4)
betaAllz = matrix(0, 1000, 4)
betaAllse = matrix(0, 1000, 4)

vij = matrix(0, nrow = n, ncol = n)
etaij = matrix(0, nrow = n, ncol = n)
v = 2

set.seed(3)
Zij = matrix(rnorm(n*n), nrow = n, ncol = n)
Zij = (Zij + t(Zij)) / 2
trail_sim = SimSetC(n, 0.5, array(Zij, c(n, n, 1)))
nn = length(trail_sim[[1]])

Nij = matrix(0, nrow = n, ncol = n)
for (i in 1:nn) {
  p1 = trail_sim[[1]][i]
  q1 = trail_sim[[2]][i]
  Nij[p1, q1] = Nij[p1, q1] + 1
}

tune = exp(seq(-3, 1, 0.25))
scadAll = matrix(0, nrow = n - 1, ncol = length(tune))

# beta = M %*% c(0, 0.3, 0.7, 1, 0.2)

for (i in 1:length(tune)) {
  test = admm_scad(beta, Nij, Zij, n, p, tune[i], 10, 2)
  scadAll[, i] = test
}

xaxis = log(tune)
yaxis = t(scadAll)
matplot(xaxis, yaxis, lty = 1, xlab = "log(tune)", ylab = "coefficients", type = "l")

plot(round(scadAll[, 4], 2))


