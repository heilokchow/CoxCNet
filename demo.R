library(Rcpp)
library(RcppEigen)
library(numDeriv)
library(expm)
library(latex2exp)

Sys.setenv("PKG_CPPFLAGS" = "-march=native")
sourceCpp(file = "Z:/network/CoxCNet/SimSet.cpp", verbose = TRUE, rebuild = TRUE)
sourceCpp(file = "Z:/network/CoxCNet/NT.cpp", verbose = TRUE, rebuild = TRUE)
sourceCpp(file = "Z:/network/CoxCNet/LikCn.cpp", verbose = TRUE, rebuild = TRUE)

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
        Pij[i, j] = Pij[i, j] + Zij[i, j, z] * y[z]
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

stmcp <- function(muij, lam, v, gam, n) {
  if (abs(muij) <= gam*lam) {
    return(st(muij, lam/v) / (1-1/(gam*v)))
  } else {
    return(muij)
  }
}

pscad <- function(x, lam, a) {
  if (abs(x) <= lam) {
    return(lam*abs(x))
  } else if (abs(x) > lam && abs(x) <= a*lam) {
    return((2*a*lam*abs(x)-x^2-lam^2)/(2*(a-1)))
  } else {
    return(lam^2*(a+1)/2)
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
        Pij[i, j] = Pij[i, j] + Zij[i, j, z] * y[z]
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

likQ <- function(beta, Nij, Zij, n, p, lam, gam) {

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
  for (i in 2:(n-1)) {
    for (j in (i+1):(n)) {
      l = l + pscad(x[i] - x[j], lam, gam) * n
    }
  }
  return(l)
}

Qpen <- function(beta, n, lam, gam) {
  l = 0
  x = c(0, beta[1:(n-1)])
  for (i in 2:(n-1)) {
    for (j in (i+1):(n)) {
      l = l + pscad(x[i] - x[j], lam, gam)
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

  while (cont1 < 2000 && diff1 > 1e-3) {

    diff = 10
    cont = 1
    while (cont < 100 && diff > 1e-3) {
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

      step = solve(h13, g13)
      l1 = LikCN(beta0, Nij, zij, n, p, vij, etaij, v, 1)[[1]]
      bet = 1
      for (i in 1:3) {
        l2 = LikCN(beta0 - bet * step, Nij, zij, n, p, vij, etaij, v, 1)[[1]]
        if (l2 > l1) {
          bet = bet / 2
        } else {
          break
        }
      }
      beta1 = beta0 - bet*step

      diff = sum(abs(beta1 - beta0))
      cont = cont + 1

      beta0 = beta1
    }

    beta = beta0[1:(n-1)]
    sum_beta = sum(beta)
    x = c(0, beta[1:(n-1)])

    for (i in 2:(n-1)) {
      for (j in (i+1):n) {
        dij[i, j] = x[i] - x[j] + vij[i, j]/v
        # etaij[i, j] = stmcp(dij[i, j], lam, v, gam, n)
        etaij[i, j] = stscad(dij[i, j], lam, v, gam, n)
        # etaij[i, j] = st(dij[i, j], lam / gam)
        vij[i, j] = vij[i, j] + v*(x[i] - x[j] - etaij[i, j])
      }
    }

    diff1 = sum(abs(beta00 - beta0)) + sum(abs(etaij0 - etaij))
    cont1 = cont1 + 1

    beta00 = beta0
    etaij0 = etaij

    cat(diff1, "\n")
  }

  return(list(beta = beta0, vij = vij, etaij = etaij))
}

ngroup <- function(beta) {
  l = length(beta)
  s = beta[1]
  cont = 1
  for (i in 2:l) {
    flag = 1
    for (j in 1:cont) {
      if (abs(beta[i] - s[j]) < 0.0001) {
        flag = 0
        break
      }
    }
    if (flag == 1) {
      s = c(s, beta[i])
      cont = cont + 1
    }
  }
  return(s)
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


vij = matrix(0, nrow = n, ncol = n)
etaij = matrix(0, nrow = n, ncol = n)
v = 2


# START HERE
sourceCpp(file = "Z:/network/CoxCNet/SimSet.cpp", verbose = TRUE, rebuild = TRUE)

betamin = matrix(0, 100, n)
bicmin = rep(0, 100)
gmin = rep(0, 100)

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

  # MCP / SCAD
  tune = exp(seq(-3.5, -0.5, 0.25))

  # L1
  # tune = exp(seq(-5, -0.1, 0.25))
  scadAll = matrix(0, nrow = n - 1, ncol = length(tune))
  scadAll1 = matrix(0, nrow = n, ncol = length(tune))
  bicall = c()
  kgall = c()
  likall = c()
  n0 = sum(Nij)

  for (i in 1:length(tune)) {
    test = admm_scad(beta, Nij, Zij, n, p, tune[i], 10, 3.7)
    scadAll[, i] = test$beta[1:(n-1)]
    scadAll1[, i] = test$beta

    kg = length(ngroup(test$beta[1:(n-1)]))
    kgall = c(kgall, kg)
    likall = c(likall, 2*lik(test$beta, Nij, Zij, n, p))
    bicall = c(bicall, 2*lik(test$beta, Nij, Zij, n, p) + log(log(n))*(kg + p)*log(n0))
  }

  cs = which.min(bicall)
  betamin[z, ] = scadAll1[, cs]
  bicmin[z] = bicall[cs]
  gmin[z] = kgall[cs]

  cat(z, gmin[z], "\n")
}


apply(betamin, 2, mean)[c(1, 16, 31, 45)]
apply(betamin, 2, sd)[c(1, 16, 31, 45)]

mean(bicmin)
mean(gmin)
sd(gmin)
sum(gmin == 3)

apply(betaAll, 2, sd)
apply(betaAllse, 2, mean)
1-apply(abs(betaAllz) > 1.96, 2, sum)/1000


betax = c(sum(test$beta[1:14])/14, sum(test$beta[15:29])/15, sum(test$beta[30:44])/15, test$beta[45])
betaT = M %*% betax
betaT = M %*% beta0

lik(test$beta, Nij, Zij, n, p)
lik(betaT, Nij, Zij, n, p)

likQ(test$beta, Nij, Zij, n, p, tune[i], 3.7)
likQ(betaT, Nij, Zij, n, p, tune[i], 3.7)
likQ(beta, Nij, Zij, n, p, tune[i], 3.7)

r = rep(0, n + p - 1)
r[1] = 0.1
likQ(betaT, Nij, Zij, n, p, 1, 3.7)
likQ(betaT - r, Nij, Zij, n, p, 1, 3.7)

Qpen(beta, n, tune[i], 3.7)
Qpen(betaT, n, tune[i], 3.7)
Qpen(test$beta, n, tune[i], 3.7)

glik1 = grad(lik, betaT, Nij = Nij, Zij = Zij, n = n, p = p)

glikQ1 = grad(likQ, betaT, Nij = Nij, Zij = Zij, n = n, p = p, lam = 1, gam = 3.7, side = rep(1, n + p - 1))
glikQ2 = grad(likQ, betaT, Nij = Nij, Zij = Zij, n = n, p = p, lam = 1, gam = 1.2, side = rep(-1, n + p - 1))

(glikQ1 * glikQ2) < 0

gQpen1 = n*grad(Qpen, betaT, n = n, lam = 1, gam = 3.7, side = rep(1, n + p - 1))
gQpen2 = n*grad(Qpen, betaT, n = n, lam = 1, gam = 3.7, side = rep(-1, n + p - 1))

likP(test$beta, Nij, Zij, n, p, vij, etaij, v)
likP(betaT, Nij, Zij, n, p, vij, etaij, v)

xaxis = log(tune)
yaxis = t(scadAll)
matplot(xaxis, yaxis, lty = 1, xlab = "log(tune)", ylab = "coefficients", type = "l")

plot(bicall)
plot(kgall)


# Censoring ---------------------------------------------------------------

sourceCpp(file = "Z:/network/CoxCNet/LikCen.cpp", verbose = TRUE, rebuild = TRUE)
sourceCpp(file = "Z:/network/CoxCNet/SimSet.cpp", verbose = TRUE, rebuild = TRUE)

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
betax = rep(0, 5)

betaAll = matrix(0, 1000, 5)
betaAllz = matrix(0, 1000, 5)
betaAllse = matrix(0, 1000, 5)

vij = matrix(0, nrow = n, ncol = n)
etaij = matrix(0, nrow = n, ncol = n)
v = 2

ceninfo = rep(0, n)
for (i in 2:n) {
  k = runif(1, 0, 1)
  if (k < 0.8) {
    ceninfo[i] = 1 - k
  }
}
ceninfo[1] = 0.7

cenord = order(ceninfo, decreasing = TRUE)
n2 = sum(ceninfo > 0)

c1 = ceninfo[cenord][n2:1]
c2 = seq(0, n-1, 1)[cenord][n2:1]
cenall = matrix(0, n2, 2)
cenall[, 1] = c1
cenall[, 2] = c2

set.seed(z)
Zij = matrix(rnorm(n*n), nrow = n, ncol = n)
Zij = (Zij + t(Zij)) / 2
trail_sim = SimSetC(n, 0.5, array(Zij, c(n, n, 1)))

nn = length(trail_sim[[1]])
trail_all = data.frame(x = trail_sim[[1]], y = trail_sim[[2]], t = trail_sim[[3]])
tsort = order(trail_all[[3]])
trail_all = trail_all[tsort,]
trail_all[,1:2] = trail_all[,1:2] - 1

for (i in 1:n2) {
  t1 = which(trail_all[, 1] == cenall[i, 2])
  t11 = which(trail_all[t1, 3] > cenall[i, 1])
  if (length(t11) > 0) {
    trail_all = trail_all[-t1[t11], ]
  }

  t2 = which(trail_all[, 2] == cenall[i, 2])
  t22 = which(trail_all[t2, 3] > cenall[i, 1])
  if (length(t22) > 0) {
    trail_all = trail_all[-t2[t22], ]
  }
}

nn = nrow(trail_all)


LikCEN(beta01, as.matrix(trail_all), cenall, Nij, Zij, nn, n, p, vij, etaij, 2, 0)

beta0 = rep(0, 4)

diff = 10
cont = 1
while (cont < 10 && diff > 1e-3) {

  beta01 = M %*% beta0

  test = GradCN(beta01, Nij, Zij, n, p, vij, etaij, v, 0)

  g2 = grad(LikCEN, beta01, Trail = as.matrix(trail_all), cen = cenall, Nij = Nij,
            Zij = Zij, nn = nn, n = n, p = p, vij = vij, etaij = etaij, v = 2, scad = 0)
  h2 = hessian(LikCEN, beta01, Trail = as.matrix(trail_all), cen = cenall, Nij = Nij,
            Zij = Zij, nn = nn, n = n, p = p, vij = vij, etaij = etaij, v = 2, scad = 0)

  g22 = t(M) %*% g2
  h22 = t(M) %*% h2 %*% M

  beta1 = beta0 - solve(h22, g22)
  cont = cont + 1
  diff = sum(abs(beta1 - beta0))

  beta0 = beta1
  cat(cont)
}

beta02 = beta0
beta02[1:3] = beta02[1:3] - 0.11
beta12 = M %*% beta02
LikCEN(beta01, as.matrix(trail_all), cenall, Nij, Zij, nn, n, p, vij, etaij, 2, 0)
LikCEN(beta12, as.matrix(trail_all), cenall, Nij, Zij, nn, n, p, vij, etaij, 2, 0)



# Enron email -------------------------------------------------------------

library(readxl)
enronSummary <- read_excel("Z:/dynamic network/dynamic/enron/enronSummary.xlsx")
employess <- read_excel("Z:/dynamic network/dynamic/enron/employess.xlsx")

enronSummary <- as.data.frame(enronSummary)
employess <- as.data.frame(employess)

n = max(enronSummary[,1])

p = 4

zij = array(0, c(n, n, p))


# Define Homo -------------------------------------------------------------


# 20 homophily Effects

# 1

for (i in 1:n) {
  for (j in 1:n) {
    if (i != j) {
      zij[i, j, 1] = (employess[i, 3] == "Legal" && employess[j, 3] == "Legal")
    }
  }
}

for (i in 1:n) {
  for (j in 1:n) {
    if (i != j) {
      zij[i, j, 2] = (employess[i, 3] == "Legal" && employess[j, 3] == "Trading") || (employess[j, 3] == "Legal" && employess[i, 3] == "Trading")
    }
  }
}

for (i in 1:n) {
  for (j in 1:n) {
    if (i != j) {
      zij[i, j, 3] = (employess[i, 3] == "Trading" && employess[j, 3] == "Trading")
    }
  }
}

for (i in 1:n) {
  for (j in 1:n) {
    if (i != j) {
      zij[i, j, 4] = (employess[i, 7] == "Junior" && employess[j, 7] == "Junior")
    }
  }
}

trail = enronSummary[,c(1,2,9)]
nn = nrow(trail)
trail = trail[-nn, ]
nn = nrow(trail)

flag = 1
start = 1
excl = c()
for (i in 2:nn) {
  if (trail[i, 1] == trail[i-1, 1] && trail[i, 3] == trail[i-1, 3]) {
    flag = flag + 1
  } else {
    if (flag > 5) {
      excl = c(excl, start:(start+flag-1))
    }
    start = i
    flag = 1
  }
}

trail = trail[-excl,]
nn = nrow(trail)

# Remove No communication individuals
n = length(unique(trail[,1]))
nid = sort(unique(trail[,1]))

map_trail = rep(0, 156)
for (i in 1:length(nid)) {
  map_trail[nid[i]] = nid[i]
}

map_trail1 = map_trail
count = 0
for (i in 2:length(map_trail1)) {
  if (map_trail1[i-1] == 0) {
    count = count + 1
  }
  if (map_trail1[i] != 0) {
    map_trail1[i] = map_trail1[i] - count
  }
}

del = which(map_trail1 == 0)

del_trail = c()
for (i in 1:4) {
  test = which(trail[, 2] == del[i])
  if (length(test) != 0) {
    del_trail = c(del_trail, test)
  }
}

trail = trail[-del_trail, ]
nn = nrow(trail)

for (i in 1:nn) {
  trail[i, 1] = map_trail1[trail[i, 1]]
  trail[i, 2] = map_trail1[trail[i, 2]]
}

zij = zij[-del, -del, ]

sum(sort(unique(trail[,1])) == sort(unique(trail[,2])))

Nij = matrix(0, n, n)
for (i in 1:nn) {
  p1 = trail[i, 1]
  q1 = trail[i, 2]
  if (p1 < q1) {
    Nij[p1, q1] = Nij[p1, q1] + 1
  } else {
    Nij[q1, p1] = Nij[q1, p1] + 1
  }
}

n0 = sum(Nij)

beta = rep(0, n + p - 1)
tune = exp(seq(-6, -0.5, 0.5))
scadAll = matrix(0, nrow = n - 1, ncol = length(tune))
scadAll1 = matrix(0, nrow = n + p - 1, ncol = length(tune))
bicall = c()
kgall = c()
likall = c()

for (i in 1:length(tune)) {
  test = admm_scad(beta, Nij, zij, n, p, tune[i], 10, 3.7)
  scadAll[, i] = test$beta[1:(n-1)]
  scadAll1[, i] = test$beta

  kg = length(ngroup(test$beta[1:(n-1)]))
  kgall = c(kgall, kg)
  likall = c(likall, 2*lik(test$beta, Nij, zij, n, p))
  bicall = c(bicall, 2*lik(test$beta, Nij, zij, n, p) + log(log(n))*(kg + p)*log(n0))
}


bicall = c()
kgall = c()
l0ikall = c()

for (i in 1:length(tune)) {
  kg = length(ngroup(scadAll[, i]))
  kgall = c(kgall, kg)
  likall = c(likall, 2*lik(scadAll1[, i], Nij, zij, n, p))
  bicall = c(bicall, 2*lik(scadAll1[, i], Nij, zij, n, p) + log(log(n))*(kg + p)*log(n0))
}

xaxis = log(tune)
yaxis = t(scadAll)
matplot(xaxis, yaxis, lty = 1, xlab = "log(tune)", ylab = "coefficients", type = "l")
abline(v = log(tune[3]), col="red", lwd=3, lty=2)

beta0 = rep(0, n+p-1)
beta = beta0[1:(n-1)]
sum_beta = sum(beta)

beta00 = beta0
etaij0 = etaij

diff = 10
cont = 1
while (cont < 100 && diff > 1e-3) {

  test = GradCN(beta0, Nij, zij, n, p, vij, etaij, v, 0)

  g13 = test$g
  h13 = test$h

  step = solve(h13, g13)
  l1 = LikCN(beta0, Nij, zij, n, p, vij, etaij, v, 0)[[1]]
  bet = 1
  for (i in 1:3) {
    l2 = LikCN(beta0 - bet * step, Nij, zij, n, p, vij, etaij, v, 0)[[1]]
    if (l2 > l1) {
      bet = bet / 2
    } else {
      break
    }
  }
  beta1 = beta0 - bet*step


  diff = sum(abs(beta1 - beta0))
  cont = cont + 1

  beta0 = beta1
}

sqrt(diag(solve(h13)))

plot((colSums(Nij) + rowSums(Nij))[2:n], exp(test$beta[1:(n-1)]))

enron_oracle = c(0, beta0[1:(n-1)])
enron_scad = c(0, scadAll[1:(n-1), 3])

ord_enron = order(enron_scad, decreasing = TRUE)
enron_oracle[ord_enron]
enron_scad[ord_enron]

(colSums(Nij) + rowSums(Nij))[1:n][ord_enron]

out_data = data.frame(rank = seq(1, n, 1), id = seq(1, n, 1)[ord_enron], edges = (colSums(Nij) + rowSums(Nij))[1:n][ord_enron],
                      oracle = enron_oracle[ord_enron],
                      sd = sqrt(diag(solve(h13)))[ord_enron], scad = enron_scad[ord_enron])

# Baseline

basehaz <- function(beta0, trailt, n) {

  Pij = matrix(0, n, n)
  x = c(0, beta0[1:(n-1)])
  y = beta0[n:(n+p-1)]
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      Pij[i, j] = x[i] + x[j]
      for (z in 1:p) {
        Pij[i, j] = Pij[i, j] + zij[i, j, z] * y[z]
      }
      Pij[i, j] = exp(Pij[i, j])
      Pij[j, i] = Pij[i, j]
    }
  }

  P0 = sum(Pij) / 2

  return(trailt / P0)

}

final_grad = grad(basehaz, beta0, trailt = sum(trail[, 3] < 1), n = n)
sd_base = final_grad %*% solve(h13) %*% final_grad


hazardtime = seq(0.01, 1, 0.01)
hazardenron = rep(0, 100)

for (i in 1:100) {
  t1 = hazardtime[i]
  hazardenron[i] = sum(trail[, 3] < t1) / P0
}

plot_base_enron = data.frame(x = hazardtime, y = hazardenron,
                             yl = hazardenron - 1.96 * hazardenron/hazardenron[100] * sd_base[[1]],
                             yu = hazardenron + 1.96 * hazardenron/hazardenron[100] * sd_base[[1]])


ggplot(plot_base_enron, aes(x = x, y = y)) +
  geom_line(size = 0.7) +
  geom_ribbon(aes(ymin = yl, ymax = yu), alpha = 0.2, linetype = "dashed", size = 0.7) +
  xlab(expression(italic("t"))) +
  ylab("Cumulative Baseline Hazard") +
  theme_classic() +
  theme(legend.position = "none")



# High School -------------------------------------------------------------


library(readxl)

highschool <- as.data.frame(highschool)

nn = nrow(highschool)
n = max(highschool[,2])
students <- seq(1, n, 1)
studentsid1 <- rep(0, n)
studentsid2 <- rep(0, n)

Nij = matrix(0, n, n)

for (i in 1:nn) {
  p1 = highschool[i, 2]
  q1 = highschool[i, 3]
  if (p1 < q1) {
    Nij[p1, q1] = Nij[p1, q1] + 1
  } else {
    Nij[q1, p1] = Nij[q1, p1] + 1
  }
}

for (i in 1:n) {
  temp = which(highschool[, 2] == students[i])
  if (length(temp) > 0) {
    studentsid1[i] = temp[1]
  } else {
    temp = which(highschool[, 3] == students[i])
    studentsid2[i] = temp[1]
  }
}

employess = rep(0, n)

for (i in 1:n) {
  if (studentsid1[i] > 0) {
    employess[i] = highschool[studentsid1[i], 4]
  } else {
    employess[i] = highschool[studentsid2[i], 5]
  }
}


p = 6

zij = array(0, c(n, n, p))


# Define Homo -------------------------------------------------------------


# 20 homophily Effects

# 1

for (i in 1:n) {
  for (j in 1:n) {
    if (i != j) {
      zij[i, j, 1] = (employess[i] == "PC" && employess[j] == "PC")
    }
  }
}

for (i in 1:n) {
  for (j in 1:n) {
    if (i != j) {
      zij[i, j, 2] = (employess[i] == "PC" && employess[j] == "PC*") || (employess[j] == "PC" && employess[i] == "PC*")
    }
  }
}

for (i in 1:n) {
  for (j in 1:n) {
    if (i != j) {
      zij[i, j, 3] = (employess[i] == "PC" && employess[j] == "PSI*") || (employess[j] == "PC" && employess[i] == "PSI*")
    }
  }
}

for (i in 1:n) {
  for (j in 1:n) {
    if (i != j) {
      zij[i, j, 4] = (employess[i] == "PC*" && employess[j] == "PC*")
    }
  }
}

for (i in 1:n) {
  for (j in 1:n) {
    if (i != j) {
      zij[i, j, 5] = (employess[i] == "PC*" && employess[j] == "PSI*") || (employess[j] == "PC*" && employess[i] == "PSI*")
    }
  }
}

for (i in 1:n) {
  for (j in 1:n) {
    if (i != j) {
      zij[i, j, 6] = (employess[i] == "PSI*" && employess[j] == "PSI*")
    }
  }
}

trail = highschool[,c(2,3,1)]
trail[,3] = (trail[,3] - trail[1,3]) / (trail[nn,3] - trail[1,3])

n0 = sum(Nij)

beta = rep(0, n + p - 1)
tune = exp(seq(-6, -0.5, 0.5))
scadAll = matrix(0, nrow = n - 1, ncol = length(tune))
scadAll1 = matrix(0, nrow = n + p - 1, ncol = length(tune))
bicall = c()
kgall = c()
likall = c()

for (i in 1:length(tune)) {
  test = admm_scad(beta, Nij, zij, n, p, tune[i], 10, 3.7)
  scadAll[, i] = test$beta[1:(n-1)]
  scadAll1[, i] = test$beta

  kg = length(ngroup(test$beta[1:(n-1)]))
  kgall = c(kgall, kg)
  likall = c(likall, 2*lik(test$beta, Nij, zij, n, p))
  bicall = c(bicall, 2*lik(test$beta, Nij, zij, n, p) + log(log(n))*(kg + p)*log(n0))
}


bicall = c()
kgall = c()
l0ikall = c()

for (i in 1:length(tune)) {
  kg = length(ngroup(scadAll[, i]))
  kgall = c(kgall, kg)
  likall = c(likall, 2*lik(scadAll1[, i], Nij, zij, n, p))
  bicall = c(bicall, 2*lik(scadAll1[, i], Nij, zij, n, p) + log(log(n))*(kg + p)*log(n0))
}

xaxis = log(tune)
yaxis = t(scadAll)
matplot(xaxis, yaxis, lty = 1, xlab = "log(tune)", ylab = "coefficients", type = "l")
abline(v = log(tune[5]), col="red", lwd=3, lty=2)

plot(xaxis, bicall, type = "b", xlab = "log(tune)", ylab = "BIC")
abline(v = log(tune[5]), col="red", lwd=3, lty=2)


beta0 = rep(0, n+p-1)
beta = beta0[1:(n-1)]
sum_beta = sum(beta)

beta00 = beta0
etaij0 = etaij

diff = 10
cont = 1
while (cont < 100 && diff > 1e-3) {

  test = GradCN(beta0, Nij, zij, n, p, vij, etaij, v, 0)

  g13 = test$g
  h13 = test$h

  step = solve(h13, g13)
  l1 = LikCN(beta0, Nij, zij, n, p, vij, etaij, v, 0)[[1]]
  bet = 1
  for (i in 1:3) {
    l2 = LikCN(beta0 - bet * step, Nij, zij, n, p, vij, etaij, v, 0)[[1]]
    if (l2 > l1) {
      bet = bet / 2
    } else {
      break
    }
  }
  beta1 = beta0 - bet*step


  diff = sum(abs(beta1 - beta0))
  cont = cont + 1

  beta0 = beta1
}

sqrt(diag(solve(h13)))

plot((colSums(Nij) + rowSums(Nij))[2:n], exp(test$beta[1:(n-1)]))

student_oracle = c(0, beta0[1:(n-1)])
student_scad = c(0, scadAll[1:(n-1), 5])

ord_student = order(student_scad, decreasing = TRUE)
student_oracle[ord_student]
student_scad[ord_student]

(colSums(Nij) + rowSums(Nij))[1:n][ord_student]

out_data = data.frame(rank = seq(1, n, 1), id = seq(1, n, 1)[ord_student], edges = (colSums(Nij) + rowSums(Nij))[1:n][ord_student],
                      oracle = student_oracle[ord_student],
                      sd = sqrt(diag(solve(h13)))[ord_student], scad = student_scad[ord_student])

scadAll1[n:(n+5), 5]

# Baseline

Pij = matrix(0, n, n)
x = c(0, beta0[1:(n-1)])
y = beta0[n:(n+p-1)]
for (i in 1:(n-1)) {
  for (j in (i+1):n) {
    Pij[i, j] = x[i] + x[j]
    for (z in 1:p) {
      Pij[i, j] = Pij[i, j] + zij[i, j, z] * y[z]
    }
    Pij[i, j] = exp(Pij[i, j])
    Pij[j, i] = Pij[i, j]
  }
}

P0 = sum(Pij) / 2

basehaz <- function(beta0, trailt, n) {

  Pij = matrix(0, n, n)
  x = c(0, beta0[1:(n-1)])
  y = beta0[n:(n+p-1)]
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      Pij[i, j] = x[i] + x[j]
      for (z in 1:p) {
        Pij[i, j] = Pij[i, j] + zij[i, j, z] * y[z]
      }
      Pij[i, j] = exp(Pij[i, j])
      Pij[j, i] = Pij[i, j]
    }
  }

  P0 = sum(Pij) / 2

  return(trailt / P0)

}

final_grad = grad(basehaz, beta0, trailt = sum(trail[, 3] < 1), n = n)
sd_base = final_grad %*% solve(h13) %*% final_grad
basehaz(beta0, sum(trail[, 3] < 1), n)

hazardtime = seq(0.01, 1, 0.01)
hazardenron = rep(0, 100)

for (i in 1:100) {
  t1 = hazardtime[i]
  hazardenron[i] = sum(trail[, 3] < t1) / P0
}

plot_base_enron = data.frame(x = hazardtime, y = hazardenron,
                             yl = hazardenron - 1.96 * hazardenron/hazardenron[100] * sd_base[[1]],
                             yu = hazardenron + 1.96 * hazardenron/hazardenron[100] * sd_base[[1]])


ggplot(plot_base_enron, aes(x = x, y = y)) +
  geom_line(size = 0.7) +
  geom_ribbon(aes(ymin = yl, ymax = yu), alpha = 0.2, linetype = "dashed", size = 0.7) +
  xlab(expression(italic("t"))) +
  ylab("Cumulative Baseline Hazard") +
  theme_classic() +
  theme(legend.position = "none")

