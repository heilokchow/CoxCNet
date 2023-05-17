// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <cmath>

#define EIGEN_NO_DEBUG
#define EIGEN_DONT_PARALLELIZE

using namespace Rcpp;
using Rcpp::NumericMatrix;
using Rcpp::NumericVector;
using Rcpp::List;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;

// [[Rcpp::export]]
double LikCEN(NumericVector beta, NumericMatrix Trail, NumericMatrix cen, NumericVector Nij, NumericVector Zij, int nn, int n, int p, NumericVector vij, NumericVector etaij, int v, bool scad) {
  VectorXd x = VectorXd::Zero(n);
  VectorXd y = VectorXd::Zero(p);

  VectorXd Ni = VectorXd::Zero(n);
  MatrixXd Eij = MatrixXd::Zero(n, n);
  VectorXd Ei = VectorXd::Zero(n);

  VectorXd zN = VectorXd::Zero(p);
  VectorXd zE = VectorXd::Zero(p);
  MatrixXd zEi = MatrixXd::Zero(p, n);

  VectorXd Excl = VectorXd::Zero(n);

  double N0(0), E0(0);

  VectorXd Gi = VectorXd::Zero(n + p - 1);
  MatrixXd Hij = MatrixXd::Zero(n + p - 1,  n + p - 1);

  for (int i = 0; i < (n - 1); i++) {
    x[i+1] = beta[i];
  }

  for (int z = 0; z < p; z++) {
    y[z] = beta[n + z - 1];
  }

  for (int i = 0; i < (n - 1); i++) {
    for (int j = (i+1); j < n; j++) {
      Eij(i, j) = x[i] + x[j];
      for (int z = 0; z < p; z++) {
        Eij(i, j) = Eij(i, j) + Zij[z*n*n + i*n + j] * y[z];
      }
      Eij(i, j) = exp(Eij(i, j));
    }
  }

  int cont = 0;
  int contall = cen.nrow();
  double l = 0;

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      E0 += Eij(i, j) + Eij(j, i);
    }
  }

  for (int i = 0; i < nn; i++) {

    int p1 = Trail(i, 0);
    int q1 = Trail(i, 1);
    if (p1 < q1) {
      l -= log(Eij(p1, q1));
    }
    if (p1 >= q1) {
      l -= log(Eij(q1, p1));
    }
    l += log(E0);

    if (cont < contall && Trail(i, 2) > cen(cont, 0)) {
      Excl[cen(cont, 1)] = 1;
      cont++;
      E0 = 0;
      for (int z1 = 0; z1 < n; z1++) {
        for (int z2 = 0; z2 < n; z2++) {
          if (Excl[z1] == 0 && Excl[z2] == 0) {
            E0 += Eij(z1, z2) + Eij(z2, z1);
          }
        }
      }
    }
  }

  return l;
}
