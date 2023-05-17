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
List LikCN(NumericVector beta, NumericVector Nij, NumericVector Zij, int n, int p, NumericVector vij, NumericVector etaij, int v, bool scad) {
  VectorXd x = VectorXd::Zero(n);
  VectorXd y = VectorXd::Zero(p);

  VectorXd Ni = VectorXd::Zero(n);
  MatrixXd Eij = MatrixXd::Zero(n, n);
  VectorXd Ei = VectorXd::Zero(n);

  VectorXd zN = VectorXd::Zero(p);
  VectorXd zE = VectorXd::Zero(p);
  MatrixXd zEi = MatrixXd::Zero(p, n);


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

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      Ei[i] += Eij(i, j) + Eij(j, i);
      Ni[i] += Nij[i*n + j] + Nij[j*n + i];
    }
    E0 += Ei[i];
    N0 += Ni[i];
  }

  E0 = E0 / 2;
  N0 = N0 / 2;

  double l = 0;
  for (int i = 0; i < (n - 1); i++) {
    for (int j = (i+1); j < n; j++) {
      l += Nij[j*n + i] * log(Eij(i, j));
    }
  }
  l -= N0 * log(E0);

  l = -l;
  for (int i = 0; i < (n - 1); i++) {
    for (int j = (i+1); j < n; j++) {
      l = l + vij[j *n + i] * (x[i] - x[j] - etaij[j *n + i]) + v/2 * (x[i] - x[j] - etaij[j *n + i]) * (x[i] - x[j] - etaij[j *n + i]);
    }
  }

  return List::create(_["Lik"] = l);
}
