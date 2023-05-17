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
List GradCN(NumericVector beta, NumericVector Nij, NumericVector Zij, int n, int p, NumericVector vij, NumericVector etaij, int v, bool scad) {
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

  for (int i = 0; i < (n - 1); i++) {
    Gi[i] = N0*Ei[i+1]/E0 - Ni[i+1];
  }

  for (int z = 0; z < p; z++) {
    for (int i = 0; i < (n-1); i++) {
      for (int j = (i+1); j < n; j++) {
        zN[z] += Zij[z*n*n + i*n + j] * Nij[j*n + i];
        zE[z] += Zij[z*n*n + i*n + j] * Eij(i, j);
      }
    }
    Gi[n - 1 + z] = N0*zE[z]/E0 - zN[z];
  }

  if (scad) {
    // for (int i = 0; i < (n - 1); i++) {
    //   int i1 = i + 1;
    //   Gi[i] -=  vij[i1 * n] + v * (0 - x[i1] - etaij[i1 * n]);
    // }

    for (int i = 0; i < (n - 2); i++) {
      for (int j = (i+1); j < (n - 1); j++) {
        int i1 = i + 1;
        int j1 = j + 1;
        Gi[i] +=  n * vij[j1 * n + i1] + n * v * (x[i1] - x[j1] - etaij[j1 * n + i1]);
        Gi[j] -=  n * vij[j1 * n + i1] + n * v * (x[i1] - x[j1] - etaij[j1 * n + i1]);
      }
    }
  }


  for (int z = 0; z < p; z++) {
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        zEi(z, i) +=  Zij[z*n*n + i*n + j] * (Eij(i, j) + Eij(j, i));
      }
    }
  }

  for (int i = 0; i < (n - 1); i++) {
    for (int j = i; j < (n - 1); j++) {
      int i1 = i + 1;
      int j1 = j + 1;
      if (i == j) {
        Hij(i, j) = N0 * (Ei[i1] * E0 - Ei[i1] * Ei[i1]) / (E0 * E0);
      }

      if (i != j) {
        Hij(i, j) = N0 * (Eij(i1, j1) * E0 - Ei[i1] * Ei[j1]) / (E0 * E0);
        Hij(j, i) = Hij(i, j);
      }
    }
  }

  if (scad) {
    // for (int i = 0; i < (n - 1); i++) {
    //   Hij(i, i) += v;
    // }

    for (int i = 0; i < (n - 2); i++) {
      for (int j = (i+1); j < (n - 1); j++) {
        Hij(i, i) += n * v;
        Hij(i, j) -= n * v;
        Hij(j, i) -= n * v;
        Hij(j, j) += n * v;
      }
    }
  }

  for (int i = 0; i < (n - 1); i++) {
    int i1 = i + 1;
    for (int z = 0; z < p; z++) {
      Hij(i, n - 1 + z) = N0 * (E0 * zEi(z, i1) - Ei[i1] * zE[z]) / (E0 * E0);
      Hij(n - 1 + z, i) = Hij(i, n - 1 + z);
    }
  }

  for (int z1 = 0; z1 < p; z1++) {
    for (int z2 = 0; z2 < p; z2++) {
      double z12E(0);
      for (int i = 0; i < (n-1); i++) {
        for (int j = (i+1); j < n; j++) {
          z12E += Zij[z1*n*n + i*n + j] * Zij[z2*n*n + i*n + j] * Eij(i, j);
        }
      }
      Hij(n - 1 + z1, n - 1 + z2) = N0 * (z12E * E0 - zE[z1] * zE[z2]) / (E0 * E0);
      Hij(n - 1 + z2, n - 1 + z1) = Hij(n - 1 + z1, n - 1 + z2);
    }
  }

  return List::create(_["g"] = Gi, _["e"] = Ei,  _["h"] = Hij,  _["zE"] = zE,  _["zN"] = zN);
}
