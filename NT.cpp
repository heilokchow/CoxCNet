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
List GradCN(NumericVector beta, NumericVector Nij, NumericVector Zij, int n, int p, NumericVector vij, NumericVector etaij, int v) {
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

  for (int i = 0; i < (n - 1); i++) {
    int i1 = i + 1;
    Gi[i] -=  vij[i1 * n] + v * (0 - x[i1] - etaij[i1 * n]);
  }

  for (int i = 0; i < (n - 2); i++) {
    for (int j = (i+1); j < (n - 1); j++) {
      int i1 = i + 1;
      int j1 = j + 1;
      Gi[i] +=  vij[j1 * n + i1] + v * (x[i1] - x[j1] - etaij[j1 * n + i1]);
      Gi[j] -=  vij[j1 * n + i1] + v * (x[i1] - x[j1] - etaij[j1 * n + i1]);
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

  for (int i = 0; i < (n - 1); i++) {
    Hij(i, i) += v;
  }

  for (int i = 0; i < (n - 2); i++) {
    for (int j = (i+1); j < (n - 1); j++) {
      Hij(i, i) += v;
      Hij(i, j) -= v;
      Hij(j, i) -= v;
      Hij(j, j) += v;
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

  return List::create(_["g"] = Gi, _["e"] = Ei,  _["h"] = Hij);
}

NumericVector NewtonMC(NumericMatrix Trail, NumericVector Zij, double t, double h1, double h2, int n, int nn, int p) {

    if (Rf_isNull(Zij.attr("dim"))) {
        throw std::runtime_error("'x' does not have 'dim' attibute.");
    }
    Rcpp::Dimension d = Zij.attr("dim");
    if (d.size() != 3) {
        throw std::runtime_error("'x' must have 3 dimensions.");
    }

    if (d[0] != n || d[1] != n || d[2] != p)
        return 0;

    MatrixXd NijT1 = MatrixXd::Zero(n, n);
    MatrixXd NijT2 = MatrixXd::Zero(n, n);
    VectorXd f1 = VectorXd::Zero(2*n + p - 1);
    VectorXd xk = VectorXd::Zero(2*n + p - 1);
    Map<VectorXd> a(xk.data(), n);
    Map<VectorXd> b(xk.data() + n, n - 1);
    Map<VectorXd> c(xk.data() + 2*n - 1, p);

    VectorXd dd = VectorXd::Zero(2*n + p - 1);
    VectorXd f = VectorXd::Zero(2*n + p - 1);
    MatrixXd D = MatrixXd::Zero(2*n + p - 1, 2*n + p - 1);
    MatrixXd Eij = MatrixXd::Zero(n, n);

    MatrixXd TrailC = MatrixXd::Map(Trail.begin(), nn, 3);

    int i(0), j(0);
    double conv(1.0), mf1(10.0), temp(0.0);
    double bw1 = t - 10*h1 > 0 ? t - 10*h1 : 0;
    double bw2 = t + 10*h1 < 1 ? t + 10*h1 : 1;
    double bw3 = t - 10*h2 > 0 ? t - 10*h2 : 0;
    double bw4 = t + 10*h2 < 1 ? t + 10*h2 : 1;

    static const double inv_sqrt_2pi = 0.3989422804014327;
    double diff(0.0);

    for (int z = 0; z < nn; z++) {
        if (TrailC(z, 2) > bw1 && TrailC(z, 2) < bw2) {
            i = TrailC(z, 0);
            j = TrailC(z, 1);
            diff = (TrailC(z, 2) - t) / h1;
            NijT1(i - 1, j - 1) += inv_sqrt_2pi / h1 * std::exp(-0.5 * diff * diff);
        }
        if (TrailC(z, 2) > bw3 && TrailC(z, 2) < bw4) {
            i = TrailC(z, 0);
            j = TrailC(z, 1);
            diff = (TrailC(z, 2) - t) / h2;
            NijT2(i - 1, j - 1) += inv_sqrt_2pi / h2 * std::exp(-0.5 * diff * diff);
        }
    }


    for (int j = 0; j < n; j++) {
        for (int i = 0; i < n; i++){
            if (j != i) {
                f1(i) += NijT1(i, j);
            }
        }
    }

    for (int j = 0; j < n - 1; j++) {
        for (int i = 0; i < n; i++) {
            if (i != j) {
                f1(n + j) += NijT1(i, j);
            }
        }
    }

    for (int k = 0; k < p; k++) {
        for (int j = 0; j < n; j++) {
            for (int i = 0; i < n; i++) {
                if (i != j) {
                    f1(2*n - 1 + k) += Zij[i + j*n + k*n*n] * NijT2(i, j);
                }
            }
        }
    }


    for (int i = 0; i < 2*n - 1; i++) {
        if (f1(i) != 0 && mf1 > f1(i))
            mf1 = f1(i);
    }

    for (int i = 0; i < 2*n - 1; i++) {
        if (f1(i) == 0)
            f1(i) = mf1;
    }

    while (conv > 0.0001) {

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n - 1; j++) {
                if (i != j) {
                    temp = a(i) + b(j);
                    for (int k = 0; k < p; k++) {
                        temp += c(k) * Zij[i + j*n + k*n*n];
                    }
                    Eij(i, j) = exp(temp);
                }
            }
            temp = a(i);
            for (int k = 0; k < p; k++) {
                temp += c(k) * Zij[i + (n-1)*n + k*n*n];
            }
            Eij(i, n - 1) = exp(temp);
        }

        for (int i = 0; i < n; i++) {
            f(i) = 0;
            for (int j = 0; j < n; j++) {
                if (j != i) {
                    f(i) -= Eij(i, j);
                }
            }
            f(i) += f1(i);
        }

        for (int j = 0; j < n - 1; j++) {
            f(n + j) = 0;
            for (int i = 0; i < n; i++) {
                if (i != j) {
                    f(n + j) -= Eij(i, j);
                }
            }
            f(n + j) += f1(n + j);
        }

        for (int k = 0; k < p; k++) {
            f(2*n - 1 + k) = 0;
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    if (i != j) {
                        f(2*n - 1 + k) -= Zij[i + j*n + k*n*n] * Eij(i, j);
                    }
                }
            }
            f(2*n - 1 + k) += f1(2*n - 1 + k);
        }


        for (int i = 0; i < n; i++) {
            D(i, i) = 0;
            for (int k = 0; k < p; k++)
                D(i, 2*n - 1 + k) = 0;
            for (int j = 0; j < n; j++) {
                if (i != j) {
                    D(i, i) += Eij(i, j);
                    for (int k = 0; k < p; k++)
                        D(i, 2*n - 1 + k) += Zij[i + j*n + k*n*n] * Eij(i, j);
                }
            }
            for (int k = 0; k < p; k++)
                D(2*n - 1 + k, i) = D(i, 2*n - 1 + k);
        }

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n - 1; j++) {
                if (i != j) {
                    D(i, j + n) = Eij(i, j);
                }
            }
        }

        for (int j = 0; j < n - 1; j++) {
            for (int i = 0; i < n; i++) {
                if (i != j) {
                    D(j + n, i) = Eij(i, j);
                }
            }
        }

        for (int j = 0; j < n - 1; j++) {
            D(j + n, j + n) = 0;
            for (int k = 0; k < p; k++)
                D(j + n, 2*n - 1 + k) = 0;
            for (int i = 0; i < n; i++) {
                if (i != j) {
                    D(j + n, j + n) += Eij(i, j);
                    for (int k = 0; k < p; k++)
                        D(j + n, 2*n - 1 + k) += Zij[i + j*n + k*n*n] * Eij(i, j);
                }
            }
            for (int k = 0; k < p; k++)
                D(2*n - 1 + k, j + n) = D(j + n, 2*n - 1 + k);
        }

        for (int k1 = 0; k1 < p; k1++) {
            for (int k2 = 0; k2 < p; k2++)
                D(2*n - 1 + k1, 2*n - 1 + k2) = 0;
        }

        for (int k1 = 0; k1 < p; k1++) {
            for (int k2 = 0; k2 < p; k2++) {
                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < n; j++) {
                        if (i != j) {
                            D(2*n - 1 + k1, 2*n - 1 + k2) += Zij[i + j*n + k1*n*n] * Zij[i + j*n + k2*n*n] * Eij(i, j);
                        }
                    }
                }
            }
        }

        dd = D.partialPivLu().solve(f);
        xk += dd;
        conv = dd.norm();
    }


    return Rcpp::wrap(xk);
}
