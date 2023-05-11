#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <random>
#include <chrono>

using namespace Rcpp;

#define PI 3.14159265

double fs(double t, const double& shift, int i, const int &n);
double fr(double t, const double& shift, int i, const int &n);
double fg(double t);

std::mt19937 gen(123);

// [[Rcpp::export]]
List SimSetC(int n, double Csparse, NumericVector Zij) {

  if (Rf_isNull(Zij.attr("dim"))) {
    throw std::runtime_error("'x' does not have 'dim' attibute.");
  }
  Rcpp::Dimension d = Zij.attr("dim");
  if (d.size() != 3) {
    throw std::runtime_error("'x' must have 3 dimensions.");
  }

  std::size_t p = d[2];

  if (d[0] != n || d[1] != n)
    return 0;

  double maxit(50.0), temp(0.0), tp(0.0), shift(Csparse*log(n)), rej(0.0);
  std::vector<int> se;
  std::vector<int> re;
  std::vector<double> te;

  int kp(0);
  std::poisson_distribution<int> Pdis(maxit);
  std::uniform_real_distribution<double> Udis(0.0,1.0);

  for (int i = 1; i <= n - 1; i++) {
    for (int j = i + 1; j <= n; j++ ) {
      if (i != j) {
        kp = Pdis(gen);

        if (kp != 0) {
          for (int z = 0; z < kp; z++) {
            tp = Udis(gen);
            temp = fs(tp, shift, i, n) + fr(tp, shift, j, n);
            for (int di = 0; di < p; di++)
              temp += fg(tp) * Zij[i-1+n*(j-1)+n*n*di];
            temp = exp(temp);

            rej = Udis(gen);
            if (rej < temp / maxit) {
              se.push_back(i);
              re.push_back(j);
              te.push_back(tp);
            }

          }
        }

      }
    }
  }

  return List::create(se, re, te);
}

double fs(double t, const double& shift, int i, const int &n) {
  if (i <= n / 3) return 0;
  if (i > n / 3 && i <= 2*n / 3) return -1;
  return 1;
}

double fr(double t, const double& shift, int i, const int &n) {
  if (i <= n / 3) return 0;
  if (i > n / 3 && i <= 2*n / 3) return -1;
  return 1;
}

//
// double fs(double t, const double& shift, int i, const int &n) {
//   if (i <= 25) return 0;
//   if (i > 25 && i <= 50) return -0.5;
//   if (i > 50 && i <= 75) return 0.5;
//   return 0;
// }
//
// double fr(double t, const double& shift, int i, const int &n) {
//   if (i <= 25) return 0;
//   if (i > 25 && i <= 50) return -0.5;
//   if (i > 50 && i <= 75) return 0.5;
//   return 0;
// }
//
// double fs(double t, const double& shift, int i, const int &n) {
//   if (i <= 40) return 0;
//   if (i > 40 && i <= 80) return 0.2;
//   if (i > 80 && i <= 120) return 0.4;
//   if (i > 120 && i <= 160) return 0.6;
//   if (i > 160 && i <= 200) return 0.8;
//   return 1;
// }
//
// double fr(double t, const double& shift, int i, const int &n) {
//   if (i <= 40) return 0;
//   if (i > 40 && i <= 80) return 0.2;
//   if (i > 80 && i <= 120) return 0.4;
//   if (i > 120 && i <= 160) return 0.6;
//   if (i > 160 && i <= 200) return 0.8;
//   return 1;
// }

double fg(double t) {
  return 0.2;
}

