// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;


#include <vector>
#include <numeric>
#include <random>
#include <ctime>
#include <chrono>
#include <algorithm>
#include <cmath>


// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// simple example of creating two matrices and
// returning the result of an operatioon on them
//
// via the exports attribute we tell Rcpp to make this function
// available from R
//

// [[Rcpp::export]]
arma::mat cp(arma::mat & W) {
  arma::mat WW = W.t()*W;
  return WW;
}


// [[Rcpp::export]]
std::vector<int> psets( std::vector<int> msets,
                        std::vector<double> setstat,
                        std::vector<double> stat,
                        int np) {
  
  
  std::random_device rd;
  unsigned int local_seed;
  local_seed = rd();
  std::mt19937 gen(local_seed);
  double u;
  int k1, k2;
  
  int ns = msets.size();
  int m = stat.size();
  
  int maxm = *max_element(msets.begin(),msets.end());
  maxm = m - maxm - 2;
  std::vector<int> p(ns); 
  
  for ( int i = 0; i < ns; i++) {
    p[i] = 0; 
    for ( int j = 0; j < np; j++) {
      std::uniform_real_distribution<double> runif(0.0, 1.0);
      u = runif(gen);
      int imax = std::floor(maxm*u);
      k1 = imax;
      k2 = k1 + msets[i];
      double sumstat= 0.0;
      for ( int k = k1; k < k2; k++) { 
        sumstat = sumstat + stat[k];
      }
      //if (sumstat > setstat[i]) p[i] = p[i] + 1;
      if (std::abs(sumstat - setstat[i]) < 1e-10 ||
          sumstat > setstat[i]) {
        p[i] = p[i] + 1;
      }
    }
  }
  
  return p;
  
}


