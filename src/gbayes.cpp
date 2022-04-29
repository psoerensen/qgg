// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;


#include <vector>
#include <numeric>
#include <random>
#include <ctime>
#include <chrono>
#include <algorithm>
#include <cmath>



// [[Rcpp::export]]
std::vector<std::vector<double>>  bayes(   std::vector<double> y,
                                           std::vector<std::vector<double>> W, 
                                           std::vector<double> b, 
                                           std::vector<double> lambda, 
                                           double pi, 
                                           double vg, 
                                           double vb, 
                                           double ve,
                                           double ssb_prior,
                                           double sse_prior,
                                           double nub,
                                           double nue,
                                           bool updateB,
                                           bool updateE,
                                           bool updatePi,
                                           int nit,
                                           int method) {
  
  
  
  // Define local variables
  int n = y.size();
  int m = W.size();
  
  double rhs, lhs, bn, conv, diff, mu;
  double rhs0, rhs1, lhs0, lhs1, like0, like1, p0;
  double ssb, sse, dfb, dfe, chi2;
  double xtau, tau, lambda_tau, mu_tau, z, z2, u;
  
  std::vector<int> d(m);
  
  std::vector<double> e(n);
  std::vector<double> wy(m),ww(m),x2(m);
  std::vector<int> order(m);
  
  std::vector<double> dm(m),bm(m);
  std::vector<double> ves(nit),vbs(nit),pis(nit),mus(nit);
  
  
  // Mean adjust y and initialize e
  mu = std::accumulate(y.begin(), y.end(), 0.0)/n;
  for ( int i = 0; i < n; i++) {
    e[i] = y[i] - mu;
  }                
  
  // Initialize variables
  for (int i = 0; i < m; i++) {
    wy[i] = 0.0;
    ww[i] = 0.0;
    for (int j = 0; j < n; j++) {
      wy[i] = wy[i] + W[i][j]*e[j];
      ww[i] = ww[i] + W[i][j]*W[i][j];
    }
    dm[i] = 0.0;
    bm[i] = 0.0;
    x2[i] = (wy[i]/ww[i])*(wy[i]/ww[i]);
  }
  
  for (int i=0; i < nit; i++) {
    vbs[i] = 0.0;
    ves[i] = 0.0;
    pis[i] = 0.0;
  }
  
  // Establish order of markers as they are entered into the model
  std::iota(order.begin(), order.end(), 0);
  std::sort(  std::begin(order), 
              std::end(order),
              [&](int i1, int i2) { return x2[i1] > x2[i2]; } );
  
  
  // Start Gibbs sampler
  std::random_device rd;
  unsigned int local_seed;
  local_seed = rd();
  std::mt19937 gen(local_seed);
  
  for ( int it = 0; it < nit; it++) {
    conv = 0.0;

    // Compute marker effects (BLUP)
    if (method==0) {
      for ( int isort = 0; isort < m; isort++) {
        int i = order[isort];
        lhs = ww[i] + lambda[i];
        rhs = 0.0;
        for ( int j = 0; j < n; j++) {
          rhs = rhs + W[i][j]*e[j]; 
        }
        rhs = rhs + ww[i]*b[i];
        bn = rhs/lhs;
        diff = bn-b[i];
        for (int j = 0; j < n; j++) {
          e[j]=e[j] - W[i][j]*(diff);
        }
        conv = conv + diff*diff;
        b[i] = bn;
      }
    }
    
    // Sample marker effects (Mixed, BayesA, Lasso)
    if ( method==1 || method==2 || method==3 ) {
      for ( int isort = 0; isort < m; isort++) {
        int i = order[isort];
        lhs = ww[i] + lambda[i];
        rhs = 0.0;
        for ( int j = 0; j < n; j++) {
          rhs = rhs + W[i][j]*e[j]; 
        }
        rhs = rhs + ww[i]*b[i];
        std::normal_distribution<double> rnorm(rhs/lhs, sqrt(ve/lhs));
        bn = rnorm(gen);
        diff = bn-b[i];
        for (int j = 0; j < n; j++) {
          e[j]=e[j] - W[i][j]*(diff);
        }
        conv = conv + diff*diff;
        b[i] = bn;
      }
    }
    // Sample marker effects (BayesC)
    if (method==4) {
      for ( int isort = 0; isort < m; isort++) {
        int i = order[isort];
        lhs0 = 1.0/vb;
        lhs1 = ww[i]/ve + 1.0/vb;
        rhs0 = 0.0;
        rhs1 = 0.0;
        for ( int j = 0; j < n; j++) {
          rhs1 = rhs1 + W[i][j]*e[j]/ve; 
        }
        rhs1 = rhs1 + ww[i]*b[i]/ve;
        like0 = sqrt((1.0/lhs0))*std::exp(0.5*(1.0/lhs0)*rhs0*rhs0);
        like1 = sqrt((1.0/lhs1))*std::exp(0.5*(1.0/lhs1)*rhs1*rhs1);
        like0 = like0*(1.0-pi); 
        like1 = like1*pi;
        p0 = like0/(like0+like1); 
        d[i]=0;
        std::uniform_real_distribution<double> runif(0.0, 1.0);
        u = runif(gen);
        if(u>p0) d[i]=1;
        bn=0.0;
        if(d[i]==1) {
          std::normal_distribution<double> rnorm(rhs1/lhs1, sqrt(1.0/lhs1));
          bn = rnorm(gen);
        } 
        diff = bn-b[i];
        for (int j = 0; j < n; j++) {
          e[j]=e[j] - W[i][j]*(diff);
        }
        conv = conv + diff*diff;
        b[i] = bn;
      }
    }
    
    // Sample marker variance
    ssb = 0.0;
    dfb = 0.0;
    for ( int i = 0; i < m; i++) {
      bm[i] = bm[i] + b[i];     
      if(d[i]==1)   {
        ssb = ssb + b[i]*b[i];  
        dfb = dfb + 1.0;
        dm[i] = dm[i] + 1.0;     
      }
    }
    if(updateB) {
      std::chi_squared_distribution<double> rchisq(dfb+nub);
      chi2 = rchisq(gen);
      vb = (ssb + ssb_prior)/chi2 ;
      vbs[it] = vb; 
    }
    
    
    // Sample residual variance
    if(updateE) {
      dfe = n + nue;
      sse = 0.0;
      for ( int j = 0; j < n; j++) {
        sse = sse + e[j]*e[j];
      }
      std::chi_squared_distribution<double> rchisq(dfe);
      chi2 = rchisq(gen);
      ve = (sse + sse_prior)/chi2 ;
      ves[it] = ve;
    }
    
    // Update lambda's for BLUP/Mixed
    if ( method<2 ) {
      for ( int i = 0; i < m; i++) {
        lambda[i] = ve/vb;
      }
    }
    
    // Sample marker specific variance for Bayes A
    if (method==2) {
      dfb = 1.0 + nub;
      for ( int i = 0; i < m; i++) { 
        ssb = b[i]*b[i];
        std::chi_squared_distribution<double> rchisq(dfb);
        chi2 = rchisq(gen);
        vb = (ssb + ssb_prior)/chi2 ;
        lambda[i] = ve/vb;
      }
    }
    
    // Sample marker specific tau for Bayes lasso
    if (method==3) { 
      lambda_tau = 2.0*2.0*0.5*0.5/vb;
      for ( int i = 0; i < m; i++) { 
        ssb = b[i]*b[i];
        mu_tau = sqrt(lambda_tau/ssb);
        std::normal_distribution<double> rnorm(0.0, 1.0);
        z = rnorm(gen);
        z2=z*z;
        xtau=mu_tau+0.5*mu_tau*mu_tau*z2/lambda_tau - 0.5*(mu_tau/lambda_tau)*sqrt(4*mu_tau*lambda_tau*z2+mu_tau*mu_tau*z2*z2);
        std::uniform_real_distribution<double> runif(0.0, 1.0);
        u = runif(gen);
        tau = mu_tau*mu_tau/xtau;
        if(u <= mu_tau/(mu_tau+xtau)) tau=xtau;
        lambda[i] = ve/tau;
      }
    }
    
    // Sample pi for Bayes C
    if(updatePi) {
      double count = dfb + 1.0;
      std::gamma_distribution<double> rgamma(count,1.0);
      double rg = rgamma(gen);
      pi = rg/(double)m;
      pis[it] = pi;
    }
    
    // Update mu and adjust residuals
    mu = std::accumulate(e.begin(), e.end(), 0.0)/n;
    for ( int i = 0; i < n; i++) {
      e[i] = e[i] - mu;
    }                
    mus[it] = mu;
    
  }
  
  // Summarize results
  std::vector<std::vector<double>> result(10);
  result[0].resize(m);
  result[1].resize(m);
  result[2].resize(nit);
  result[3].resize(nit);
  result[4].resize(nit);
  result[5].resize(nit);
  result[6].resize(n);
  result[7].resize(n);
  result[8].resize(3);
  result[9].resize(m);
  
  for (int i=0; i < m; i++) {
    result[0][i] = bm[i]/nit;
    result[1][i] = dm[i]/nit;
  }
  for (int i=0; i < nit; i++) {
    result[2][i] = mus[i];
    result[3][i] = vbs[i];
    result[4][i] = ves[i];
    result[5][i] = pis[i];
  }
  for (int i=0; i < n; i++) {
    result[6][i] = y[i]- mu - e[i];
    result[7][i] = e[i];
  }
  result[8][0] = vb;
  result[8][1] = ve;
  result[8][2] = pi;
  for (int i=0; i < m; i++) {
    result[9][i] = b[i];
  }
  
  return result;
}


// [[Rcpp::export]]
std::vector<std::vector<double>>  sbayes( std::vector<double> wy,
                                          std::vector<std::vector<double>> LD,
                                          std::vector<double> b,
                                          std::vector<double> lambda,
                                          double yy,
                                          double pi,
                                          double vg,
                                          double vb,
                                          double ve,
                                          double ssb_prior,
                                          double sse_prior,
                                          double nub,
                                          double nue,
                                          bool updateB,
                                          bool updateE,
                                          bool updatePi,
                                          int n,
                                          int nit,
                                          int method) {
  
  // Define local variables
  int m = b.size();
  
  double rhs, lhs, bn, conv, diff;
  double rhs0, rhs1, lhs0, lhs1, like0, like1, p0;
  double ssb, sse, ssg, dfb, dfe, chi2;
  double xtau, tau, lambda_tau, mu_tau, z, z2, u;
  
  std::vector<int> d(m);
  
  std::vector<double> ww(m),r(m);
  
  std::vector<double> dm(m),bm(m);
  std::vector<double> ves(nit),vbs(nit),pis(nit);
  
  std::vector<double> x2(m);
  std::vector<int> order(m);
  
  // Initialize variables
  for ( int i = 0; i < m; i++) {
    dm[i] = 0.0;
    bm[i] = 0.0;
    ww[i] = LD[i][i];
    r[i] = wy[i];
    x2[i] = (wy[i]/ww[i])*(wy[i]/ww[i]);
  }
  
  for (int i=0; i < nit; i++) {
    vbs[i] = 0.0;
    ves[i] = 0.0;
    pis[i] = 0.0;
  }
  
  // Establish order of markers as they are entered into the model
  std::iota(order.begin(), order.end(), 0);
  std::sort(  std::begin(order),
              std::end(order),
              [&](int i1, int i2) { return x2[i1] > x2[i2]; } );
  
  // Wy - W'Wb
  for ( int i = 0; i < m; i++) {
    if (b[i]!= 0.0) {
      for (int j = 0; j < m; j++) {
        r[j]=r[j] - LD[i][j]*b[i];
      }
    }
  }
  
  // Start Gibbs sampler
  std::random_device rd;
  unsigned int local_seed;
  local_seed = rd();
  std::mt19937 gen(local_seed);
  
  for ( int it = 0; it < nit; it++) {
    conv = 0.0;
    
    // Compute marker effects (BLUP)
    if (method==0) {
      for ( int isort = 0; isort < m; isort++) {
        int i = order[isort];
        lhs = ww[i] + lambda[i];
        rhs = r[i] + ww[i]*b[i];
        bn = rhs/lhs;
        diff = bn-b[i];
        for (int j = 0; j < m; j++) {
          r[j]=r[j] - LD[i][j]*(diff);
        }
        conv = conv + diff*diff;
        b[i] = bn;
      }
    }
    
    // Sample marker effects (Mixed, BayesA, Lasso)
    if ( method==1 || method==2 || method==3 ) {
      for ( int isort = 0; isort < m; isort++) {
        int i = order[isort];
        lhs = ww[i] + lambda[i];
        rhs = r[i] + ww[i]*b[i];
        std::normal_distribution<double> rnorm(rhs/lhs, sqrt(ve/lhs));
        bn = rnorm(gen);
        diff = bn-b[i];
        for (int j = 0; j < m; j++) {
          r[j]=r[j] - LD[i][j]*(diff);
        }
        conv = conv + diff*diff;
        b[i] = bn;
      }
    }
    
    // Sample marker effects (BayesC)
    if (method==4) {
      for ( int isort = 0; isort < m; isort++) {
        int i = order[isort];
        lhs0 = 1.0/vb;
        lhs1 = ww[i]/ve + 1.0/vb;
        rhs0 = 0.0;
        rhs1 = 0.0;
        rhs1 = r[i]/ve + ww[i]*b[i]/ve;
        like0 = sqrt((1.0/lhs0))*std::exp(0.5*(1.0/lhs0)*rhs0*rhs0);
        like1 = sqrt((1.0/lhs1))*std::exp(0.5*(1.0/lhs1)*rhs1*rhs1);
        like0 = like0*(1.0-pi);
        like1 = like1*pi;
        p0 = like0/(like0+like1);
        d[i]=0;
        std::uniform_real_distribution<double> runif(0.0, 1.0);
        u = runif(gen);
        if(u>p0) d[i]=1;
        bn=0.0;
        if(d[i]==1) {
          std::normal_distribution<double> rnorm(rhs1/lhs1, sqrt(1.0/lhs1));
          bn = rnorm(gen);
        }
        diff = bn-b[i];
        if(diff!=0.0) {
          for (int j = 0; j < m; j++) {
            r[j]=r[j] - LD[i][j]*(diff);
          }
          conv = conv + diff*diff;
        }
        b[i] = bn;
      }
    }
    
    // Sample marker variance
    ssb = 0.0;
    dfb = 0.0;
    for ( int i = 0; i < m; i++) {
      bm[i] = bm[i] + b[i];
      if(d[i]==1)   {
        ssb = ssb + b[i]*b[i];
        dfb = dfb + 1.0;
        dm[i] = dm[i] + 1.0;
      }
    }
    if(updateB) {
      std::chi_squared_distribution<double> rchisq(dfb+nub);
      chi2 = rchisq(gen);
      vb = (ssb + ssb_prior)/chi2 ;
      vbs[it] = vb; 
    }
    
    
    // Sample residual variance
    if(updateE) {
      dfe = n + nue;
      ssg = 0.0;
      sse = 0.0;
      for ( int i = 0; i < m; i++) {
        ssg = ssg + b[i] * (wy[i] -  r[i]);
        sse = sse + b[i] * (r[i] + wy[i]);
      }
      sse = yy - sse;
      std::chi_squared_distribution<double> rchisq(dfe);
      chi2 = rchisq(gen);
      ve = (sse + sse_prior)/chi2 ;
      ves[it] = ve;
    }
    
    // Update lambda's for BLUP/Mixed
    if ( method==1) {
      for ( int i = 0; i < m; i++) {
        lambda[i] = ve/vb;
      }
    }
    
    // Sample marker specific variance for Bayes A
    if (method==2) {
      dfb = 1.0 + nub;
      for ( int i = 0; i < m; i++) {
        ssb = b[i]*b[i];
        std::chi_squared_distribution<double> rchisq(dfb);
        chi2 = rchisq(gen);
        vb = (ssb + ssb_prior)/chi2 ;
        lambda[i] = ve/vb;
      }
    }
    
    // Sample marker specific tau for Bayesian lasso
    if (method==3) {
      lambda_tau = 2.0*2.0*0.5*0.5/vb;
      for ( int i = 0; i < m; i++) {
        ssb = b[i]*b[i];
        mu_tau = sqrt(lambda_tau/ssb);
        std::normal_distribution<double> norm(0.0, 1.0);
        z = norm(gen);
        z2=z*z;
        xtau=mu_tau+0.5*mu_tau*mu_tau*z2/lambda_tau - 0.5*(mu_tau/lambda_tau)*sqrt(4*mu_tau*lambda_tau*z2+mu_tau*mu_tau*z2*z2);
        std::uniform_real_distribution<double> runif(0.0, 1.0);
        u = runif(gen);
        tau = mu_tau*mu_tau/xtau;
        if(u <= mu_tau/(mu_tau+xtau)) tau=xtau;
        lambda[i] = ve/tau;
      }
    }
    
    // Sample pi for Bayes C
    if(method==4 && updatePi) {
      double count = dfb + 1.0;
      std::gamma_distribution<double> rgamma(count,1.0);
      double rg = rgamma(gen);
      pi = rg/(double)m;
      pis[it] = pi;
    }
  }
  
  // Summarize results
  std::vector<std::vector<double>> result(10);
  result[0].resize(m);
  result[1].resize(m);
  result[2].resize(nit);
  result[3].resize(nit);
  result[4].resize(nit);
  result[5].resize(nit);
  result[6].resize(m);
  result[7].resize(m);
  result[8].resize(3);
  result[9].resize(m);
  
  for (int i=0; i < m; i++) {
    result[0][i] = bm[i]/nit;
    result[1][i] = dm[i]/nit;
  }
  for (int i=0; i < nit; i++) {
    //result[2][i] = mus[i];
    result[2][i] = 0.0;
    result[3][i] = vbs[i];
    result[4][i] = ves[i];
    result[5][i] = pis[i];
  }
  for (int i=0; i < m; i++) {
    result[6][i] = wy[i];
    result[7][i] = r[i];
  }
  result[8][0] = vb;
  result[8][1] = ve;
  result[8][2] = pi;
  for (int i=0; i < m; i++) {
    result[9][i] = b[i];
  }
  
  
  return result;
}


// [[Rcpp::export]]
std::vector<std::vector<double>>  sbayes_spa( std::vector<double> wy,
                                              std::vector<std::vector<double>> LDvalues, 
                                              std::vector<std::vector<int>> LDindices, 
                                              std::vector<double> b, 
                                              std::vector<double> lambda, 
                                              double yy, 
                                              double pi, 
                                              double vg, 
                                              double vb, 
                                              double ve,
                                              double ssb_prior,
                                              double sse_prior,
                                              double nub,
                                              double nue,
                                              bool updateB,
                                              bool updateE,
                                              bool updatePi,
                                              int n, 
                                              int nit,
                                              int method) {
  
  // Define local variables
  int m = b.size();
  
  double rhs, lhs, bn, conv, diff;
  double rhs0, rhs1, lhs0, lhs1, like0, like1, p0, v0. v1;
  double ssb, sse, ssg, dfb, dfe, chi2;
  double xtau, tau, lambda_tau, mu_tau, z, z2, u;
  
  std::vector<int> d(m);
  
  std::vector<double> ww(m),r(m);
  
  std::vector<int> mask(m);
  std::vector<double> dm(m),bm(m);
  std::vector<double> ves(nit),vbs(nit),pis(nit);
  
  std::vector<double> x2(m);
  std::vector<int> order(m);
  
  
  // Initialize variables
  for ( int i = 0; i < m; i++) {
    mask[i]=1;
    dm[i] = 0.0;
    bm[i] = 0.0;
    //ww[i] = LD[i][i];
    ww[i] = (double)n;
    r[i] = wy[i];
    x2[i] = (wy[i]/ww[i])*(wy[i]/ww[i]);
    if(wy[i]==0.0) mask[i]=0;
  }
  for (int i=0; i < nit; i++) {
    vbs[i] = 0.0;
    ves[i] = 0.0;
    pis[i] = 0.0;
  }
  
  // Establish order of markers as they are entered into the model
  std::iota(order.begin(), order.end(), 0);
  std::sort(  std::begin(order), 
              std::end(order),
              [&](int i1, int i2) { return x2[i1] > x2[i2]; } );
  
  // Wy - W'Wb
  for ( int i = 0; i < m; i++) {
    if (b[i]!= 0.0) {
      for (size_t j = 0; j < LDindices[i].size(); j++) {
        r[LDindices[i][j]]=r[LDindices[i][j]] - LDvalues[i][j]*b[i];
      }
    }
  }
  
  // Start Gibbs sampler
  std::random_device rd;
  unsigned int local_seed;
  local_seed = rd();
  std::mt19937 gen(local_seed);
  
  for ( int it = 0; it < nit; it++) {
    conv = 0.0;
    
    // Compute marker effects (BLUP)
    if (method==0) {
      for ( int isort = 0; isort < m; isort++) {
        int i = order[isort];
        lhs = ww[i] + lambda[i];
        rhs = r[i] + ww[i]*b[i];
        bn = rhs/lhs;
        diff = bn-b[i];
        for (size_t j = 0; j < LDindices[i].size(); j++) {
          r[LDindices[i][j]]=r[LDindices[i][j]] - LDvalues[i][j]*diff;
        }
        conv = conv + diff*diff;
        b[i] = bn;
      }
    }
    
    // Sample marker effects (Mixed, BayesA, Lasso)
    if ( method==1 || method==2 || method==3 ) {
      for ( int isort = 0; isort < m; isort++) {
        int i = order[isort];
        lhs = ww[i] + lambda[i];
        rhs = r[i] + ww[i]*b[i];
        std::normal_distribution<double> rnorm(rhs/lhs, sqrt(ve/lhs));
        bn=0.0;
        if(mask[i]==1) {
          bn = rnorm(gen);
          diff = bn-b[i];
          for (size_t j = 0; j < LDindices[i].size(); j++) {
            r[LDindices[i][j]]=r[LDindices[i][j]] - LDvalues[i][j]*diff;
          }
          conv = conv + diff*diff;
        }
        b[i] = bn;
      }
    }
    // Sample marker effects (BayesC)
    if (method==4) {
      for ( int isort = 0; isort < m; isort++) {
        int i = order[isort];
        //lhs0 = 1.0/vb;
        //lhs0 = ww[i]/ve + 1.0/vb;
        //lhs0 = ww[i]/ve;
        //lhs0 = ww[i];
        //lhs1 = ww[i]/ve + 1.0/vb;
        //lhs1 = ww[i] + ve/vb;
        //rhs0 = 0.0;
        //rhs1 = 0.0;
        //rhs0 = r[i] + ww[i]*b[i];
        //rhs0 = r[i]/ve + ww[i]*b[i]/ve;
        //rhs1 = r[i]/ve + ww[i]*b[i]/ve;
        //rhs1 = r[i] + ww[i]*b[i];
        lhs = ww[i] + 1/vb;
        rhs = r[i]/ve + ww[i]*b[i]/ve;
        ri =r[i] + ww[i]*b[i];
        v0 = ww[i]*ve;
        v1 = ww[i]*ve + ww[i]*ww[i]*va;
        //like0 = sqrt((1.0/lhs))*std::exp(-0.5*(1.0/lhs0)*rhs0*rhs0);
        //like1 = sqrt((1.0/lhs))*std::exp(-0.5*(1.0/lhs1)*rhs1*rhs1);
        like0 = sqrt((1.0/v0))*std::exp(-0.5*((ri*ri)/v0));
        like1 = sqrt((1.0/v1))*std::exp(-0.5*((ri*ri)/v1));
        like0 = like0*(1.0-pi); 
        like1 = like1*pi;
        p0 = like0/(like0+like1); 
        d[i]=0;
        std::uniform_real_distribution<double> runif(0.0, 1.0);
        u = runif(gen);
        if(u>p0) d[i]=1;
        bn=0.0;
        if(d[i]==1) {
          std::normal_distribution<double> rnorm(rhs1/lhs1, sqrt(1.0/lhs1));
          bn = rnorm(gen);
        } 
        diff = bn-b[i];
        if(diff!=0.0) {
          for (size_t j = 0; j < LDindices[i].size(); j++) {
            r[LDindices[i][j]]=r[LDindices[i][j]] - LDvalues[i][j]*diff;
          }
          conv = conv + diff*diff;
        }
        b[i] = bn;
      }
    }
    
    // Sample marker variance
    ssb = 0.0;
    dfb = 0.0;
    for ( int i = 0; i < m; i++) {
      bm[i] = bm[i] + b[i];
      if(d[i]==1)   {
        ssb = ssb + b[i]*b[i];
        dfb = dfb + 1.0;
        dm[i] = dm[i] + 1.0;
      }
    }
    if(updateB) {
      std::chi_squared_distribution<double> rchisq(dfb+nub);
      chi2 = rchisq(gen);
      vb = (ssb + ssb_prior)/chi2 ;
      vbs[it] = vb; 
    }
    
    // Sample residual variance
    if(updateE) {
      dfe = n + nue;
      ssg = 0.0;
      sse = 0.0;
      for ( int i = 0; i < m; i++) {
        ssg = ssg + b[i] * (wy[i] -  r[i]);
        sse = sse + b[i] * (r[i] + wy[i]);
      }
      sse = yy - sse;
      std::chi_squared_distribution<double> rchisq(dfe);
      chi2 = rchisq(gen);
      ve = (sse + sse_prior)/chi2 ;
      ves[it] = ve;
    }
    
    // Update lambda's for BLUP/Mixed
    if ( method==1 ) {
      for ( int i = 0; i < m; i++) {
        lambda[i] = ve/vb;
      }
    }
    
    // Sample marker specific variance for Bayes A
    if (method==2) {
      dfb = 1.0 + nub;
      for ( int i = 0; i < m; i++) { 
        ssb = b[i]*b[i];
        std::chi_squared_distribution<double> rchisq(dfb);
        chi2 = rchisq(gen);
        vb = (ssb + ssb_prior)/chi2 ;
        lambda[i] = ve/vb;
      }
    }
    
    // Sample marker specific tau for Bayes lasso
    if (method==3) { 
      lambda_tau = 2.0*2.0*0.5*0.5/vb;
      for ( int i = 0; i < m; i++) { 
        ssb = b[i]*b[i];
        mu_tau = sqrt(lambda_tau/ssb);
        std::normal_distribution<double> norm(0.0, 1.0);
        z = norm(gen);
        z2=z*z;
        xtau=mu_tau+0.5*mu_tau*mu_tau*z2/lambda_tau - 0.5*(mu_tau/lambda_tau)*sqrt(4*mu_tau*lambda_tau*z2+mu_tau*mu_tau*z2*z2);
        std::uniform_real_distribution<double> runif(0.0, 1.0);
        u = runif(gen);
        tau = mu_tau*mu_tau/xtau;
        if(u <= mu_tau/(mu_tau+xtau)) tau=xtau;
        lambda[i] = ve/tau;
      }
    }
    
    // Sample pi for Bayes C
    if(method==4 && updatePi) {
      double count = dfb + 1.0;
      std::gamma_distribution<double> rgamma(count,1.0);
      double rg = rgamma(gen);
      pi = rg/(double)m;
      pis[it] = pi;
    }
  }
  
  // Summarize results
  std::vector<std::vector<double>> result(10);
  result[0].resize(m);
  result[1].resize(m);
  result[2].resize(nit);
  result[3].resize(nit);
  result[4].resize(nit);
  result[5].resize(nit);
  result[6].resize(m);
  result[7].resize(m);
  result[8].resize(3);
  result[9].resize(m);
  
  for (int i=0; i < m; i++) {
    result[0][i] = bm[i]/nit;
    result[1][i] = dm[i]/nit;
  }
  for (int i=0; i < nit; i++) {
    //result[2][i] = mus[i];
    result[2][i] = 0.0;
    result[3][i] = vbs[i];
    result[4][i] = ves[i];
    result[5][i] = pis[i];
  }
  for (int i=0; i < m; i++) {
    result[6][i] = wy[i];
    result[7][i] = r[i];
  }
  result[8][0] = vb;
  result[8][1] = ve;
  result[8][2] = pi;
  for (int i=0; i < m; i++) {
    result[9][i] = b[i];
  }
  
  
  return result;
}

