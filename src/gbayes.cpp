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
                                           std::vector<double> pi, 
                                           std::vector<double> gamma, 
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
  int nc = pi.size();
  
    
  double rhs, lhs, bn, conv, diff, mu;
  double rhs0, rhs1, lhs0, lhs1, like0, like1, p0;
  double ssb, ssg, sse, dfb, dfe, dfg, chi2;
  double ssg_prior, nug;
  double xtau, tau, lambda_tau, mu_tau, z, z2, u;
  
  std::vector<int> d(m);
  
  std::vector<double> e(n),g(n);
  std::vector<double> wy(m),ww(m),x2(m);
  std::vector<int> order(m), mask(m);
  
  std::vector<double> dm(m),bm(m),vbi(m),vei(m);
  std::vector<double> ves(nit),vbs(nit),vgs(nit),pis(nit),mus(nit);

  std::vector<double> vbscale(nc), pic(nc), pim(nc), probc(nc), logLc(nc);
  double cumprobc, vbc, v0, v1, logLcAdj;
  
  for ( int i = 0; i < nc; i++) {
    vbscale[i]=gamma[i];
    pic[i]=pi[i];
  }
  
  // Mean adjust y and initialize e
  mu = std::accumulate(y.begin(), y.end(), 0.0)/n;
  for ( int i = 0; i < n; i++) {
    e[i] = y[i] - mu;
  }                

  std::fill(g.begin(), g.end(), 0.0);
  std::fill(wy.begin(), wy.end(), 0.0);
  std::fill(ww.begin(), ww.end(), 0.0);
  std::fill(dm.begin(), dm.end(), 0.0);
  std::fill(bm.begin(), bm.end(), 0.0);
  std::fill(vbs.begin(), vbs.end(), 0.0);
  std::fill(vgs.begin(), vgs.end(), 0.0);
  std::fill(ves.begin(), ves.end(), 0.0);
  std::fill(pis.begin(), pis.end(), 0.0);
  std::fill(pim.begin(), pim.end(), 0.0);
  
  
  // Initialize variables
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      wy[i] = wy[i] + W[i][j]*e[j];
      ww[i] = ww[i] + W[i][j]*W[i][j];
    }
    x2[i] = (wy[i]/ww[i])*(wy[i]/ww[i]);
    mask[i]=1;
    vbi[i]=vb;
  }
  
  nug=nub;
  ssg_prior=((nug-2.0)/nug)*vg;
  
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
        if(!mask[i])   continue;
        rhs = 0.0;
        for ( int j = 0; j < n; j++) {
          rhs = rhs + W[i][j]*e[j];
        }
        rhs = rhs/ve + ww[i]*b[i]/ve;
        lhs = ww[i]/ve + 1.0/vb;
        bn = rhs/lhs;
        diff = bn-b[i];
        for (int j = 0; j < n; j++) {
          e[j]=e[j] - W[i][j]*(diff);
        }
        b[i] = bn;
      }
    }
    
    // Compute marker effects (BayesN or BayesRR)
    if (method==1) {
      for ( int isort = 0; isort < m; isort++) {
        int i = order[isort];
        if(!mask[i])   continue;
        rhs = 0.0;
        for ( int j = 0; j < n; j++) {
          rhs = rhs + W[i][j]*e[j];
        }
        //rhs = rhs/ve + ww[i]*b[i]/ve;
        //lhs = ww[i]/ve + 1.0/vb;
        rhs = rhs + ww[i]*b[i];
        lhs = ww[i] + ve/vb;
        std::normal_distribution<double> rnorm(rhs/lhs, sqrt(ve/lhs));
        bn = rnorm(gen);
        diff = bn-b[i];
        for (int j = 0; j < n; j++) {
          e[j]=e[j] - W[i][j]*(diff);
        }
        b[i] = bn;
      }
    }

    // Compute marker effects (BayesA)
    if (method==2) {
      dfb = 1.0 + nub;
      for ( int isort = 0; isort < m; isort++) {
        int i = order[isort];
        if(!mask[i])   continue;
        ssb = b[i]*b[i];
        std::chi_squared_distribution<double> rchisq(dfb);
        chi2 = rchisq(gen);
        vbi[i] = (ssb + ssb_prior*nub)/chi2 ;
        rhs = 0.0;
        for ( int j = 0; j < n; j++) {
          rhs = rhs + W[i][j]*e[j];
        }
        rhs = rhs + ww[i]*b[i];
        lhs = ww[i] + ve/vbi[i];
        std::normal_distribution<double> rnorm(rhs/lhs, sqrt(ve/lhs));
        bn = rnorm(gen);
        diff = bn-b[i];
        for (int j = 0; j < n; j++) {
          e[j]=e[j] - W[i][j]*(diff);
        }
        b[i] = bn;
        // ssb = b[i]*b[i];
        // std::chi_squared_distribution<double> rchisq(dfb);
        // chi2 = rchisq(gen);
        // vbi[i] = (ssb + ssb_prior*nub)/chi2 ;
      }
    }
    
    // Sample marker effects (Lasso)
    if ( method==3 ) {
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
        rhs0 = 0.0;
        rhs1 = 0.0;
        for ( int j = 0; j < n; j++) {
          rhs1 = rhs1 + W[i][j]*e[j]/ve;
        }
        rhs1 = rhs1 + ww[i]*b[i]/ve;
        lhs0 = 1.0/vb;
        lhs1 = ww[i]/ve + 1.0/vb;
        like0 = std::log(1.0/std::sqrt(lhs0)) + 0.5*(rhs0*rhs0)/lhs0 + std::log(pi[0]);
        like1 = std::log(1.0/std::sqrt(lhs1)) + 0.5*(rhs1*rhs1)/lhs1 + std::log(pi[1]);
        p0 = 1.0/(std::exp(like1 - like0) + 1.0);
        d[i]=0;
        std::uniform_real_distribution<double> runif(0.0, 1.0);
        u = runif(gen);
        if(u>p0) d[i]=1;
        bn=0.0;
        if(d[i]==1) {
          //std::normal_distribution<double> rnorm(rhs1/lhs1, sqrt(ve/lhs1));
          std::normal_distribution<double> rnorm(rhs1/lhs1, sqrt(1.0/lhs1));
          bn = rnorm(gen);
        }
        diff = bn-b[i];
        if(diff!=0.0) {
          for (int j = 0; j < n; j++) {
            e[j]=e[j] - W[i][j]*(diff);
          }
        }
        conv = conv + diff*diff;
        b[i] = bn;
      }

      // Sample pi for Bayes C
      if(updatePi) {
        dfb=0.0;
        for (int i = 0; i<m ; i++) {
          if(d[i]==1)   {
            dfb = dfb + 1.0;
          }
        }
        double count = dfb + 1.0;
        std::gamma_distribution<double> rgamma(count,1.0);
        double rg = rgamma(gen);
        pi[1] = rg/(double)m;
        pi[0] = 1.0 - pi[1];
        pis[it] = pi[1];
        pim[0] = pim[0] + pi[0];
        pim[1] = pim[1] + pi[1];
      }  
      
            
    }

    // Sample marker effects (BayesR)
    if (method==5) {
      for ( int isort = 0; isort < m; isort++) {
        int i = order[isort];
        if(!mask[i])   continue;
        // variance class likelihood version 1 
        rhs = 0.0;
        for ( int j = 0; j < n; j++) {
          rhs = rhs + W[i][j]*e[j];
        }
        rhs = rhs + ww[i]*b[i];
        v0 = ww[i]*ve;
        logLc[0] = -0.5*std::log(v0) -0.5*((rhs*rhs)/v0) + std::log(pi[0]);
        for (int j = 1; j<nc ; j++) {
          vbc = vb * gamma[j];
          v1 = ww[i]*ve + ww[i]*ww[i]*vbc;
          logLc[j] = -0.5*std::log(v1) -0.5*((rhs*rhs)/v1) + std::log(pi[j]); 
        }
        
        // variance class probability 
        std::fill(probc.begin(), probc.end(), 0.0);
        for (int j = 0; j<nc ; j++) {
          logLcAdj = 0.0;
          for (int k = 0; k<nc ; k++) {
            logLcAdj += std::exp(logLc[k] - logLc[j]);
          }
          probc[j] = 1.0/logLcAdj;
        }
        // sample variance class indicator
        std::uniform_real_distribution<double> runif(0.0, 1.0);
        u = runif(gen);
        d[i]=0;
        cumprobc = 0.0;
        for (int j = 0; j<nc ; j++) {
          cumprobc += probc[j];
          if(u < cumprobc){
            d[i] = j;
            break;
          }
        }
        // sample marker effect
        bn=0.0;
        if(d[i]>0) {
          vbc = vb * gamma[d[i]];
          lhs =ww[i]+ve/vbc;
          std::normal_distribution<double> rnorm(rhs/lhs, sqrt(ve/lhs));
          //std::normal_distribution<double> rnorm(rhs/lhs, sqrt(1.0/lhs));
          bn = rnorm(gen);
        }
        diff = bn-b[i];
        if(diff!=0.0) {
          for (int j = 0; j < n; j++) {
            e[j]=e[j] - W[i][j]*(diff);
          }
        }
        conv = conv + diff*diff;
        b[i] = bn;
      }
      // Sample pi for Bayes R
      if(updatePi) {
        std::vector<double> mc(nc);
        std::fill(mc.begin(), mc.end(), 0.0);
        for (int i = 0; i<m ; i++) {
          mc[d[i]] = mc[d[i]] + 1.0;
        }
        double pisum=0.0;
        for (int j = 0; j<nc ; j++) {
          std::gamma_distribution<double> rgamma(mc[j]+1.0,1.0);
          double rg = rgamma(gen);
          pi[j] = rg/m;
          pisum = pisum + pi[j];
        }
        for (int j = 0; j<nc ; j++) {
          pi[j] = pi[j]/pisum;
          pim[j] = pim[j] + pi[j];
        }
      }
    }
    
    

    // Sample marker variance
    ssb = 0.0;
    dfb = 0.0;
    if (method<4) {
      for ( int i = 0; i < m; i++) {
        bm[i] = bm[i] + b[i];
        ssb = ssb + b[i]*b[i];
        dfb = dfb + 1.0;
        dm[i] = dm[i] + 1.0;
      }
    }
    if (method==4) {
      for ( int i = 0; i < m; i++) {
        bm[i] = bm[i] + b[i];
        if(d[i]==1)   {
          ssb = ssb + b[i]*b[i];
          dfb = dfb + 1.0;
          dm[i] = dm[i] + 1.0;
        }
      }
    }
    if (method==5) {
      for ( int i = 0; i < m; i++) {
        bm[i] = bm[i] + b[i];
        if(d[i]>0)   {
          ssb = ssb + (b[i]*b[i])/gamma[d[i]];
          dfb = dfb + 1.0;
          dm[i] = dm[i] + (double)d[i];
        }
      }
    }
    
    // marker variance
    if(updateB) {
      std::chi_squared_distribution<double> rchisq(dfb+nub);
      chi2 = rchisq(gen);
      vb = (ssb + ssb_prior*nub)/chi2 ;
      vbs[it] = (ssb + ssb_prior*nub)/chi2;
    }
    
    // Sample residual variance
    if(updateE) {
      dfe = (double)n + nue;
      sse = 0.0;
      for ( int j = 0; j < n; j++) {
        sse = sse + e[j]*e[j];
      }
      std::chi_squared_distribution<double> rchisq(dfe);
      chi2 = rchisq(gen);
      ve = (sse + sse_prior*nue)/chi2 ;
      ves[it] = ve;
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
    
    // Update mu and adjust residuals
    mu = std::accumulate(e.begin(), e.end(), 0.0)/n;
    for ( int i = 0; i < n; i++) {
      e[i] = e[i] - mu;
    }                
    mus[it] = mu;
    
    // Sample genetic variance
    ssg = 0.0;
    for ( int i = 0; i < n; i++) {
      g[i] = y[i] - e[i] - mu;
      ssg = ssg + g[i]*g[i];
    }                
    dfg = n + nug;
    std::chi_squared_distribution<double> rchisq(dfg);
    chi2 = rchisq(gen);
    vg = (ssg + ssg_prior*nug)/chi2;
    vgs[it] = vg;

    
  }
  
  // Summarize results
  std::vector<std::vector<double>> result(11);
  result[0].resize(m);
  result[1].resize(m);
  result[2].resize(nit);
  result[3].resize(nit);
  result[4].resize(nit);
  result[5].resize(nit);
  result[6].resize(nc);
  result[7].resize(n);
  result[8].resize(m);
  result[9].resize(m);
  result[10].resize(3);
  
  for (int i=0; i < m; i++) {
    result[0][i] = bm[i]/nit;
    result[1][i] = dm[i]/nit;
  }
  for (int i=0; i < nit; i++) {
    result[2][i] = mus[i];
    result[3][i] = vbs[i];
    result[4][i] = vgs[i];
    result[5][i] = ves[i];
  }
  for (int j=0; j < nc; j++) {
    result[6][j] = pim[j]/nit;
  }  
  
  for (int i=0; i < n; i++) {
    result[7][i] = y[i]- mu - e[i];
  }
  for (int i=0; i < m; i++) {
    result[8][i] = b[i];
    result[9][i] = d[i];
  }
  result[10][0] = vb;
  result[10][1] = ve;
  result[10][2] = pi[0];
  
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
      vb = (ssb + ssb_prior*nub)/chi2 ;
      vbs[it] = vb; 
    }
    
    
    // Sample residual variance
    if(updateE) {
      dfe = (double)n + nue;
      ssg = 0.0;
      sse = 0.0;
      for ( int i = 0; i < m; i++) {
        ssg = ssg + b[i] * (wy[i] -  r[i]);
        sse = sse + b[i] * (r[i] + wy[i]);
      }
      sse = yy - sse;
      std::chi_squared_distribution<double> rchisq(dfe);
      chi2 = rchisq(gen);
      ve = (sse + sse_prior*nue)/chi2 ;
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
        vb = (ssb + ssb_prior*nub)/chi2 ;
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
                                              std::vector<double> ww, 
                                              std::vector<std::vector<double>> LDvalues, 
                                              std::vector<std::vector<int>> LDindices, 
                                              std::vector<double> b, 
                                              std::vector<double> lambda, 
                                              double yy, 
                                              std::vector<double> pi, 
                                              std::vector<double> gamma, 
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
                                              bool updateG,
                                              bool adjustE,
                                              int n, 
                                              int nit,
                                              int method) {
  
  // Define local variables
  int m = b.size();
  int nc = pi.size();
  
  double rhs, lhs, bn, diff;
  double rhs1, lhs1, like0, like1, p0, v0, v1;
  double ssb, sse, ssg, dfb, dfe, dfg, chi2;
  double x_tau, tau, lambda_tau, mu_tau, z, z2, u, vbin;
  double shape, shape0, rate, rate0, lambda2;
  
  std::vector<double> vbscale(nc), pic(nc), pim(nc), probc(nc), logLc(nc);
  double cumprobc, vbc, logLcAdj;
  
  
  std::vector<int> d(m);

  std::vector<double> r(m),vei(m);
  
  std::vector<int> mask(m);
  std::vector<double> dm(m),bm(m);
  std::vector<double> ves(nit),vbs(nit),vgs(nit),pis(nit);
  
  std::vector<double> x2(m),vadj(m),vbi(m);
  std::vector<int> order(m);
  
  
  // Initialize variables
  for ( int i = 0; i < m; i++) {
    mask[i]=1;
    if(wy[i]==0) mask[i]=0;
    vbi[i]=vb;
    r[i] = wy[i];
    x2[i] = (wy[i]/ww[i])*(wy[i]/ww[i]);
    if(wy[i]==0.0) mask[i]=0;
  }
  
  std::fill(bm.begin(), bm.end(), 0.0);
  std::fill(dm.begin(), dm.end(), 0.0);
  std::fill(vbs.begin(), vbs.end(), 0.0);
  std::fill(vgs.begin(), vgs.end(), 0.0);
  std::fill(ves.begin(), ves.end(), 0.0);
  std::fill(pis.begin(), pis.end(), 0.0);
  std::fill(pim.begin(), pim.end(), 0.0);
  
  // adjust sparseld
  for ( int i = 0; i < m; i++) {
    vadj[i] = 0.0;
    if(adjustE) {
      vadj[i] = ((double)m-(double)LDindices[i].size())/(double)m;
    }  
    vei[i] = vadj[i]*vg + ve;
  }
  
  // should be added as argument to function
  //vy=yy/(n-1);
  //nug=nub;
  //ssg_prior=((nug-2.0)/nug)*vg;
  // initialize BayesL parameters
  //if (method==3) {
    dfb = (nub - 2.0)/nub;
    lambda2 = 2.0*(1.0 - dfb)/(dfb)*n;
    shape0 = 1.1;
    rate0 = (shape0 - 1.0) / lambda2;
    for ( int i = 0; i < m; i++) {
      lambda[i] = sqrt(lambda2); 
    }
  //}
  for ( int i = 0; i < nc; i++) {
    vbscale[i]=gamma[i];
    pic[i]=pi[i];
  }

  // Establish order of markers as they are entered into the model
  std::iota(order.begin(), order.end(), 0);
  std::sort(  std::begin(order), 
              std::end(order),
              [&](int i1, int i2) { return x2[i1] > x2[i2]; } );

  // Adjust LD values
  for ( int i = 0; i < m; i++) {
    for (size_t j = 0; j < LDindices[i].size(); j++) {
      LDvalues[i][j] = LDvalues[i][j]*std::sqrt(ww[i])*std::sqrt(ww[LDindices[i][j]]);
    }
  }
  
  // // Wy - W'Wb
  // for ( int i = 0; i < m; i++) {
  //   if (b[i]!= 0.0) {
  //     diff = b[i]*ww[i];
  //     for (size_t j = 0; j < LDindices[i].size(); j++) {
  //       r[LDindices[i][j]]=r[LDindices[i][j]] - LDvalues[i][j]*diff;
  //     }
  //   }
  // }
  

  // Start Gibbs sampler
  std::random_device rd;
  unsigned int local_seed;
  local_seed = rd();
  std::mt19937 gen(local_seed);
  
  for ( int it = 0; it < nit; it++) {
    
    // Compute marker effects (BLUP)
    if (method==0) {
      for ( int isort = 0; isort < m; isort++) {
        int i = order[isort];
        if(!mask[i])   continue;
        lhs = ww[i] + vei[i]/vb;
        rhs = r[i] + ww[i]*b[i];
        bn = rhs/lhs;
        diff = (bn-b[i]);
        for (size_t j = 0; j < LDindices[i].size(); j++) {
          r[LDindices[i][j]] += -LDvalues[i][j]*diff;
        }
        b[i] = bn;
      }
    }
    
    // Compute marker effects (BayesN or BayesRR)
    if (method==1) {
      for ( int isort = 0; isort < m; isort++) {
        int i = order[isort];
        if(!mask[i])   continue;
        lhs = ww[i] + vei[i]/vb;
        rhs = r[i] + ww[i]*b[i];
        std::normal_distribution<double> rnorm(rhs/lhs, sqrt(vei[i]/lhs));
        bn = rnorm(gen);
        diff = (bn-b[i]);
        for (size_t j = 0; j < LDindices[i].size(); j++) {
          r[LDindices[i][j]] += -LDvalues[i][j]*diff;
        }
        b[i] = bn;
      }
    }

    // Compute marker effects (BayesA)
    if (method==2) {
      dfb = 1.0 + nub;
      for ( int isort = 0; isort < m; isort++) {
        int i = order[isort];
        if(!mask[i])   continue;
        ssb = b[i]*b[i];
        std::chi_squared_distribution<double> rchisq(dfb);
        chi2 = rchisq(gen);
        vbi[i] = (ssb + ssb_prior*nub)/chi2 ;
        lhs = ww[i] + vei[i]/vbi[i];
        rhs = r[i] + ww[i]*b[i];
        std::normal_distribution<double> rnorm(rhs/lhs, sqrt(vei[i]/lhs));
        bn = rnorm(gen);
        diff = (bn-b[i]);
        for (size_t j = 0; j < LDindices[i].size(); j++) {
          r[LDindices[i][j]] += -LDvalues[i][j]*diff;
        }
        b[i] = bn;
      }
    }

    // Compute marker effects (BayesL)
    if (method==3) {
      dfb = 1.0 + nub;
      for ( int isort = 0; isort < m; isort++) {
        int i = order[isort];
        if(!mask[i])   continue;
        lhs = ww[i] + vei[i]/vbi[i];
        rhs = r[i] + ww[i]*b[i];
        std::normal_distribution<double> rnorm(rhs/lhs, sqrt(vei[i]/lhs));
        bn = rnorm(gen);
        diff = (bn-b[i]);
        for (size_t j = 0; j < LDindices[i].size(); j++) {
          r[LDindices[i][j]] += -LDvalues[i][j]*diff;
        }
        b[i] = bn;
        mu_tau=sqrt(vei[i])*lambda[i]/std::abs(b[i]);
        lambda_tau=lambda2;  
        std::normal_distribution<double> norm(0.0, 1.0);
        z = norm(gen);
        z2=z*z;
        x_tau=mu_tau+0.5*mu_tau*mu_tau*z2/lambda_tau - 0.5*(mu_tau/lambda_tau)*sqrt(4*mu_tau*lambda_tau*z2+mu_tau*mu_tau*z2*z2);
        std::uniform_real_distribution<double> runif(0.0, 1.0);
        u = runif(gen);
        tau = mu_tau*mu_tau/x_tau;
        if(u <= mu_tau/(mu_tau+x_tau)) tau=x_tau;
        //vbin = 1.0/tau;
        vbin = tau;
        if(vbin > 0)   vbi[i] = vbin;
        
      }
      // update hyperparameters
      ssb = 0.0;
      dfb = 0.0;
      for ( int i = 0; i < m; i++) {
        ssb = ssb + vbi[i]*vbi[i];
        dfb = dfb + 1.0;
      }
      shape = shape0 + dfb;
      rate = rate0 + ssb/ 2.0;
      std::gamma_distribution<double> rgamma(shape, 1.0/rate);
      lambda2 = rgamma(gen);
      for ( int i = 0; i < m; i++) {
        lambda[i] = sqrt(lambda2);
      }
    }
    
    // lambda_tau = 2.0*2.0*0.5*0.5/vb;
    // ssb = b[i]*b[i];
    // mu_tau = std::sqrt(lambda_tau/ssb);
    // std::normal_distribution<double> rnorm(0.0, 1.0);
    // z = rnorm(gen);
    // z2=z*z;
    // xtau=mu_tau+0.5*mu_tau*mu_tau*z2/lambda_tau - 0.5*(mu_tau/lambda_tau)*sqrt(4*mu_tau*lambda_tau*z2+mu_tau*mu_tau*z2*z2);
    // std::uniform_real_distribution<double> runif(0.0, 1.0);
    // u = runif(gen);
    // tau = mu_tau*mu_tau/xtau;
    // if(u <= mu_tau/(mu_tau+xtau)) tau=xtau;
    // lambda[i] = ve/tau;
    

    // Sample marker effects (BayesC)
    if (method==4) {
      for ( int isort = 0; isort < m; isort++) {
        int i = order[isort];
        // version 1
        //rhs0 = 0.0;
        //rhs1 = r[i] + ww[i]*b[i];
        //lhs0 = 1.0/vb;
        //lhs1 = ww[i]/vei[i] + 1.0/vb;
        //like0 = std::log(1.0/std::sqrt(lhs0)) + 0.5*(rhs0*rhs0)/lhs0 + std::log(1.0-pi); 
        //like1 = std::log(1.0/std::sqrt(lhs1)) + 0.5*(rhs1*rhs1)/lhs1 + std::log(pi); 
        //p0 = 1.0/(std::exp(like1 - like0) + 1.0);
        // version 2
        rhs = r[i] + ww[i]*b[i];
        v0 = ww[i]*vei[i];
        v1 = ww[i]*vei[i] + ww[i]*ww[i]*vb;
        like0 = -0.5*std::log(v0) -0.5*((rhs*rhs)/v0) + std::log(pi[0]);
        like1 = -0.5*std::log(v1) -0.5*((rhs*rhs)/v1) + std::log(pi[1]);
        p0 = 1.0/(std::exp(like1 - like0) + 1.0);
        d[i]=0;
        std::uniform_real_distribution<double> runif(0.0, 1.0);
        u = runif(gen);
        if(u>p0) d[i]=1;
        bn=0.0;
        if(d[i]==1) {
          rhs1 = r[i] + ww[i]*b[i];
          lhs1 = ww[i] + vei[i]/vb;
          std::normal_distribution<double> rnorm(rhs1/lhs1, sqrt(vei[i]/lhs1));
          bn = rnorm(gen);
        } 
        diff = (bn-b[i]);
        if(diff!=0.0) {
          for (size_t j = 0; j < LDindices[i].size(); j++) {
            r[LDindices[i][j]] += -LDvalues[i][j]*diff;
          }
        }
        b[i] = bn;
      }
      
      // Sample pi for Bayes C
      if(updatePi) {
        dfb=0.0;
        for (int i = 0; i<m ; i++) {
          if(d[i]==1)   {
            dfb = dfb + 1.0;
          }
        }
        double count = dfb + 1.0;
        std::gamma_distribution<double> rgamma(count,1.0);
        double rg = rgamma(gen);
        pi[1] = rg/(double)m;
        pi[0] = 1.0 - pi[1];
        pis[it] = pi[1];
        pim[0] = pim[0] + pi[0];
        pim[1] = pim[1] + pi[1];
      }
      
    }
    
    // Sample marker effects (BayesR)
    if (method==5) {
      
      for ( int isort = 0; isort < m; isort++) {
        int i = order[isort];
        if(!mask[i])   continue;
        // variance class likelihood 
        rhs = r[i] + ww[i]*b[i];
        v0 = ww[i]*vei[i];
        logLc[0] = -0.5*std::log(v0) -0.5*((rhs*rhs)/v0) + std::log(pi[0]);
        for (int j = 1; j<nc ; j++) {
          vbc = vb * gamma[j];
          v1 = ww[i]*vei[i] + ww[i]*ww[i]*vbc;
          logLc[j] = -0.5*std::log(v1) -0.5*((rhs*rhs)/v1) + std::log(pi[j]); 
        }
        // variance class probability 
        std::fill(probc.begin(), probc.end(), 0.0);
        for (int j = 0; j<nc ; j++) {
          logLcAdj = 0.0;
          for (int k = 0; k<nc ; k++) {
            logLcAdj += std::exp(logLc[k] - logLc[j]);
          }
          probc[j] = 1.0/logLcAdj;
        }
        // sample variance class indicator
        std::uniform_real_distribution<double> runif(0.0, 1.0);
        u = runif(gen);
        d[i]=0;
        cumprobc = 0.0;
        for (int j = 0; j<nc ; j++) {
          cumprobc += probc[j];
          if(u < cumprobc){
            d[i] = j;
            break;
          }
        }
        // sample marker effect
        bn=0.0;
        if(d[i]>0) {
          vbc = vb * gamma[d[i]];
          lhs =ww[i]+vei[i]/vbc;
          std::normal_distribution<double> rnorm(rhs/lhs, sqrt(vei[i]/lhs));
          bn = rnorm(gen);
        }
        diff = (bn-b[i]);
        if(diff!=0.0) {
          for (size_t j = 0; j < LDindices[i].size(); j++) {
            r[LDindices[i][j]] += -LDvalues[i][j]*diff;
          }
        }
        b[i] = bn;
      }
      // Sample pi for Bayes R
      if(updatePi) {
        std::vector<double> mc(nc);
        std::fill(mc.begin(), mc.end(), 0.0);
        for (int i = 0; i<m ; i++) {
          mc[d[i]] = mc[d[i]] + 1.0;
        }
        double pisum=0.0;
        for (int j = 0; j<nc ; j++) {
          std::gamma_distribution<double> rgamma(mc[j]+1.0,1.0);
          double rg = rgamma(gen);
          pi[j] = rg/m;
          pisum = pisum + pi[j];
        }
        for (int j = 0; j<nc ; j++) {
          pi[j] = pi[j]/pisum;
          pim[j] = pim[j] + pi[j];
        }
      }
    }
    
    
    // Sample marker variance
    ssb = 0.0;
    dfb = 0.0;
    if (method<4) {
      for ( int i = 0; i < m; i++) {
        bm[i] = bm[i] + b[i];
        ssb = ssb + b[i]*b[i];
        dfb = dfb + 1.0;
        dm[i] = dm[i] + 1.0;
      }
    }
    if (method==4) {
      for ( int i = 0; i < m; i++) {
        bm[i] = bm[i] + b[i];
        if(d[i]==1)   {
          ssb = ssb + b[i]*b[i];
          dfb = dfb + 1.0;
          dm[i] = dm[i] + 1.0;
        }
      }
    }
    if (method==5) {
      for ( int i = 0; i < m; i++) {
        bm[i] = bm[i] + b[i];
        if(d[i]>0)   {
          ssb = ssb + (b[i]*b[i])/gamma[d[i]];
          dfb = dfb + 1.0;
          dm[i] = dm[i] + (double)d[i];
        }
      }
    }
    
    // marker variance
    if(updateB) {
      std::chi_squared_distribution<double> rchisq(dfb+nub);
      chi2 = rchisq(gen);
      vb = (ssb + ssb_prior*nub)/chi2 ;
      vbs[it] = vb;
    }

    // Sample residual variance
    if(updateE) {
      sse = 0.0;
      for ( int i = 0; i < m; i++) {
        sse = sse + b[i] * (r[i] + wy[i]);
      }
      dfe = (double)n + nue;
      sse = yy - sse;
      std::chi_squared_distribution<double> rchisq(dfe);
      chi2 = rchisq(gen);
      ve = (sse + sse_prior*nue)/chi2 ;
      for ( int i = 0; i < m; i++) {
        vei[i] = vadj[i]*vg + ve;
      }
      ves[it] = ve;
    }

    
    // Sample genetic variance
    ssg = 0.0;
    for ( int i = 0; i < m; i++) {
      ssg = ssg + b[i] * (wy[i] -  r[i]);
    }
    //dfg = (double)n + nug;
    //std::chi_squared_distribution<double> rchisq(dfg);
    //chi2 = rchisq(gen);
    //vg = (ssg + ssg_prior*nug)/chi2;
    dfg = (double)n;
    vgs[it] = ssg/dfg;
    if(updateG) {
      vg = ssg/dfg;
    }
    if(adjustE) {
      for ( int i = 0; i < m; i++) {
        vei[i] = vadj[i]*vg + ve;
      }
    }
  }
  
  // Summarize results
  std::vector<std::vector<double>> result(11);
  result[0].resize(m);
  result[1].resize(m);
  result[2].resize(nit);
  result[3].resize(nit);
  result[4].resize(nit);
  result[5].resize(nit);
  result[6].resize(nc);
  result[7].resize(m);
  result[8].resize(m);
  result[9].resize(m);
  result[10].resize(3);

  for (int i=0; i < m; i++) {
    result[0][i] = bm[i]/nit;
    result[1][i] = dm[i]/nit;
  }
  for (int i=0; i < nit; i++) {
    //result[2][i] = mus[i];
    result[2][i] = 0.0;
    result[3][i] = vbs[i];
    result[4][i] = vgs[i];
    result[5][i] = ves[i];
  }
  for (int j=0; j < nc; j++) {
    result[6][j] = pim[j]/nit;
  }  
  
  for (int i=0; i < m; i++) {
    result[7][i] = r[i];
    result[8][i] = b[i];
    result[9][i] = d[i];
  }
  result[10][0] = vb;
  result[10][1] = ve;
  result[10][2] = pi[0];

  
  return result;
}

