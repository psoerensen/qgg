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

#include <iostream>
#include <random>
#include <cmath>

// Function to sample from an inverse Gaussian distribution
double rinvgauss(double mu_tau, double lambda_tau, std::mt19937 &gen) {
  // Step 1: Generate a random sample from the normal distribution
  std::normal_distribution<double> rnorm(0.0, 1.0);
  double z = rnorm(gen);
  double z2 = z * z;
  
  // Step 2: Calculate xtau using the given formula
  double xtau = mu_tau + 0.5 * mu_tau * mu_tau * z2 / lambda_tau -
    0.5 * (mu_tau / lambda_tau) * std::sqrt(4 * mu_tau * lambda_tau * z2 + mu_tau * mu_tau * z2 * z2);
  
  // Step 3: Generate a random sample from the uniform distribution
  std::uniform_real_distribution<double> runif(0.0, 1.0);
  double u = runif(gen);
  
  // Step 4: Calculate tau based on the comparison
  double tau = mu_tau * mu_tau / xtau;
  if (u <= mu_tau / (mu_tau + xtau)) {
    tau = xtau;
  }
  
  return tau;
}



// [[Rcpp::export]]
std::vector<std::vector<double>>  bayes(   std::vector<double> y,
                                           std::vector<std::vector<double>> W, 
                                           std::vector<double> b, 
                                           std::vector<double> lambda, 
                                           std::vector<double> pi, 
                                           std::vector<double> gamma, 
                                           double vb, 
                                           double vg, 
                                           double ve,
                                           double ssb_prior,
                                           double ssg_prior,
                                           double sse_prior,
                                           double nub,
                                           double nug,
                                           double nue,
                                           bool updateB,
                                           bool updateG,
                                           bool updateE,
                                           bool updatePi,
                                           int nit,
                                           int nburn,
                                           int nthin,
                                           int method,
                                           int seed) {
  
  // Define local variables
  int n = y.size();
  int m = W.size();
  int nc = pi.size();
  double nsamples=0.0;
  
  double rhs, lhs, bn, diff, mu;
  double rhs0, rhs1, lhs0, lhs1, like0, like1, p0;
  double ssb, ssg, sse, dfb, dfe, dfg, chi2;
  double lambda_tau, mu_tau, u, lambda2;

  std::vector<int> d(m);
  
  std::vector<double> e(n),g(n);
  std::vector<double> wy(m),ww(m),x2(m);
  std::vector<int> order(m), mask(m);
  
  std::vector<double> dm(m),bm(m),vbi(m),vei(m);
  std::vector<double> ves(nit+nburn),vbs(nit+nburn),vgs(nit+nburn),pis(nit+nburn),mus(nit+nburn);
  
  std::vector<double> pim(nc), probc(nc), logLc(nc);
  double cumprobc, vbc, v0, v1, logLcAdj;
  
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
  

  // Mean adjust y and initialize e
  mu = std::accumulate(y.begin(), y.end(), 0.0)/n;
  for ( int i = 0; i < n; i++) {
    e[i] = y[i] - mu;
  }                
  
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
  
  lambda2 = m;
  
  // Establish order of markers as they are entered into the model
  std::iota(order.begin(), order.end(), 0);
  std::sort(  std::begin(order), 
              std::end(order),
              [&](int i1, int i2) { return x2[i1] > x2[i2]; } );
  
  
  // Start Gibbs sampler
  std::random_device rd;
  std::mt19937 gen(seed);

  for ( int it = 0; it < nit+nburn; it++) {
    
    if ( (it > nburn) && (it % nthin == 0) ) {
      nsamples = nsamples + 1.0;
    }
    
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
      }
    }
    
    // Sample marker effects (Lasso)
    if ( method==3 ) {
      for ( int isort = 0; isort < m; isort++) {
        int i = order[isort];
        // Legarra et al 2011
        lhs = ww[i]/ve + lambda[i];
        // Park & Casella - Campos et al.
        //lhs = ww[i]/ve + lambda[i]/ve;
        rhs = 0.0;
        for ( int j = 0; j < n; j++) {
          rhs = rhs + W[i][j]*e[j];
        }
        rhs = rhs + ww[i]*b[i];
        std::normal_distribution<double> rnorm((rhs/ve)/lhs, sqrt(1.0/lhs));
        bn = rnorm(gen);
        diff = bn-b[i];
        for (int j = 0; j < n; j++) {
          e[j]=e[j] - W[i][j]*(diff);
        }
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
        b[i] = bn;
      }

      // Sample pi for Bayes C
      if(updatePi) {
        std::vector<double> mc(2);
        std::fill(mc.begin(), mc.end(), 0.0);
        for (int i = 0; i<m ; i++) {
          mc[d[i]] = mc[d[i]] + 1.0;
        }
        double pisum=0.0;
        for (int j = 0; j<2 ; j++) {
          std::gamma_distribution<double> rgamma(mc[j]+1.0,1.0);
          double rg = rgamma(gen);
          pi[j] = rg/m;
          pisum = pisum + pi[j];
        }
        for (int j = 0; j<2 ; j++) {
          pi[j] = pi[j]/pisum;
          if(it>nburn) pim[j] = pim[j] + pi[j];
        }
        pis[it] = pi[1];
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
          if(it>nburn) pim[j] = pim[j] + pi[j];
        }
      }
    }
    
    
    // Sample marker variance
    ssb = 0.0;
    dfb = 0.0;
    if (method<4) {
      for ( int i = 0; i < m; i++) {
        if(it>nburn && (it % nthin == 0)) bm[i] = bm[i] + b[i];
        ssb = ssb + b[i]*b[i];
        dfb = dfb + 1.0;
        if(it>nburn && (it % nthin == 0)) dm[i] = dm[i] + 1.0;
      }
    }
    if (method==4) {
      for ( int i = 0; i < m; i++) {
        if(it>nburn && (it % nthin == 0)) bm[i] = bm[i] + b[i];
        if(d[i]==1)   {
          ssb = ssb + b[i]*b[i];
          dfb = dfb + 1.0;
          if(it>nburn && (it % nthin == 0)) dm[i] = dm[i] + 1.0;
        }
      }
    }
    if (method==5) {
      for ( int i = 0; i < m; i++) {
        if(it>nburn && (it % nthin == 0)) bm[i] = bm[i] + b[i];
        if(d[i]>0)   {
          ssb = ssb + (b[i]*b[i])/gamma[d[i]];
          dfb = dfb + 1.0;
          if(it>nburn && (it % nthin == 0)) dm[i] = dm[i] + (double)d[i];
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
      // Hyper parameters for the inverse Gaussian distribution
      //double lrate0 = 0.1;  // lambda0 hyperparameter rate
      double lshape0 = 1.0; // lambda0 hyperparameter shape
      
      // Set a small tolerance value (epsilon)
      const double epsilon = 1e-8;
      
      lambda_tau = lambda2;
      lambda_tau = std::max(lambda_tau, epsilon);
      
      // Compute mutau for each element in vector b
      double sum_tau = 0.0;
      for (size_t i = 0; i < b.size(); ++i) {
        // If b[i] is too small (close to zero), use epsilon instead
        double abs_b = std::abs(b[i]) < epsilon ? epsilon : std::abs(b[i]);
        // Calculate mu_tau with the tolerance for very small b[i]
        // Legarra et al 2011
        mu_tau = std::sqrt(lambda_tau) / abs_b;
        // Park & Casella - Campos et al.
        //mu_tau = std::sqrt(lambda_tau*ve) / abs_b;
        mu_tau = std::max(mu_tau, epsilon);
        double invtau = rinvgauss(mu_tau, lambda_tau, gen);
        lambda[i] = invtau;
        sum_tau += 1.0/invtau;
      }
      
      // Calculate ratel2 and shl2
      double shl2 = m + lshape0;    
      //double ratel2 = sum_tau / 2.0 + lrate0;
      //std::gamma_distribution<double> rgamma(shl2, 1.0/ratel2);
      double scl2 = 2.0/sum_tau;
      std::gamma_distribution<double> rgamma(shl2, scl2);
      lambda2 = rgamma(gen);
      lambda2 = std::max(1.0, std::min(lambda2, 1000000.0));
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
  std::vector<std::vector<double>> result(12);
  result[0].resize(m);
  result[1].resize(m);
  result[2].resize(nit+nburn);
  result[3].resize(nit+nburn);
  result[4].resize(nit+nburn);
  result[5].resize(nit+nburn);
  result[6].resize(nit+nburn);
  result[7].resize(nc);
  result[8].resize(n);
  result[9].resize(m);
  result[10].resize(m);
  result[11].resize(3);
  
  for (int i=0; i < m; i++) {
    //result[0][i] = bm[i]/nit;
    //result[1][i] = dm[i]/nit;
    result[0][i] = bm[i]/nsamples;
    result[1][i] = dm[i]/nsamples;
  }
  for (int i=0; i < nit+nburn; i++) {
    result[2][i] = mus[i];
    result[3][i] = vbs[i];
    result[4][i] = vgs[i];
    result[5][i] = ves[i];
    result[6][i] = pis[i];
  }
  for (int j=0; j < nc; j++) {
    result[7][j] = pim[j]/nit;
  }  
  for (int i=0; i < n; i++) {
    result[8][i] = y[i]- mu - e[i];
  }
  for (int i=0; i < m; i++) {
    result[9][i] = b[i];
    result[10][i] = d[i];
  }
  result[11][0] = vb;
  result[11][1] = ve;
  result[11][2] = pi[0];
  return result;
}




// [[Rcpp::export]]
std::vector<std::vector<double>>  sbayes( double yy,
                                          std::vector<double>& wy,
                                          std::vector<double>& ww, 
                                          std::vector<std::vector<double>>& LDvalues, 
                                          std::vector<std::vector<int>>& LDindices, 
                                          std::vector<double> b, 
                                          std::vector<double> lambda,
                                          std::vector<bool> mask, 
                                          std::vector<double> pi, 
                                          std::vector<double> gamma, 
                                          double vg, 
                                          double vb, 
                                          double ve,
                                          double ssb_prior,
                                          double ssg_prior,
                                          double sse_prior,
                                          double nub,
                                          double nug,
                                          double nue,
                                          bool updateB,
                                          bool updateG,
                                          bool updateE,
                                          bool updatePi,
                                          bool adjustE,
                                          int n, 
                                          int nit,
                                          int nburn,
                                          int nthin,
                                          int method,
                                          int algo,
                                          int seed) {
  
  
    
  // Define local variables
  int m = b.size();
  int nc = pi.size();
  double nsamples=0.0;
  
  
  double rhs, lhs, bn, bj, diff;
  double rhs1, lhs1, like0, like1, p0, v0, v1;
  //double rhs0, lhs0;
  double ssb, sse, ssg, dfb, dfe, dfg, chi2;
  double x_tau, tau, lambda_tau, mu_tau, z, z2, u, vbin;
  double shape, shape0, rate, rate0, lambda2;
  
  std::vector<double> vbscale(nc), pic(nc), pim(nc), probc(nc), logLc(nc);
  double cumprobc, vbc, logLcAdj;
  
  
  std::vector<int> d(m);

  std::vector<double> r(m),vei(m);
  
  std::vector<double> dm(m),bm(m);
  std::vector<double> ves(nit+nburn),vbs(nit+nburn),vgs(nit+nburn),pis(nit+nburn);
  
  std::vector<double> x2(m),vadj(m),vbi(m);
  std::vector<int> order(m);
  
  
  // Initialize variables
  for ( int i = 0; i < m; i++) {
    vbi[i]=vb;
    r[i] = wy[i];
    x2[i] = (wy[i]/ww[i])*(wy[i]/ww[i]);
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
  std::mt19937 gen(seed);
  
  for ( int it = 0; it < nit+nburn; it++) {
  
  if ( (it > nburn) && (it % nthin == 0) ) {
    nsamples = nsamples + 1.0;
  }
  
  
    // Compute marker effects (BLUP)
    if (method==0) {
      for ( int isort = 0; isort < m; isort++) {
        int i = order[isort];
        if(mask[i])   continue;
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
        if(mask[i])   continue;
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
        if(mask[i])   continue;
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
        if(mask[i])   continue;
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
        if(mask[i])   continue;
        // version 1
        //rhs0 = 0.0;
        //rhs1 = r[i] + ww[i]*b[i];
        //lhs0 = 1.0/vb;
        //lhs1 = ww[i]/vei[i] + 1.0/vb;
        //like0 = std::log(1.0/std::sqrt(lhs0)) + 0.5*(rhs0*rhs0)/lhs0 + std::log(pi[0]); 
        //like1 = std::log(1.0/std::sqrt(lhs1)) + 0.5*(rhs1*rhs1)/lhs1 + std::log(pi[1]); 
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
          if(algo==2 && it<nburn) {
            bn = (1.0-p0)*bn;
          }
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
        std::vector<double> mc(2);
        std::fill(mc.begin(), mc.end(), 0.0);
        for (int i = 0; i<m ; i++) {
          mc[d[i]] = mc[d[i]] + 1.0;
        }
        double pisum=0.0;
        for (int j = 0; j<2 ; j++) {
          std::gamma_distribution<double> rgamma(mc[j]+1.0,1.0);
          double rg = rgamma(gen);
          pi[j] = rg/m;
          pisum = pisum + pi[j];
        }
        for (int j = 0; j<2 ; j++) {
          pi[j] = pi[j]/pisum;
          if(it>nburn) pim[j] = pim[j] + pi[j];
        }
        pis[it] = pi[1];
        // dfb=0.0;
        // for (int i = 0; i<m ; i++) {
        //   if(d[i]==1)   {
        //     dfb = dfb + 1.0;
        //   }
        // }
        // double count = dfb + 1.0;
        // std::gamma_distribution<double> rgamma(count,1.0);
        // double rg = rgamma(gen);
        // double pisum=0.0;
        // pi[1] = rg/(double)m;
        // pi[0] = 1.0 - pi[1];
        // pisum = pi[0] + pi[1];
        // pi[0] = pi[0]/pisum;
        // pi[1] = pi[1]/pisum;
        // pis[it] = pi[1];
        // if(it>nburn) pim[0] = pim[0] + pi[0];
        // if(it>nburn) pim[1] = pim[1] + pi[1];
      }
      
    }
    
    // Sample marker effects (BayesR)
    if (method==5) {
      
      for ( int isort = 0; isort < m; isort++) {
        int i = order[isort];
        if(mask[i])   continue;
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
          if(algo==2 && it<nburn) {
          bn=0.0;
          for (size_t j = 1; j < gamma.size(); j++) {
            vbc = vb * gamma[j];
            lhs =ww[i]+vei[i]/vbc;
            std::normal_distribution<double> rnorm(rhs/lhs, sqrt(vei[i]/lhs));
            bj = rnorm(gen);
            bn += probc[j]*bj;
          }
          }
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
          if(it>nburn) pim[j] = pim[j] + pi[j];
        }
        pis[it] = 1.0 - pi[0];
      }
    }
    
    
    // Sample marker variance
    ssb = 0.0;
    dfb = 0.0;
    if (method<4) {
      for ( int i = 0; i < m; i++) {
        if(it>nburn && (it % nthin == 0)) bm[i] = bm[i] + b[i];
        ssb = ssb + b[i]*b[i];
        dfb = dfb + 1.0;
        if(it>nburn && (it % nthin == 0)) dm[i] = dm[i] + 1.0;
      }
    }
    if (method==4) {
      for ( int i = 0; i < m; i++) {
        if(it>nburn && (it % nthin == 0)) bm[i] = bm[i] + b[i];
        if(d[i]==1)   {
          ssb = ssb + b[i]*b[i];
          dfb = dfb + 1.0;
          if(it>nburn && (it % nthin == 0)) dm[i] = dm[i] + 1.0;
        }
      }
    }
    if (method==5) {
      for ( int i = 0; i < m; i++) {
        if(it>nburn && (it % nthin == 0)) bm[i] = bm[i] + b[i];
        if(d[i]>0)   {
          ssb = ssb + (b[i]*b[i])/gamma[d[i]];
          dfb = dfb + 1.0;
          if(it>nburn && (it % nthin == 0)) dm[i] = dm[i] + 1.0;
          //if(it>nburn) dm[i] = dm[i] + (double)d[i];
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
        if(mask[i])   continue;
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
      if(mask[i])   continue;
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
  result[2].resize(nit+nburn);
  result[3].resize(nit+nburn);
  result[4].resize(nit+nburn);
  result[5].resize(nit+nburn);
  result[6].resize(nit+nburn);
  result[7].resize(nc);
  result[8].resize(m);
  result[9].resize(m);
  result[10].resize(3);

  for (int i=0; i < m; i++) {
    //result[0][i] = bm[i]/nit;
    //result[1][i] = dm[i]/nit;
    result[0][i] = bm[i]/nsamples;
    result[1][i] = dm[i]/nsamples;
  }
  for (int i=0; i < nit+nburn; i++) {
    //result[2][i] = mus[i];
    result[2][i] = 0.0;
    result[3][i] = vbs[i];
    result[4][i] = vgs[i];
    result[5][i] = ves[i];
    result[6][i] = pis[i];
  }
  for (int j=0; j < nc; j++) {
    result[7][j] = pim[j]/nit;
  }  
  
  for (int i=0; i < m; i++) {
    result[8][i] = r[i];
    result[9][i] = b[i];
    //result[9][i] = d[i];
  }
  result[10][0] = vb;
  result[10][1] = ve;
  result[10][2] = pi[0];

  
  return result;
}


// [[Rcpp::export]]
std::vector<std::vector<double>>  sbayes_reg( double yy,
                                              std::vector<double>& wy,
                                              std::vector<double>& ww, 
                                              std::vector<std::vector<double>>& LDvalues, 
                                              std::vector<std::vector<int>>& LDindices, 
                                              std::vector<double> b, 
                                              std::vector<double> lambda, 
                                              std::vector<bool> mask, 
                                              std::vector<double> pi, 
                                              std::vector<double> gamma, 
                                              double vb, 
                                              double vg, 
                                              double ve,
                                              double ssb_prior,
                                              double ssg_prior,
                                              double sse_prior,
                                              double nub,
                                              double nug,
                                              double nue,
                                              bool updateB,
                                              bool updateG,
                                              bool updateE,
                                              bool updatePi,
                                              int n, 
                                              int nit,
                                              int nburn,
                                              int nthin,
                                              int method,
                                              int algo,
                                              int seed) {
  
  // Define local variables
  int m = b.size();
  int nc = pi.size();
  double nsamples=0.0;
  
  double rhs, lhs, bn, diff, logcpo;
  double rhs1, lhs1, like0, like1, p0, v0, v1;
  double ssb, sse, ssg, dfb, dfe, dfg, chi2;
  double u;
  
  std::vector<double> vbscale(nc), probc(nc), logLc(nc), pim(nc);
  double cumprobc, vbc, logLcAdj;
  
  std::vector<int> d(m),order(m);
  std::vector<double> r(m),vei(m),x2(m),vadj(m);
  std::vector<double> dm(m),bm(m);
  std::vector<double> logsum(m);
  
  std::vector<double> ves(nit+nburn),vbs(nit+nburn),vgs(nit+nburn),pis(nit+nburn);
  std::vector<std::vector<double>> bs(nit+nburn, std::vector<double>(m, 0.0));  
  std::vector<std::vector<int>> ds(nit+nburn, std::vector<int>(m, 0));  
  std::vector<std::vector<double>> prob(nit+nburn, std::vector<double>(m, 0.0));  
  
  // Initialize variables
  for ( int i = 0; i < m; i++) {
    r[i] = wy[i];
    x2[i] = (wy[i]/ww[i])*(wy[i]/ww[i]);
  }
  
  std::fill(bm.begin(), bm.end(), 0.0);
  std::fill(dm.begin(), dm.end(), 0.0);
  std::fill(vbs.begin(), vbs.end(), 0.0);
  std::fill(vgs.begin(), vgs.end(), 0.0);
  std::fill(ves.begin(), ves.end(), 0.0);
  std::fill(pis.begin(), pis.end(), 0.0);
  std::fill(pim.begin(), pim.end(), 0.0);
  std::fill(logsum.begin(), logsum.end(), 0.0);
  
  // Adjust residual
  for ( int i = 0; i < m; i++) {
    vadj[i] = 0.0;
    vadj[i] = ((double)m-(double)LDindices[i].size())/(double)m;
    vei[i] = vadj[i]*vg + ve;
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
  
  // Start Gibbs sampler
  std::random_device rd;
  std::mt19937 gen(seed);
  
  for ( int it = 0; it < nit+nburn; it++) {
    
    if ( (it > nburn) && (it % nthin == 0) ) {
      nsamples = nsamples + 1.0;
    }
    
    // Sample marker effects (BayesC)
    if (method==4) {
      for ( int isort = 0; isort < m; isort++) {
        int i = order[isort];
        if(mask[i])   continue;
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
        std::vector<double> mc(2);
        std::fill(mc.begin(), mc.end(), 0.0);
        for (int i = 0; i<m ; i++) {
          mc[d[i]] = mc[d[i]] + 1.0;
        }
        double pisum=0.0;
        for (int j = 0; j<2 ; j++) {
          std::gamma_distribution<double> rgamma(mc[j]+1.0,1.0);
          double rg = rgamma(gen);
          pi[j] = rg/m;
          pisum = pisum + pi[j];
        }
        for (int j = 0; j<2 ; j++) {
          pi[j] = pi[j]/pisum;
          if(it>nburn && (it % nthin == 0)) pim[j] = pim[j] + pi[j];
        }
        pis[it] = pi[1];
      }
      
    }
    
    // Sample marker effects (BayesR)
    if (method==5) {
      for ( int isort = 0; isort < m; isort++) {
        int i = order[isort];
        if(mask[i])   continue;
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
          if(it>nburn && (it % nthin == 0)) pim[j] = pim[j] + pi[j];
        }
        pis[it] = 1.0 - pi[0];
      }
    }
    
    // Sample marker variance
    ssb = 0.0;
    dfb = 0.0;
    if (method==4) {
      for ( int i = 0; i < m; i++) {
        if(it>nburn && (it % nthin == 0)) bm[i] = bm[i] + b[i];
        if(d[i]==1)   {
          ssb = ssb + b[i]*b[i];
          dfb = dfb + 1.0;
          if(it>nburn && (it % nthin == 0)) dm[i] = dm[i] + 1.0;
        }
      }
    }
    if (method==5) {
      for ( int i = 0; i < m; i++) {
        if(it>nburn && (it % nthin == 0)) bm[i] = bm[i] + b[i];
        if(d[i]>0)   {
          ssb = ssb + (b[i]*b[i])/gamma[d[i]];
          dfb = dfb + 1.0;
          if(it>nburn && (it % nthin == 0)) dm[i] = dm[i] + 1.0;
        }
      }
    }
    
    // marker variance
    std::chi_squared_distribution<double> rchisq(dfb+nub);
    chi2 = rchisq(gen);
    vbs[it] = (ssb + ssb_prior*nub)/chi2 ;
    if(updateB) {
      vb = vbs[it]; 
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
    
    
    // Update genetic variance
    ssg = 0.0;
    for ( int i = 0; i < m; i++) {
      ssg = ssg + b[i] * (wy[i] -  r[i]);
    }
    dfg = (double)n;
    vgs[it] = ssg/dfg;
    if(updateG) {
      vg = ssg/dfg;
    }
    for ( int i = 0; i < m; i++) {
      vei[i] = vadj[i]*vg + ve;
    }
    
    // Compute logcpo
    if(it>nburn && (it % nthin == 0)) {
      double constant = 1.0 / std::sqrt(2 * M_PI);  
      for (int i = 0; i < m; i++) {
        logsum[i] += 1.0/(constant * std::exp(-0.5 * r[i] * r[i]));
      }
    }
  
    for ( int i = 0; i < m; i++) {
      bs[it][i] = b[i];
      ds[it][i] = d[i];
    }
  }

  // Compute logcpo
  logcpo=0.0;
  for (int i=0; i < m; i++) {
    logcpo += std::log(nsamples*(1/logsum[i]));
  }
  
  // Summarize results
  std::vector<std::vector<double>> result(14);
  
  result[0].resize(m);
  result[1].resize(m);
  result[2].resize(nit+nburn);
  result[3].resize(nit+nburn);
  result[4].resize(nit+nburn);
  result[5].resize(nit+nburn);
  result[6].resize(nit+nburn);
  result[7].resize(nc);
  result[8].resize(m);
  result[9].resize(m);
  result[10].resize(4);
  result[11].resize((nit+nburn)*m);
  result[12].resize((nit+nburn)*m);
  result[13].resize((nit+nburn)*m);
  for (int i=0; i < m; i++) {
    result[0][i] = bm[i]/nsamples;
    result[1][i] = dm[i]/nsamples;
  }
  for (int i=0; i < nit+nburn; i++) {
    result[2][i] = 0.0;
    result[3][i] = vbs[i];
    result[4][i] = vgs[i];
    result[5][i] = ves[i];
    result[6][i] = pis[i];
  }
  for (int j=0; j < nc; j++) {
    result[7][j] = pim[j]/nsamples;
  }  
  
  for (int i=0; i < m; i++) {
    result[8][i] = r[i];
    result[9][i] = b[i];
  }
  result[10][0] = vb;
  result[10][1] = ve;
  result[10][2] = pi[0];
  result[10][3] = logcpo;
  
  for ( int it = 0; it < nit+nburn; it++) {
    for ( int i = 0; i < m; i++) {
      result[11][it*m + i] = bs[it][i];
      result[12][it*m + i] = (double)ds[it][i];
      result[13][it*m + i] = prob[it][i];
    }
  }
  return result;
}



// [[Rcpp::export]]
std::vector<std::vector<double>>  sbayes_reg_eigen( std::vector<double>& wy,
                                                    std::vector<double>& ww, 
                                                    std::vector<std::vector<double>>& LDvalues, 
                                                    std::vector<std::vector<int>>& LDindices, 
                                                    std::vector<double> b, 
                                                    std::vector<double> lambda, 
                                                    std::vector<bool> mask, 
                                                    std::vector<double> pi, 
                                                    std::vector<double> gamma, 
                                                    double vb, 
                                                    double vg, 
                                                    double ve,
                                                    double ssb_prior,
                                                    double ssg_prior,
                                                    double sse_prior,
                                                    double nub,
                                                    double nug,
                                                    double nue,
                                                    bool updateB,
                                                    bool updateG,
                                                    bool updateE,
                                                    bool updatePi,
                                                    int n, 
                                                    int nit,
                                                    int nburn,
                                                    int nthin,
                                                    int method,
                                                    int algo,
                                                    int seed) {
  
  // Define local variables
  int m = b.size();
  int nc = pi.size();
  double nsamples=0.0;
  int q=wy.size();
  
  double rhs, lhs, bn, diff;
  double ssb, sse, ssg, dfb, dfe, chi2;
  double u;
  
  std::vector<double> probc(nc), logLc(nc), pim(nc);
  double cumprobc, vbc, logLcAdj, logcpo;
  
  std::vector<int> d(m);
  std::vector<double> r(m);
  std::vector<double> dm(m),bm(m);
  std::vector<double> logsum(q);
  
  std::vector<double> ves(nit+nburn),vbs(nit+nburn),vgs(nit+nburn),pis(nit+nburn);
  std::vector<std::vector<double>> bs(nit+nburn, std::vector<double>(m, 0.0));  
  std::vector<std::vector<int>> ds(nit+nburn, std::vector<int>(m, 0));  
  std::vector<std::vector<double>> prob(nit+nburn, std::vector<double>(m, 0.0));  
  
  // Initialize variables
  
  std::fill(bm.begin(), bm.end(), 0.0);
  std::fill(dm.begin(), dm.end(), 0.0);
  std::fill(vbs.begin(), vbs.end(), 0.0);
  std::fill(vgs.begin(), vgs.end(), 0.0);
  std::fill(ves.begin(), ves.end(), 0.0);
  std::fill(pis.begin(), pis.end(), 0.0);
  std::fill(pim.begin(), pim.end(), 0.0);
  std::fill(logsum.begin(), logsum.end(), 0.0);
  
  for ( int i = 0; i < q; i++) {
    r[i] = wy[i];
  }
  
  
  // Start Gibbs sampler
  std::random_device rd;
  std::mt19937 gen(seed);
  
  
  for ( int it = 0; it < nit+nburn; it++) {
    
    if ( (it > nburn) && (it % nthin == 0) ) {
      nsamples = nsamples + 1.0;
    }
    
    
    // Sample marker effects (BayesR)
    if (method == 4 || method == 5) {
      for ( int i = 0; i < m; i++) {
        if(mask[i])   continue;
        
        // variance class likelihood 
        rhs=0.0;
        for (size_t j = 0; j < LDvalues[i].size(); j++) {
          rhs +=r[LDindices[i][j]]*LDvalues[i][LDindices[i][j]];
        }
        rhs = (rhs + b[i])*(ww[i]/ve);
        logLc[0] = std::log(pi[0]);
        for (int j = 1; j<nc ; j++) {
          vbc = vb*gamma[j];
          lhs = (ww[i] / ve) + (1 / vbc);
          logLc[j] = 0.5*( (std::log(1.0/lhs)-std::log(vbc)) + (rhs*rhs)/lhs) + std::log(pi[j]); 
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
        prob[it][i] = 1.0 - probc[0];
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
          lhs = (ww[i] / ve) + (1 / vbc);
          std::normal_distribution<double> rnorm(rhs/lhs, sqrt(1.0/lhs));
          bn = rnorm(gen);
        }
        diff = (bn-b[i]);
        if(diff!=0.0) {
          for (size_t j = 0; j < LDindices[i].size(); j++) {
            r[LDindices[i][j]] += -LDvalues[i][LDindices[i][j]]*diff;
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
          if(it>nburn && (it % nthin == 0)) pim[j] = pim[j] + pi[j];
        }
        pis[it] = 1.0 - pi[0];
      }
    }
    
    // Sample marker variance
    ssb = 0.0;
    dfb = 0.0;
    if (method == 4 || method == 5) {
      for ( int i = 0; i < m; i++) {
        if(it>nburn && (it % nthin == 0)) bm[i] = bm[i] + b[i];
        if(d[i]>0)   {
          ssb += (b[i]*b[i])/gamma[d[i]];
          dfb += 1.0;
          if(it>nburn && (it % nthin == 0)) dm[i] = dm[i] + 1.0;
        }
      }
    }
    std::chi_squared_distribution<double> rchisq_vb(dfb+nub);
    chi2 = rchisq_vb(gen);
    vbs[it] = (ssb + ssb_prior*nub)/chi2 ;
    if(updateB) vb = vbs[it];
    
    // Sample residual variance
    sse = 0.0;
    for ( int i = 0; i < q; i++) {
      sse += r[i]*r[i];
    }
    dfe = (double)q + nue;
    sse= sse*dfe;
    std::chi_squared_distribution<double> rchisq_ve(dfe);
    chi2 = rchisq_ve(gen);
    ves[it] = (sse + sse_prior*nue)/chi2;
    if(updateE) ve = ves[it];
    
    // Sample genetic variance
    ssg = 0.0;
    for (int i = 0; i < q; i++) {
      ssg += (wy[i] - r[i]) * (wy[i] - r[i]);
    }
    vgs[it] = ssg;
    if(updateG) vg = ssg;
    
    // Compute logcpo
    if(it>nburn && (it % nthin == 0)) {
      double constant = 1.0 / std::sqrt(2 * M_PI);  
      for (int i = 0; i < q; i++) {
        logsum[i] += 1.0/(constant * std::exp(-0.5 * r[i] * r[i]));
      }
    }

    for ( int i = 0; i < m; i++) {
      bs[it][i] = b[i];
      ds[it][i] = d[i];
    }
  }
  
  // Compute logcpo
  logcpo=0.0;
  for (int i=0; i < q; i++) {
    logcpo += std::log(nsamples*(1/logsum[i]));
  }
  
  
  // Summarize results
  std::vector<std::vector<double>> result(14);
  
  result[0].resize(m);
  result[1].resize(m);
  result[2].resize(nit+nburn);
  result[3].resize(nit+nburn);
  result[4].resize(nit+nburn);
  result[5].resize(nit+nburn);
  result[6].resize(nit+nburn);
  result[7].resize(nc);
  result[8].resize(q);
  result[9].resize(m);
  result[10].resize(4);
  result[11].resize((nit+nburn)*m);
  result[12].resize((nit+nburn)*m);
  result[13].resize((nit+nburn)*m);
  for (int i=0; i < m; i++) {
    result[0][i] = bm[i]/nsamples;
    result[1][i] = dm[i]/nsamples;
  }
  for (int i=0; i < nit+nburn; i++) {
    result[2][i] = 0.0;
    result[3][i] = vbs[i];
    result[4][i] = vgs[i];
    result[5][i] = ves[i];
    result[6][i] = pis[i];
  }
  for (int j=0; j < nc; j++) {
    result[7][j] = pim[j]/nsamples;
  }  
  
  for (int i=0; i < q; i++) {
    result[8][i] = r[i];
  }
  
  for (int i=0; i < m; i++) {
    result[9][i] = b[i];
  }
  result[10][0] = vb;
  result[10][1] = ve;
  result[10][2] = pi[0];
  result[10][3] = logcpo;
  
  for ( int it = 0; it < nit+nburn; it++) {
    for ( int i = 0; i < m; i++) {
      result[11][it*m + i] = bs[it][i];
      result[12][it*m + i] = (double)ds[it][i];
      result[13][it*m + i] = prob[it][i];
    }
  }
  return result;
}

// // [[Rcpp::export]]
// std::vector<std::vector<double>>  sbayes( std::vector<double> wy,
//                                           std::vector<std::vector<double>> LD,
//                                           std::vector<double> b,
//                                           std::vector<double> lambda,
//                                           double yy,
//                                           double pi,
//                                           double vg,
//                                           double vb,
//                                           double ve,
//                                           double ssb_prior,
//                                           double sse_prior,
//                                           double nub,
//                                           double nue,
//                                           bool updateB,
//                                           bool updateE,
//                                           bool updatePi,
//                                           int n,
//                                           int nit,
//                                           int method) {
//   
//   // Define local variables
//   int m = b.size();
//   
//   double rhs, lhs, bn, conv, diff;
//   double rhs0, rhs1, lhs0, lhs1, like0, like1, p0;
//   double ssb, sse, ssg, dfb, dfe, chi2;
//   double xtau, tau, lambda_tau, mu_tau, z, z2, u;
//   
//   std::vector<int> d(m);
//   
//   std::vector<double> ww(m),r(m);
//   
//   std::vector<double> dm(m),bm(m);
//   std::vector<double> ves(nit),vbs(nit),pis(nit);
//   
//   std::vector<double> x2(m);
//   std::vector<int> order(m);
//   
//   // Initialize variables
//   for ( int i = 0; i < m; i++) {
//     dm[i] = 0.0;
//     bm[i] = 0.0;
//     ww[i] = LD[i][i];
//     r[i] = wy[i];
//     x2[i] = (wy[i]/ww[i])*(wy[i]/ww[i]);
//   }
//   
//   for (int i=0; i < nit; i++) {
//     vbs[i] = 0.0;
//     ves[i] = 0.0;
//     pis[i] = 0.0;
//   }
//   
//   // Establish order of markers as they are entered into the model
//   std::iota(order.begin(), order.end(), 0);
//   std::sort(  std::begin(order),
//               std::end(order),
//               [&](int i1, int i2) { return x2[i1] > x2[i2]; } );
//   
//   // Wy - W'Wb
//   for ( int i = 0; i < m; i++) {
//     if (b[i]!= 0.0) {
//       for (int j = 0; j < m; j++) {
//         r[j]=r[j] - LD[i][j]*b[i];
//       }
//     }
//   }
//   
//   // Start Gibbs sampler
//   std::random_device rd;
//   unsigned int local_seed;
//   local_seed = rd();
//   std::mt19937 gen(local_seed);
//   
//   for ( int it = 0; it < nit; it++) {
//     conv = 0.0;
//     
//     // Compute marker effects (BLUP)
//     if (method==0) {
//       for ( int isort = 0; isort < m; isort++) {
//         int i = order[isort];
//         lhs = ww[i] + lambda[i];
//         rhs = r[i] + ww[i]*b[i];
//         bn = rhs/lhs;
//         diff = bn-b[i];
//         for (int j = 0; j < m; j++) {
//           r[j]=r[j] - LD[i][j]*(diff);
//         }
//         conv = conv + diff*diff;
//         b[i] = bn;
//       }
//     }
//     
//     // Sample marker effects (Mixed, BayesA, Lasso)
//     if ( method==1 || method==2 || method==3 ) {
//       for ( int isort = 0; isort < m; isort++) {
//         int i = order[isort];
//         lhs = ww[i] + lambda[i];
//         rhs = r[i] + ww[i]*b[i];
//         std::normal_distribution<double> rnorm(rhs/lhs, sqrt(ve/lhs));
//         bn = rnorm(gen);
//         diff = bn-b[i];
//         for (int j = 0; j < m; j++) {
//           r[j]=r[j] - LD[i][j]*(diff);
//         }
//         conv = conv + diff*diff;
//         b[i] = bn;
//       }
//     }
//     
//     // Sample marker effects (BayesC)
//     if (method==4) {
//       for ( int isort = 0; isort < m; isort++) {
//         int i = order[isort];
//         lhs0 = 1.0/vb;
//         lhs1 = ww[i]/ve + 1.0/vb;
//         rhs0 = 0.0;
//         rhs1 = 0.0;
//         rhs1 = r[i]/ve + ww[i]*b[i]/ve;
//         like0 = sqrt((1.0/lhs0))*std::exp(0.5*(1.0/lhs0)*rhs0*rhs0);
//         like1 = sqrt((1.0/lhs1))*std::exp(0.5*(1.0/lhs1)*rhs1*rhs1);
//         like0 = like0*(1.0-pi);
//         like1 = like1*pi;
//         p0 = like0/(like0+like1);
//         d[i]=0;
//         std::uniform_real_distribution<double> runif(0.0, 1.0);
//         u = runif(gen);
//         if(u>p0) d[i]=1;
//         bn=0.0;
//         if(d[i]==1) {
//           std::normal_distribution<double> rnorm(rhs1/lhs1, sqrt(1.0/lhs1));
//           bn = rnorm(gen);
//         }
//         diff = bn-b[i];
//         if(diff!=0.0) {
//           for (int j = 0; j < m; j++) {
//             r[j]=r[j] - LD[i][j]*(diff);
//           }
//           conv = conv + diff*diff;
//         }
//         b[i] = bn;
//       }
//     }
//     
//     // Sample marker variance
//     ssb = 0.0;
//     dfb = 0.0;
//     for ( int i = 0; i < m; i++) {
//       bm[i] = bm[i] + b[i];
//       if(d[i]==1)   {
//         ssb = ssb + b[i]*b[i];
//         dfb = dfb + 1.0;
//         dm[i] = dm[i] + 1.0;
//       }
//     }
//     if(updateB) {
//       std::chi_squared_distribution<double> rchisq(dfb+nub);
//       chi2 = rchisq(gen);
//       vb = (ssb + ssb_prior*nub)/chi2 ;
//       vbs[it] = vb; 
//     }
//     
//     
//     // Sample residual variance
//     if(updateE) {
//       dfe = (double)n + nue;
//       ssg = 0.0;
//       sse = 0.0;
//       for ( int i = 0; i < m; i++) {
//         ssg = ssg + b[i] * (wy[i] -  r[i]);
//         sse = sse + b[i] * (r[i] + wy[i]);
//       }
//       sse = yy - sse;
//       std::chi_squared_distribution<double> rchisq(dfe);
//       chi2 = rchisq(gen);
//       ve = (sse + sse_prior*nue)/chi2 ;
//       ves[it] = ve;
//     }
//     
//     // Update lambda's for BLUP/Mixed
//     if ( method==1) {
//       for ( int i = 0; i < m; i++) {
//         lambda[i] = ve/vb;
//       }
//     }
//     
//     // Sample marker specific variance for Bayes A
//     if (method==2) {
//       dfb = 1.0 + nub;
//       for ( int i = 0; i < m; i++) {
//         ssb = b[i]*b[i];
//         std::chi_squared_distribution<double> rchisq(dfb);
//         chi2 = rchisq(gen);
//         vb = (ssb + ssb_prior*nub)/chi2 ;
//         lambda[i] = ve/vb;
//       }
//     }
//     
//     // Sample marker specific tau for Bayesian lasso
//     if (method==3) {
//       lambda_tau = 2.0*2.0*0.5*0.5/vb;
//       for ( int i = 0; i < m; i++) {
//         ssb = b[i]*b[i];
//         mu_tau = sqrt(lambda_tau/ssb);
//         std::normal_distribution<double> norm(0.0, 1.0);
//         z = norm(gen);
//         z2=z*z;
//         xtau=mu_tau+0.5*mu_tau*mu_tau*z2/lambda_tau - 0.5*(mu_tau/lambda_tau)*sqrt(4*mu_tau*lambda_tau*z2+mu_tau*mu_tau*z2*z2);
//         std::uniform_real_distribution<double> runif(0.0, 1.0);
//         u = runif(gen);
//         tau = mu_tau*mu_tau/xtau;
//         if(u <= mu_tau/(mu_tau+xtau)) tau=xtau;
//         lambda[i] = ve/tau;
//       }
//     }
//     
//     // Sample pi for Bayes C
//     if(method==4 && updatePi) {
//       double count = dfb + 1.0;
//       std::gamma_distribution<double> rgamma(count,1.0);
//       double rg = rgamma(gen);
//       pi = rg/(double)m;
//       pis[it] = pi;
//     }
//   }
//   
//   // Summarize results
//   std::vector<std::vector<double>> result(10);
//   result[0].resize(m);
//   result[1].resize(m);
//   result[2].resize(nit);
//   result[3].resize(nit);
//   result[4].resize(nit);
//   result[5].resize(nit);
//   result[6].resize(m);
//   result[7].resize(m);
//   result[8].resize(3);
//   result[9].resize(m);
//   
//   for (int i=0; i < m; i++) {
//     result[0][i] = bm[i]/nit;
//     result[1][i] = dm[i]/nit;
//   }
//   for (int i=0; i < nit; i++) {
//     //result[2][i] = mus[i];
//     result[2][i] = 0.0;
//     result[3][i] = vbs[i];
//     result[4][i] = ves[i];
//     result[5][i] = pis[i];
//   }
//   for (int i=0; i < m; i++) {
//     result[6][i] = wy[i];
//     result[7][i] = r[i];
//   }
//   result[8][0] = vb;
//   result[8][1] = ve;
//   result[8][2] = pi;
//   for (int i=0; i < m; i++) {
//     result[9][i] = b[i];
//   }
//   
//   
//   return result;
// }
