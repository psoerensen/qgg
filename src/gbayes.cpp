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
                                          double vara, 
                                          double varb, 
                                          double vare,
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
  
  double pi0=pi;
  //double vara0=vara;
  
  double rhs, lhs, bnew, conv, diff, mu;
  double rhs0, rhs1, lhs0, lhs1, like0, like1, p0;
  double ssb, sse, dfb, dfe, chi2;
  double xtau, tau, lambda_tau, mu_tau, z, z2, u;
  
  std::vector<int> d(m);
  
  std::vector<double> e(n);
  std::vector<double> wy(m),ww(m),x2(m);
  std::vector<int> order(m);
  
  std::vector<double> d_post_mean(m),b_post_mean(m);
  std::vector<double> vare_post(nit),varb_post(nit),pi_post(nit),mu_post(nit);
  
  // Prior variance and degrees of freedom
  dfe = n + nue;
  //ssb_prior = nub*varb;
  //ssb_prior =  (nub-2.0)/nub * (vara0/(pi0*m*0.5));
  //sse_prior = nue*vare;
  
  // Mean adjust y
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
      d[i] = 1;
      d_post_mean[i] = 0.0;
      b_post_mean[i] = 0.0;
    }
    x2[i] = (wy[i]/ww[i])*(wy[i]/ww[i]);
  }
  
  std::iota(order.begin(), order.end(), 0);
  std::sort(  std::begin(order), 
              std::end(order),
              [&](int i1, int i2) { return x2[i1] > x2[i2]; } );
  
  for (int i=0; i < nit; i++) {
    varb_post[i] = 0.0;
    vare_post[i] = 0.0;
    pi_post[i] = 0.0;
  }
  

  


  // Start Gibbs sampler
  
  std::cout << "  " << "\n";
  std::cout << "Starting Gibbs sampler" << "\n";
  std::cout << "  " << "\n";
  

  std::random_device rd;
  unsigned int local_seed;
  local_seed = rd();
  std::mt19937 gen(local_seed);
  
  for ( int it = 0; it < nit; it++) {
    conv = 0.0;
    
    // Sample marker effects (Mixed, BayesA, Lasso)
    if (method<=3) {
      for ( int i0 = 0; i0 < m; i0++) {
        int i = order[i0];
        lhs = ww[i] + lambda[i];
        rhs = 0.0;
        for ( int j = 0; j < n; j++) {
          rhs = rhs + W[i][j]*e[j]; 
        }
        rhs = rhs + ww[i]*b[i];
        std::normal_distribution<double> rnorm_b(rhs/lhs, sqrt(vare/lhs));
        bnew = rnorm_b(gen);
        diff = bnew-b[i];
        for (int j = 0; j < n; j++) {
          e[j]=e[j] - W[i][j]*(diff);
        }
        conv = conv + diff*diff;
        b[i] = bnew;
      }
    }
    // Sample marker effects (BayesC)
    if (method==4) {
      for ( int i0 = 0; i0 < m; i0++) {
        int i = order[i0];
        lhs0 = 1.0/varb;
        lhs1 = ww[i]/vare + 1.0/varb;
        rhs0 = 0.0;
        rhs1 = 0.0;
        for ( int j = 0; j < n; j++) {
          rhs1 = rhs1 + W[i][j]*e[j]/vare; 
        }
        rhs1 = rhs1 + ww[i]*b[i]/vare;
        like0 = sqrt((1.0/lhs0))*std::exp(0.5*(1.0/lhs0)*rhs0*rhs0);
        like1 = sqrt((1.0/lhs1))*std::exp(0.5*(1.0/lhs1)*rhs1*rhs1);
        like0 = like0*(1.0-pi); 
        like1 = like1*pi;
        p0 = like0/(like0+like1); 
        //p1 = like1/(like0+like1);
        d[i]=0;
        std::uniform_real_distribution<double> runif(0.0, 1.0);
        u = runif(gen);
        if(u>p0) d[i]=1;
        bnew=0.0;
        if(d[i]==1) {
          std::normal_distribution<double> rnorm_b(rhs1/lhs1, sqrt(1.0/lhs1));
          //std::normal_distribution<double> rnorm_b(rhs1/lhs1, sqrt(vare/lhs1));
          bnew = rnorm_b(gen);
        } 
        diff = bnew-b[i];
        for (int j = 0; j < n; j++) {
          e[j]=e[j] - W[i][j]*(diff);
        }
        conv = conv + diff*diff;
        b[i] = bnew;
      }
    }
    
    //std::cout << "Sampling variance components: " << it + 1 << "\n";
    
    // Sample marker variance
    ssb = 0.0;
    dfb = 0.0;
    for ( int i = 0; i < m; i++) {
      b_post_mean[i] = b_post_mean[i] + b[i];     
      if(d[i]==1)   {
        ssb = ssb + b[i]*b[i];  
        dfb = dfb + 1.0;
        d_post_mean[i] = d_post_mean[i] + 1.0;     
      }
    }
    if(updateB) {
      dfb = dfb + nub;
      std::chi_squared_distribution<double> rchisq_b(dfb);
      chi2 = rchisq_b(gen);
      varb = (ssb + ssb_prior)/chi2 ;
      varb_post[it] = varb;
    }
    

    // Sample residual variance
    if(updateE) {
      sse = 0.0;
      for ( int j = 0; j < n; j++) {
        sse = sse + e[j]*e[j];
      }
      std::chi_squared_distribution<double> rchisq_e(dfe);
      chi2 = rchisq_e(gen);
      vare = (sse + sse_prior)/chi2 ;
      vare_post[it] = vare;
    }
    
    // Update lambda's
    if ( method<2 || method==4 ) {
      for ( int i = 0; i < m; i++) {
        lambda[i] = vare/varb;
      }
    }
    
    // Bayes A
    if (method==2) {
      dfb = 1.0 + nub;
      for ( int i = 0; i < m; i++) { 
        ssb = b[i]*b[i];
        std::chi_squared_distribution<double> rchisqb(dfb);
        chi2 = rchisqb(gen);
        varb = (ssb + ssb_prior)/chi2 ;
        lambda[i] = vare/varb;
      }
    }
    
    // Bayesian lasso
    if (method==3) { 
      lambda_tau = 2.0*2.0*0.5*0.5/varb;
      for ( int i = 0; i < m; i++) { 
        ssb = b[i]*b[i];
        mu_tau = sqrt(lambda_tau/ssb);
        // sample from inverse gauss
        std::normal_distribution<double> dist(0.0, 1.0);
        z = dist(gen);
        z2=z*z;
        xtau=mu_tau+0.5*mu_tau*mu_tau*z2/lambda_tau - 0.5*(mu_tau/lambda_tau)*sqrt(4*mu_tau*lambda_tau*z2+mu_tau*mu_tau*z2*z2);
        std::uniform_real_distribution<double> runif(0.0, 1.0);
        u = runif(gen);
        tau = mu_tau*mu_tau/xtau;
        //if(u <= mu/(mu+xtau)) tau=xtau;
        if(u <= mu_tau/(mu_tau+xtau)) tau=xtau;
        lambda[i] = vare/tau;
      }
    }

    // Sample pi
    if(updatePi) {
      //double beta_a = 1.0 + (double)dfb;
      //double beta_b = (double)m - (double)dfb; 
      //pi = (double) R::rbeta(beta_a,beta_b);
      //double count0 = (double) m - dfb + 1.0;
      double count1 = dfb + 1.0;
      //std::gamma_distribution<double> rgamma0(count0,1.0);
      //double rg0 = rgamma0(gen);
      std::gamma_distribution<double> rgamma1(count1,1.0);
      double rg1 = rgamma1(gen);
      pi = rg1/(double)m;
      pi_post[it] = pi;
      //std::cout << "Pi: " << pi << "\n";
      //ssb_prior =  (nub-2.0)/nub * vara0/(pi0*m*0.5);
    }
    
    
    // Update mu and adjust residuals
    mu = std::accumulate(e.begin(), e.end(), 0.0)/n;
    for ( int i = 0; i < n; i++) {
      e[i] = e[i] - mu;
    }                
    mu_post[it] = mu;
    
    //std::cout << "Finished iteration: " << it + 1 << "\n";
    //std::cout << "Convergence: " << conv << "\n";
    
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
    result[0][i] = b_post_mean[i]/nit;
    result[1][i] = d_post_mean[i]/nit;
  }
  for (int i=0; i < nit; i++) {
    result[2][i] = mu_post[i];
    result[3][i] = varb_post[i];
    result[4][i] = vare_post[i];
    result[5][i] = pi_post[i];
  }
  for (int i=0; i < n; i++) {
    result[6][i] = y[i]- mu - e[i];
    result[7][i] = e[i];
  }
  result[8][0] = pi;
  result[8][1] = varb;
  result[8][2] = vare;
  for (int i=0; i < m; i++) {
    result[9][i] = b[i];
  }
  
  return result;
}


// [[Rcpp::export]]
std::vector<std::vector<double>>  fbayes(  std::vector<double> y,
                                          std::vector<std::vector<double>> W, 
                                          std::vector<std::vector<double>> LD, 
                                          std::vector<double> b, 
                                          std::vector<double> lambda, 
                                          double pi, 
                                          double vara, 
                                          double varb, 
                                          double vare,
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
  
  double pi0=pi;
  double vara0=vara;
  
  double rhs, lhs, bnew, conv, diff, mu;
  double rhs0, rhs1, lhs0, lhs1, like0, like1, p0;
  double yy, ssb, sse, ssg, ssb_prior, sse_prior, dfb, dfe, chi2;
  double xtau, tau, lambda_tau, mu_tau, z, z2, u;
  
  std::vector<int> d(m);
  
  std::vector<double> e(n);
  std::vector<double> ww(m),r(m),wy(m);
  
  std::vector<double> d_post_mean(m),b_post_mean(m);
  std::vector<double> vare_post(nit),varb_post(nit),pi_post(nit),mu_post(nit);
  
  std::vector<double> x2(m);
  std::vector<int> order(m);
  
  //std::vector<std::vector<double>> LD(m, std::vector<double>(m, 0.0));
  
  
  // Prior variance and degrees of freedom
  dfe = n + nue;
  //ssb_prior = nub*varb;
  ssb_prior =  (nub-2.0)/nub * (vara0/(pi0*m*0.5));
  sse_prior = nue*vare;
  

  //std::clock_t startcputime = std::clock();
  
  // Initialize variables
  for (int i = 0; i < m; i++) {
    ww[i] = 0.0;
    for (int j = 0; j < n; j++) {
      ww[i] = ww[i] + W[i][j]*W[i][j];
      d[i] = 1;
      d_post_mean[i] = 0.0;
      b_post_mean[i] = 0.0;
    }
  }
  
  
  // Mean adjust y
  mu = std::accumulate(y.begin(), y.end(), 0.0)/n;
  yy = 0.0;
  for ( int i = 0; i < n; i++) {
    e[i] = y[i] - mu;
    yy = yy + e[i]*e[i];
  }                

  for ( int i = 0; i < m; i++) {
    r[i] = 0.0;
    for (int j = 0; j < n; j++) {
      r[i] = r[i] + W[i][j]*e[j];
    }
    wy[i] = r[i];
    x2[i] = (wy[i]/ww[i])*(wy[i]/ww[i]);
  }

  std::iota(order.begin(), order.end(), 0);
  std::sort(  std::begin(order), 
              std::end(order),
              [&](int i1, int i2) { return x2[i1] > x2[i2]; } );
  
  

  //for ( int i = 0; i < m; i++) {
  //  for (int j = i; j < m; j++) {
  //    LD[i][j] = 0.0;
  //    for (int k = 0; k < n; k++) {
  //      LD[i][j] = LD[i][j] + W[i][k]*W[j][k];
  //    }
  //    LD[j][i] = LD[i][j];
  //  }
  //}                
  
  //double cpu_duration = (std::clock() - startcputime)/ (double)CLOCKS_PER_SEC;
  //Rcout << "Finished in " << cpu_duration << " seconds [CPU Clock] " << "\n";

    
  // Start Gibbs sampler
  
  //std::cout << "  " << "\n";
  //std::cout << "Starting solver" << "\n";
  //std::cout << "  " << "\n";
  Rcout << "  " << "\n";
  Rcout << "Starting Gibbs sampler" << "\n";
  Rcout << "  " << "\n";

  //startcputime = std::clock();
  
  std::random_device rd;
  unsigned int local_seed;
  local_seed = rd();
  std::mt19937 gen(local_seed);
  
  for ( int it = 0; it < nit; it++) {
    conv = 0.0;
    
    // Sample marker effects (Mixed, BayesA, Lasso)
    if (method<=3) {
      for ( int i0 = 0; i0 < m; i0++) {
        int i = order[i0];
        lhs = ww[i] + lambda[i];
        rhs = r[i] + ww[i]*b[i];
        std::normal_distribution<double> rnorm_b(rhs/lhs, sqrt(vare/lhs));
        bnew = rnorm_b(gen);
        diff = bnew-b[i];
        for (int j = 0; j < m; j++) {
          r[j]=r[j] - LD[i][j]*(diff);
        }
        conv = conv + diff*diff;
        b[i] = bnew;
      }
    }
    // Sample marker effects (BayesC)
    if (method==4) {
      for ( int i0 = 0; i0 < m; i0++) {
        int i = order[i0];
        lhs0 = 1.0/varb;
        lhs1 = ww[i]/vare + 1.0/varb;
        rhs0 = 0.0;
        rhs1 = 0.0;
        rhs1 = r[i]/vare + ww[i]*b[i]/vare;
        
        like0 = sqrt((1.0/lhs0))*std::exp(0.5*(1.0/lhs0)*rhs0*rhs0);
        like1 = sqrt((1.0/lhs1))*std::exp(0.5*(1.0/lhs1)*rhs1*rhs1);
        like0 = like0*(1.0-pi); 
        like1 = like1*pi;
        p0 = like0/(like0+like1); 
        //p1 = like1/(like0+like1);
        d[i]=0;
        std::uniform_real_distribution<double> runif(0.0, 1.0);
        u = runif(gen);
        if(u>p0) d[i]=1;
        bnew=0.0;
        if(d[i]==1) {
          std::normal_distribution<double> rnorm_b(rhs1/lhs1, sqrt(1.0/lhs1));
          bnew = rnorm_b(gen);
        } 
        diff = bnew-b[i];
        if(diff!=0.0) {
          for (int j = 0; j < m; j++) {
            r[j]=r[j] - LD[i][j]*(diff);
          }
          conv = conv + diff*diff;
        }
        b[i] = bnew;
      }
    }
    
    //std::cout << "Sampling variance components: " << it + 1 << "\n";
    
    // Sample marker variance
    ssb = 0.0;
    dfb = 0.0;
    for ( int i = 0; i < m; i++) {
      b_post_mean[i] = b_post_mean[i] + b[i];     
      if(d[i]==1)   {
        ssb = ssb + b[i]*b[i];  
        dfb = dfb + 1.0;
        d_post_mean[i] = d_post_mean[i] + 1.0;     
      }
    }
    dfb = dfb + nub;
    std::chi_squared_distribution<double> rchisq_b(dfb);
    chi2 = rchisq_b(gen);
    varb = (ssb + ssb_prior)/chi2 ;
    varb_post[it] = varb;

    // Sample residual variance
    ssg = 0.0;
    sse = 0.0;
    //for ( int j = 0; j < n; j++) {
    //  sse = sse + e[j]*e[j];
    //}
    for ( int i = 0; i < m; i++) {
      ssg = ssg + b[i] * (wy[i] -  r[i]);
      sse = sse + b[i] * (r[i] + wy[i]);
    }
    sse = yy - sse;
    
    std::chi_squared_distribution<double> rchisq_e(dfe);
    chi2 = rchisq_e(gen);
    vare = (sse + sse_prior)/chi2 ;
    vare_post[it] = vare;
    
    
    // Update lambda's
    if ( method<2 || method==4 ) {
      for ( int i = 0; i < m; i++) {
        lambda[i] = vare/varb;
      }
    }
    
    // Bayes A
    if (method==2) {
      dfb = 1.0 + nub;
      for ( int i = 0; i < m; i++) { 
        ssb = b[i]*b[i];
        std::chi_squared_distribution<double> rchisqb(dfb);
        chi2 = rchisqb(gen);
        varb = (ssb + ssb_prior)/chi2 ;
        lambda[i] = vare/varb;
      }
    }
    
    // Bayesian lasso
    if (method==3) { 
      lambda_tau = 2.0*2.0*0.5*0.5/varb;
      for ( int i = 0; i < m; i++) { 
        ssb = b[i]*b[i];
        mu_tau = sqrt(lambda_tau/ssb);
        // sample from inverse gauss
        std::normal_distribution<double> dist(0.0, 1.0);
        z = dist(gen);
        z2=z*z;
        xtau=mu_tau+0.5*mu_tau*mu_tau*z2/lambda_tau - 0.5*(mu_tau/lambda_tau)*sqrt(4*mu_tau*lambda_tau*z2+mu_tau*mu_tau*z2*z2);
        std::uniform_real_distribution<double> runif(0.0, 1.0);
        u = runif(gen);
        tau = mu_tau*mu_tau/xtau;
        //if(u <= mu/(mu+xtau)) tau=xtau;
        if(u <= mu_tau/(mu_tau+xtau)) tau=xtau;
        lambda[i] = vare/tau;
      }
    }

    // Sample pi
    if(updatePi) {
      //double beta_a = 1.0 + (double)dfb;
      //double beta_b = (double)m - (double)dfb; 
      //pi = (double) R::rbeta(beta_a,beta_b);
      //double count0 = (double) m - dfb + 1.0;
      double count1 = dfb + 1.0;
      //std::gamma_distribution<double> rgamma0(count0,1.0);
      //double rg0 = rgamma0(gen);
      std::gamma_distribution<double> rgamma1(count1,1.0);
      double rg1 = rgamma1(gen);
      pi = rg1/(double)m;
      //std::cout << "Pi: " << pi << "\n";
      pi_post[it] = pi;
    }
    
    // Update mu and adjust residuals
    mu = std::accumulate(e.begin(), e.end(), 0.0)/n;
    for ( int i = 0; i < n; i++) {
      e[i] = e[i] - mu;
    }                
    mu_post[it] = mu;
    
    //std::cout << "Finished iteration: " << it + 1 << "\n";
    //std::cout << "Convergence: " << conv << "\n";
    
  }
  
  //cpu_duration = (std::clock() - startcputime)/ (double)CLOCKS_PER_SEC;
  //Rcout << "Finished in " << cpu_duration << " seconds [CPU Clock] " << "\n";
  

  // Summarize results
  std::vector<std::vector<double>> result(6);
  result[0].resize(m);
  result[1].resize(m);
  result[2].resize(nit);
  result[3].resize(nit);
  result[4].resize(nit);
  result[5].resize(nit);
  
  for (int i=0; i < m; i++) {
    result[0][i] = b_post_mean[i]/nit;
    result[1][i] = d_post_mean[i]/nit;
  }
  for (int i=0; i < nit; i++) {
    result[2][i] = mu_post[i];
    result[3][i] = varb_post[i];
    result[4][i] = vare_post[i];
    result[5][i] = pi_post[i];
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
                                          double vara, 
                                          double varb, 
                                          double vare,
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
  //int n = y.size();
  int m = b.size();
  
  double pi0=pi;
  //double vara0=vara;
  
  double rhs, lhs, bnew, conv, diff;
  double rhs0, rhs1, lhs0, lhs1, like0, like1, p0;
  double ssb, sse, ssg, dfb, dfe, chi2;
  double xtau, tau, lambda_tau, mu_tau, z, z2, u;
  
  std::vector<int> d(m);
  
  std::vector<double> ww(m),r(m);
  
  std::vector<double> d_post_mean(m),b_post_mean(m);
  std::vector<double> vare_post(nit),varb_post(nit),pi_post(nit);
  
  std::vector<double> x2(m);
  std::vector<int> order(m);
  
  //std::vector<std::vector<double>> LD(m, std::vector<double>(m, 0.0));
  
  
  // Prior variance and degrees of freedom
  dfe = n + nue;

  
  //std::clock_t startcputime = std::clock();
  
  // Initialize variables
  //yy = 0.0;
  for ( int i = 0; i < m; i++) {
    d[i] = 1;
    d_post_mean[i] = 0.0;
    b_post_mean[i] = 0.0;
    //ww[i] = (double)n;
    ww[i] = LD[i][i];
    r[i] = wy[i];
    x2[i] = (wy[i]/ww[i])*(wy[i]/ww[i]);
    //yy = yy + seb2[i]*ww[i]*(n-2.0) + b2[i]*ww[i];
  }
  //yy = yy/m;

  
  // Wy - W'Wb
  for ( int i = 0; i < m; i++) {
    if (b[i]!= 0.0) {
      for (int j = 0; j < m; j++) {
        r[j]=r[j] - LD[i][j]*b[i];
      }
    }
  }
  
  for (int i=0; i < nit; i++) {
    varb_post[i] = 0.0;
    vare_post[i] = 0.0;
    pi_post[i] = 0.0;
  }
  
  std::iota(order.begin(), order.end(), 0);
  std::sort(  std::begin(order), 
              std::end(order),
              [&](int i1, int i2) { return x2[i1] > x2[i2]; } );
  
  
  // Start Gibbs sampler
  
  //std::cout << "  " << "\n";
  //std::cout << "Starting solver" << "\n";
  //std::cout << "  " << "\n";
  //Rcout << "  " << "\n";
  //Rcout << "Starting Gibbs sampler" << "\n";
  //Rcout << "  " << "\n";
  
  //startcputime = std::clock();
  
  std::random_device rd;
  unsigned int local_seed;
  local_seed = rd();
  std::mt19937 gen(local_seed);
  
  for ( int it = 0; it < nit; it++) {
    conv = 0.0;

    // Compute marker effects (BLUP)
    if (method==0) {
      for ( int i0 = 0; i0 < m; i0++) {
        int i = order[i0];
        lhs = ww[i] + lambda[i];
        rhs = r[i] + ww[i]*b[i];
        bnew = rhs/lhs;
        diff = bnew-b[i];
        for (int j = 0; j < m; j++) {
          r[j]=r[j] - LD[i][j]*(diff);
        }
        conv = conv + diff*diff;
        b[i] = bnew;
      }
    }
    
    // Sample marker effects (Mixed, BayesA, Lasso)
    if ( method==1 || method==2 || method==3 ) {
      for ( int i0 = 0; i0 < m; i0++) {
        int i = order[i0];
        lhs = ww[i] + lambda[i];
        rhs = r[i] + ww[i]*b[i];
        std::normal_distribution<double> rnorm_b(rhs/lhs, sqrt(vare/lhs));
        bnew = rnorm_b(gen);
        diff = bnew-b[i];
        for (int j = 0; j < m; j++) {
          r[j]=r[j] - LD[i][j]*(diff);
        }
        conv = conv + diff*diff;
        b[i] = bnew;
      }
    }
    // Sample marker effects (BayesC)
    if (method==4) {
      for ( int i0 = 0; i0 < m; i0++) {
        int i = order[i0];
        lhs0 = 1.0/varb;
        lhs1 = ww[i]/vare + 1.0/varb;
        rhs0 = 0.0;
        rhs1 = 0.0;
        rhs1 = r[i]/vare + ww[i]*b[i]/vare;
        
        like0 = sqrt((1.0/lhs0))*std::exp(0.5*(1.0/lhs0)*rhs0*rhs0);
        like1 = sqrt((1.0/lhs1))*std::exp(0.5*(1.0/lhs1)*rhs1*rhs1);
        like0 = like0*(1.0-pi); 
        like1 = like1*pi;
        p0 = like0/(like0+like1); 
        //p1 = like1/(like0+like1);
        d[i]=0;
        std::uniform_real_distribution<double> runif(0.0, 1.0);
        u = runif(gen);
        if(u>p0) d[i]=1;
        bnew=0.0;
        if(d[i]==1) {
          std::normal_distribution<double> rnorm_b(rhs1/lhs1, sqrt(1.0/lhs1));
          bnew = rnorm_b(gen);
        } 
        diff = bnew-b[i];
        if(diff!=0.0) {
          for (int j = 0; j < m; j++) {
            r[j]=r[j] - LD[i][j]*(diff);
          }
          conv = conv + diff*diff;
        }
        b[i] = bnew;
      }
    }
    
    //std::cout << "Sampling variance components: " << it + 1 << "\n";
    
    // Sample marker variance
    ssb = 0.0;
    dfb = 0.0;
    for ( int i = 0; i < m; i++) {
      b_post_mean[i] = b_post_mean[i] + b[i];     
      if(d[i]==1)   {
        ssb = ssb + b[i]*b[i];  
        dfb = dfb + 1.0;
        d_post_mean[i] = d_post_mean[i] + 1.0;     
      }
    }
    dfb = dfb + nub;
    std::chi_squared_distribution<double> rchisq_b(dfb);
    chi2 = rchisq_b(gen);
    if(updateB) {
      varb = (ssb + ssb_prior)/chi2 ;
    }
    varb_post[it] = varb;
    
    
    
    // Sample residual variance
    if(updateE) {
      ssg = 0.0;
      sse = 0.0;
      for ( int i = 0; i < m; i++) {
        ssg = ssg + b[i] * (wy[i] -  r[i]);
        sse = sse + b[i] * (r[i] + wy[i]);
      }
      sse = yy - sse;
      
      std::chi_squared_distribution<double> rchisq_e(dfe);
      chi2 = rchisq_e(gen);
      vare = (sse + sse_prior)/chi2 ;
      vare_post[it] = vare;
    }
    

    
    // Update lambda's
    if ( method==1 || method==4 ) {
      for ( int i = 0; i < m; i++) {
        lambda[i] = vare/varb;
      }
    }
    
    // Bayes A
    if (method==2 && updateB) {
      dfb = 1.0 + nub;
      for ( int i = 0; i < m; i++) { 
        ssb = b[i]*b[i];
        std::chi_squared_distribution<double> rchisqb(dfb);
        chi2 = rchisqb(gen);
        varb = (ssb + ssb_prior)/chi2 ;
        lambda[i] = vare/varb;
      }
    }
    
    // Bayesian lasso
    if (method==3 && updateB) { 
      lambda_tau = 2.0*2.0*0.5*0.5/varb;
      for ( int i = 0; i < m; i++) { 
        ssb = b[i]*b[i];
        mu_tau = sqrt(lambda_tau/ssb);
        // sample from inverse gauss
        std::normal_distribution<double> dist(0.0, 1.0);
        z = dist(gen);
        z2=z*z;
        xtau=mu_tau+0.5*mu_tau*mu_tau*z2/lambda_tau - 0.5*(mu_tau/lambda_tau)*sqrt(4*mu_tau*lambda_tau*z2+mu_tau*mu_tau*z2*z2);
        std::uniform_real_distribution<double> runif(0.0, 1.0);
        u = runif(gen);
        tau = mu_tau*mu_tau/xtau;
        //if(u <= mu/(mu+xtau)) tau=xtau;
        if(u <= mu_tau/(mu_tau+xtau)) tau=xtau;
        lambda[i] = vare/tau;
      }
    }
    
    // Sample pi
    if(updatePi) {
      double count1 = dfb + 1.0;
      std::gamma_distribution<double> rgamma1(count1,1.0);
      double rg1 = rgamma1(gen);
      pi = rg1/(double)m;
      pi_post[it] = pi;
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
    result[0][i] = b_post_mean[i]/nit;
    result[1][i] = d_post_mean[i]/nit;
  }
  for (int i=0; i < nit; i++) {
    //result[2][i] = mu_post[i];
    result[2][i] = 0.0;
    result[3][i] = varb_post[i];
    result[4][i] = vare_post[i];
    result[5][i] = pi_post[i];
  }
  for (int i=0; i < m; i++) {
    result[6][i] = wy[i];
    result[7][i] = r[i];
  }
  result[8][0] = pi;
  result[8][1] = varb;
  result[8][2] = vare;
  for (int i=0; i < m; i++) {
    result[9][i] = b[i];
  }
  
  // // Summarize results
  // std::vector<std::vector<double>> result(6);
  // result[0].resize(m);
  // result[1].resize(m);
  // result[2].resize(nit);
  // result[3].resize(nit);
  // result[4].resize(nit);
  // result[5].resize(nit);
  // 
  // for (int i=0; i < m; i++) {
  //   result[0][i] = b_post_mean[i]/nit;
  //   result[1][i] = d_post_mean[i]/nit;
  // }
  // for (int i=0; i < nit; i++) {
  //   result[2][i] = 0.0;
  //   result[3][i] = varb_post[i];
  //   result[4][i] = vare_post[i];
  //   result[5][i] = pi_post[i];
  // }
  
  return result;
}
