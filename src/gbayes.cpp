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



// // [[Rcpp::export]]
// std::vector<std::vector<double>>  bayes(   std::vector<double> y,
//                                            std::vector<std::vector<double>> W, 
//                                            std::vector<double> b, 
//                                            std::vector<double> lambda, 
//                                            std::vector<double> pi, 
//                                            std::vector<double> gamma, 
//                                            double vb, 
//                                            double vg, 
//                                            double ve,
//                                            double ssb_prior,
//                                            double ssg_prior,
//                                            double sse_prior,
//                                            double nub,
//                                            double nug,
//                                            double nue,
//                                            bool updateB,
//                                            bool updateG,
//                                            bool updateE,
//                                            bool updatePi,
//                                            int nit,
//                                            int nburn,
//                                            int nthin,
//                                            int method,
//                                            int seed) {
//   
//   // Define local variables
//   int n = y.size();
//   int m = W.size();
//   int nc = pi.size();
//   double nsamples=0.0;
//   
//   double rhs, lhs, bn, diff, mu;
//   double rhs0, rhs1, lhs0, lhs1, like0, like1, p0;
//   double ssb, ssg, sse, dfb, dfe, dfg, chi2;
//   double lambda_tau, mu_tau, u, lambda2;
// 
//   std::vector<int> d(m);
//   
//   std::vector<double> e(n),g(n);
//   std::vector<double> wy(m),ww(m),x2(m);
//   std::vector<int> order(m), mask(m);
//   
//   std::vector<double> dm(m),bm(m),vbi(m),vei(m);
//   std::vector<double> ves(nit+nburn),vbs(nit+nburn),vgs(nit+nburn),pis(nit+nburn),mus(nit+nburn);
//   
//   std::vector<double> pim(nc), probc(nc), logLc(nc);
//   double cumprobc, vbc, v0, v1, logLcAdj;
//   
//   std::fill(g.begin(), g.end(), 0.0);
//   std::fill(wy.begin(), wy.end(), 0.0);
//   std::fill(ww.begin(), ww.end(), 0.0);
//   std::fill(dm.begin(), dm.end(), 0.0);
//   std::fill(bm.begin(), bm.end(), 0.0);
//   std::fill(vbs.begin(), vbs.end(), 0.0);
//   std::fill(vgs.begin(), vgs.end(), 0.0);
//   std::fill(ves.begin(), ves.end(), 0.0);
//   std::fill(pis.begin(), pis.end(), 0.0);
//   std::fill(pim.begin(), pim.end(), 0.0);
//   
// 
//   // Mean adjust y and initialize e
//   mu = std::accumulate(y.begin(), y.end(), 0.0)/n;
//   for ( int i = 0; i < n; i++) {
//     e[i] = y[i] - mu;
//   }                
//   
//   // Initialize variables
//   for (int i = 0; i < m; i++) {
//     for (int j = 0; j < n; j++) {
//       wy[i] = wy[i] + W[i][j]*e[j];
//       ww[i] = ww[i] + W[i][j]*W[i][j];
//     }
//     x2[i] = (wy[i]/ww[i])*(wy[i]/ww[i]);
//     mask[i]=1;
//     vbi[i]=vb;
//   }
//   
//   lambda2 = m;
//   
//   // Establish order of markers as they are entered into the model
//   std::iota(order.begin(), order.end(), 0);
//   std::sort(  std::begin(order), 
//               std::end(order),
//               [&](int i1, int i2) { return x2[i1] > x2[i2]; } );
//   
//   
//   // Start Gibbs sampler
//   std::random_device rd;
//   std::mt19937 gen(seed);
// 
//   for ( int it = 0; it < nit+nburn; it++) {
//     
//     if ( (it > nburn) && (it % nthin == 0) ) {
//       nsamples = nsamples + 1.0;
//     }
//     
//     // Compute marker effects (BLUP)
//     if (method==0) {
//       for ( int isort = 0; isort < m; isort++) {
//         int i = order[isort];
//         if(!mask[i])   continue;
//         rhs = 0.0;
//         for ( int j = 0; j < n; j++) {
//           rhs = rhs + W[i][j]*e[j];
//         }
//         rhs = rhs/ve + ww[i]*b[i]/ve;
//         lhs = ww[i]/ve + 1.0/vb;
//         bn = rhs/lhs;
//         diff = bn-b[i];
//         for (int j = 0; j < n; j++) {
//           e[j]=e[j] - W[i][j]*(diff);
//         }
//         b[i] = bn;
//       }
//     }
//     
//     // Compute marker effects (BayesN or BayesRR)
//     if (method==1) {
//       for ( int isort = 0; isort < m; isort++) {
//         int i = order[isort];
//         if(!mask[i])   continue;
//         rhs = 0.0;
//         for ( int j = 0; j < n; j++) {
//           rhs = rhs + W[i][j]*e[j];
//         }
//         //rhs = rhs/ve + ww[i]*b[i]/ve;
//         //lhs = ww[i]/ve + 1.0/vb;
//         rhs = rhs + ww[i]*b[i];
//         lhs = ww[i] + ve/vb;
//         std::normal_distribution<double> rnorm(rhs/lhs, sqrt(ve/lhs));
//         bn = rnorm(gen);
//         diff = bn-b[i];
//         for (int j = 0; j < n; j++) {
//           e[j]=e[j] - W[i][j]*(diff);
//         }
//         b[i] = bn;
//       }
//     }
// 
//     // Compute marker effects (BayesA)
//     if (method==2) {
//       dfb = 1.0 + nub;
//       for ( int isort = 0; isort < m; isort++) {
//         int i = order[isort];
//         if(!mask[i])   continue;
//         ssb = b[i]*b[i];
//         std::chi_squared_distribution<double> rchisq(dfb);
//         chi2 = rchisq(gen);
//         vbi[i] = (ssb + ssb_prior*nub)/chi2 ;
//         rhs = 0.0;
//         for ( int j = 0; j < n; j++) {
//           rhs = rhs + W[i][j]*e[j];
//         }
//         rhs = rhs + ww[i]*b[i];
//         lhs = ww[i] + ve/vbi[i];
//         std::normal_distribution<double> rnorm(rhs/lhs, sqrt(ve/lhs));
//         bn = rnorm(gen);
//         diff = bn-b[i];
//         for (int j = 0; j < n; j++) {
//           e[j]=e[j] - W[i][j]*(diff);
//         }
//         b[i] = bn;
//       }
//     }
//     
//     // Sample marker effects (Lasso)
//     if ( method==3 ) {
//       for ( int isort = 0; isort < m; isort++) {
//         int i = order[isort];
//         // Legarra et al 2011
//         lhs = ww[i]/ve + lambda[i];
//         // Park & Casella - Campos et al.
//         //lhs = ww[i]/ve + lambda[i]/ve;
//         rhs = 0.0;
//         for ( int j = 0; j < n; j++) {
//           rhs = rhs + W[i][j]*e[j];
//         }
//         rhs = rhs + ww[i]*b[i];
//         std::normal_distribution<double> rnorm((rhs/ve)/lhs, sqrt(1.0/lhs));
//         bn = rnorm(gen);
//         diff = bn-b[i];
//         for (int j = 0; j < n; j++) {
//           e[j]=e[j] - W[i][j]*(diff);
//         }
//         b[i] = bn;
//       }
//     }
// 
// 
//    // Sample marker effects (BayesC)
//     if (method==4) {
//       for ( int isort = 0; isort < m; isort++) {
//         int i = order[isort];
//         rhs0 = 0.0;
//         rhs1 = 0.0;
//         for ( int j = 0; j < n; j++) {
//           rhs1 = rhs1 + W[i][j]*e[j]/ve;
//         }
//         rhs1 = rhs1 + ww[i]*b[i]/ve;
//         lhs0 = 1.0/vb;
//         lhs1 = ww[i]/ve + 1.0/vb;
//         like0 = std::log(1.0/std::sqrt(lhs0)) + 0.5*(rhs0*rhs0)/lhs0 + std::log(pi[0]);
//         like1 = std::log(1.0/std::sqrt(lhs1)) + 0.5*(rhs1*rhs1)/lhs1 + std::log(pi[1]);
//         p0 = 1.0/(std::exp(like1 - like0) + 1.0);
//         d[i]=0;
//         std::uniform_real_distribution<double> runif(0.0, 1.0);
//         u = runif(gen);
//         if(u>p0) d[i]=1;
//         bn=0.0;
//         if(d[i]==1) {
//           //std::normal_distribution<double> rnorm(rhs1/lhs1, sqrt(ve/lhs1));
//           std::normal_distribution<double> rnorm(rhs1/lhs1, sqrt(1.0/lhs1));
//           bn = rnorm(gen);
//         }
//         diff = bn-b[i];
//         if(diff!=0.0) {
//           for (int j = 0; j < n; j++) {
//             e[j]=e[j] - W[i][j]*(diff);
//           }
//         }
//         b[i] = bn;
//       }
// 
//       // Sample pi for Bayes C
//       if(updatePi) {
//         std::vector<double> mc(2);
//         std::fill(mc.begin(), mc.end(), 0.0);
//         for (int i = 0; i<m ; i++) {
//           mc[d[i]] = mc[d[i]] + 1.0;
//         }
//         double pisum=0.0;
//         for (int j = 0; j<2 ; j++) {
//           std::gamma_distribution<double> rgamma(mc[j]+1.0,1.0);
//           double rg = rgamma(gen);
//           pi[j] = rg/m;
//           pisum = pisum + pi[j];
//         }
//         for (int j = 0; j<2 ; j++) {
//           pi[j] = pi[j]/pisum;
//           if(it>nburn) pim[j] = pim[j] + pi[j];
//         }
//         pis[it] = pi[1];
//       }  
//     }
// 
//     // Sample marker effects (BayesR)
//     if (method==5) {
//       for ( int isort = 0; isort < m; isort++) {
//         int i = order[isort];
//         if(!mask[i])   continue;
//         
//         // variance class likelihood version 1 
//         rhs = 0.0;
//         for ( int j = 0; j < n; j++) {
//           rhs = rhs + W[i][j]*e[j];
//         }
//         rhs = rhs + ww[i]*b[i];
//         v0 = ww[i]*ve;
//         logLc[0] = -0.5*std::log(v0) -0.5*((rhs*rhs)/v0) + std::log(pi[0]);
//         for (int j = 1; j<nc ; j++) {
//           vbc = vb * gamma[j];
//           v1 = ww[i]*ve + ww[i]*ww[i]*vbc;
//           logLc[j] = -0.5*std::log(v1) -0.5*((rhs*rhs)/v1) + std::log(pi[j]); 
//         }
//         
//         // variance class probability 
//         std::fill(probc.begin(), probc.end(), 0.0);
//         for (int j = 0; j<nc ; j++) {
//           logLcAdj = 0.0;
//           for (int k = 0; k<nc ; k++) {
//             logLcAdj += std::exp(logLc[k] - logLc[j]);
//           }
//           probc[j] = 1.0/logLcAdj;
//         }
//         // variance class probability (potential improved version)
//         //double m = *std::max_element(logLc.begin(), logLc.end());
//         //double denom = 0.0;
//         //for (int k = 0; k < nc; ++k) denom += std::exp(logLc[k] - m);
//         //for (int j = 0; j < nc; ++j) probc[j] = std::exp(logLc[j] - m) / denom;
//         
//         
//         // sample variance class indicator
//         std::uniform_real_distribution<double> runif(0.0, 1.0);
//         u = runif(gen);
//         d[i]=0;
//         cumprobc = 0.0;
//         for (int j = 0; j<nc ; j++) {
//           cumprobc += probc[j];
//           if(u < cumprobc){
//             d[i] = j;
//             break;
//           }
//         }
//         // sample marker effect
//         bn=0.0;
//         if(d[i]>0) {
//           vbc = vb * gamma[d[i]];
//           lhs =ww[i]+ve/vbc;
//           std::normal_distribution<double> rnorm(rhs/lhs, sqrt(ve/lhs));
//           //std::normal_distribution<double> rnorm(rhs/lhs, sqrt(1.0/lhs));
//           bn = rnorm(gen);
//         }
//         diff = bn-b[i];
//         if(diff!=0.0) {
//           for (int j = 0; j < n; j++) {
//             e[j]=e[j] - W[i][j]*(diff);
//           }
//         }
//         b[i] = bn;
//       }
//       // Sample pi for Bayes R
//       if(updatePi) {
//         std::vector<double> mc(nc);
//         std::fill(mc.begin(), mc.end(), 0.0);
//         for (int i = 0; i<m ; i++) {
//           mc[d[i]] = mc[d[i]] + 1.0;
//         }
//         double pisum=0.0;
//         for (int j = 0; j<nc ; j++) {
//           std::gamma_distribution<double> rgamma(mc[j]+1.0,1.0);
//           double rg = rgamma(gen);
//           pi[j] = rg/m;
//           pisum = pisum + pi[j];
//         }
//         for (int j = 0; j<nc ; j++) {
//           pi[j] = pi[j]/pisum;
//           if(it>nburn) pim[j] = pim[j] + pi[j];
//         }
//       }
//     }
//     
//     
//     // Sample marker variance
//     ssb = 0.0;
//     dfb = 0.0;
//     if (method<4) {
//       for ( int i = 0; i < m; i++) {
//         if(it>nburn && (it % nthin == 0)) bm[i] = bm[i] + b[i];
//         ssb = ssb + b[i]*b[i];
//         dfb = dfb + 1.0;
//         if(it>nburn && (it % nthin == 0)) dm[i] = dm[i] + 1.0;
//       }
//     }
//     if (method==4) {
//       for ( int i = 0; i < m; i++) {
//         if(it>nburn && (it % nthin == 0)) bm[i] = bm[i] + b[i];
//         if(d[i]==1)   {
//           ssb = ssb + b[i]*b[i];
//           dfb = dfb + 1.0;
//           if(it>nburn && (it % nthin == 0)) dm[i] = dm[i] + 1.0;
//         }
//       }
//     }
//     if (method==5) {
//       for ( int i = 0; i < m; i++) {
//         if(it>nburn && (it % nthin == 0)) bm[i] = bm[i] + b[i];
//         if(d[i]>0)   {
//           ssb = ssb + (b[i]*b[i])/gamma[d[i]];
//           dfb = dfb + 1.0;
//           if(it>nburn && (it % nthin == 0)) dm[i] = dm[i] + (double)d[i];
//         }
//       }
//     }
//     
//     // marker variance
//     if(updateB) {
//       std::chi_squared_distribution<double> rchisq(dfb+nub);
//       chi2 = rchisq(gen);
//       vb = (ssb + ssb_prior*nub)/chi2 ;
//       vbs[it] = (ssb + ssb_prior*nub)/chi2;
//     }
//     
//     // Sample residual variance
//     if(updateE) {
//       dfe = (double)n + nue;
//       sse = 0.0;
//       for ( int j = 0; j < n; j++) {
//         sse = sse + e[j]*e[j];
//       }
//       std::chi_squared_distribution<double> rchisq(dfe);
//       chi2 = rchisq(gen);
//       ve = (sse + sse_prior*nue)/chi2 ;
//       ves[it] = ve;
//     }
//     
//     // Sample marker specific tau for Bayes lasso
//     if (method==3) { 
//       // Hyper parameters for the inverse Gaussian distribution
//       //double lrate0 = 0.1;  // lambda0 hyperparameter rate
//       double lshape0 = 1.0; // lambda0 hyperparameter shape
//       
//       // Set a small tolerance value (epsilon)
//       const double epsilon = 1e-8;
//       
//       lambda_tau = lambda2;
//       lambda_tau = std::max(lambda_tau, epsilon);
//       
//       // Compute mutau for each element in vector b
//       double sum_tau = 0.0;
//       for (size_t i = 0; i < b.size(); ++i) {
//         // If b[i] is too small (close to zero), use epsilon instead
//         double abs_b = std::abs(b[i]) < epsilon ? epsilon : std::abs(b[i]);
//         // Calculate mu_tau with the tolerance for very small b[i]
//         // Legarra et al 2011
//         mu_tau = std::sqrt(lambda_tau) / abs_b;
//         // Park & Casella - Campos et al.
//         //mu_tau = std::sqrt(lambda_tau*ve) / abs_b;
//         mu_tau = std::max(mu_tau, epsilon);
//         double invtau = rinvgauss(mu_tau, lambda_tau, gen);
//         lambda[i] = invtau;
//         sum_tau += 1.0/invtau;
//       }
//       
//       // Calculate ratel2 and shl2
//       double shl2 = m + lshape0;    
//       //double ratel2 = sum_tau / 2.0 + lrate0;
//       //std::gamma_distribution<double> rgamma(shl2, 1.0/ratel2);
//       double scl2 = 2.0/sum_tau;
//       std::gamma_distribution<double> rgamma(shl2, scl2);
//       lambda2 = rgamma(gen);
//       lambda2 = std::max(1.0, std::min(lambda2, 1000000.0));
//     }
//     
//     // Update mu and adjust residuals
//     mu = std::accumulate(e.begin(), e.end(), 0.0)/n;
//     for ( int i = 0; i < n; i++) {
//       e[i] = e[i] - mu;
//     }                
//     mus[it] = mu;
//     
//     // Sample genetic variance
//     ssg = 0.0;
//     for ( int i = 0; i < n; i++) {
//       g[i] = y[i] - e[i] - mu;
//       ssg = ssg + g[i]*g[i];
//     }                
//     dfg = n + nug;
//     std::chi_squared_distribution<double> rchisq(dfg);
//     chi2 = rchisq(gen);
//     vg = (ssg + ssg_prior*nug)/chi2;
//     vgs[it] = vg;
// 
//     
//   }
// 
// 
//   // Summarize results
//   std::vector<std::vector<double>> result(12);
//   result[0].resize(m);
//   result[1].resize(m);
//   result[2].resize(nit+nburn);
//   result[3].resize(nit+nburn);
//   result[4].resize(nit+nburn);
//   result[5].resize(nit+nburn);
//   result[6].resize(nit+nburn);
//   result[7].resize(nc);
//   result[8].resize(n);
//   result[9].resize(m);
//   result[10].resize(m);
//   result[11].resize(3);
//   
//   for (int i=0; i < m; i++) {
//     //result[0][i] = bm[i]/nit;
//     //result[1][i] = dm[i]/nit;
//     result[0][i] = bm[i]/nsamples;
//     result[1][i] = dm[i]/nsamples;
//   }
//   for (int i=0; i < nit+nburn; i++) {
//     result[2][i] = mus[i];
//     result[3][i] = vbs[i];
//     result[4][i] = vgs[i];
//     result[5][i] = ves[i];
//     result[6][i] = pis[i];
//   }
//   for (int j=0; j < nc; j++) {
//     result[7][j] = pim[j]/nit;
//   }  
//   for (int i=0; i < n; i++) {
//     result[8][i] = y[i]- mu - e[i];
//   }
//   for (int i=0; i < m; i++) {
//     result[9][i] = b[i];
//     result[10][i] = d[i];
//   }
//   result[11][0] = vb;
//   result[11][1] = ve;
//   result[11][2] = pi[0];
//   return result;
// }

// [[Rcpp::export]]
std::vector<std::vector<double>> bayes(
    std::vector<double> y,
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
  
  // -----------------------------
  // Helpers
  // -----------------------------
  auto log_sum_exp = [](const std::vector<double>& x) -> double {
    double mx = *std::max_element(x.begin(), x.end());
    double s = 0.0;
    for (double v : x) s += std::exp(v - mx);
    return mx + std::log(s);
  };
  
  auto logistic_from_logs = [](double logp1, double logp0) -> double {
    double d = logp0 - logp1;
    if (d > 35.0) return 0.0;
    if (d < -35.0) return 1.0;
    return 1.0 / (1.0 + std::exp(d));
  };
  
  auto is_kept_sample = [&](int it) -> bool {
    if (it < nburn) return false;
    return ((it - nburn) % nthin) == 0;
  };
  
  // Inverse-Gaussian sampler IG(mu, lambda) using Michael–Schucany–Haas
  // Returns a draw from IG(mean=mu, shape=lambda)
  auto rinvgauss = [&](double mu, double shape, std::mt19937& gen,
                       std::normal_distribution<double>& norm01,
                       std::uniform_real_distribution<double>& runif) -> double {
                         const double eps = 1e-12;
                         mu = std::max(mu, eps);
                         shape = std::max(shape, eps);
                         
                         double z = norm01(gen);
                         double z2 = z * z;
                         double x = mu + (mu * mu * z2) / (2.0 * shape)
                           - (mu / (2.0 * shape)) *
                             std::sqrt(std::max(4.0 * mu * shape * z2 + mu * mu * z2 * z2, 0.0));
                         x = std::max(x, eps);
                         
                         double u = runif(gen);
                         double out = (u <= mu / (mu + x)) ? x : (mu * mu) / x;
                         return std::max(out, eps);
                       };
  
  // -----------------------------
  // Dimensions / RNG
  // -----------------------------
  const int n  = static_cast<int>(y.size());
  const int m  = static_cast<int>(W.size()); // W is m x n in your code
  const int nc = static_cast<int>(pi.size());
  const int total_it = nit + nburn;
  
  std::mt19937 gen(static_cast<uint32_t>(seed));
  std::normal_distribution<double> norm01(0.0, 1.0);
  std::uniform_real_distribution<double> runif(0.0, 1.0);
  
  // -----------------------------
  // State / storage
  // -----------------------------
  std::vector<int> d(m, 0);
  
  std::vector<double> e(n, 0.0);
  std::vector<double> g(n, 0.0);
  
  std::vector<double> wy(m, 0.0), ww(m, 0.0), x2(m, 0.0);
  std::vector<int> order(m, 0), mask(m, 1); // mask==1 means included/active (your original)
  
  std::vector<double> bm(m, 0.0), dm(m, 0.0);
  
  std::vector<double> ves(total_it, 0.0), vbs(total_it, 0.0), vgs(total_it, 0.0),
  pis(total_it, 0.0), mus(total_it, 0.0);
  
  std::vector<double> pim(nc, 0.0), probc(nc, 0.0), logLc(nc, 0.0);
  
  // ---- Posterior storage (post-burn only) ----
  int nstore = 0;
  for (int it = 0; it < total_it; ++it) {
    if (it >= nburn && ((it - nburn) % nthin == 0)) {
      nstore++;
    }
  }
  std::vector<std::vector<double>> bs(nstore, std::vector<double>(m));

  std::vector<double> gs(n, 0.0);
  std::vector<double> ys(n, 0.0);
  
  std::vector<double> vqs(total_it, 0.0);
  std::vector<double> pmse2(total_it, 0.0);
  std::vector<double> pmse3(total_it, 0.0);
  
  // log-CPO accumulator
  std::vector<double> sumpyinvt(n, 0.0);
  
  // BayesA full (global scale + local scales)
  double s2 = 1.0;
  double a0 = 3.0, b0 = 2.0;               // matches your R defaults
  std::vector<double> vba(m, vb);          // local scales (tau_i), init to vb like R's vba <- rep(vb,m)

  
  
    
  // -----------------------------
  // Initialize mu and residuals e
  // (R starts with mu=0 and e=y; your observed-level code centers y once.
  // We'll follow your current behavior: start from centered y, but then sample mu each iter as in R.)
  // -----------------------------
  double mu = std::accumulate(y.begin(), y.end(), 0.0) / std::max(n, 1);
  for (int i = 0; i < n; i++) e[i] = y[i] - mu;

  std::vector<double> varw(m, 0.0);
  
  // Precompute wy and ww based on centered residuals
  for (int i = 0; i < m; i++) {
    double sum_wye = 0.0;
    double sum_ww  = 0.0;
    for (int j = 0; j < n; j++) {
      sum_wye += W[i][j] * e[j];
      sum_ww  += W[i][j] * W[i][j];
    }
    wy[i] = sum_wye;
    ww[i] = std::max(sum_ww, 1e-300);
    varw[i] = ww[i] / std::max(n - 1.0, 1.0);
    double bhat = wy[i] / ww[i];
    x2[i] = bhat * bhat;
  }

  // Sort markers by descending x2
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(),
            [&](int i1, int i2) { return x2[i1] > x2[i2]; });
  
  // -----------------------------
  // BayesL init (R-identical intent)
  // -----------------------------
  double lshape0 = 1.0;
  double lambda2 = std::max(static_cast<double>(m), 10.0);
  {
    double init_lambda = std::sqrt(std::max(ve / vb, 1e-12));
    for (int i = 0; i < m; i++) lambda[i] = init_lambda;
  }
  
  std::vector<double> yhat_mean(n, 0.0);
  std::vector<double> yhat_M2(n, 0.0);
  int Kpost = 0;  // number of post-burn iterations actually accumulated
  
  
  // -----------------------------
  // Sampling counters
  // -----------------------------
  double nsamples = 0.0;
  int store_idx = 0;
  
  // -----------------------------
  // Gibbs sampler
  // -----------------------------
  for (int it = 0; it < total_it; it++) {
    
    if (is_kept_sample(it)) nsamples += 1.0;
    
    // ---- Sample intercept mu ----
    {
      // add back previous mu
      for (int j = 0; j < n; j++) e[j] += mu;
      
      double mean_e = std::accumulate(e.begin(), e.end(), 0.0) / std::max(n, 1);
      double sd_mu = std::sqrt(std::max(ve, 1e-300) / std::max(n, 1));
      
      std::normal_distribution<double> rnorm_mu(mean_e, sd_mu);
      mu = rnorm_mu(gen);
      mus[it] = mu;
      
      // remove new mu
      for (int j = 0; j < n; j++) e[j] -= mu;
    }
    
    // -------------------------
    // Marker updates
    // -------------------------
    
    // method==0: BLUP (deterministic; consistent with RR form)
    if (method == 0) {
      for (int isort = 0; isort < m; isort++) {
        int i = order[isort];
        if (!mask[i]) continue;
        
        double rhs = 0.0;
        for (int j = 0; j < n; j++) rhs += W[i][j] * e[j];
        rhs += ww[i] * b[i];
        
        double lhs = ww[i] + ve / vb;
        double bn  = rhs / lhs;
        
        double diff = bn - b[i];
        if (diff != 0.0) {
          for (int j = 0; j < n; j++) e[j] -= W[i][j] * diff;
        }
        b[i] = bn;
      }
    }
    
    // method==1: BayesN / BayesRR (matches R)
    if (method == 1) {
      for (int isort = 0; isort < m; isort++) {
        int i = order[isort];
        if (!mask[i]) continue;
        
        double rhs = 0.0;
        for (int j = 0; j < n; j++) rhs += W[i][j] * e[j];
        rhs += ww[i] * b[i];
        
        double lhs = ww[i] + ve / vb;
        
        std::normal_distribution<double> rnorm(rhs / lhs, std::sqrt(ve / lhs));
        double bn = rnorm(gen);
        
        double diff = bn - b[i];
        if (diff != 0.0) {
          for (int j = 0; j < n; j++) e[j] -= W[i][j] * diff;
        }
        b[i] = bn;
      }
    }
    
    // method==2: BayesA FULL (matches your R BayesA: local vba + global s2)
    if (method == 2) {
      const double df_post = nub + 1.0;
      
      for (int isort = 0; isort < m; isort++) {
        int i = order[isort];
        if (!mask[i]) continue;
        
        // σ_i^2 = s2 * vba[i]
        double vbi = std::max(s2 * vba[i], 1e-12);
        
        // sample b_i | vbi
        double rhs = 0.0;
        for (int j = 0; j < n; j++) rhs += W[i][j] * e[j];
        rhs += ww[i] * b[i];
        
        double lhs = ww[i] + ve / vbi;
        std::normal_distribution<double> rnorm(rhs / lhs, std::sqrt(ve / lhs));
        double bn = rnorm(gen);
        
        double diff = bn - b[i];
        if (diff != 0.0) {
          for (int j = 0; j < n; j++) e[j] -= W[i][j] * diff;
        }
        b[i] = bn;
        
        // sample vba[i] | b_i, s2  (Scale-inv-chi^2)
        double ssb = (b[i] * b[i]) / std::max(s2, 1e-12);
        std::chi_squared_distribution<double> rchisq(df_post);
        double chi2 = std::max(rchisq(gen), 1e-300);
        vba[i] = (ssb + nub * ssb_prior) / chi2;
      }
      
      // sample global s2 | b, vba  (Inv-Gamma via 1/Gamma)
      double shape = a0 + 0.5 * static_cast<double>(m);
      double rate  = b0;
      for (int i = 0; i < m; i++) {
        rate += 0.5 * (b[i] * b[i]) / std::max(vba[i], 1e-12);
      }
      std::gamma_distribution<double> rgamma(shape, 1.0 / std::max(rate, 1e-300));
      s2 = 1.0 / std::max(rgamma(gen), 1e-300);
    }
    
    // method==3: BayesL (matches your R BayesL exactly)
    if (method == 3) {
      const double eps = 1e-12;
      
      // 1) sample b | lambda
      for (int isort = 0; isort < m; isort++) {
        int i = order[isort];
        if (!mask[i]) continue;
        
        double rhs = 0.0;
        for (int j = 0; j < n; j++) rhs += W[i][j] * e[j];
        rhs += ww[i] * b[i];
        
        double tau_i = std::max(lambda[i], eps);            // local precision
        double lhs   = ww[i] / std::max(ve, eps) + tau_i;
        
        double mean_b = (rhs / std::max(ve, eps)) / lhs;
        double sd_b   = std::sqrt(1.0 / lhs);
        
        std::normal_distribution<double> rnorm(mean_b, sd_b);
        double bn = rnorm(gen);
        
        double diff = bn - b[i];
        if (diff != 0.0) {
          for (int j = 0; j < n; j++) e[j] -= W[i][j] * diff;
        }
        b[i] = bn;
      }
      
      // 2) sample local precisions lambda[i] via IG(mean=sqrt(lambda2)/|b|, shape=lambda2)
      lambda2 = std::max(lambda2, eps);
      double sum_tau = 0.0;
      for (int i = 0; i < m; i++) {
        if (!mask[i]) continue;
        double absb = std::max(std::abs(b[i]), eps);
        double mu_tau = std::sqrt(lambda2) / absb;
        double psi = rinvgauss(mu_tau, lambda2, gen, norm01, runif);
        lambda[i] = std::max(psi, eps);
        sum_tau += 1.0 / lambda[i];
      }
      
      // 3) sample global lambda2
      double shl2 = static_cast<double>(m) + lshape0;
      double scale = 2.0 / std::max(sum_tau, eps);
      std::gamma_distribution<double> rgamma(shl2, scale);
      lambda2 = rgamma(gen);
      lambda2 = std::min(std::max(lambda2, 1e-6), 1e6);
    }
    
    // method==4: BayesC (log Bayes factor form; matches your R BayesC)
    if (method == 4) {
      // normalize pi
      pi[0] = std::max(pi[0], 1e-300);
      pi[1] = std::max(pi[1], 1e-300);
      double ps = pi[0] + pi[1];
      pi[0] /= ps; pi[1] /= ps;
      
      for (int isort = 0; isort < m; isort++) {
        int i = order[isort];
        if (!mask[i]) continue;
        
        // r_i = e + W_i * b_i
        double rdot = 0.0;
        for (int j = 0; j < n; j++) rdot += W[i][j] * (e[j] + W[i][j] * b[i]);
        
        double denom = ve + ww[i] * vb;
        denom = std::max(denom, 1e-300);
        
        double logBF = 0.5 * std::log(std::max(ve, 1e-300) / denom)
          + 0.5 * (rdot * rdot) * vb / (std::max(ve, 1e-300) * denom);
        
        double logp1 = std::log(pi[1]) + logBF;
        double logp0 = std::log(pi[0]);
        double p1 = logistic_from_logs(logp1, logp0);
        
        d[i] = (runif(gen) < p1) ? 1 : 0;
        
        double bn = 0.0;
        if (d[i] == 1) {
          double lhs = ww[i] + ve / vb;
          double mean_b = rdot / lhs;
          double sd_b   = std::sqrt(ve / lhs);
          std::normal_distribution<double> rnorm(mean_b, sd_b);
          bn = rnorm(gen);
        }
        
        double diff = bn - b[i];
        if (diff != 0.0) {
          for (int j = 0; j < n; j++) e[j] -= W[i][j] * diff;
        }
        b[i] = bn;
      }
      
      // update pi ~ Dirichlet(1,1)
      if (updatePi) {
        double mc0 = 0.0, mc1 = 0.0;
        for (int i = 0; i < m; i++) (d[i] == 0) ? (mc0 += 1.0) : (mc1 += 1.0);
        
        std::gamma_distribution<double> rg0(mc0 + 1.0, 1.0);
        std::gamma_distribution<double> rg1(mc1 + 1.0, 1.0);
        double g0 = std::max(rg0(gen), 1e-300);
        double g1 = std::max(rg1(gen), 1e-300);
        double s = g0 + g1;
        
        pi[0] = g0 / s;
        pi[1] = g1 / s;
        
        pis[it] = pi[1];
        if (it >= nburn) { pim[0] += pi[0]; pim[1] += pi[1]; }
      } else {
        pis[it] = pi[1];
      }
    }
    
    // method==5: BayesR (log Bayes factor form; matches your R BayesR)
    if (method == 5) {
      // normalize pi
      double ps = 0.0;
      for (int k = 0; k < nc; k++) { pi[k] = std::max(pi[k], 1e-300); ps += pi[k]; }
      for (int k = 0; k < nc; k++) pi[k] /= ps;
      
      for (int isort = 0; isort < m; isort++) {
        int i = order[isort];
        if (!mask[i]) continue;
        
        // r_i = e + W_i*b_i, rdot = w_i^T r_i
        double rdot = 0.0;
        for (int j = 0; j < n; j++) rdot += W[i][j] * (e[j] + W[i][j] * b[i]);
        
        // spike
        logLc[0] = std::log(pi[0]);
        
        for (int k = 1; k < nc; k++) {
          double v_k = vb * gamma[k];
          if (v_k <= 0.0) { logLc[k] = -INFINITY; continue; }
          
          double denom = ve + ww[i] * v_k;
          denom = std::max(denom, 1e-300);
          
          double logBF = 0.5 * std::log(std::max(ve, 1e-300) / denom)
            + 0.5 * (rdot * rdot) * v_k / (std::max(ve, 1e-300) * denom);
          
          logLc[k] = std::log(pi[k]) + logBF;
        }
        
        double lse = log_sum_exp(logLc);
        for (int k = 0; k < nc; k++) probc[k] = std::exp(logLc[k] - lse);
        
        double u = runif(gen);
        double cum = 0.0;
        int cls = 0;
        for (int k = 0; k < nc; k++) {
          cum += probc[k];
          if (u <= cum) { cls = k; break; }
        }
        d[i] = cls;
        
        double bn = 0.0;
        if (cls > 0) {
          double v_k = vb * gamma[cls];
          double lhs = ww[i] + ve / v_k;
          double mean_b = rdot / lhs;
          double sd_b   = std::sqrt(ve / lhs);
          std::normal_distribution<double> rnorm(mean_b, sd_b);
          bn = rnorm(gen);
        }
        
        double diff = bn - b[i];
        if (diff != 0.0) {
          for (int j = 0; j < n; j++) e[j] -= W[i][j] * diff;
        }
        b[i] = bn;
      }
      
      // update pi ~ Dirichlet(1,...,1)
      if (updatePi) {
        std::vector<double> mc(nc, 0.0);
        for (int i = 0; i < m; i++) mc[d[i]] += 1.0;
        
        double s = 0.0;
        for (int k = 0; k < nc; k++) {
          std::gamma_distribution<double> rg(mc[k] + 1.0, 1.0);
          double gk = std::max(rg(gen), 1e-300);
          pi[k] = gk;
          s += gk;
        }
        for (int k = 0; k < nc; k++) {
          pi[k] /= s;
          if (it >= nburn) pim[k] += pi[k];
        }
        pis[it] = 1.0 - pi[0];
      } else {
        pis[it] = 1.0 - pi[0];
      }
    }
    
    // -------------------------
    // Accumulate bm/dm on thinned post-burn samples
    // -------------------------
    if (is_kept_sample(it)) {
      for (int i = 0; i < m; i++) {
        if (!mask[i]) continue;
        bm[i] += b[i];

        if (method < 4) dm[i] += 1.0;
        else if (method == 4) dm[i] += (d[i] == 1) ? 1.0 : 0.0;
        else if (method == 5) dm[i] += (d[i] > 0) ? 1.0 : 0.0;
      }
      bs[store_idx] = b;   // copy current marker effects
      store_idx++;
    }

    // -------------------------
    // Sample vb (matches your R "common vb update" logic)
    // -------------------------
    if (updateB) {
      double ssb = 0.0;
      double dfb = 0.0;
      
      if (method == 5) {
        for (int i = 0; i < m; i++) {
          if (!mask[i]) continue;
          if (d[i] > 0) {
            ssb += (b[i] * b[i]) / std::max(gamma[d[i]], 1e-300);
            dfb += 1.0;
          }
        }
      } else if (method == 4) {
        for (int i = 0; i < m; i++) {
          if (!mask[i]) continue;
          if (d[i] == 1) {
            ssb += b[i] * b[i];
            dfb += 1.0;
          }
        }
      } else {
        for (int i = 0; i < m; i++) {
          if (!mask[i]) continue;
          ssb += b[i] * b[i];
          dfb += 1.0;
        }
      }
      
      std::chi_squared_distribution<double> rchisq(dfb + nub);
      double chi2 = std::max(rchisq(gen), 1e-300);
      vb = (ssb + ssb_prior * nub) / chi2;
      vbs[it] = vb;
    } else {
      vbs[it] = vb;
    }
    
    // -------------------------
    // Sample ve (matches your R residual variance update)
    // -------------------------
    if (updateE) {
      double sse = 0.0;
      for (int j = 0; j < n; j++) sse += e[j] * e[j];
      
      double dfe = static_cast<double>(n) + nue;
      std::chi_squared_distribution<double> rchisq(dfe);
      double chi2 = std::max(rchisq(gen), 1e-300);
      ve = (sse + sse_prior * nue) / chi2;
      ves[it] = ve;
    } else {
      ves[it] = ve;
    }
    
    // -------------------------
    // Sample vg (should in principal be a derived parameter)
    // -------------------------
    double ssg = 0.0;
    for ( int i = 0; i < n; i++) {
      g[i] = y[i] - e[i] - mu;
      ssg = ssg + g[i]*g[i];
    }
    double dfg = static_cast<double>(n) + nug;
    std::chi_squared_distribution<double> rchisq(dfg);
    double chi2 = std::max(rchisq(gen), 1e-300);
    vg = (ssg + ssg_prior * nug) / chi2;
    vgs[it] = vg;

    // ---- Posterior diagnostics (post-burn only; NOT thinned) ----
    
    // g = y - e - mu
    double mean_g = 0.0;
    for (int j = 0; j < n; j++) {
      g[j] = y[j] - e[j] - mu;
      mean_g += g[j];
    }
    mean_g /= n;
    
    // center g
    for (int j = 0; j < n; j++) {
      g[j] -= mean_g;
    }
    
    // vqs[k] = sum(varw * b^2)
    double vq = 0.0;
    for (int i = 0; i < m; i++) {
      vq += varw[i] * b[i] * b[i];
    }
    vqs[it] = vq;
    
    
    // fitted values: yhat = y - e
    for (int j = 0; j < n; j++) {
      gs[j] = y[j] - e[j];
    }
    
    Kpost++;
    
    for (int j = 0; j < n; j++) {
      double x = gs[j];                 // yhat_j at this iteration
      double delta = x - yhat_mean[j];
      yhat_mean[j] += delta / Kpost;
      double delta2 = x - yhat_mean[j];
      yhat_M2[j] += delta * delta2;
    }
    
    // posterior predictive: yhat + N(0, ve)
    std::normal_distribution<double> rnorm(0.0, std::sqrt(ve));
    for (int j = 0; j < n; j++) {
      ys[j] = gs[j] + rnorm(gen);
    }
    
    
    
    // PMSE2 / PMSE3
    double mse2 = 0.0, mse3 = 0.0;
    for (int j = 0; j < n; j++) {
      double d2 = y[j] - gs[j];
      double d3 = y[j] - ys[j];
      mse2 += d2 * d2;
      mse3 += d3 * d3;
    }
    pmse2[it] = mse2 / static_cast<double>(n);
    pmse3[it] = mse3 / static_cast<double>(n);
    
    if (it >= nburn) {
      const double log_norm = std::log(2.0 * M_PI * ve);
      for (int j = 0; j < n; j++) {
        double diff = y[j] - gs[j];   // yhat
        double logdens = -0.5 * (log_norm + (diff * diff) / ve);
        sumpyinvt[j] += std::exp(-logdens);  // 1 / dnorm
      }
    }
    
    
  }

  // ---- Posterior summaries ----
  int K = nit;   // number of post-burn diagnostic samples
  
  // posterior mean residual variance
  double mve = 0.0;
  for (int k = 0; k < K; k++) {
    mve += ves[nburn + k];
  }
  mve /= static_cast<double>(K);
  
  // posterior mean variance of fitted values
  double mvxb = 0.0;
  if (Kpost > 1) {
    for (int j = 0; j < n; j++) {
      mvxb += yhat_M2[j] / (Kpost - 1.0);   // Var(yhat_j) across draws
    }
    mvxb /= static_cast<double>(n);         // mean across j
  } else {
    mvxb = 0.0;
  }
  
    
  // log-CPO
  double logcpo = 0.0;
  for (int j = 0; j < n; j++) {
    double denom = std::max(sumpyinvt[j], 1e-300);
    double phatyt_j = static_cast<double>(K) / denom;
    logcpo += std::log(phatyt_j);
  }
  
  // EpMSE
  double epmse1 = mve + mvxb;
  double epmse2 = mve + 2.0 * mvxb;
  double epmse3 = 2.0 * mve + 2.0 * mvxb;
  
    
  // -----------------------------
  // Summarize results (same layout as your current function)
  // -----------------------------
  if (nsamples <= 0.0) nsamples = 1.0;
  
  std::vector<std::vector<double>> result(13);
  result[0].resize(m);        // bm
  result[1].resize(m);        // dm
  result[2] = vqs;            // mus
  result[3] = vbs;            // vbs
  result[4] = vgs;            // vgs
  result[5] = ves;            // ves
  result[6] = pis;            // pis
  result[7].resize(nc);       // pim
  result[8].resize(n);        // fitted = y - mu - e (same as before)
  result[9].resize(m);        // b
  result[10].resize(m);       // d
  result[11].resize(3);       // param
  result[12].resize(4);       // diagnostics
  
  for (int i = 0; i < m; i++) {
    result[0][i] = bm[i] / nsamples;
    result[1][i] = dm[i] / nsamples;
  }
  
  // pim averaged over post-burn iterations (unthinned), matching accumulation (it>=nburn)
  double npi = static_cast<double>(std::max(0, total_it - nburn));
  if (npi <= 0.0) npi = 1.0;
  for (int k = 0; k < nc; k++) result[7][k] = pim[k] / npi;
  
  // fitted values
  for (int j = 0; j < n; j++) result[8][j] = y[j] - mu - e[j];
  
  // final b and d
  for (int i = 0; i < m; i++) {
    result[9][i]  = b[i];
    result[10][i] = static_cast<double>(d[i]);
  }
  
  result[11][0] = vb;
  result[11][1] = ve;
  result[11][2] = (nc > 0 ? pi[0] : 0.0);
  
  // diagnostics
  result[12][0] = logcpo;
  result[12][1] = epmse1;
  result[12][2] = epmse2;
  result[12][3] = epmse3;
  
  return result;
}


// [[Rcpp::export]]
std::vector<std::vector<double>> sbayes(
    double yy,
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
    int algo,   // kept for interface compatibility; not used (algo==2 skipped)
    int seed) {
  
  (void)algo; // silence unused warning
  
  // -----------------------------
  // Helpers
  // -----------------------------
  auto log_sum_exp = [](const std::vector<double>& x) -> double {
    double mx = *std::max_element(x.begin(), x.end());
    double s = 0.0;
    for (double v : x) s += std::exp(v - mx);
    return mx + std::log(s);
  };
  
  auto logistic_from_logs = [](double logp1, double logp0) -> double {
    // p1 = 1 / (1 + exp(logp0 - logp1))
    double d = logp0 - logp1;
    if (d > 35.0) return 0.0;
    if (d < -35.0) return 1.0;
    return 1.0 / (1.0 + std::exp(d));
  };
  
  auto update_r_from_diff = [&](std::vector<double>& r, int i, double diff) {
    // r[j] -= LD_ji * diff (LDvalues already scaled to match ww)
    for (size_t k = 0; k < LDindices[i].size(); k++) {
      r[LDindices[i][k]] -= LDvalues[i][k] * diff;
    }
  };
  
  // -----------------------------
  // Dimensions / RNG
  // -----------------------------
  const int m  = static_cast<int>(b.size());
  const int nc = static_cast<int>(pi.size());
  const int total_it = nit + nburn;
  
  std::mt19937 gen(static_cast<uint32_t>(seed));
  std::uniform_real_distribution<double> runif(0.0, 1.0);
  
  // -----------------------------
  // State / storage
  // -----------------------------
  std::vector<int> d(m, 0);
  std::vector<double> r(m, 0.0);      // summary-stat "residual": wy - (W'W)b
  std::vector<double> vba(m, 1.0);  // local scales τ_i
  std::vector<double> vbi(m, vb);   // σ_i^2 = s2 * vba[i]
  
  std::vector<double> vadj(m, 0.0);
  std::vector<double> vei(m, ve);
  
  std::vector<double> bm(m, 0.0), dm(m, 0.0);
  
  std::vector<double> ves(total_it, 0.0), vbs(total_it, 0.0), vgs(total_it, 0.0),
  pis(total_it, 0.0), vqs(total_it, 0.0);

  //std::vector<std::vector<double>> bs(nit+nburn, std::vector<double>(m, 0.0));  

  // ---- Posterior storage (post-burn only) ----
  int nstore = 0;
  for (int it = 0; it < total_it; ++it) {
    if (it >= nburn && ((it - nburn) % nthin == 0)) {
      nstore++;
    }
  }
  std::vector<std::vector<double>> bs(nstore, std::vector<double>(m));
  
    
  std::vector<double> pim(nc, 0.0);
  std::vector<double> x2(m, 0.0);
  std::vector<int> order(m, 0);
  
  std::vector<double> probc(nc, 0.0), logLc(nc, 0.0);
  
  // ---- BayesA global scale ----
  double s2 = 1.0;        // global variance scale
  double a0 = 3.0;        // matches your R defaults
  double b0 = 2.0;
  
  
  // -----------------------------
  // Initialize r, order, adjustments
  // -----------------------------
  for (int i = 0; i < m; i++) {
    r[i]  = wy[i]; // assumes initial b is zero in R interface (it is)
    // ranking statistic
    if (ww[i] > 0.0) {
      double bhat = wy[i] / ww[i];
      x2[i] = bhat * bhat;
    } else {
      x2[i] = 0.0;
    }
  }
  
  // adjustE: per-marker residual variance component via vg + ve
  for (int i = 0; i < m; i++) {
    vadj[i] = 0.0;
    if (adjustE) {
      // fraction of markers NOT in LD neighborhood (your original heuristic)
      vadj[i] = (static_cast<double>(m) - static_cast<double>(LDindices[i].size())) /
        static_cast<double>(m);
    }
    vei[i] = vadj[i] * vg + ve;
  }
  
  // order by descending x2
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(), [&](int i1, int i2) { return x2[i1] > x2[i2]; });
  
  // Scale LD values once (as you did): LDvalues *= sqrt(ww_i)*sqrt(ww_j)
  for (int i = 0; i < m; i++) {
    double si = std::sqrt(std::max(ww[i], 0.0));
    for (size_t k = 0; k < LDindices[i].size(); k++) {
      int j = LDindices[i][k];
      double sj = std::sqrt(std::max(ww[j], 0.0));
      LDvalues[i][k] *= (si * sj);
    }
  }
  
  // -----------------------------
  // BayesL initialization (R-identical intent)
  // -----------------------------
  
  // Matches R defaults:
  //   lambda  <- rep(sqrt(ve / vb), m)
  //   lambda2 <- max(m, 10)
  
  double shape0 = 1.0;                 // R: lshape0 <- 1.0
  double lambda2 = std::max(static_cast<double>(m), 10.0);
  
  // Local precisions
  double init_lambda = std::sqrt(std::max(ve / vb, 1e-12));
  for (int i = 0; i < m; i++) {
    lambda[i] = init_lambda;
  }
  
  // -----------------------------
  // Sampling counters
  // -----------------------------
  double nsamples = 0.0;
  
  auto is_kept_sample = [&](int it) -> bool {
    // burn-in: it in [0, nburn-1]
    if (it < nburn) return false;
    // thin in the post-burn segment: keep every nthin
    return ((it - nburn) % nthin) == 0;
  };
  
  // -----------------------------
  // Gibbs sampler
  // -----------------------------
  
  int store_idx = 0;
  
  
  for (int it = 0; it < total_it; it++) {
    
    if (is_kept_sample(it)) nsamples += 1.0;
    
    // -------------------------
    // Marker updates by method
    // -------------------------
    
    // method==0: BLUP (deterministic)
    if (method == 0) {
      for (int isort = 0; isort < m; isort++) {
        int i = order[isort];
        if (mask[i]) continue;
        
        double lhs = ww[i] + vei[i] / vb;
        double rhs = r[i] + ww[i] * b[i];
        double bn  = rhs / lhs;
        
        double diff = bn - b[i];
        if (diff != 0.0) update_r_from_diff(r, i, diff);
        b[i] = bn;
      }
    }
    
    // method==1: BayesN / BayesRR
    if (method == 1) {
      for (int isort = 0; isort < m; isort++) {
        int i = order[isort];
        if (mask[i]) continue;
        
        double lhs = ww[i] + vei[i] / vb;
        double rhs = r[i] + ww[i] * b[i];
        
        std::normal_distribution<double> rnorm(rhs / lhs, std::sqrt(vei[i] / lhs));
        double bn = rnorm(gen);
        
        double diff = bn - b[i];
        if (diff != 0.0) update_r_from_diff(r, i, diff);
        b[i] = bn;
      }
    }
    
    // method == 2 : FULL BayesA (local + global scale)
    if (method == 2) {
      
      const double df_post = nub + 1.0;
      
      // ---- 1. Update marker effects b_i | σ_i^2 ----
      for (int isort = 0; isort < m; isort++) {
        int i = order[isort];
        if (mask[i]) continue;
        
        // σ_i^2 = s2 * τ_i
        vbi[i] = std::max(s2 * vba[i], 1e-12);
        
        double lhs = ww[i] + vei[i] / vbi[i];
        double rhs = r[i] + ww[i] * b[i];
        
        std::normal_distribution<double> rnorm(rhs / lhs,
                                               std::sqrt(vei[i] / lhs));
        double bn = rnorm(gen);
        
        double diff = bn - b[i];
        if (diff != 0.0) update_r_from_diff(r, i, diff);
        b[i] = bn;
        
        // ---- 2. Update local scale τ_i | b_i, s2 ----
        double ssb = (b[i] * b[i]) / s2;
        std::chi_squared_distribution<double> rchisq(df_post);
        double chi2 = std::max(rchisq(gen), 1e-300);
        
        vba[i] = (ssb + nub * ssb_prior) / chi2;
      }
      
      // ---- 3. Update global scale s2 | {b_i, τ_i} ----
      double shape = a0 + 0.5 * static_cast<double>(m);
      double rate  = b0;
      
      for (int i = 0; i < m; i++) {
        rate += 0.5 * (b[i] * b[i]) / std::max(vba[i], 1e-12);
      }
      
      std::gamma_distribution<double> rgamma(shape,
                                             1.0 / std::max(rate, 1e-300));
      s2 = 1.0 / std::max(rgamma(gen), 1e-300);
    }
    // method==3: BayesL (Bayesian Lasso) -- matches your R formulation, adapted to summary stats
    if (method == 3) {
      
      const double eps = 1e-12;
      const double lshape0 = 1.0;  // matches your R: lshape0 <- 1.0
      
      // Guard lambda2
      lambda2 = std::max(lambda2, eps);
      
      // --- 1) Sample regression effects b_i | lambda[i], vei[i] ---
      for (int isort = 0; isort < m; isort++) {
        int i = order[isort];
        if (mask[i]) continue;
        
        // r_i_dot = w_i^T r_i, where r_i excludes marker i (in summary-stat bookkeeping)
        double r_i_dot = r[i] + ww[i] * b[i];
        
        // R: lhs <- ww[i]/ve + lambda[i]
        //    bn  <- rnorm(mean = (rhs/ve)/lhs, sd = sqrt(1/lhs))
        double tau_i = std::max(lambda[i], eps);              // local precision
        double ve_i  = std::max(vei[i], eps);
        
        double lhs = ww[i] / ve_i + tau_i;
        double mean_b = (r_i_dot / ve_i) / lhs;
        double sd_b   = std::sqrt(1.0 / lhs);
        
        std::normal_distribution<double> rnorm(mean_b, sd_b);
        double bn = rnorm(gen);
        
        double diff = bn - b[i];
        if (diff != 0.0) update_r_from_diff(r, i, diff);
        b[i] = bn;
      }
      
      // --- 2) Sample local precisions lambda[i] via Inv-Gaussian (R: rinvgauss(mean=mu_tau, shape=lambda2)) ---
      double sum_tau = 0.0;  // R accumulates sum(1/lambda[i]) into sum_tau
      
      // Note: In your existing code, you already have the Michael–Schucany–Haas IG sampler.
      // We'll reuse that sampler but with parameters matching R:
      //   mu_tau = sqrt(lambda2) / |b_i|
      //   shape  = lambda2
      std::normal_distribution<double> norm01(0.0, 1.0);
      
      for (int i = 0; i < m; i++) {
        if (mask[i]) continue;
        
        double absb  = std::max(std::abs(b[i]), eps);
        double mu_tau = std::sqrt(lambda2) / absb;            // R: mu_tau <- sqrt(lambda2) / abs(b[i])
        mu_tau = std::max(mu_tau, eps);
        
        double shape_ig = lambda2;                            // R: shape=lambda2
        shape_ig = std::max(shape_ig, eps);
        
        // Michael–Schucany–Haas sampler for IG(mu, lambda)
        double z  = norm01(gen);
        double z2 = z * z;
        
        double x = mu_tau
        + (mu_tau * mu_tau * z2) / (2.0 * shape_ig)
          - (mu_tau / (2.0 * shape_ig)) *
          std::sqrt(std::max(4.0 * mu_tau * shape_ig * z2 + mu_tau * mu_tau * z2 * z2, 0.0));
        
        x = std::max(x, eps);
        
        double u = runif(gen);
        double psi;
        if (u <= mu_tau / (mu_tau + x)) psi = x;
        else                            psi = (mu_tau * mu_tau) / x;
        
        // In R: lambda[i] <- psi_i  (lambda[i] is the *local precision*)
        lambda[i] = std::max(psi, eps);
        
        // R: sum_tau <- sum_tau + 1 / lambda[i]
        sum_tau += 1.0 / lambda[i];
      }
      
      // --- 3) Sample global lambda2 ---
      // R:
      // shl2    <- m + lshape0
      // lambda2 <- rgamma(shape=shl2, scale = 2 / sum_tau)
      // then clamp (optional)
      double shl2 = static_cast<double>(m) + lshape0;
      
      // Gamma(shape, scale) in R is Gamma(shape, scale). In C++ std::gamma_distribution uses shape, scale.
      double scale = 2.0 / std::max(sum_tau, eps);
      std::gamma_distribution<double> rgamma(shl2, scale);
      lambda2 = rgamma(gen);
      
      // same clamps you used in R
      lambda2 = std::min(std::max(lambda2, 1e-6), 1e6);
    }
    
    // method==3: BayesL (kept compatible with your current implementation)
    // if (method == 3) {
    //   const double df_post = 1.0 + nub;
    //   
    //   for (int isort = 0; isort < m; isort++) {
    //     int i = order[isort];
    //     if (mask[i]) continue;
    //     
    //     double vloc = std::max(vbi[i], 1e-12);
    //     double lhs  = ww[i] + vei[i] / vloc;
    //     double rhs  = r[i] + ww[i] * b[i];
    //     
    //     std::normal_distribution<double> rnorm(rhs / lhs, std::sqrt(vei[i] / lhs));
    //     double bn = rnorm(gen);
    //     
    //     double diff = bn - b[i];
    //     if (diff != 0.0) update_r_from_diff(r, i, diff);
    //     b[i] = bn;
    //     
    //     // Local IG update (as in your existing code)
    //     double eps = 1e-12;
    //     double absb = std::max(std::abs(b[i]), eps);
    //     double mu_tau = std::sqrt(std::max(vei[i], eps)) * lambda[i] / absb;
    //     double lambda_tau = std::max(lambda2, eps);
    //     
    //     std::normal_distribution<double> norm01(0.0, 1.0);
    //     double z = norm01(gen);
    //     double z2 = z * z;
    //     
    //     double x_tau = mu_tau
    //     + 0.5 * mu_tau * mu_tau * z2 / lambda_tau
    //     - 0.5 * (mu_tau / lambda_tau) *
    //     std::sqrt(std::max(4.0 * mu_tau * lambda_tau * z2 + mu_tau * mu_tau * z2 * z2, 0.0));
    //     
    //     double u = runif(gen);
    //     double tau = mu_tau * mu_tau / std::max(x_tau, eps);
    //     if (u <= mu_tau / (mu_tau + x_tau)) tau = x_tau;
    //     
    //     double vbin = std::max(tau, eps);
    //     vbi[i] = vbin;
    //   }
    //   
    //   // update lambda2 hyperparameter (as in your existing code)
    //   double ss = 0.0;
    //   double df = 0.0;
    //   for (int i = 0; i < m; i++) {
    //     ss += vbi[i] * vbi[i];
    //     df += 1.0;
    //   }
    //   double shape = shape0 + df;
    //   double rate  = rate0 + ss / 2.0;
    //   std::gamma_distribution<double> rgamma(shape, 1.0 / std::max(rate, 1e-300));
    //   lambda2 = rgamma(gen);
    //   lambda2 = std::max(lambda2, 1e-12);
    //   for (int i = 0; i < m; i++) lambda[i] = std::sqrt(lambda2);
    // }
    
    // method==4: BayesC (log Bayes factor; matches your R form with ve -> vei[i])
    if (method == 4) {
      // normalize pi
      pi[0] = std::max(pi[0], 1e-300);
      pi[1] = std::max(pi[1], 1e-300);
      double ps = pi[0] + pi[1];
      pi[0] /= ps; pi[1] /= ps;
      
      for (int isort = 0; isort < m; isort++) {
        int i = order[isort];
        if (mask[i]) continue;
        
        // r_i_dot = w_i^T r_i  (r_i = residual excluding marker i)
        double r_i_dot = r[i] + ww[i] * b[i];
        
        double denom = vei[i] + ww[i] * vb;
        denom = std::max(denom, 1e-300);
        
        double logBF = 0.5 * std::log(std::max(vei[i], 1e-300) / denom)
          + 0.5 * (r_i_dot * r_i_dot) * vb / (std::max(vei[i], 1e-300) * denom);
        
        double logp1 = std::log(pi[1]) + logBF;
        double logp0 = std::log(pi[0]);
        double p1 = logistic_from_logs(logp1, logp0);
        
        int di = (runif(gen) < p1) ? 1 : 0;
        d[i] = di;
        
        double bn = 0.0;
        if (di == 1) {
          double lhs = ww[i] + vei[i] / vb;
          double mean_b = r_i_dot / lhs;
          double sd_b   = std::sqrt(vei[i] / lhs);
          std::normal_distribution<double> rnorm(mean_b, sd_b);
          bn = rnorm(gen);
        }
        
        double diff = bn - b[i];
        if (diff != 0.0) update_r_from_diff(r, i, diff);
        b[i] = bn;
      }
      
      // pi update (Dirichlet(1,1))
      if (updatePi) {
        double mc0 = 0.0, mc1 = 0.0;
        for (int i = 0; i < m; i++) (d[i] == 0) ? (mc0 += 1.0) : (mc1 += 1.0);
        
        std::gamma_distribution<double> rg0(mc0 + 1.0, 1.0);
        std::gamma_distribution<double> rg1(mc1 + 1.0, 1.0);
        double g0 = std::max(rg0(gen), 1e-300);
        double g1 = std::max(rg1(gen), 1e-300);
        double s = g0 + g1;
        
        pi[0] = g0 / s;
        pi[1] = g1 / s;
        
        pis[it] = pi[1];
        if (it >= nburn) { pim[0] += pi[0]; pim[1] += pi[1]; }
      }
    }
    
    // method==5: BayesR (log Bayes factor; matches your R form with ve -> vei[i])
    if (method == 5) {
      // normalize pi
      double ps = 0.0;
      for (int k = 0; k < nc; k++) {
        pi[k] = std::max(pi[k], 1e-300);
        ps += pi[k];
      }
      for (int k = 0; k < nc; k++) pi[k] /= ps;
      
      for (int isort = 0; isort < m; isort++) {
        int i = order[isort];
        if (mask[i]) continue;
        
        double r_i_dot = r[i] + ww[i] * b[i];
        
        // spike class
        logLc[0] = std::log(pi[0]);
        
        // slab classes
        for (int k = 1; k < nc; k++) {
          double v_k = vb * gamma[k];
          if (v_k <= 0.0) {
            logLc[k] = -INFINITY;
            continue;
          }
          double denom = vei[i] + ww[i] * v_k;
          denom = std::max(denom, 1e-300);
          
          double logBF = 0.5 * std::log(std::max(vei[i], 1e-300) / denom)
            + 0.5 * (r_i_dot * r_i_dot) * v_k / (std::max(vei[i], 1e-300) * denom);
          
          logLc[k] = std::log(pi[k]) + logBF;
        }
        
        // softmax -> probc
        double lse = log_sum_exp(logLc);
        for (int k = 0; k < nc; k++) probc[k] = std::exp(logLc[k] - lse);
        
        // sample class
        double u = runif(gen);
        double cum = 0.0;
        int cls = 0;
        for (int k = 0; k < nc; k++) {
          cum += probc[k];
          if (u <= cum) { cls = k; break; }
        }
        d[i] = cls;
        
        // sample effect
        double bn = 0.0;
        if (cls > 0) {
          double v_k = vb * gamma[cls];
          double lhs = ww[i] + vei[i] / v_k;
          double mean_b = r_i_dot / lhs;
          double sd_b   = std::sqrt(vei[i] / lhs);
          std::normal_distribution<double> rnorm(mean_b, sd_b);
          bn = rnorm(gen);
        }
        
        double diff = bn - b[i];
        if (diff != 0.0) update_r_from_diff(r, i, diff);
        b[i] = bn;
      }
      
      // pi update (Dirichlet(1,...,1))
      if (updatePi) {
        std::vector<double> mc(nc, 0.0);
        for (int i = 0; i < m; i++) mc[d[i]] += 1.0;
        
        double s = 0.0;
        for (int k = 0; k < nc; k++) {
          std::gamma_distribution<double> rg(mc[k] + 1.0, 1.0);
          double g = std::max(rg(gen), 1e-300);
          pi[k] = g;
          s += g;
        }
        for (int k = 0; k < nc; k++) {
          pi[k] /= s;
          if (it >= nburn) pim[k] += pi[k];
        }
        pis[it] = 1.0 - pi[0];
      }
    }
    
    // -------------------------
    // Accumulate posterior means (bm, dm) on kept samples
    // -------------------------
    if (is_kept_sample(it)) {
      for (int i = 0; i < m; i++) {
        if (mask[i]) continue;
        bm[i] += b[i];
        
        if (method < 4) {
          dm[i] += 1.0;
        } else if (method == 4) {
          dm[i] += (d[i] == 1) ? 1.0 : 0.0;
        } else if (method == 5) {
          dm[i] += (d[i] > 0) ? 1.0 : 0.0;
        }
      }
      bs[store_idx] = b;   // copy current marker effects
      store_idx++;
    }
    
    

    
    
    // -------------------------
    // Update vb (marker variance)
    // -------------------------
    if (updateB) {
      double ssb = 0.0;
      double ssbr = 0.0;
      double dfb = 0.0;
      
      if (method < 4) {
        for (int i = 0; i < m; i++) {
          if (mask[i]) continue;
          ssb += b[i] * b[i];
          dfb += 1.0;
        }
        vqs[it] = ssb;
      } else if (method == 4) {
        for (int i = 0; i < m; i++) {
          if (mask[i]) continue;
          if (d[i] == 1) {
            ssb += b[i] * b[i];
            dfb += 1.0;
          }
        }
        vqs[it] = ssb;
      } else if (method == 5) {
        for (int i = 0; i < m; i++) {
          if (mask[i]) continue;
          if (d[i] > 0) {
            ssb += (b[i] * b[i]) / std::max(gamma[d[i]], 1e-300);
            ssbr += (b[i]* b[i]);
            dfb += 1.0;
          }
        }
        vqs[it] = ssbr;
      }
      
      std::chi_squared_distribution<double> rchisq(dfb + nub);
      double chi2 = std::max(rchisq(gen), 1e-300);
      vb = (ssb + ssb_prior * nub) / chi2;
      
      vbs[it]  = vb;
    } else {
      vbs[it] = vb;
      // ssbs[it] left as 0 unless you want it always recorded; keep current behavior
    }
    
    // -------------------------
    // Update ve (residual variance) using your summary-stat SSE formula
    // -------------------------
    if (updateE) {
      double sse = 0.0;
      for (int i = 0; i < m; i++) {
        if (mask[i]) continue;
        sse += b[i] * (r[i] + wy[i]);
      }
      double dfe = static_cast<double>(n) + nue;
      sse = yy - sse;
      
      std::chi_squared_distribution<double> rchisq(dfe);
      double chi2 = std::max(rchisq(gen), 1e-300);
      ve = (sse + sse_prior * nue) / chi2;
      
      // refresh per-marker vei if adjustE
      for (int i = 0; i < m; i++) {
        vei[i] = vadj[i] * vg + ve;
      }
      
      ves[it] = ve;
    } else {
      ves[it] = ve;
    }
    
    // -------------------------
    // Update vg (genetic variance) using your summary-stat estimator
    // -------------------------
    {
      double ssg = 0.0;
      for (int i = 0; i < m; i++) {
        if (mask[i]) continue;
        ssg += b[i] * (wy[i] - r[i]);
      }
      double dfg = static_cast<double>(n);
      double vg_new = ssg / std::max(dfg, 1.0);
      
      vgs[it] = vg_new;
      if (updateG) vg = vg_new;
      
      if (adjustE) {
        for (int i = 0; i < m; i++) {
          vei[i] = vadj[i] * vg + ve;
        }
      }
    }
    
    // If pi not updated in this iteration, still store something sensible for pis[it]
    if (!updatePi) {
      if (method == 4 && nc >= 2) pis[it] = pi[1];
      if (method == 5 && nc >= 1) pis[it] = 1.0 - pi[0];
    }
  }
  
  // -----------------------------
  // Finalize summaries / outputs
  // -----------------------------
  if (nsamples <= 0.0) nsamples = 1.0;
  
  std::vector<std::vector<double>> result(12);
  result[0].assign(m, 0.0);          // bm
  result[1].assign(m, 0.0);          // dm
  result[2] = vqs;                  // ssbs
  result[3] = vbs;                   // vbs
  result[4] = vgs;                   // vgs
  result[5] = ves;                   // ves
  result[6] = pis;                   // pis
  result[7].assign(nc, 0.0);         // pim
  result[8] = r;                     // r (final)
  result[9] = b;                     // b (final)
  result[10].assign(3, 0.0);         // param
  result[11].resize(nstore*m);
  
  for (int i = 0; i < m; i++) {
    result[0][i] = bm[i] / nsamples;
    result[1][i] = dm[i] / nsamples;
  }
  
  // pim averaged over post-burn iterations (unthinned), matching accumulation condition (it >= nburn)
  double npi = static_cast<double>(std::max(0, total_it - nburn));
  if (npi <= 0.0) npi = 1.0;
  for (int k = 0; k < nc; k++) {
    result[7][k] = pim[k] / npi;
  }
  
  result[10][0] = vb;
  result[10][1] = ve;
  result[10][2] = (nc > 0 ? pi[0] : 0.0);

  for ( int it = 0; it < nstore; it++) {
    for ( int i = 0; i < m; i++) {
      result[11][it*m + i] = bs[it][i];
    }
  }
  
  return result;
}




// // [[Rcpp::export]]
// std::vector<std::vector<double>>  sbayes( double yy,
//                                           std::vector<double>& wy,
//                                           std::vector<double>& ww, 
//                                           std::vector<std::vector<double>>& LDvalues, 
//                                           std::vector<std::vector<int>>& LDindices, 
//                                           std::vector<double> b, 
//                                           std::vector<double> lambda,
//                                           std::vector<bool> mask, 
//                                           std::vector<double> pi, 
//                                           std::vector<double> gamma, 
//                                           double vg, 
//                                           double vb, 
//                                           double ve,
//                                           double ssb_prior,
//                                           double ssg_prior,
//                                           double sse_prior,
//                                           double nub,
//                                           double nug,
//                                           double nue,
//                                           bool updateB,
//                                           bool updateG,
//                                           bool updateE,
//                                           bool updatePi,
//                                           bool adjustE,
//                                           int n, 
//                                           int nit,
//                                           int nburn,
//                                           int nthin,
//                                           int method,
//                                           int algo,
//                                           int seed) {
//   
//   
//     
//   // Define local variables
//   int m = b.size();
//   int nc = pi.size();
//   double nsamples=0.0;
//   
//   
//   double rhs, lhs, bn, bj, diff;
//   double rhs1, lhs1, like0, like1, p0, v0, v1;
//   //double rhs0, lhs0;
//   double ssb, sse, ssg, dfb, dfe, dfg, chi2;
//   double x_tau, tau, lambda_tau, mu_tau, z, z2, u, vbin;
//   double shape, shape0, rate, rate0, lambda2;
//   
//   std::vector<double> vbscale(nc), pic(nc), pim(nc), probc(nc), logLc(nc);
//   double cumprobc, vbc, logLcAdj;
//   
//   
//   std::vector<int> d(m);
// 
//   std::vector<double> r(m),vei(m);
//   
//   std::vector<double> dm(m),bm(m);
//   std::vector<double> ves(nit+nburn),vbs(nit+nburn),vgs(nit+nburn),pis(nit+nburn),ssbs(nit+nburn);
//   
//   std::vector<double> x2(m),vadj(m),vbi(m);
//   std::vector<int> order(m);
//   
//   
//   // Initialize variables
//   for ( int i = 0; i < m; i++) {
//     vbi[i]=vb;
//     r[i] = wy[i];
//     x2[i] = (wy[i]/ww[i])*(wy[i]/ww[i]);
//   }
//   
//   std::fill(bm.begin(), bm.end(), 0.0);
//   std::fill(dm.begin(), dm.end(), 0.0);
//   std::fill(vbs.begin(), vbs.end(), 0.0);
//   std::fill(vgs.begin(), vgs.end(), 0.0);
//   std::fill(ves.begin(), ves.end(), 0.0);
//   std::fill(pis.begin(), pis.end(), 0.0);
//   std::fill(ssbs.begin(), ssbs.end(), 0.0);
//   std::fill(pim.begin(), pim.end(), 0.0);
//   
//   // adjust sparseld
//   for ( int i = 0; i < m; i++) {
//     vadj[i] = 0.0;
//     if(adjustE) {
//       vadj[i] = ((double)m-(double)LDindices[i].size())/(double)m;
//     }  
//     vei[i] = vadj[i]*vg + ve;
//   }
//   
//   // initialize BayesL parameters
//   //if (method==3) {
//   dfb = (nub - 2.0)/nub;
//   lambda2 = 2.0*(1.0 - dfb)/(dfb)*n;
//   shape0 = 1.1;
//   rate0 = (shape0 - 1.0) / lambda2;
//   for ( int i = 0; i < m; i++) {
//     lambda[i] = sqrt(lambda2); 
//   }
//   //}
//   for ( int i = 0; i < nc; i++) {
//     vbscale[i]=gamma[i];
//     pic[i]=pi[i];
//   }
// 
//   // Establish order of markers as they are entered into the model
//   std::iota(order.begin(), order.end(), 0);
//   std::sort(  std::begin(order), 
//               std::end(order),
//               [&](int i1, int i2) { return x2[i1] > x2[i2]; } );
// 
//   // Adjust LD values
//   for ( int i = 0; i < m; i++) {
//     for (size_t j = 0; j < LDindices[i].size(); j++) {
//       LDvalues[i][j] = LDvalues[i][j]*std::sqrt(ww[i])*std::sqrt(ww[LDindices[i][j]]);
//     }
//   }
//   
//   // // Wy - W'Wb
//   // for ( int i = 0; i < m; i++) {
//   //   if (b[i]!= 0.0) {
//   //     diff = b[i]*ww[i];
//   //     for (size_t j = 0; j < LDindices[i].size(); j++) {
//   //       r[LDindices[i][j]]=r[LDindices[i][j]] - LDvalues[i][j]*diff;
//   //     }
//   //   }
//   // }
//   
// 
//   // Start Gibbs sampler
//   std::random_device rd;
//   std::mt19937 gen(seed);
//   
//   for ( int it = 0; it < nit+nburn; it++) {
//   
//   if ( (it > nburn) && (it % nthin == 0) ) {
//     nsamples = nsamples + 1.0;
//   }
//   
//   
//     // Compute marker effects (BLUP)
//     if (method==0) {
//       for ( int isort = 0; isort < m; isort++) {
//         int i = order[isort];
//         if(mask[i])   continue;
//         lhs = ww[i] + vei[i]/vb;
//         rhs = r[i] + ww[i]*b[i];
//         bn = rhs/lhs;
//         diff = (bn-b[i]);
//         for (size_t j = 0; j < LDindices[i].size(); j++) {
//           r[LDindices[i][j]] += -LDvalues[i][j]*diff;
//         }
//         b[i] = bn;
//       }
//     }
//     
//     // Compute marker effects (BayesN or BayesRR)
//     if (method==1) {
//       for ( int isort = 0; isort < m; isort++) {
//         int i = order[isort];
//         if(mask[i])   continue;
//         lhs = ww[i] + vei[i]/vb;
//         rhs = r[i] + ww[i]*b[i];
//         std::normal_distribution<double> rnorm(rhs/lhs, sqrt(vei[i]/lhs));
//         bn = rnorm(gen);
//         diff = (bn-b[i]);
//         for (size_t j = 0; j < LDindices[i].size(); j++) {
//           r[LDindices[i][j]] += -LDvalues[i][j]*diff;
//         }
//         b[i] = bn;
//       }
//     }
// 
//     // Compute marker effects (BayesA)
//     if (method==2) {
//       dfb = 1.0 + nub;
//       for ( int isort = 0; isort < m; isort++) {
//         int i = order[isort];
//         if(mask[i])   continue;
//         ssb = b[i]*b[i];
//         std::chi_squared_distribution<double> rchisq(dfb);
//         chi2 = rchisq(gen);
//         vbi[i] = (ssb + ssb_prior*nub)/chi2 ;
//         lhs = ww[i] + vei[i]/vbi[i];
//         rhs = r[i] + ww[i]*b[i];
//         std::normal_distribution<double> rnorm(rhs/lhs, sqrt(vei[i]/lhs));
//         bn = rnorm(gen);
//         diff = (bn-b[i]);
//         for (size_t j = 0; j < LDindices[i].size(); j++) {
//           r[LDindices[i][j]] += -LDvalues[i][j]*diff;
//         }
//         b[i] = bn;
//       }
//     }
// 
//     // Compute marker effects (BayesL)
//     if (method==3) {
//       dfb = 1.0 + nub;
//       for ( int isort = 0; isort < m; isort++) {
//         int i = order[isort];
//         if(mask[i])   continue;
//         lhs = ww[i] + vei[i]/vbi[i];
//         rhs = r[i] + ww[i]*b[i];
//         std::normal_distribution<double> rnorm(rhs/lhs, sqrt(vei[i]/lhs));
//         bn = rnorm(gen);
//         diff = (bn-b[i]);
//         for (size_t j = 0; j < LDindices[i].size(); j++) {
//           r[LDindices[i][j]] += -LDvalues[i][j]*diff;
//         }
//         b[i] = bn;
//         mu_tau=sqrt(vei[i])*lambda[i]/std::abs(b[i]);
//         lambda_tau=lambda2;  
//         std::normal_distribution<double> norm(0.0, 1.0);
//         z = norm(gen);
//         z2=z*z;
//         x_tau=mu_tau+0.5*mu_tau*mu_tau*z2/lambda_tau - 0.5*(mu_tau/lambda_tau)*sqrt(4*mu_tau*lambda_tau*z2+mu_tau*mu_tau*z2*z2);
//         std::uniform_real_distribution<double> runif(0.0, 1.0);
//         u = runif(gen);
//         tau = mu_tau*mu_tau/x_tau;
//         if(u <= mu_tau/(mu_tau+x_tau)) tau=x_tau;
//         //vbin = 1.0/tau;
//         vbin = tau;
//         if(vbin > 0)   vbi[i] = vbin;
//         
//       }
//       // update hyperparameters
//       ssb = 0.0;
//       dfb = 0.0;
//       for ( int i = 0; i < m; i++) {
//         ssb = ssb + vbi[i]*vbi[i];
//         dfb = dfb + 1.0;
//       }
//       shape = shape0 + dfb;
//       rate = rate0 + ssb/ 2.0;
//       std::gamma_distribution<double> rgamma(shape, 1.0/rate);
//       lambda2 = rgamma(gen);
//       for ( int i = 0; i < m; i++) {
//         lambda[i] = sqrt(lambda2);
//       }
//     }
//     
//     // lambda_tau = 2.0*2.0*0.5*0.5/vb;
//     // ssb = b[i]*b[i];
//     // mu_tau = std::sqrt(lambda_tau/ssb);
//     // std::normal_distribution<double> rnorm(0.0, 1.0);
//     // z = rnorm(gen);
//     // z2=z*z;
//     // xtau=mu_tau+0.5*mu_tau*mu_tau*z2/lambda_tau - 0.5*(mu_tau/lambda_tau)*sqrt(4*mu_tau*lambda_tau*z2+mu_tau*mu_tau*z2*z2);
//     // std::uniform_real_distribution<double> runif(0.0, 1.0);
//     // u = runif(gen);
//     // tau = mu_tau*mu_tau/xtau;
//     // if(u <= mu_tau/(mu_tau+xtau)) tau=xtau;
//     // lambda[i] = ve/tau;
//     
// 
//     // Sample marker effects (BayesC)
//     if (method==4) {
//       for ( int isort = 0; isort < m; isort++) {
//         int i = order[isort];
//         if(mask[i])   continue;
//         // version 1
//         //rhs0 = 0.0;
//         //rhs1 = r[i] + ww[i]*b[i];
//         //lhs0 = 1.0/vb;
//         //lhs1 = ww[i]/vei[i] + 1.0/vb;
//         //like0 = std::log(1.0/std::sqrt(lhs0)) + 0.5*(rhs0*rhs0)/lhs0 + std::log(pi[0]); 
//         //like1 = std::log(1.0/std::sqrt(lhs1)) + 0.5*(rhs1*rhs1)/lhs1 + std::log(pi[1]); 
//         //p0 = 1.0/(std::exp(like1 - like0) + 1.0);
//         // version 2
//         rhs = r[i] + ww[i]*b[i];
//         v0 = ww[i]*vei[i];
//         v1 = ww[i]*vei[i] + ww[i]*ww[i]*vb;
//         like0 = -0.5*std::log(v0) -0.5*((rhs*rhs)/v0) + std::log(pi[0]);
//         like1 = -0.5*std::log(v1) -0.5*((rhs*rhs)/v1) + std::log(pi[1]);
//         p0 = 1.0/(std::exp(like1 - like0) + 1.0);
//         d[i]=0;
//         std::uniform_real_distribution<double> runif(0.0, 1.0);
//         u = runif(gen);
//         if(u>p0) d[i]=1;
//         bn=0.0;
//         if(d[i]==1) {
//           rhs1 = r[i] + ww[i]*b[i];
//           lhs1 = ww[i] + vei[i]/vb;
//           std::normal_distribution<double> rnorm(rhs1/lhs1, sqrt(vei[i]/lhs1));
//           bn = rnorm(gen);
//           if(algo==2 && it<nburn) {
//             bn = (1.0-p0)*bn;
//           }
//         } 
//         diff = (bn-b[i]);
//         if(diff!=0.0) {
//           for (size_t j = 0; j < LDindices[i].size(); j++) {
//             r[LDindices[i][j]] += -LDvalues[i][j]*diff;
//           }
//         }
//         b[i] = bn;
//       }
//       
//       // Sample pi for Bayes C
//       if(updatePi) {
//         std::vector<double> mc(2);
//         std::fill(mc.begin(), mc.end(), 0.0);
//         for (int i = 0; i<m ; i++) {
//           mc[d[i]] = mc[d[i]] + 1.0;
//         }
//         double pisum=0.0;
//         for (int j = 0; j<2 ; j++) {
//           std::gamma_distribution<double> rgamma(mc[j]+1.0,1.0);
//           double rg = rgamma(gen);
//           pi[j] = rg/m;
//           pisum = pisum + pi[j];
//         }
//         for (int j = 0; j<2 ; j++) {
//           pi[j] = pi[j]/pisum;
//           if(it>nburn) pim[j] = pim[j] + pi[j];
//         }
//         pis[it] = pi[1];
//         // dfb=0.0;
//         // for (int i = 0; i<m ; i++) {
//         //   if(d[i]==1)   {
//         //     dfb = dfb + 1.0;
//         //   }
//         // }
//         // double count = dfb + 1.0;
//         // std::gamma_distribution<double> rgamma(count,1.0);
//         // double rg = rgamma(gen);
//         // double pisum=0.0;
//         // pi[1] = rg/(double)m;
//         // pi[0] = 1.0 - pi[1];
//         // pisum = pi[0] + pi[1];
//         // pi[0] = pi[0]/pisum;
//         // pi[1] = pi[1]/pisum;
//         // pis[it] = pi[1];
//         // if(it>nburn) pim[0] = pim[0] + pi[0];
//         // if(it>nburn) pim[1] = pim[1] + pi[1];
//       }
//       
//     }
//     
//     // Sample marker effects (BayesR)
//     if (method==5) {
//       
//       for ( int isort = 0; isort < m; isort++) {
//         int i = order[isort];
//         if(mask[i])   continue;
//         // variance class likelihood 
//         rhs = r[i] + ww[i]*b[i];
//         v0 = ww[i]*vei[i];
//         logLc[0] = -0.5*std::log(v0) -0.5*((rhs*rhs)/v0) + std::log(pi[0]);
//         for (int j = 1; j<nc ; j++) {
//           vbc = vb * gamma[j];
//           v1 = ww[i]*vei[i] + ww[i]*ww[i]*vbc;
//           logLc[j] = -0.5*std::log(v1) -0.5*((rhs*rhs)/v1) + std::log(pi[j]); 
//         }
//         // variance class probability 
//         std::fill(probc.begin(), probc.end(), 0.0);
//         for (int j = 0; j<nc ; j++) {
//           logLcAdj = 0.0;
//           for (int k = 0; k<nc ; k++) {
//             logLcAdj += std::exp(logLc[k] - logLc[j]);
//           }
//           probc[j] = 1.0/logLcAdj;
//         }
//         // sample variance class indicator
//         std::uniform_real_distribution<double> runif(0.0, 1.0);
//         u = runif(gen);
//         d[i]=0;
//         cumprobc = 0.0;
//         for (int j = 0; j<nc ; j++) {
//           cumprobc += probc[j];
//           if(u < cumprobc){
//             d[i] = j;
//             break;
//           }
//         }
//         // sample marker effect
//         bn=0.0;
//         if(d[i]>0) {
//           vbc = vb * gamma[d[i]];
//           lhs =ww[i]+vei[i]/vbc;
//           std::normal_distribution<double> rnorm(rhs/lhs, sqrt(vei[i]/lhs));
//           bn = rnorm(gen);
//           if(algo==2 && it<nburn) {
//           bn=0.0;
//           for (size_t j = 1; j < gamma.size(); j++) {
//             vbc = vb * gamma[j];
//             lhs =ww[i]+vei[i]/vbc;
//             std::normal_distribution<double> rnorm(rhs/lhs, sqrt(vei[i]/lhs));
//             bj = rnorm(gen);
//             bn += probc[j]*bj;
//           }
//           }
//         }
//         diff = (bn-b[i]);
//         if(diff!=0.0) {
//           for (size_t j = 0; j < LDindices[i].size(); j++) {
//             r[LDindices[i][j]] += -LDvalues[i][j]*diff;
//           }
//         }
//         b[i] = bn;
//       }
//       // Sample pi for Bayes R
//       if(updatePi) {
//         std::vector<double> mc(nc);
//         std::fill(mc.begin(), mc.end(), 0.0);
//         for (int i = 0; i<m ; i++) {
//           mc[d[i]] = mc[d[i]] + 1.0;
//         }
//         double pisum=0.0;
//         for (int j = 0; j<nc ; j++) {
//           std::gamma_distribution<double> rgamma(mc[j]+1.0,1.0);
//           double rg = rgamma(gen);
//           pi[j] = rg/m;
//           pisum = pisum + pi[j];
//         }
//         for (int j = 0; j<nc ; j++) {
//           pi[j] = pi[j]/pisum;
//           if(it>nburn) pim[j] = pim[j] + pi[j];
//         }
//         pis[it] = 1.0 - pi[0];
//       }
//     }
//     
//     
//     // Sample marker variance
//     ssb = 0.0;
//     dfb = 0.0;
//     if (method<4) {
//       for ( int i = 0; i < m; i++) {
//         if(it>nburn && (it % nthin == 0)) bm[i] = bm[i] + b[i];
//         ssb = ssb + b[i]*b[i];
//         dfb = dfb + 1.0;
//         if(it>nburn && (it % nthin == 0)) dm[i] = dm[i] + 1.0;
//       }
//     }
//     if (method==4) {
//       for ( int i = 0; i < m; i++) {
//         if(it>nburn && (it % nthin == 0)) bm[i] = bm[i] + b[i];
//         if(d[i]==1)   {
//           ssb = ssb + b[i]*b[i];
//           dfb = dfb + 1.0;
//           if(it>nburn && (it % nthin == 0)) dm[i] = dm[i] + 1.0;
//         }
//       }
//     }
//     if (method==5) {
//       for ( int i = 0; i < m; i++) {
//         if(it>nburn && (it % nthin == 0)) bm[i] = bm[i] + b[i];
//         if(d[i]>0)   {
//           ssb = ssb + (b[i]*b[i])/gamma[d[i]];
//           dfb = dfb + 1.0;
//           if(it>nburn && (it % nthin == 0)) dm[i] = dm[i] + 1.0;
//           //if(it>nburn) dm[i] = dm[i] + (double)d[i];
//         }
//       }
//     }
//     
//     // marker variance
//     if(updateB) {
//       std::chi_squared_distribution<double> rchisq(dfb+nub);
//       chi2 = rchisq(gen);
//       vb = (ssb + ssb_prior*nub)/chi2 ;
//       vbs[it] = vb;
//       ssbs[it] = ssb;
//     }
// 
//     // Sample residual variance
//     if(updateE) {
//       sse = 0.0;
//       for ( int i = 0; i < m; i++) {
//         if(mask[i])   continue;
//         sse = sse + b[i] * (r[i] + wy[i]);
//       }
//       dfe = (double)n + nue;
//       sse = yy - sse;
//       std::chi_squared_distribution<double> rchisq(dfe);
//       chi2 = rchisq(gen);
//       ve = (sse + sse_prior*nue)/chi2 ;
//       for ( int i = 0; i < m; i++) {
//         vei[i] = vadj[i]*vg + ve;
//       }
//       ves[it] = ve;
//     }
// 
//     
//     // Sample genetic variance
//     ssg = 0.0;
//     for ( int i = 0; i < m; i++) {
//       if(mask[i])   continue;
//       ssg = ssg + b[i] * (wy[i] -  r[i]);
//     }
//     //dfg = (double)n + nug;
//     //std::chi_squared_distribution<double> rchisq(dfg);
//     //chi2 = rchisq(gen);
//     //vg = (ssg + ssg_prior*nug)/chi2;
//     dfg = (double)n;
//     vgs[it] = ssg/dfg;
//     if(updateG) {
//       vg = ssg/dfg;
//     }
//     if(adjustE) {
//       for ( int i = 0; i < m; i++) {
//         vei[i] = vadj[i]*vg + ve;
//       }
//     }
//   }
//   
//   // Summarize results
//   std::vector<std::vector<double>> result(11);
//   result[0].resize(m);
//   result[1].resize(m);
//   result[2].resize(nit+nburn);
//   result[3].resize(nit+nburn);
//   result[4].resize(nit+nburn);
//   result[5].resize(nit+nburn);
//   result[6].resize(nit+nburn);
//   result[7].resize(nc);
//   result[8].resize(m);
//   result[9].resize(m);
//   result[10].resize(3);
// 
//   for (int i=0; i < m; i++) {
//     //result[0][i] = bm[i]/nit;
//     //result[1][i] = dm[i]/nit;
//     result[0][i] = bm[i]/nsamples;
//     result[1][i] = dm[i]/nsamples;
//   }
//   for (int i=0; i < nit+nburn; i++) {
//     //result[2][i] = mus[i];
//     result[2][i] = ssbs[i];
//     result[3][i] = vbs[i];
//     result[4][i] = vgs[i];
//     result[5][i] = ves[i];
//     result[6][i] = pis[i];
//   }
//   for (int j=0; j < nc; j++) {
//     result[7][j] = pim[j]/nit;
//   }  
//   
//   for (int i=0; i < m; i++) {
//     result[8][i] = r[i];
//     result[9][i] = b[i];
//     //result[9][i] = d[i];
//   }
//   result[10][0] = vb;
//   result[10][1] = ve;
//   result[10][2] = pi[0];
// 
//   
//   return result;
// }
// 

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
                                              double mh_p,
                                              double mh_r2,
                                              bool updateB,
                                              bool updateG,
                                              bool updateE,
                                              bool updatePi,
                                              bool updateMH,
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
  
  // Create a deep copy of LDvalues
  std::vector<std::vector<double>> XXvalues = LDvalues;
  
  
  // Adjust LD values
  for ( int i = 0; i < m; i++) {
    for (size_t j = 0; j < LDindices[i].size(); j++) {
      XXvalues[i][j] = XXvalues[i][j]*std::sqrt(ww[i])*std::sqrt(ww[LDindices[i][j]]);
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
            r[LDindices[i][j]] += -XXvalues[i][j]*diff;
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
            r[LDindices[i][j]] += -XXvalues[i][j]*diff;
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
    
    // MH step
    if (updateMH) {
      for (int i = 0; i < m; i++) {
        
        if (mask[i]) continue;     // Skip masked SNPs
        if (d[i] == 0) continue;   // Only consider SNPs currently included
        
        int nc = LDindices[i].size();  
        if (nc == 0) continue;
        
        // Sample whether to attempt a jump
        std::uniform_real_distribution<double> runif(0.0, 1.0);
        double mh_u = runif(gen);
        if (mh_u > mh_p) continue;         // Skip with 1 - mh_p probability
        
        // Copy LD values and indices
        std::vector<double> pld = LDvalues[i];
        std::vector<int> indices = LDindices[i];
        
        // Compute squared LD, apply r² threshold, and normalize
        double sumsq = 0.0;
        for (int j = 0; j < nc; j++) {
          pld[j] *= pld[j];  // r²
          if (pld[j] < mh_r2) pld[j] = 0.0;
          sumsq += pld[j];
        }
        
        if (sumsq == 0.0) continue;  // No strong LD neighbors, skip
        
        // Normalize LD-based probabilities
        for (int j = 0; j < nc; j++) {
          pld[j] /= sumsq;
        }
        
        // Build cumulative distribution
        std::vector<double> cumpld(nc);
        cumpld[0] = pld[0];
        for (int j = 1; j < nc; j++) {
          cumpld[j] = cumpld[j - 1] + pld[j];
        }
        
        // Sample target SNP j from LD cluster
        double u = runif(gen);
        int tj = 0;
        while (tj < nc && u > cumpld[tj]) {
          tj++;
        }
        if (tj >= nc) continue;  // Safety
        int j = indices[tj];
        
        if (j == i) continue;  // No-op (stay), skip
        
        // Effect transfer
        double diff = -b[i];  // Remove i's effect
        for (size_t k = 0; k < LDindices[i].size(); k++) {
          r[LDindices[i][k]] += XXvalues[i][k] * diff;
        }
        
        // Transferring scaled effect
        double bn = b[i] * LDvalues[i][tj]; 
        
        // Update residuals for j
        diff = bn - b[j];
        for (size_t k = 0; k < LDindices[j].size(); k++) {
          r[LDindices[j][k]] += XXvalues[j][k] * diff;
        }
        
        // Update effects and inclusion indicators
        b[j] = bn;
        b[i] = 0.0;
        d[j] = d[i];
        d[i] = 0;
      }
    }
    // End MH step
    
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
                                                    double mh_p,
                                                    double mh_r2,
                                                    bool updateB,
                                                    bool updateG,
                                                    bool updateE,
                                                    bool updatePi,
                                                    bool updateMH,
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
          //rhs +=r[LDindices[i][j]]*LDvalues[i][LDindices[i][j]];
          rhs +=r[LDindices[i][j]]*LDvalues[i][j];
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
            //r[LDindices[i][j]] += -LDvalues[i][LDindices[i][j]]*diff;
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
    
    // MH step
    if (updateMH) {
      for (int i = 0; i < m; i++) {
        
        if (mask[i]) continue;     // Skip masked SNPs
        if (d[i] == 0) continue;   // Only consider SNPs currently included
        
        int nc = LDindices[i].size();  
        if (nc == 0) continue;
        
        // Sample whether to attempt a jump
        std::uniform_real_distribution<double> runif(0.0, 1.0);
        double mh_u = runif(gen);
        if (mh_u > mh_p) continue;         // Skip with 1 - mh_p probability
        
        // Copy LD values and indices
        std::vector<double> pld = LDvalues[i];
        std::vector<int> indices = LDindices[i];
        
        // Compute squared LD, apply r² threshold, and normalize
        double sumsq = 0.0;
        for (int j = 0; j < nc; j++) {
          pld[j] *= pld[j];  // r²
          if (pld[j] < mh_r2) pld[j] = 0.0;
          sumsq += pld[j];
        }
        
        if (sumsq == 0.0) continue;  // No strong LD neighbors, skip
        
        // Normalize LD-based probabilities
        for (int j = 0; j < nc; j++) {
          pld[j] /= sumsq;
        }
        
        // Build cumulative distribution
        std::vector<double> cumpld(nc);
        cumpld[0] = pld[0];
        for (int j = 1; j < nc; j++) {
          cumpld[j] = cumpld[j - 1] + pld[j];
        }
        
        // Sample target SNP j from LD cluster
        double u = runif(gen);
        int tj = 0;
        while (tj < nc && u > cumpld[tj]) {
          tj++;
        }
        if (tj >= nc) continue;  // Safety
        int j = indices[tj];
        
        if (j == i) continue;  // No-op (stay), skip
        
        // Effect transfer
        double diff = -b[i];  // Remove i's effect
        for (size_t k = 0; k < LDindices[i].size(); k++) {
          r[LDindices[i][k]] += LDvalues[i][k] * diff;
        }
        
        // Transferring scaled effect
        double bn = b[i] * LDvalues[i][tj]; 
        
        // Update residuals for j
        diff = bn - b[j];
        for (size_t k = 0; k < LDindices[j].size(); k++) {
          r[LDindices[j][k]] += LDvalues[j][k] * diff;
        }
        
        // Update effects and inclusion indicators
        b[j] = bn;
        b[i] = 0.0;
        d[j] = d[i];
        d[i] = 0;
      }
    }
    // End MH step
    
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
