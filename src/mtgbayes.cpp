// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;
using namespace arma;


#include <vector>
#include <numeric>
#include <random>

arma::mat mmult(arma::mat A, arma::mat B) {
  return A * B;
}

// [[Rcpp::export]]
arma::mat mvrnormARMA(arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(1, ncols);
  return (Y * arma::chol(sigma));
}


// [[Rcpp::export]]
arma::mat rwishart(unsigned int df, const arma::mat& S){
  // Dimension of returned wishart
  unsigned int m = S.n_rows;
  
  // Z composition:
  // sqrt chisqs on diagonal
  // random normals below diagonal
  // misc above diagonal
  arma::mat Z(m,m);
  
  // Fill the diagonal
  for(unsigned int i = 0; i < m; i++){
    Z(i,i) = sqrt(R::rchisq(df-i));
  }
  
  // Fill the lower matrix with random guesses
  for(unsigned int j = 0; j < m; j++){  
    for(unsigned int i = j+1; i < m; i++){    
      Z(i,j) = R::rnorm(0,1);
    }
  }
  
  // Lower triangle * chol decomp
  arma::mat C = arma::trimatl(Z).t() * arma::chol(S);
  
  // Return random wishart
  return C.t()*C;
}

// [[Rcpp::export]]
arma::mat riwishart(unsigned int df, const arma::mat& S){
  return rwishart(df,S.i()).i();
}

#include <iostream>
#include <random>
#include <vector>
#include <algorithm>
#include <numeric>


// Sample pi
void samplePi(std::vector<double>& cmodel, 
              std::vector<double>& pi, 
              std::mt19937& gen) {
  // Iterate over elements of cmodel
  for (size_t k = 0; k < cmodel.size(); k++) {
    // Create a gamma distribution with shape parameter cmodel[k] and scale parameter 1.0
    std::gamma_distribution<double> rgamma(cmodel[k], 1.0);
    // Generate a random gamma-distributed value using the provided random number generator gen
    double rg = rgamma(gen);
    // Store the generated value in the pi vector
    pi[k] = rg;
  }
  // Calculate the sum of all values in the pi vector
  double psum = std::accumulate(pi.begin(), pi.end(), 0.0);
  // Normalize the pi vector by dividing each element by the sum
  for (size_t k = 0; k < cmodel.size(); k++) {
    pi[k] = pi[k] / psum;
  }
  // Reset all elements in the cmodel vector to 1.0
  std::fill(cmodel.begin(), cmodel.end(), 1.0);
}

// Sample marker effect covariance matrix B
void sampleB(int nt,
             int m,
             int nub,
             arma::mat& B,
             const std::vector<std::vector<int>>& d,
             const std::vector<std::vector<double>>& b,
             const std::vector<std::vector<double>>& ssb_prior,
             std::mt19937& gen) {
  // Initialize matrices and vectors for intermediate calculations
  arma::mat Sb(nt, nt, arma::fill::zeros);   // Sum of squared b values
  arma::mat corb(nt, nt, arma::fill::zeros); // Correlation matrix
  arma::mat dfB(nt, nt, arma::fill::zeros);  // Degrees of freedom matrix
  
  // Calculate Sb, dfB, and stdv based on input data
  for (int t1 = 0; t1 < nt; t1++) {
    double ssb = 0.0; // Initialize sum of squared b values for t1
    double dfb = 0.0; // Initialize degrees of freedom for t1
    
    // Calculate ssb and dfb for t1 based on d and b matrices
    for (int i = 0; i < m; i++) {
      if (d[t1][i] == 1) {
        ssb += b[t1][i] * b[t1][i];
        dfb += 1.0;
      }
    }
    
    // Store ssb and dfb in the corresponding matrices
    Sb(t1, t1) = ssb;
    dfB(t1, t1) = dfb;
    
    // Calculate correlations and off-diagonal elements of Sb and dfB
    for (int t2 = t1; t2 < nt; t2++) {
      ssb = 0.0; // Reset ssb for t2
      dfb = 0.0; // Reset dfb for t2
      
      if (t1 != t2) {
        for (int i = 0; i < m; i++) {
          if (d[t1][i] == 1 && d[t2][i] == 1) {
            ssb += b[t1][i] * b[t2][i];
            dfb += 1.0;
          }
        }
        
        // Populate the lower and upper triangles of Sb and dfB
        dfB(t1, t2) = dfb;
        dfB(t2, t1) = dfb;
        Sb(t1, t2) = ssb;
        Sb(t2, t1) = ssb;
      }
    }
  }
  
  // Calculate standard deviations and correlations (corb)
  std::vector<double> stdv(nt);
  for (int t = 0; t < nt; t++) {
    stdv[t] = sqrt(Sb(t,t));
  }
  
  for (int t1 = 0; t1 < nt; t1++) {
    for (int t2 = 0; t2 < nt; t2++) {
      if (t1 == t2) {
        corb(t1, t2) = 1.0; // Diagonal elements of corb are set to 1.0
      }
      
      if (t1 != t2 && stdv[t1] != 0.0 &&  stdv[t2] != 0.0) {
        //corb(t1, t2) = 0.0;
        corb(t1, t2) = Sb(t1, t2) / (stdv[t1] * stdv[t2]); // Calculate correlations
        corb(t2, t1) = corb(t1, t2);                      // Correlation matrix is symmetric
      }
    }
  }
  
  // Add a small constant to the diagonal elements of corb to ensure PD matrix
  for (int t1 = 0; t1 < nt; t1++) {
    corb(t1, t1) += 0.0001;
  }
  
  // Calculate diagonal elements of B using chi-squared distribution
  for (int t = 0; t < nt; t++) {
    std::chi_squared_distribution<double> rchisq(dfB(t, t) + nub);
    double chi2 = rchisq(gen);
    B(t, t) = (Sb(t, t) + nub * ssb_prior[t][t]) / chi2;
  }
  
  // Calculate off-diagonal elements of B based on correlations and standard deviations
  for (int t1 = 0; t1 < nt; t1++) {
    for (int t2 = 0; t2 < nt; t2++) {
      if (t1 != t2) {
        B(t1, t2) = corb(t1, t2) * std::sqrt(B(t1, t1)) * std::sqrt(B(t2, t2));
      }
    }
  }
  
  // Check if B is symmetric; if not, make it symmetric
  bool issym = B.is_symmetric();
  if (!issym) {
    B = 0.5 * (B + B.t());
  }
  
  // Calculate the inverse of B and check for success
  arma::mat Bi(nt, nt, arma::fill::zeros);
  //bool success;
  //success = arma::inv_sympd(Bi, B, arma::inv_opts::allow_approx);
  
  // If the inverse calculation fails, print an error message
  //if (!success) {
  //  std::cerr << "Error: Condition is false." << std::endl;
  //}
}

// Sample residual covariance matrix E
void sampleE(
    int nt,
    int m,
    int nue,
    arma::mat& E, 
    const std::vector<std::vector<double>>& b,
    const std::vector<std::vector<double>>& wy,
    const std::vector<std::vector<double>>& r,
    const std::vector<std::vector<double>>& sse_prior,
    const std::vector<double>& yy,
    const std::vector<int>& n,
    std::mt19937& gen
) {
  // Initialize the Se matrix to store the sum of squared residuals
  arma::mat Se(nt, nt, arma::fill::zeros);
  
  // Calculate the sum of squared residuals for each time point
  for (int t1 = 0; t1 < nt; t1++) {
    double sse = 0.0; // Initialize the sum of squared residuals for t1
    
    // Calculate the sum of squared residuals using b, wy, and r vectors
    for (int i = 0; i < m; i++) {
      sse += b[t1][i] * (r[t1][i] + wy[t1][i]);
    }
    
    // Calculate the adjusted sum of squared residuals
    sse = yy[t1] - sse;
    
    // Store the adjusted sum of squared residuals in the Se matrix
    Se(t1, t1) = sse + nue * sse_prior[t1][t1];
  }
  
  // Generate chi-squared random variables and update the E matrix
  for (int t = 0; t < nt; t++) {
    std::chi_squared_distribution<double> rchisq(n[t] + nue);
    double chi2 = rchisq(gen);
    E(t, t) = Se(t, t) / chi2;
  }
  
  // Calculate the inverse of the E matrix
  arma::mat Ei = arma::inv(E);
}


// Compute genetic covariance matrix G
void computeG(
    int nt,
    int m,
    const std::vector<std::vector<double>>& b,
    const std::vector<std::vector<double>>& wy,
    const std::vector<std::vector<double>>& r,
    const std::vector<int>& n,
    arma::mat& G
) {
  for (int t1 = 0; t1 < nt; t1++) {
    for (int t2 = t1; t2 < nt; t2++) {
      double ssg = 0.0; // Initialize the sum of squared G values
      
      if (t1 == t2) {
        // Calculate the sum of squared G values for the diagonal elements
        for (int i = 0; i < m; i++) {
          ssg += b[t1][i] * (wy[t1][i] - r[t1][i]);
        }
        
        // Store the G value in the matrix, normalized by the square root of n[t1] and n[t2]
        G(t1, t2) = ssg / (std::sqrt(static_cast<double>(n[t1])) * std::sqrt(static_cast<double>(n[t2])));
      }
      
      if (t1 != t2) {
        ssg = 0.0; // Reset the sum of squared G values
        
        // Calculate the sum of squared G values for off-diagonal elements
        for (int i = 0; i < m; i++) {
          ssg += b[t1][i] * (wy[t1][i] - r[t1][i]);
          ssg += b[t2][i] * (wy[t2][i] - r[t2][i]);
        }
        
        // Adjust ssg and store the G value in the matrix, normalized by the square root of n[t1] and n[t2]
        ssg = ssg / 2.0;
        G(t1, t2) = ssg / (std::sqrt(static_cast<double>(n[t1])) * std::sqrt(static_cast<double>(n[t2])));
        
        // Since G is symmetric, set the corresponding off-diagonal element
        G(t2, t1) = G(t1, t2);
      }
    }
  }
}


// Sample marker effects using multiple trait BayesC
void sampleBetaCMt(int i, 
                 int nt, 
                 int nmodels, 
                 const std::vector<std::vector<int>>& models, 
                 std::vector<double>& cmodel,
                 const std::vector<double>& pi,
                 const arma::mat& Ei, 
                 const arma::mat& Bi,
                 const std::vector<std::vector<double>>& ww,
                 std::vector<std::vector<double>>& r, 
                 std::vector<std::vector<double>>& b,
                 std::vector<std::vector<int>>& d,
                 const std::vector<std::vector<int>>& XXindices, 
                 const std::vector<std::vector<std::vector<double>>>& XXvalues,
                 std::mt19937& gen) {
  std::vector<double> rhs(nt);
  double detC;
  std::vector<double> loglik(nmodels);
  std::vector<double> pmodel(nmodels);
  std::vector<double> logliksum(nmodels, 0.0);
  double u;
  int mselect = 0;
  double cumprobc = 0.0;
  
  // Compute rhs
  for (int t = 0; t < nt; t++) {
    rhs[t] = Ei(t, t) * r[t][i] + Ei(t, t) * ww[t][i] * b[t][i];
  }
  
  for (int k = 0; k < nmodels; k++) {
    arma::mat C = Bi;
    for (int t1 = 0; t1 < nt; t1++) {
      if (models[k][t1] == 1) {
        C(t1, t1) = C(t1, t1) + ww[t1][i] * Ei(t1, t1);
      }
    }
    arma::mat Ci = arma::inv(C);
    detC = arma::det(Ci);
    loglik[k] = 0.5 * std::log(detC) + std::log(pi[k]);
    for (int t1 = 0; t1 < nt; t1++) {
      for (int t2 = t1; t2 < nt; t2++) {
        if (models[k][t1] == 1 && models[k][t2] == 1) {
          loglik[k] = loglik[k] + 0.5 * rhs[t1] * rhs[t2] * Ci(t1, t2);
          if (t1 != t2) {
            loglik[k] = loglik[k] + 0.5 * rhs[t2] * rhs[t1] * Ci(t2, t1);
          }
        }
      }
    }
  }
  
  std::fill(pmodel.begin(), pmodel.end(), 0.0);
  for (int k = 0; k < nmodels; k++) {
    logliksum[k] = 0.0;
    for (int l = 0; l < nmodels; l++) {
      logliksum[k] += std::exp(loglik[l] - loglik[k]);
    }
    pmodel[k] = 1.0 / logliksum[k];
  }
  
  // Sample variance class indicator
  std::uniform_real_distribution<double> runif(0.0, 1.0);
  u = runif(gen);
  mselect=0;
  cumprobc = 0.0;
  
  for (int k = 0; k<nmodels ; k++) {
    cumprobc += pmodel[k];
    if(u < cumprobc){
      mselect = k;
      break;
    }
  }
  cmodel[mselect] = cmodel[mselect] + 1.0; 
  for ( int t = 0; t < nt; t++) { 
    d[t][i] = models[mselect][t];
  }
  
  // Sample marker effect conditional on variance class indicator
  arma::mat C = Bi;
  for (int t1 = 0; t1 < nt; t1++) {
    if (models[mselect][t1] == 1) {
      C(t1, t1) = C(t1, t1) + ww[t1][i] * Ei(t1, t1);
    }
  }
  arma::mat Ci = arma::inv(C);
  arma::mat mub = mvrnormARMA(Ci);
  for (int t1 = 0; t1 < nt; t1++) {
    if (models[mselect][t1] == 0) {
      mub(0, t1) = 0.0;
    }
    if (models[mselect][t1] == 1) {
      mub(0, t1) = mub(0, t1) + Ci(t1, t1) * rhs[t1];
      for (int t2 = 0; t2 < nt; t2++) {
        if (t1 != t2 && models[mselect][t2] == 1) {
          mub(0, t1) = mub(0, t1) + Ci(t1, t2) * rhs[t2];
        }
      }
    }
  }
  
  // Adjust residuals based on sample marker effects
  double diff;
  for (int t = 0; t < nt; t++) {
    diff = (mub(0, t) - b[t][i]);
    if (diff != 0.0) {
      for (size_t j = 0; j < XXindices[i].size(); j++) {
        r[t][XXindices[i][j]] = r[t][XXindices[i][j]] - XXvalues[t][i][j] * diff;
      }
    }
    b[t][i] = mub(0, t);
  }
}


// Sample pi across trait
void samplePiMt(int nt,
               std::vector<double>& pimarker,
               const std::vector<int>& dmarker,
               std::mt19937& gen) {
  
  std::vector<double> mc(2);
  
  // Sample marker specific pi0
  std::fill(mc.begin(), mc.end(), 0.0);
  for (size_t i = 0; i<dmarker.size() ; i++) {
    mc[dmarker[i]] = mc[dmarker[i]] + 1.0;
  }
  double pisum=0.0;
  for (int j = 0; j<2 ; j++) {
    std::gamma_distribution<double> rgamma(mc[j]+1.0,1.0);
    double rg = rgamma(gen);
    pimarker[j] = rg/dmarker.size();
    pisum = pisum + pimarker[j];
  }
  for (int j = 0; j<2 ; j++) {
    pimarker[j] = pimarker[j]/pisum;
  }
  
}

// Sample pi within trait BayesC
void samplePiC(int nt,
               std::vector<std::vector<double>>& pitrait, 
               const std::vector<std::vector<int>>& d,
               std::mt19937& gen) {

  std::vector<double> mc(2);
  
  // Sample trait specific pi
  for (int t = 0; t < nt; t++) {
    std::fill(mc.begin(), mc.end(), 0.0);
    for (size_t i = 0; i<d[t].size() ; i++) {
      mc[d[t][i]] = mc[d[t][i]] + 1.0;
    }
    double pisum=0.0;
    for (int j = 0; j<2 ; j++) {
      std::gamma_distribution<double> rgamma(mc[j]+1.0,1.0);
      double rg = rgamma(gen);
      pitrait[t][j] = rg/d[t].size();
      pisum = pisum + pitrait[t][j];
    }
    for (int j = 0; j<2 ; j++) {
      pitrait[t][j] = pitrait[t][j]/pisum;
    }
  }
}


// Sample marker effects based on within and across pi BayesC 
void sampleBetaCSt(int i, 
                  int nt, 
                  std::vector<int>& dmarker,
                  std::vector<double>& pimarker,
                  std::vector<std::vector<double>>& pitrait, 
                  const arma::mat& E, 
                  const arma::mat& B,
                  const std::vector<std::vector<double>>& ww,
                  std::vector<std::vector<double>>& r, 
                  std::vector<std::vector<double>>& b,
                  std::vector<std::vector<int>>& d,
                  const std::vector<std::vector<int>>& XXindices, 
                  const std::vector<std::vector<std::vector<double>>>& XXvalues,
                  std::mt19937& gen) {
  
  std::vector<double> rhs(nt), loglik0(nt), loglik1(nt), bn(nt);
  double lik0t, lik1t, v0, v1, p0, rhs1, lhs1, diff;
  double u;

  // Compute rhs
  lik0t = 0.0;
  lik1t = 0.0;
  for (int t = 0; t < nt; t++) {
    rhs[t] = r[t][i] + ww[t][i] * b[t][i];
    v0 = ww[t][i]*E(t,t);
    v1 = ww[t][i]*E(t,t) + ww[t][i]*ww[t][i]*B(t,t);
    loglik0[t] = -0.5*std::log(v0) -0.5*((rhs[t]*rhs[t])/v0) + std::log(pitrait[t][0]);
    loglik1[t] = -0.5*std::log(v1) -0.5*((rhs[t]*rhs[t])/v1) + std::log(pitrait[t][1]);
    //lik0t = lik0t + loglik0[t];
    lik0t = lik0t - 0.5*std::log(v0) - 0.5*((rhs[t]*rhs[t])/v0);
    lik1t = lik1t + std::log(std::exp(loglik0[t])+std::exp(loglik1[t]));
  }
  lik0t = lik0t + std::log(pimarker[0]);
  lik1t = lik1t + std::log(pimarker[1]);
  p0 = 1.0/(std::exp(lik1t - lik0t) + 1.0);
  dmarker[i]=0;
  std::uniform_real_distribution<double> runif(0.0, 1.0);
  u = runif(gen);
  if(u>p0) dmarker[i]=1;
  
  for (int t = 0; t<nt ; t++) {
    d[t][i]=0;
    bn[t]=0.0;
  }
  if(dmarker[i]==1) {
    for (int t = 0; t<nt ; t++) {
      p0 = 1.0/(std::exp(loglik1[t] - loglik0[t]) + 1.0);
      d[t][i]=0;
      std::uniform_real_distribution<double> runif(0.0, 1.0);
      u = runif(gen);
      if(u>p0) d[t][i]=1;
      if(d[t][i]==1) {
        rhs1 = r[t][i] + ww[t][i]*b[t][i];
        lhs1 = ww[t][i] + E(t,t)/B(t,t);
        std::normal_distribution<double> rnorm(rhs1/lhs1, sqrt(E(t,t)/lhs1));
        bn[t] = rnorm(gen);
      } 
    }
  }
  // Adjust residuals based on sample marker effects
  for (int t = 0; t < nt; t++) {
    diff = (bn[t] - b[t][i]);
    if (diff != 0.0) {
      for (size_t j = 0; j < XXindices[i].size(); j++) {
        r[t][XXindices[i][j]] = r[t][XXindices[i][j]] - XXvalues[t][i][j] * diff;
      }
    }
    b[t][i] = bn[t];
  }
  
  // // Adjust residuals based on sampled marker effects
  // for ( int t = 0; t < nt; t++) {
  //   diff = (bn[t]-b[t][i]);
  //   if(diff!=0.0) {
  //     diff = (bn[t]-b[t][i])*std::sqrt(ww[t][i]);
  //     for (size_t j = 0; j < LDindices[i].size(); j++) {
  //       r[t][LDindices[i][j]] = r[t][LDindices[i][j]] - LDvalues[i][j]*diff*std::sqrt(ww[t][LDindices[i][j]]);
  //     }
  //   }
  //   b[t][i] = bn[t];
  // }
  
}

// Sample pi within trait BayesR
void samplePiR(int nt,
                 std::vector<std::vector<double>>& pitrait, 
                 const std::vector<std::vector<int>>& d,
                 std::mt19937& gen) {
  
  std::vector<double> mc(4);
  
  // Sample trait specific pi
  for (int t = 0; t < nt; t++) {
    std::fill(mc.begin(), mc.end(), 0.0);
    for (size_t i = 0; i<d[t].size() ; i++) {
      mc[d[t][i]] = mc[d[t][i]] + 1.0;
    }
    double pisum=0.0;
    for (int j = 0; j<4 ; j++) {
      std::gamma_distribution<double> rgamma(mc[j]+1.0,1.0);
      double rg = rgamma(gen);
      pitrait[t][j] = rg/d[t].size();
      pisum = pisum + pitrait[t][j];
    }
    for (int j = 0; j<4 ; j++) {
      pitrait[t][j] = pitrait[t][j]/pisum;
    }
  }
}

// Sample marker effects based on within and across pi BayesC
void sampleBetaR(int i,
                   int nt,
                   const std::vector<double>& gamma,
                   std::vector<int>& dmarker,
                   const std::vector<double>& pimarker,
                   const std::vector<std::vector<double>>& pitrait,
                   const arma::mat& E,
                   const arma::mat& B,
                   const std::vector<std::vector<double>>& ww,
                   std::vector<std::vector<double>>& r,
                   std::vector<std::vector<double>>& b,
                   std::vector<std::vector<int>>& d,
                   const std::vector<std::vector<int>>& XXindices,
                   const std::vector<std::vector<std::vector<double>>>& XXvalues,
                   std::mt19937& gen) {

  std::vector<double> rhs(nt), loglik0(nt), loglik1(nt), loglik2(nt), loglik3(nt), bn(nt);
  double lik0t, lik1t, v0, v1, v2, v3, p0, rhs1, lhs1, diff;
  double u;

  int nc = gamma.size();

  std::vector<std::vector<double>> probc(nt, std::vector<double>(nc, 0.0));
  std::vector<std::vector<double>> logLc(nt, std::vector<double>(nc, 0.0));

  double cumprobc, vbc, logLcAdj;


  // Compute rhs
  lik0t = 0.0;
  lik1t = 0.0;
  for (int t = 0; t < nt; t++) {
    rhs[t] = r[t][i] + ww[t][i] * b[t][i];
    v0 = ww[t][i]*E(t,t);
    v1 = ww[t][i]*E(t,t) + ww[t][i]*ww[t][i]*B(t,t)*0.01;
    v2 = ww[t][i]*E(t,t) + ww[t][i]*ww[t][i]*B(t,t)*0.1;
    v3 = ww[t][i]*E(t,t) + ww[t][i]*ww[t][i]*B(t,t);
    logLc[t][0] = -0.5*std::log(v0) - 0.5*((rhs[t]*rhs[t])/v0) + std::log(pitrait[t][0]);
    logLc[t][1] = -0.5*std::log(v1) - 0.5*((rhs[t]*rhs[t])/v1) + std::log(pitrait[t][1]);
    logLc[t][2] = -0.5*std::log(v2) - 0.5*((rhs[t]*rhs[t])/v2) + std::log(pitrait[t][2]);
    logLc[t][3] = -0.5*std::log(v3) - 0.5*((rhs[t]*rhs[t])/v3) + std::log(pitrait[t][3]);
    lik0t = lik0t  - 0.5*std::log(v0) - 0.5*((rhs[t]*rhs[t])/v0);
    lik1t = lik1t + std::log( std::exp(logLc[t][0])
                                + std::exp(logLc[t][1])
                                + std::exp(logLc[t][2])
                                + std::exp(logLc[t][3]));
  }
  // for (int t = 0; t < nt; t++) {
  //   rhs[t] = r[t][i] + ww[t][i] * b[t][i];
  //   v0 = E(t,t)/ww[t][i];
  //   v1 = E(t,t)/ww[t][i] + B(t,t)*gamma[1];
  //   v2 = E(t,t)/ww[t][i] + B(t,t)*gamma[2];
  //   v3 = E(t,t)/ww[t][i] + B(t,t)*gamma[3];
  //   b2 = (rhs[t]*rhs[t])/(ww[t][i]*ww[t][i]);
  //   logLc[t][0] = -0.5*std::log(v0) - 0.5*(b2/v0) + std::log(pitrait[t][0]);
  //   logLc[t][1] = -0.5*std::log(v1) - 0.5*(b2/v1) + std::log(pitrait[t][1]);
  //   logLc[t][2] = -0.5*std::log(v2) - 0.5*(b2/v2) + std::log(pitrait[t][2]);
  //   logLc[t][3] = -0.5*std::log(v3) - 0.5*(b2/v3) + std::log(pitrait[t][3]);
  //   lik0t = lik0t  - 0.5*std::log(v0) - 0.5*(b2/v0);
  //   lik1t = lik1t + std::log( std::exp(logLc[t][0]) 
  //                               + std::exp(logLc[t][1]) 
  //                               + std::exp(logLc[t][2])
  //                               + std::exp(logLc[t][3]));
  // }
  lik0t = lik0t + std::log(pimarker[0]);
  lik1t = lik1t + std::log(pimarker[1]);
  p0 = 1.0/(std::exp(lik1t - lik0t) + 1.0);
  dmarker[i]=0;
  std::uniform_real_distribution<double> runif(0.0, 1.0);
  u = runif(gen);
  if(u>p0) dmarker[i]=1;

  for (int t = 0; t<nt ; t++) {
    // rhs[t] = r[t][i] + ww[t][i] * b[t][i];
    // v0 = E(t,t)/ww[t][i];
    // v1 = E(t,t)/ww[t][i] + B(t,t)*0.01;
    // v2 = E(t,t)/ww[t][i] + B(t,t)*0.1;
    // v3 = E(t,t)/ww[t][i] + B(t,t);
    // b2 = (rhs[t]*rhs[t])/(ww[t][i]*ww[t][i]);
    // logLc[t][0] = -0.5*std::log(v0) - 0.5*(b2/v0) + std::log(pitrait[t][0]);
    // logLc[t][1] = -0.5*std::log(v1) - 0.5*(b2/v1) + std::log(pitrait[t][1]);
    // logLc[t][2] = -0.5*std::log(v2) - 0.5*(b2/v2) + std::log(pitrait[t][2]);
    // logLc[t][3] = -0.5*std::log(v3) - 0.5*(b2/v3) + std::log(pitrait[t][3]);
    d[t][i]=0;
    bn[t]=0.0;
  }
  if(dmarker[i]==1) {
    for (int t = 0; t<nt ; t++) {
      // variance class probability
      for (int j = 0; j<nc ; j++) {
        logLcAdj = 0.0;
        for (int k = 0; k<nc ; k++) {
          logLcAdj += std::exp(logLc[t][k] - logLc[t][j]);
        }
        probc[t][j] = 1.0/logLcAdj;
      }
      // sample variance class indicator
      std::uniform_real_distribution<double> runif(0.0, 1.0);
      u = runif(gen);
      d[t][i]=0;
      cumprobc = 0.0;
      for (int j = 0; j<nc ; j++) {
        cumprobc += probc[t][j];
        if(u < cumprobc){
          d[t][i] = j;
          break;
        }
      }
      // sample marker effect
      bn[t]=0.0;
      if(d[t][i]>0) {
        vbc = B(t,t)*gamma[d[t][i]];
        rhs1 = r[t][i] + ww[t][i]*b[t][i];
        lhs1 = ww[t][i] + E(t,t)/vbc;
        std::normal_distribution<double> rnorm(rhs1/lhs1, sqrt(E(t,t)/lhs1));
        bn[t] = rnorm(gen);
      }
    }
  }
  // Adjust residuals based on sample marker effects
  for (int t = 0; t < nt; t++) {
    diff = (bn[t] - b[t][i]);
    if (diff != 0.0) {
      for (size_t j = 0; j < XXindices[i].size(); j++) {
        r[t][XXindices[i][j]] = r[t][XXindices[i][j]] - XXvalues[t][i][j] * diff;
      }
    }
    b[t][i] = bn[t];
  }
}

// Sample marker effects based on within and across pi BayesC
void sampleBetaRS(int i,
                 int nt,
                 const std::vector<double>& gamma,
                 std::vector<int>& dmarker,
                 const std::vector<double>& pimarker,
                 const std::vector<std::vector<double>>& pitrait,
                 const arma::mat& E,
                 const arma::mat& B,
                 const std::vector<std::vector<double>>& ww,
                 std::vector<std::vector<double>>& r,
                 std::vector<std::vector<double>>& b,
                 std::vector<std::vector<int>>& d,
                 const std::vector<std::vector<int>>& LDindices,
                 const std::vector<std::vector<double>>& LDvalues,
                 std::mt19937& gen) {
  
  std::vector<double> rhs(nt), loglik0(nt), loglik1(nt), loglik2(nt), loglik3(nt), bn(nt);
  double lik0t, lik1t, v0, v1, v2, v3, p0, rhs1, lhs1, diff;
  double u;
  
  int nc = gamma.size();
  
  std::vector<std::vector<double>> probc(nt, std::vector<double>(nc, 0.0));
  std::vector<std::vector<double>> logLc(nt, std::vector<double>(nc, 0.0));
  
  double cumprobc, vbc, logLcAdj;
  
  
  // Compute rhs
  lik0t = 0.0;
  lik1t = 0.0;
  for (int t = 0; t < nt; t++) {
    rhs[t] = r[t][i] + ww[t][i] * b[t][i];
    v0 = ww[t][i]*E(t,t);
    v1 = ww[t][i]*E(t,t) + ww[t][i]*ww[t][i]*B(t,t)*0.01;
    v2 = ww[t][i]*E(t,t) + ww[t][i]*ww[t][i]*B(t,t)*0.1;
    v3 = ww[t][i]*E(t,t) + ww[t][i]*ww[t][i]*B(t,t);
    logLc[t][0] = -0.5*std::log(v0) - 0.5*((rhs[t]*rhs[t])/v0) + std::log(pitrait[t][0]);
    logLc[t][1] = -0.5*std::log(v1) - 0.5*((rhs[t]*rhs[t])/v1) + std::log(pitrait[t][1]);
    logLc[t][2] = -0.5*std::log(v2) - 0.5*((rhs[t]*rhs[t])/v2) + std::log(pitrait[t][2]);
    logLc[t][3] = -0.5*std::log(v3) - 0.5*((rhs[t]*rhs[t])/v3) + std::log(pitrait[t][3]);
    lik0t = lik0t  - 0.5*std::log(v0) - 0.5*((rhs[t]*rhs[t])/v0);
    lik1t = lik1t + std::log( std::exp(logLc[t][0])
                                + std::exp(logLc[t][1])
                                + std::exp(logLc[t][2])
                                + std::exp(logLc[t][3]));
  }
  lik0t = lik0t + std::log(pimarker[0]);
  lik1t = lik1t + std::log(pimarker[1]);
  p0 = 1.0/(std::exp(lik1t - lik0t) + 1.0);
  dmarker[i]=0;
  std::uniform_real_distribution<double> runif(0.0, 1.0);
  u = runif(gen);
  if(u>p0) dmarker[i]=1;
  
  for (int t = 0; t<nt ; t++) {
    d[t][i]=0;
    bn[t]=0.0;
  }
  if(dmarker[i]==1) {
    for (int t = 0; t<nt ; t++) {
      // variance class probability
      for (int j = 0; j<nc ; j++) {
        logLcAdj = 0.0;
        for (int k = 0; k<nc ; k++) {
          logLcAdj += std::exp(logLc[t][k] - logLc[t][j]);
        }
        probc[t][j] = 1.0/logLcAdj;
      }
      // sample variance class indicator
      std::uniform_real_distribution<double> runif(0.0, 1.0);
      u = runif(gen);
      d[t][i]=0;
      cumprobc = 0.0;
      for (int j = 0; j<nc ; j++) {
        cumprobc += probc[t][j];
        if(u < cumprobc){
          d[t][i] = j;
          break;
        }
      }
      // sample marker effect
      bn[t]=0.0;
      if(d[t][i]>0) {
        vbc = B(t,t)*gamma[d[t][i]];
        rhs1 = r[t][i] + ww[t][i]*b[t][i];
        lhs1 = ww[t][i] + E(t,t)/vbc;
        std::normal_distribution<double> rnorm(rhs1/lhs1, sqrt(E(t,t)/lhs1));
        bn[t] = rnorm(gen);
      }
      
    }
  }
  // Adjust residuals based on sampled marker effects
  for ( int t = 0; t < nt; t++) {
    diff = (bn[t]-b[t][i]);
    if(diff!=0.0) {
      diff = (bn[t]-b[t][i])*std::sqrt(ww[t][i]);
      for (size_t j = 0; j < LDindices[i].size(); j++) {
        r[t][LDindices[i][j]] = r[t][LDindices[i][j]] - LDvalues[i][j]*diff*std::sqrt(ww[t][LDindices[i][j]]);
      }
    }
    b[t][i] = bn[t];
  }
  
}


// Sample marker effect covariance matrix B
void sampleBR(int nt,
             int m,
             int nub,
             std::vector<double>& gamma,
             arma::mat& B,
             const std::vector<std::vector<int>>& d,
             const std::vector<std::vector<double>>& b,
             const std::vector<std::vector<double>>& ssb_prior,
             std::mt19937& gen) {
  // Initialize matrices and vectors for intermediate calculations
  arma::mat Sb(nt, nt, arma::fill::zeros);   // Sum of squared b values
  arma::mat corb(nt, nt, arma::fill::zeros); // Correlation matrix
  arma::mat dfB(nt, nt, arma::fill::zeros);  // Degrees of freedom matrix
  
  // Calculate Sb, dfB, and stdv based on input data
  for (int t1 = 0; t1 < nt; t1++) {
    double ssb = 0.0; // Initialize sum of squared b values for t1
    double dfb = 0.0; // Initialize degrees of freedom for t1
    
    // Calculate ssb and dfb for t1 based on d and b matrices
    for (int i = 0; i < m; i++) {
      if (d[t1][i]>0) {
        ssb += b[t1][i]*b[t1][i]/gamma[d[t1][i]];
        //ssb += b[t1][i]*b[t1][i];
        dfb += 1.0;
      }
    }
    
    // Store ssb and dfb in the corresponding matrices
    Sb(t1, t1) = ssb;
    dfB(t1, t1) = dfb;
    
    // Calculate correlations and off-diagonal elements of Sb and dfB
    for (int t2 = t1; t2 < nt; t2++) {
      ssb = 0.0; // Reset ssb for t2
      dfb = 0.0; // Reset dfb for t2
      
      if (t1 != t2) {
        for (int i = 0; i < m; i++) {
          if (d[t1][i] > 0 && d[t2][i] > 0) {
            //ssb += b[t1][i] * b[t2][i]/(gamma[d[t1][i]]*gamma[d[t2][i]]);
            ssb += b[t1][i] * b[t2][i];
            dfb += 1.0;
          }
        }
        
        // Populate the lower and upper triangles of Sb and dfB
        dfB(t1, t2) = dfb;
        dfB(t2, t1) = dfb;
        Sb(t1, t2) = ssb;
        Sb(t2, t1) = ssb;
      }
    }
  }
  
  // Calculate standard deviations and correlations (corb)
  std::vector<double> stdv(nt);
  for (int t = 0; t < nt; t++) {
    stdv[t] = sqrt(Sb(t,t));
  }
  
  for (int t1 = 0; t1 < nt; t1++) {
    for (int t2 = 0; t2 < nt; t2++) {
      if (t1 == t2) {
        corb(t1, t2) = 1.0; // Diagonal elements of corb are set to 1.0
      }
      
      if (t1 != t2 && stdv[t1] != 0.0 &&  stdv[t2] != 0.0) {
        //corb(t1, t2) = 0.0;
        corb(t1, t2) = Sb(t1, t2) / (stdv[t1] * stdv[t2]); // Calculate correlations
        corb(t2, t1) = corb(t1, t2);                      // Correlation matrix is symmetric
      }
    }
  }
  
  // Add a small constant to the diagonal elements of corb to ensure PD matrix
  for (int t1 = 0; t1 < nt; t1++) {
    corb(t1, t1) += 0.0001;
  }
  
  // Calculate diagonal elements of B using chi-squared distribution
  for (int t = 0; t < nt; t++) {
    std::chi_squared_distribution<double> rchisq(dfB(t, t) + nub);
    double chi2 = rchisq(gen);
    B(t, t) = (dfB(t, t)*Sb(t, t) + nub * ssb_prior[t][t]) / chi2;
  }
  
  // Calculate off-diagonal elements of B based on correlations and standard deviations
  for (int t1 = 0; t1 < nt; t1++) {
    for (int t2 = 0; t2 < nt; t2++) {
      if (t1 != t2) {
        B(t1, t2) = corb(t1, t2) * std::sqrt(B(t1, t1)) * std::sqrt(B(t2, t2));
      }
    }
  }
  
  // Check if B is symmetric; if not, make it symmetric
  bool issym = B.is_symmetric();
  if (!issym) {
    B = 0.5 * (B + B.t());
  }
  
  // Calculate the inverse of B and check for success
  arma::mat Bi(nt, nt, arma::fill::zeros);
  //bool success;
  //success = arma::inv_sympd(Bi, B, arma::inv_opts::allow_approx);
  
  // If the inverse calculation fails, print an error message
  //if (!success) {
  // std::cerr << "Error: Condition is false." << std::endl;
  //}
}


// [[Rcpp::export]]
std::vector<std::vector<std::vector<double>>>  mtbayes(   std::vector<std::vector<double>> y,
                                                          std::vector<std::vector<double>> W,
                                                          std::vector<std::vector<double>> b,
                                                          arma::mat B,
                                                          arma::mat E,
                                                          std::vector<std::vector<double>> ssb_prior,
                                                          std::vector<std::vector<double>> sse_prior,
                                                          std::vector<std::vector<int>> models,
                                                          std::vector<double> pi,
                                                          double nub,
                                                          double nue,
                                                          bool updateB,
                                                          bool updateE,
                                                          bool updatePi,
                                                          int nit,
                                                          int method) {
  
  
  
  // Define local variables
  int nt = y.size();
  int n = y[0].size();
  int m = W.size();
  int nmodels = models.size();
  
  double ssb, sse, dfb, u, logliksum, detC, diff, cumprobc;
  int mselect;
  
  std::vector<std::vector<double>> e(nt, std::vector<double>(n, 0.0));
  std::vector<std::vector<double>> g(nt, std::vector<double>(n, 0.0));
  std::vector<std::vector<double>> ww(nt, std::vector<double>(m, 0.0));
  std::vector<std::vector<int>> d(nt, std::vector<int>(m, 0));
  
  std::vector<double> mu(nt), rhs(nt), conv(nt);
  std::vector<double> pmodel(nmodels), pcum(nmodels), loglik(nmodels), cmodel(nmodels);
  std::vector<double> pis(nmodels);
  
  std::vector<std::vector<double>> bm(nt, std::vector<double>(m, 0.0));
  std::vector<std::vector<double>> dm(nt, std::vector<double>(m, 0.0));
  std::vector<std::vector<double>> ves(nt, std::vector<double>(nit, 0.0));
  std::vector<std::vector<double>> vbs(nt, std::vector<double>(nit, 0.0));
  std::vector<std::vector<double>> cvbm(nt, std::vector<double>(nt, 0.0));
  std::vector<std::vector<double>> cvem(nt, std::vector<double>(nt, 0.0));
  std::vector<std::vector<double>> mus(nt, std::vector<double>(nit, 0.0));
  
  std::vector<double> x2t(m);
  std::vector<std::vector<double>> x2(nt, std::vector<double>(m, 0.0));
  std::vector<std::vector<double>> wy(nt, std::vector<double>(m, 0.0));
  std::vector<std::vector<double>> r(nt, std::vector<double>(m, 0.0));
  std::vector<int> order(m);
  
  
  // Prior variance and degrees of freedom
  // fix this dfe and n should be vectors of length nt
  // dfe = n;
  
  // Mean adjust y and initialize e
  for ( int t = 0; t < nt; t++) {
    mu[t] = std::accumulate(y[t].begin(), y[t].end(), 0.0)/n;
    for ( int i = 0; i < n; i++) {
      e[t][i] = y[t][i] - mu[t];
    }                
  }
  
  // Initialize variables
  for (int i = 0; i < m; i++) {
    for (int t = 0; t < nt; t++) {
      ww[t][i] = 0.0;
      for (int j = 0; j < n; j++) {
        ww[t][i] = ww[t][i] + W[i][j]*W[i][j];
      }
    }
  }
  
  std::fill(cmodel.begin(), cmodel.end(), 1.0);
  std::fill(pis.begin(), pis.end(), 0.0);
  
  for ( int t = 0; t < nt; t++) {
    for ( int i = 0; i < m; i++) {
      r[t][i] = 0.0;
      for (int j = 0; j < n; j++) {
        r[t][i] = r[t][i] + W[i][j]*e[t][j];
      }
      r[t][i] = ww[t][i]*b[t][i];
      wy[t][i] = r[t][i];
      x2[t][i] = (wy[t][i]/ww[t][i])*(wy[t][i]/ww[t][i]);
    }
  }
  for ( int i = 0; i < m; i++) {
    x2t[i] = 0.0;
    for ( int t = 0; t < nt; t++) { 
      x2t[i] = x2t[i] + x2[t][i];
    }
  }

  // Establish order of markers as they are entered into the model
  std::iota(order.begin(), order.end(), 0);
  std::sort(  std::begin(order), 
              std::end(order),
              [&](int i1, int i2) { return x2t[i1] > x2t[i2]; } );
  
  
  // Initialize (co)variance matrices
  arma::mat C(nt,nt, fill::zeros);
  arma::mat Bi = arma::inv(B);
  arma::mat Ei = arma::inv(E);
  
  
  // Start Gibbs sampler
  std::random_device rd;
  unsigned int seed;
  seed = rd();
  std::mt19937 gen(seed);
  
  for ( int it = 0; it < nit; it++) {
    for ( int t = 0; t < nt; t++) {
      conv[t] = 0.0;
    }
    
    // Sample marker effects (BayesN)
    if (method<2) {
      for ( int isort = 0; isort < m; isort++) {
        int i = order[isort];
        for ( int t = 0; t < nt; t++) {
          rhs[t] = 0.0;
          for ( int j = 0; j < n; j++) {
            rhs[t] = rhs[t] + Ei(t,t)*W[i][j]*e[t][j]; 
          }
          //rhs[t] = Ei(t,t)*rhs[t] + Ei(t,t)*ww[t][i]*b[t][i];
          rhs[t] = rhs[t] + Ei(t,t)*ww[t][i]*b[t][i];
        }
        
        for ( int t = 0; t < nt; t++) { 
          d[t][i] = 1;
        }
        arma::mat C = Bi;
        for ( int t = 0; t < nt; t++) {
          C(t,t) = C(t,t) + ww[t][i]*Ei(t,t);       
        }
        arma::mat Ci = arma::inv(C);
        arma::mat mub = mvrnormARMA(Ci);
        
        for ( int t1 = 0; t1 < nt; t1++) {
          if(method==0) {
            mub(0,t1) = 0.0;
          }  
          for ( int t2 = 0; t2 < nt; t2++) {
            mub(0,t1) = mub(0,t1) + Ci(t1,t2)*rhs[t2];
          }
        } 
        for ( int t = 0; t < nt; t++) {
          diff = mub(0,t)-b[t][i];
          for (int j = 0; j < n; j++) {
            e[t][j]=e[t][j] - W[i][j]*(diff);
          }
          conv[t] = conv[t] + diff*diff;
          b[t][i] = mub(0,t);
          dm[t][i] = dm[t][i] + 1.0;
          bm[t][i] = bm[t][i] + b[t][i];     
        }
      }
    }
    
    // Sample marker effects (BayesC)
    if (method==4) {
      for ( int isort = 0; isort < m; isort++) {
        int i = order[isort];
        
        // Compute rhs
        for ( int t = 0; t < nt; t++) {
          rhs[t] = 0.0;
          for ( int j = 0; j < n; j++) {
            rhs[t] = rhs[t] + Ei(t,t)*W[i][j]*e[t][j]; 
          }
          rhs[t] = rhs[t] + Ei(t,t)*ww[t][i]*b[t][i];
        }
        
        // Sample variance class indicator
        for ( int k = 0; k < nmodels; k++) {
          arma::mat C = Bi;
          for ( int t1 = 0; t1 < nt; t1++) {
            if(models[k][t1]==1) {
              C(t1,t1) = C(t1,t1) + ww[t1][i]*Ei(t1,t1);       
            } 
          }
          arma::mat Ci = arma::inv(C);
          detC = arma::det(Ci);
          loglik[k] = 0.5*std::log(detC) + std::log(pi[k]);
          for ( int t1 = 0; t1 < nt; t1++) {
            for ( int t2 = t1; t2 < nt; t2++) {
              if(models[k][t1]==1 && models[k][t2]==1) {
                loglik[k] = loglik[k] + 0.5*rhs[t1]*rhs[t2]*Ci(t1,t2);
                if(t1!=t2) {
                  loglik[k] = loglik[k] + 0.5*rhs[t2]*rhs[t1]*Ci(t2,t1);
                  //loglik[k] = loglik[k] + 0.5*rhs[t1]*rhs[t2]*Ci(t1,t2);
                }
              }
            }
          }
        }
        
        // Compute marker variance class probability
        std::fill(pmodel.begin(), pmodel.end(), 0.0);
        for (int k = 0; k<nmodels ; k++) {
          logliksum = 0.0;
          for (int l = 0; l<nmodels ; l++) {
            logliksum += std::exp(loglik[l] - loglik[k]);
          }
          pmodel[k] = 1.0/logliksum;
        }
        
        // Sample variance class indicator
        std::uniform_real_distribution<double> runif(0.0, 1.0);
        u = runif(gen);
        mselect=0;
        cumprobc = 0.0;
        for (int k = 0; k<nmodels ; k++) {
          cumprobc += pmodel[k];
          if(u < cumprobc){
            mselect = k;
            break;
          }
        }
        cmodel[mselect] = cmodel[mselect] + 1.0; 
        for ( int t = 0; t < nt; t++) { 
          d[t][i] = models[mselect][t];
        }
        
        // Sample marker effect conditional on variance class indicator
        arma::mat C = Bi;
        for ( int t = 0; t < nt; t++) {
          if(models[mselect][t]==1) {
            C(t,t) = C(t,t) + ww[t][i]*Ei(t,t);       
          } 
        }
        arma::mat Ci = arma::inv(C);
        arma::mat mub = mvrnormARMA(Ci);
        for ( int t1 = 0; t1 < nt; t1++) {
          if(models[mselect][t1]!=1) {
            mub(0,t1) = 0.0;  
          } 
          for ( int t2 = 0; t2 < nt; t2++) {
            if(models[mselect][t2]==1) {
              mub(0,t1) = mub(0,t1) + Ci(t1,t2)*rhs[t2];
            }
          }
        } 
        
        // Adjust residuals based on sampled marker effects
        for ( int t = 0; t < nt; t++) {
          diff = mub(0,t)-b[t][i];
          for (int j = 0; j < n; j++) {
            e[t][j]=e[t][j] - W[i][j]*(diff);
          }
          conv[t] = conv[t] + diff*diff;
          b[t][i] = mub(0,t);
          if(d[t][i]==1) {
            dm[t][i] = dm[t][i] + 1.0;
            bm[t][i] = bm[t][i] + b[t][i];     
          }
        }
      }
      
      // Sample pi for Bayes C
      if(updatePi) {
        for (int k = 0; k<nmodels ; k++) {
          std::gamma_distribution<double> rgamma(cmodel[k],1.0);
          double rg = rgamma(gen);
          pi[k] = rg;
        }
        double psum = std::accumulate(pi.begin(), pi.end(), 0.0);
        for (int k = 0; k<nmodels ; k++) {
          pi[k] = pi[k]/psum;
          pis[k] = pis[k] + pi[k];
        }
        std::fill(cmodel.begin(), cmodel.end(), 1.0);
      }
    }
    
    
    // Sample marker variance
    if(updateB) {
      arma::mat Sb(nt,nt, fill::zeros);
      dfb = 0.0;
      for (int t1 = 0; t1 < nt; t1++) {
        for (int t2 = t1; t2 < nt; t2++) {
          ssb = 0.0;
          if (t1==t2) {
            for (int i=0; i < m; i++) {
              if( d[t1][i]==1 ) {
                ssb = ssb + b[t1][i]*b[t1][i];  
                dfb = dfb + 1.0;
              }
            }
            Sb(t1,t1) = ssb + nub*ssb_prior[t1][t1];
          } 
          if (t1!=t2) {
            for (int i=0; i < m; i++) {
              if( d[t1][i]==1 && d[t2][i]==1 ) {
                ssb = ssb + b[t1][i]*b[t2][i];  
              }
            }
            Sb(t1,t2) = ssb + nub*ssb_prior[t1][t2];
            Sb(t2,t1) = ssb + nub*ssb_prior[t2][t1];
          } 
        }
      }
      int dfSb = dfb/nt + nub;
      arma::mat B = riwishart(dfSb, Sb);
      for (int t = 0; t < nt; t++) {
        vbs[t][it] = B(t,t);
      }
      for (int t1 = 0; t1 < nt; t1++) {
        for (int t2 = 0; t2 < nt; t2++) {
          cvbm[t1][t2] = cvbm[t1][t2] + B(t1,t2);
        } 
      } 
      arma::mat Bi = arma::inv(B);
    }
    
    
    // Sample residual variance
    if(updateE) {
      arma::mat Se(nt,nt, fill::zeros);
      for ( int t1 = 0; t1 < nt; t1++) {
        for ( int t2 = t1; t2 < nt; t2++) {
          sse = 0.0;
          for ( int j = 0; j < n; j++) {
            sse = sse + e[t1][j]*e[t2][j];
          }
          Se(t1,t2) = sse + nue*sse_prior[t1][t2];
          if (t1!=t2) {
            Se(t1,t2) = Se(t2,t1);
          }
        }
      }
      int dfSe = n + nue;
      arma::mat E = riwishart(dfSe, Se);
      arma::mat Ei = arma::inv(E);
      for (int t = 0; t < nt; t++) {
        ves[t][it] = E(t,t);
      }
      for (int t1 = 0; t1 < nt; t1++) {
        for (int t2 = 0; t2 < nt; t2++) {
          cvem[t1][t2] = cvem[t1][t2] + E(t1,t2);
        } 
      } 
    }
    
    // Update mu and adjust residuals
    for ( int t = 0; t < nt; t++) {
      mu[t] = std::accumulate(e[t].begin(), e[t].end(), 0.0)/n;
      for ( int i = 0; i < n; i++) {
        e[t][i] = e[t][i] - mu[t];
      }                
      mus[t][it] = mu[t];
    }
  }
  
  // Summarize results
  std::vector<std::vector<std::vector<double>>> result;
  result.resize(15);
  
  result[0].resize(nt);
  result[1].resize(nt);
  result[2].resize(nt);
  result[3].resize(nt);
  result[4].resize(nt);
  result[5].resize(nt);
  result[6].resize(nt);
  result[7].resize(nt);
  result[8].resize(nt);
  result[9].resize(nt);
  result[10].resize(nt);
  result[11].resize(nt);
  result[12].resize(nt);
  result[13].resize(nt);
  result[14].resize(nt);
  
  for (int t=0; t < nt; t++) {
    result[0][t].resize(m);
    result[1][t].resize(m);
    result[2][t].resize(nit);
    result[3][t].resize(nit);
    result[4][t].resize(nit);
    result[5][t].resize(nt);
    result[6][t].resize(nt);
    result[7][t].resize(n);
    result[8][t].resize(n);
    result[9][t].resize(m);
    result[10][t].resize(nt);
    result[11][t].resize(nt);
    result[12][t].resize(nmodels);
    result[13][t].resize(nmodels);
    result[14][t].resize(m);
  }
  
  for (int t=0; t < nt; t++) {
    for (int i=0; i < m; i++) {
      result[0][t][i] = bm[t][i]/nit;
      result[1][t][i] = dm[t][i]/nit;
      result[9][t][i] = b[t][i];
      result[14][t][i] = order[i];
    }
    for (int i=0; i < n; i++) {
      result[7][t][i] = y[t][i] - mu[t] - e[t][i];
      result[8][t][i] = e[t][i];
    }
  }
  
  for (int t=0; t < nt; t++) {
    for (int i=0; i < nit; i++) {
      result[2][t][i] = mus[t][i];
      result[3][t][i] = vbs[t][i];
      result[4][t][i] = ves[t][i];
    }
  }
  for (int t1=0; t1 < nt; t1++) {
    for (int t2=0; t2 < nt; t2++) {
      result[5][t1][t2] = cvbm[t1][t2]/nit;
      result[6][t1][t2] = cvem[t1][t2]/nit;
    }
  }
  for (int t1=0; t1 < nt; t1++) {
    for (int t2=0; t2 < nt; t2++) {
      result[5][t1][t2] = cvbm[t1][t2]/nit;
      result[6][t1][t2] = cvem[t1][t2]/nit;
      result[10][t1][t2] = B(t1,t2);
      //result[10][t2][t1] = B(t2,t1);
      result[11][t1][t2] = E(t1,t2);
      //result[11][t2][t1] = E(t2,t1);
    }
  }
  for (int t=0; t < nt; t++) {
    for (int i=0; i < nmodels; i++) {
      result[12][t][i] = pi[i];
      result[13][t][i] = pis[i]/nit;
    }
  }  
  return result;
  
}


// [[Rcpp::export]]
std::vector<std::vector<std::vector<double>>>  mtsbayes(   std::vector<std::vector<double>> wy,
                                                           std::vector<std::vector<double>> ww,
                                                           std::vector<double> yy,
                                                           std::vector<std::vector<double>> b,
                                                           std::vector<std::vector<double>> LDvalues, 
                                                           std::vector<std::vector<int>> LDindices, 
                                                           arma::mat B,
                                                           arma::mat E,
                                                           std::vector<std::vector<double>> ssb_prior,
                                                           std::vector<std::vector<double>> sse_prior,
                                                           std::vector<std::vector<int>> models,
                                                           std::vector<double> pi,
                                                           double nub,
                                                           double nue,
                                                           bool updateB,
                                                           bool updateE,
                                                           bool updatePi,
                                                           std::vector<int> n,
                                                           int nit,
                                                           int nburn,
                                                           int nthin,
                                                           int seed,
                                                           int method) {
  
  // Define local variables
  int nt = wy.size();
  int m = wy[0].size();
  int nmodels = models.size();
  double nsamples=0.0;
  
  double ssb, sse,ssg, dfb, u, logliksum, detC, diff, cumprobc;
  int mselect;
  
  std::vector<std::vector<int>> d(nt, std::vector<int>(m, 0));
  
  std::vector<double> mu(nt), rhs(nt);
  std::vector<double> pmodel(nmodels), pcum(nmodels), loglik(nmodels), cmodel(nmodels);
  std::vector<double> pis(nmodels);
  
  std::vector<std::vector<double>> bm(nt, std::vector<double>(m, 0.0));
  std::vector<std::vector<double>> dm(nt, std::vector<double>(m, 0.0));
  std::vector<std::vector<double>> ves(nt, std::vector<double>(nit+nburn, 0.0));
  std::vector<std::vector<double>> vbs(nt, std::vector<double>(nit+nburn, 0.0));
  std::vector<std::vector<double>> vgs(nt, std::vector<double>(nit+nburn, 0.0));
  std::vector<std::vector<double>> cvbm(nt, std::vector<double>(nt, 0.0));
  std::vector<std::vector<double>> cvem(nt, std::vector<double>(nt, 0.0));
  std::vector<std::vector<double>> cvgm(nt, std::vector<double>(nt, 0.0));
  std::vector<std::vector<double>> mus(nt, std::vector<double>(nit+nburn, 0.0));
  
  std::vector<double> x2t(m);
  std::vector<std::vector<double>> x2(nt, std::vector<double>(m, 0.0));
  std::vector<std::vector<double>> r(nt, std::vector<double>(m, 0.0));
  std::vector<int> order(m);
  std::vector<double> vei(m),vadj(m);

  std::vector<double> gamma(4, 0.0);
  
  std::vector<std::vector<double>> pitrait(nt, std::vector<double>(4, 0.0)); 
  std::vector<std::vector<double>> pistrait(nt, std::vector<double>(4, 0.0)); 
  std::vector<double> pimarker(2, 0.0);
  std::vector<double> pismarker(2, 0.0);
  std::vector<int> dmarker(m, 0);
  
  gamma[0] = 0.0;
  gamma[1] = 0.01;
  gamma[2] = 0.1;
  gamma[3] = 1.0;
  
  for (int t = 0; t < nt; t++) {
    pitrait[t][0] = 1.0-pi[0];
    pitrait[t][1] = pi[0];
  }
  if(method==5){
    for (int t = 0; t < nt; t++) {
      pitrait[t][0] = 0.95;
      pitrait[t][1] = 0.02;
      pitrait[t][2] = 0.02;
      pitrait[t][3] = 0.01;
    }
  }
  
  pimarker[0] = 1.0-pi[0];
  pimarker[1] = pi[0];
  
  // Initialize variables
  for (int i = 0; i < m; i++) {
    for (int t = 0; t < nt; t++) {
      r[t][i] = wy[t][i];
      x2[t][i] = (wy[t][i]/ww[t][i])*(wy[t][i]/ww[t][i]);
    }
  }

  // Wy - W'Wb
  for ( int i = 0; i < m; i++) {
    for (int t = 0; t < nt; t++) {
      if (b[t][i]!= 0.0) {
        diff = b[t][i]*ww[t][i];
        for (size_t j = 0; j < LDindices[i].size(); j++) {
          r[t][LDindices[i][j]]=r[t][LDindices[i][j]] - LDvalues[i][j]*diff;
        }
      }
    }
  }
  
  // for ( int i = 0; i < m; i++) {
  //   x2t[i] = 0.0;
  //   for ( int t = 0; t < nt; t++) { 
  //     x2t[i] = x2t[i] + x2[t][i];
  //   }
  // }
  for ( int i = 0; i < m; i++) {
    x2t[i] = 0.0;
    for ( int t = 0; t < nt; t++) { 
      if(x2[t][i]>x2t[i]) {
        x2t[i] = x2[t][i];
      }
    }
  }
  
  
  std::fill(cmodel.begin(), cmodel.end(), 1.0);
  std::fill(pis.begin(), pis.end(), 0.0);
  
  
  // Establish order of markers as they are entered into the model
  std::iota(order.begin(), order.end(), 0);
  std::sort(  std::begin(order), 
              std::end(order),
              [&](int i1, int i2) { return x2t[i1] > x2t[i2]; } );
  
  
  // Initialize (co)variance matrices
  arma::mat C(nt,nt, fill::zeros);
  arma::mat Bi = arma::inv(B);
  arma::mat Ei = arma::inv(E);
  arma::mat G(nt,nt, fill::zeros);
  
  arma::rowvec probs(nmodels, fill::zeros);
  
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
        
        // Compute rhs
        for ( int t = 0; t < nt; t++) {
          rhs[t] = Ei(t,t)*r[t][i] + Ei(t,t)*ww[t][i]*b[t][i];
        }
        
        // Compute log likelihood for each marker variance class  
        for ( int k = 0; k < nmodels; k++) {
          arma::mat C = Bi;
          for ( int t1 = 0; t1 < nt; t1++) {
            if(models[k][t1]==1) {
              C(t1,t1) = C(t1,t1) + ww[t1][i]*Ei(t1,t1);
            }
          }
          arma::mat Ci = arma::inv(C);
          detC = arma::det(Ci);
          loglik[k] = 0.5*std::log(detC) + std::log(pi[k]);
          for ( int t1 = 0; t1 < nt; t1++) {
            for ( int t2 = t1; t2 < nt; t2++) {
              if(models[k][t1]==1 && models[k][t2]==1) {
                loglik[k] = loglik[k] + 0.5*rhs[t1]*rhs[t2]*Ci(t1,t2);
                if(t1!=t2) {
                  loglik[k] = loglik[k] + 0.5*rhs[t2]*rhs[t1]*Ci(t2,t1);
                }
              }
            }
          }
        }

        // Compute marker variance class probability
        std::fill(pmodel.begin(), pmodel.end(), 0.0);
        for (int k = 0; k<nmodels ; k++) {
          logliksum = 0.0;
          for (int l = 0; l<nmodels ; l++) {
            logliksum += std::exp(loglik[l] - loglik[k]);
          }
          pmodel[k] = 1.0/logliksum;
        }

        // Sample variance class indicator
        std::uniform_real_distribution<double> runif(0.0, 1.0);
        u = runif(gen);
        mselect=0;
        cumprobc = 0.0;
        
        for (int k = 0; k<nmodels ; k++) {
          cumprobc += pmodel[k];
          if(u < cumprobc){
            mselect = k;
            break;
          }
        }
        cmodel[mselect] = cmodel[mselect] + 1.0; 
        for ( int t = 0; t < nt; t++) { 
          d[t][i] = models[mselect][t];
        }

        // Sample marker effect conditional on variance class indicator
        arma::mat C = Bi;
        for ( int t1 = 0; t1 < nt; t1++) {
          if(models[mselect][t1]==1) {
            C(t1,t1) = C(t1,t1) + ww[t1][i]*Ei(t1,t1);       
          } 
        }
        arma::mat Ci = arma::inv(C);
        arma::mat mub = mvrnormARMA(Ci);
        for ( int t1 = 0; t1 < nt; t1++) {
          if(models[mselect][t1]==0) {
            mub(0,t1) = 0.0;
          }
          if(models[mselect][t1]==1) {
            mub(0,t1) = mub(0,t1) + Ci(t1,t1)*rhs[t1];
            for ( int t2 = 0; t2 < nt; t2++) {
              if(t1!=t2 && models[mselect][t2]==1) {
                mub(0,t1) = mub(0,t1) + Ci(t1,t2)*rhs[t2];
              }
            }
          }
        }

        // Adjust residuals based on sampled marker effects
        for ( int t = 0; t < nt; t++) {
          diff = (mub(0,t)-b[t][i]);
          if(diff!=0.0) {
            diff = (mub(0,t)-b[t][i])*std::sqrt(ww[t][i]);
            for (size_t j = 0; j < LDindices[i].size(); j++) {
              r[t][LDindices[i][j]] = r[t][LDindices[i][j]] - LDvalues[i][j]*diff*std::sqrt(ww[t][LDindices[i][j]]);
            }
          }
          b[t][i] = mub(0,t);
          if ( d[t][i]==1 && (it > nburn) && (it % nthin == 0) ) {
            dm[t][i] = dm[t][i] + 1.0;
            bm[t][i] = bm[t][i] + b[t][i]; 
          }
        }
      }
    }
    
    // Sample marker effects (BayesR)
    if (method==5) {
      for ( int isort = 0; isort < m; isort++) {
        int i = order[isort];
        sampleBetaRS(i, nt, gamma,
                    dmarker, pimarker, pitrait,
                    E, B,
                    ww, r, b, d,
                    LDindices, LDvalues,
                    gen);
      }
      // Store values
      for (int t = 0; t < nt; t++) {
        for ( int i = 0; i < m; i++) {
          if (d[t][i] > 0) {
            if ((it > nburn) && (it % nthin == 0)) {
              dm[t][i] = dm[t][i] + 1.0;
              bm[t][i] = bm[t][i] + b[t][i];
            }
          }
        }
      }
      
    }
    
    // Sample pi for Bayes C
    if(updatePi && method==4) {
      for (int k = 0; k<nmodels ; k++) {
        std::gamma_distribution<double> rgamma(cmodel[k],1.0);
        double rg = rgamma(gen);
        pi[k] = rg;
      }
      double psum = std::accumulate(pi.begin(), pi.end(), 0.0);
      for (int k = 0; k<nmodels ; k++) {
        pi[k] = pi[k]/psum;
        if(it>nburn) pis[k] = pis[k] + pi[k];
      }
      std::fill(cmodel.begin(), cmodel.end(), 1.0);
    }
    
    // Sample pi for Bayes R
    if(updatePi && method==5) {
      samplePiMt(nt, pimarker, dmarker, gen);
      if(pimarker[0]<0.9) {
        pimarker[0]=0.9;
        pimarker[1]=1.0-pimarker[0];
      }
      //samplePiC(nt, pitrait, d, gen);
      samplePiR(nt, pitrait, d, gen);
      //Rcpp::Rcout << "pimarker[0]: " << pimarker[0] << std::endl;
      //Rcpp::Rcout << "pimarker[1]: " << pimarker[1] << std::endl;
      if(it>nburn) {
        for (int t = 0; t<nt ; t++) {
          for (size_t k = 0; k<gamma.size() ; k++) {
            pistrait[t][k] = pistrait[t][k] + pitrait[t][k];
          }
        }
        for (int k = 0; k<2 ; k++) {
          pismarker[k] = pismarker[k] + pimarker[k];
        }
      }
    }
    
    
    // Sample marker variance
    // if(updateB) {
    //   arma::mat Sb(nt,nt, fill::zeros);
    //   dfb = 0.0;
    //   for (int t1 = 0; t1 < nt; t1++) {
    //     for (int t2 = t1; t2 < nt; t2++) {
    //       ssb = 0.0;
    //       if (t1==t2) {
    //         for (int i=0; i < m; i++) {
    //           if( d[t1][i]==1 ) {
    //             ssb = ssb + b[t1][i]*b[t1][i];
    //             dfb = dfb + 1.0;
    //           }
    //         }
    //         Sb(t1,t1) = ssb + nub*ssb_prior[t1][t1];
    //       }
    //       if (t1!=t2) {
    //         for (int i=0; i < m; i++) {
    //           if( d[t1][i]==1 && d[t2][i]==1 ) {
    //             ssb = ssb + b[t1][i]*b[t2][i];
    //           }
    //         }
    //         Sb(t1,t2) = ssb + nub*ssb_prior[t1][t2];
    //         Sb(t2,t1) = ssb + nub*ssb_prior[t2][t1];
    //       }
    //     }
    //   }
    //   dfb = 0.0;
    //   for (int i=0; i < m; i++) {
    //     double dfb0=0.0;
    //     for (int t1 = 0; t1 < nt; t1++) {
    //       if( d[t1][i]==1 ) {
    //         dfb0=1.0;
    //       }
    //     }
    //     dfb = dfb + dfb0;
    //   }
    //   int dfSb = dfb + nub;
    //   arma::mat B = riwishart(dfSb, Sb);
    //   for (int t = 0; t < nt; t++) {
    //     vbs[t][it] = B(t,t);
    //   }
    //   for (int t1 = 0; t1 < nt; t1++) {
    //     for (int t2 = 0; t2 < nt; t2++) {
    //       if(it>nburn) cvbm[t1][t2] = cvbm[t1][t2] + B(t1,t2);
    //     }
    //   }
    //   arma::mat Bi = arma::inv(B);
    // }
    if(updateB && method==4) {
      arma::mat Sb(nt,nt, fill::zeros);
      arma::mat corb(nt,nt, fill::zeros);
      arma::mat dfB(nt,nt, fill::zeros);
      for (int t1 = 0; t1 < nt; t1++) {
        ssb = 0.0;
        dfb = 0.0;
        for (int i=0; i < m; i++) {
          if( d[t1][i]==1 ) {
            ssb = ssb + b[t1][i]*b[t1][i];
            dfb = dfb + 1.0;
          }
        }
        Sb(t1,t1) = ssb;
        dfB(t1,t1) = dfb;
        for (int t2 = t1; t2 < nt; t2++) {
          ssb = 0.0;
          dfb = 0.0;
          if (t1!=t2) {
            for (int i=0; i < m; i++) {
              if( d[t1][i]==1 && d[t2][i]==1 ) {
                ssb = ssb + b[t1][i]*b[t2][i];
                dfb = dfb + 1.0;
              }
            }
            dfB(t1,t2) = dfb;
            dfB(t2,t1) = dfb;
            Sb(t1,t2) = ssb;
            Sb(t2,t1) = ssb;
          }
        }
      }
      std::vector<double> stdv(nt);
      for (int t = 0; t < nt; t++) {
        stdv[t] = sqrt(Sb(t,t));
      }
      for (int t1 = 0; t1 < nt; t1++) {
        for (int t2 = 0; t2 < nt; t2++) {
          if (t1==t2) {
            corb(t1,t2) = 1.0;
          }
          if (t1!=t2 && stdv[t1]!=0.0 && stdv[t2]!=0.0 ) {
            //corb(t1,t2) = 0.0;
            corb(t1,t2) = Sb(t1,t2)/(stdv[t1] * stdv[t2]);
            corb(t2,t1) = corb(t1,t2);
          }
        }
      }
      
      arma::mat B(nt,nt, fill::zeros);
      for (int t = 0; t < nt; t++) {
        std::chi_squared_distribution<double> rchisq(dfB(t,t) + nub);
        double chi2 = rchisq(gen);
        B(t,t) = (Sb(t,t) + nub*ssb_prior[t][t])/chi2;
      }
      for (int t1 = 0; t1 < nt; t1++) {
        for (int t2 = 0; t2 < nt; t2++) {
          if (t1!=t2) {
            B(t1,t2) = corb(t1,t2)*sqrt(B(t1,t1))*sqrt(B(t2,t2));
          }
        }
      }
      for (int t1 = 0; t1 < nt; t1++) {
        B(t1,t1) = B(t1,t1) + 0.0001;
      }
      for (int t = 0; t < nt; t++) {
        vbs[t][it] = B(t,t);
      }
      for (int t1 = 0; t1 < nt; t1++) {
        for (int t2 = 0; t2 < nt; t2++) {
          if(it>nburn) cvbm[t1][t2] = cvbm[t1][t2] + B(t1,t2);
        }
      }
      bool issym = B.is_symmetric();
      if (!issym){
        B = 0.5*(B + B.t());
      }
      arma::mat Bi(nt,nt, fill::zeros);
      //bool success;
      //success = inv_sympd(Bi,B,inv_opts::allow_approx);
      //if(!success) {
      //  std::cerr << "Error: Condition is false." << std::endl;
      //}
      //arma::mat Bi = arma::inv(B);
    }

    // Sample marker variance
    if(updateB && method==5) {
      sampleBR(nt, m, nub, gamma, B, d, b, ssb_prior, gen);
      for (int t = 0; t < nt; t++) {
        vbs[t][it] = B(t,t);
      }
      for (int t1 = 0; t1 < nt; t1++) {
        for (int t2 = 0; t2 < nt; t2++) {
          if(it>nburn) cvbm[t1][t2] = cvbm[t1][t2] + B(t1,t2);
        }
      }
    }
    
    //Update genetic variance
    arma::mat G(nt,nt, fill::zeros);
    for ( int t1 = 0; t1 < nt; t1++) {
      for ( int t2 = t1; t2 < nt; t2++) {
        if (t1==t2) {
          ssg = 0.0;
          for ( int i = 0; i < m; i++) {
            ssg = ssg + b[t1][i]*(wy[t1][i]-r[t1][i]);
          }
          G(t1,t2) = ssg/(std::sqrt((double)n[t1])*std::sqrt((double)n[t2]));
        }
        if (t1!=t2) {
          ssg = 0.0;
          for ( int i = 0; i < m; i++) {
            ssg = ssg + b[t1][i]*(wy[t1][i]-r[t1][i]);
            ssg = ssg + b[t2][i]*(wy[t2][i]-r[t2][i]);
          }
          ssg = ssg/2.0;
          G(t1,t2) = ssg/(std::sqrt((double)n[t1])*std::sqrt((double)n[t2]));
          G(t2,t1) = G(t1,t2);
        }
      }
    }
    for (int t = 0; t < nt; t++) {
      vgs[t][it] = G(t,t);
    }
    for (int t1 = 0; t1 < nt; t1++) {
      for (int t2 = 0; t2 < nt; t2++) {
        if(it>nburn) cvgm[t1][t2] = cvgm[t1][t2] + G(t1,t2);
      } 
    } 
    
    
    // Sample residual variance 
    // (current version assume uncorrelated residuals among traits traits)
    if(updateE) {
      arma::mat Se(nt,nt, fill::zeros);
      for ( int t1 = 0; t1 < nt; t1++) {
        for ( int t2 = t1; t2 < nt; t2++) {
          if (t1==t2) {
            sse = 0.0;
            for ( int i = 0; i < m; i++) {
              sse = sse + b[t1][i] * (r[t1][i] + wy[t1][i]);
            }
            sse = yy[t1] - sse;
            Se(t1,t1) = sse + nue*sse_prior[t1][t1];
          }
        }
      }
      arma::mat E(nt,nt, fill::zeros);
      for (int t = 0; t < nt; t++) {
        std::chi_squared_distribution<double> rchisq(n[t] + nue);
        double chi2 = rchisq(gen);
        E(t,t) = (Se(t,t))/chi2;
      }
      //int dfSe = n[0] + nue;
      //arma::mat E = riwishart(dfSe, Se);
      arma::mat Ei = arma::inv(E);
      for (int t = 0; t < nt; t++) {
        ves[t][it] = E(t,t);
      }
      for (int t1 = 0; t1 < nt; t1++) {
        for (int t2 = 0; t2 < nt; t2++) {
          if(it>nburn) cvem[t1][t2] = cvem[t1][t2] + E(t1,t2);
        } 
      } 
    }
    
  
  }
  
  // // Summarize results
  // std::vector<std::vector<std::vector<double>>> result;
  // result.resize(15);
  // 
  // result[0].resize(nt);
  // result[1].resize(nt);
  // result[2].resize(nt);
  // result[3].resize(nt);
  // result[4].resize(nt);
  // result[5].resize(nt);
  // result[6].resize(nt);
  // result[7].resize(nt);
  // result[8].resize(nt);
  // result[9].resize(nt);
  // result[10].resize(nt);
  // result[11].resize(nt);
  // result[12].resize(nt);
  // result[13].resize(nt);
  // result[14].resize(nt);
  // 
  // for (int t=0; t < nt; t++) {
  //   result[0][t].resize(m);
  //   result[1][t].resize(m);
  //   result[2][t].resize(nit+nburn);
  //   result[3][t].resize(nit+nburn);
  //   result[4][t].resize(nit+nburn);
  //   result[5][t].resize(nt);
  //   result[6][t].resize(nt);
  //   result[7][t].resize(m);
  //   result[8][t].resize(m);
  //   result[9][t].resize(m);
  //   result[10][t].resize(nt);
  //   result[11][t].resize(nt);
  //   result[12][t].resize(nmodels);
  //   result[13][t].resize(nmodels);
  //   result[14][t].resize(m);
  // }
  // 
  // for (int t=0; t < nt; t++) {
  //   for (int i=0; i < m; i++) {
  //     //result[0][t][i] = bm[t][i]/nit;
  //     //result[1][t][i] = dm[t][i]/nit;
  //     result[0][t][i] = bm[t][i]/nsamples;
  //     result[1][t][i] = dm[t][i]/nsamples;
  //     result[9][t][i] = b[t][i];
  //     result[14][t][i] = order[i];
  //   }
  //   for (int i=0; i < m; i++) {
  //     result[7][t][i] = wy[t][i];
  //     result[8][t][i] = r[t][i];
  //   }
  // }
  // 
  // for (int t=0; t < nt; t++) {
  //   for (int i=0; i < nit+nburn; i++) {
  //     result[2][t][i] = vgs[t][i];
  //     result[3][t][i] = vbs[t][i];
  //     result[4][t][i] = ves[t][i];
  //   }
  // }
  // for (int t1=0; t1 < nt; t1++) {
  //   for (int t2=0; t2 < nt; t2++) {
  //     result[5][t1][t2] = cvbm[t1][t2]/nit;
  //     result[6][t1][t2] = cvem[t1][t2]/nit;
  //     result[10][t1][t2] = cvgm[t1][t2]/nit;
  //     result[11][t1][t2] = E(t1,t2);
  //   }
  // }
  // for (int t=0; t < nt; t++) {
  //   for (int i=0; i < nmodels; i++) {
  //     result[12][t][i] = pi[i];
  //     result[13][t][i] = pis[i]/nit;
  //   }
  // }
  // Summarize results
  std::vector<std::vector<std::vector<double>>> result;
  result.resize(20);
  
  result[0].resize(nt);
  result[1].resize(nt);
  result[2].resize(nt);
  result[3].resize(nt);
  result[4].resize(nt);
  result[5].resize(nt);
  result[6].resize(nt);
  result[7].resize(nt);
  result[8].resize(nt);
  result[9].resize(nt);
  result[10].resize(nt);
  result[11].resize(nt);
  result[12].resize(nt);
  result[13].resize(nt);
  result[14].resize(nt);
  result[15].resize(nt);
  result[16].resize(nt);
  result[17].resize(nt);
  result[18].resize(nt);
  result[19].resize(nt);
  
  for (int t=0; t < nt; t++) {
    result[0][t].resize(m);
    result[1][t].resize(m);
    result[2][t].resize(m);
    result[3][t].resize(m);
    result[4][t].resize(m);
    result[5][t].resize(m);
    result[6][t].resize(m);
    result[7][t].resize(nit+nburn);
    result[8][t].resize(nit+nburn);
    result[9][t].resize(nit+nburn);
    result[10][t].resize(nt);
    result[11][t].resize(nt);
    result[12][t].resize(nt);
    result[13][t].resize(nt);
    result[14][t].resize(nt);
    result[15][t].resize(nt);
    result[16][t].resize(nmodels);
    result[17][t].resize(nmodels);
    result[18][t].resize(4);
    result[19][t].resize(2);
  }
  
  for (int t=0; t < nt; t++) {
    for (int i=0; i < m; i++) {
      result[0][t][i] = bm[t][i]/nsamples;
      result[1][t][i] = dm[t][i]/nsamples;
      result[2][t][i] = wy[t][i];
      result[3][t][i] = r[t][i];
      result[4][t][i] = b[t][i];
      result[5][t][i] = d[t][i];
      result[6][t][i] = order[i];
    }
  }
  
  for (int t=0; t < nt; t++) {
    for (int i=0; i < nit+nburn; i++) {
      result[7][t][i] = vbs[t][i];
      result[8][t][i] = vgs[t][i];
      result[9][t][i] = ves[t][i];
    }
  }
  for (int t1=0; t1 < nt; t1++) {
    for (int t2=0; t2 < nt; t2++) {
      result[10][t1][t2] = cvbm[t1][t2]/nit;
      result[11][t1][t2] = cvgm[t1][t2]/nit;
      result[12][t1][t2] = cvem[t1][t2]/nit;
      result[13][t1][t2] = B(t1,t2);
      result[14][t1][t2] = G(t1,t2);
      result[15][t1][t2] = E(t1,t2);
    }
  }
  for (int t=0; t < nt; t++) {
    for (int i=0; i < nmodels; i++) {
      result[16][t][i] = pi[i];
      result[17][t][i] = pis[i]/nit;
    }
  }
  for (int t=0; t < nt; t++) {
    for (int i=0; i < 4; i++) {
      result[18][t][i] = pistrait[t][i]/nit;
    }
    for (int i=0; i < 2; i++) {
      result[19][t][i] = pismarker[i]/nit;
    }
  }
  return result;
}

// [[Rcpp::export]]
std::vector<std::vector<std::vector<double>>>  mtblr(   std::vector<std::vector<double>> wy,
                                                        std::vector<std::vector<double>> ww,
                                                        std::vector<double> yy,
                                                        std::vector<std::vector<double>> b,
                                                        const std::vector<std::vector<std::vector<double>>>& XXvalues,
                                                        const std::vector<std::vector<int>>& XXindices,
                                                        arma::mat B,
                                                        arma::mat E,
                                                        std::vector<std::vector<double>> ssb_prior,
                                                        std::vector<std::vector<double>> sse_prior,
                                                        std::vector<std::vector<int>> models,
                                                        std::vector<double> pi,
                                                        double nub,
                                                        double nue,
                                                        bool updateB,
                                                        bool updateE,
                                                        bool updatePi,
                                                        std::vector<int> n,
                                                        int nit,
                                                        int nburn,
                                                        int nthin,
                                                        int seed,
                                                        int method) {
  
  // Define local variables
  int nt = wy.size();
  int m = wy[0].size();
  int nmodels = models.size();
  double nsamples=0.0;
  
  //double logliksum, detC, diff, cumprobc;
  //int mselect;
  
  std::vector<std::vector<int>> d(nt, std::vector<int>(m, 0));
  
  std::vector<double> mu(nt), rhs(nt), conv(nt);
  std::vector<double> pmodel(nmodels), pcum(nmodels), loglik(nmodels), cmodel(nmodels);
  std::vector<double> pis(nmodels, 0.0);
  
  std::vector<std::vector<double>> bm(nt, std::vector<double>(m, 0.0));
  std::vector<std::vector<double>> dm(nt, std::vector<double>(m, 0.0));
  std::vector<std::vector<double>> ves(nt, std::vector<double>(nit+nburn, 0.0));
  std::vector<std::vector<double>> vbs(nt, std::vector<double>(nit+nburn, 0.0));
  std::vector<std::vector<double>> vgs(nt, std::vector<double>(nit+nburn, 0.0));
  std::vector<std::vector<double>> cvbm(nt, std::vector<double>(nt, 0.0));
  std::vector<std::vector<double>> cvem(nt, std::vector<double>(nt, 0.0));
  std::vector<std::vector<double>> cvgm(nt, std::vector<double>(nt, 0.0));
  std::vector<std::vector<double>> mus(nt, std::vector<double>(nit+nburn, 0.0));
  
  std::vector<double> x2t(m);
  std::vector<std::vector<double>> x2(nt, std::vector<double>(m, 0.0));
  std::vector<std::vector<double>> r(nt, std::vector<double>(m, 0.0));
  std::vector<int> order(m);
  std::vector<double> gamma(4, 0.0);

  std::vector<std::vector<double>> pitrait(nt, std::vector<double>(4, 0.0)); 
  std::vector<std::vector<double>> pistrait(nt, std::vector<double>(4, 0.0)); 
  std::vector<double> pimarker(2, 0.0);
  std::vector<double> pismarker(2, 0.0);
  std::vector<int> dmarker(m, 0);

  gamma[0] = 0.0;
  gamma[1] = 0.01;
  gamma[2] = 0.1;
  gamma[3] = 1.0;
  
  for (int t = 0; t < nt; t++) {
    pitrait[t][0] = 1.0-pi[0];
    pitrait[t][1] = pi[0];
  }
  if(method==5){
    for (int t = 0; t < nt; t++) {
      pitrait[t][0] = 0.95;
      pitrait[t][1] = 0.02;
      pitrait[t][2] = 0.02;
      pitrait[t][3] = 0.01;
    }
  }
  
  pimarker[0] = 1.0-pi[0];
  pimarker[1] = pi[0];
  
  // Initialize variables
  for (int i = 0; i < m; i++) {
    for (int t = 0; t < nt; t++) {
      r[t][i] = wy[t][i];
      x2[t][i] = (wy[t][i]/ww[t][i])*(wy[t][i]/ww[t][i]);
    }
  }
  
  // Wy - W'Wb
  for ( int i = 0; i < m; i++) {
    for (int t = 0; t < nt; t++) {
      if (b[t][i]!= 0.0) {
        for (size_t j = 0; j < XXindices[i].size(); j++) {
          r[t][XXindices[i][j]]=r[t][XXindices[i][j]] - XXvalues[t][i][j]*b[t][i];
        }
      }
    }
  }
  
  for ( int i = 0; i < m; i++) {
    x2t[i] = 0.0;
    for ( int t = 0; t < nt; t++) { 
      if(x2[t][i]>x2t[i]) {
        x2t[i] = x2[t][i];
      }
    }
  }
  
  
  // Establish order of markers as they are entered into the model
  std::iota(order.begin(), order.end(), 0);
  std::sort(  std::begin(order),
              std::end(order),
              [&](int i1, int i2) { return x2t[i1] > x2t[i2]; } );
  
  
  // Initialize (co)variance matrices
  arma::mat C(nt,nt, fill::zeros);
  arma::mat Bi = arma::inv(B);
  arma::mat Ei = arma::inv(E);
  arma::mat G(nt,nt, fill::zeros);
  arma::rowvec probs(nmodels, fill::zeros);
  
  // Start Gibbs sampler
  std::random_device rd;
  std::mt19937 gen(seed);
  
  for ( int it = 0; it < nit+nburn; it++) {
    
    std::fill(cmodel.begin(), cmodel.end(), 1.0);
    
    // Sample marker effects (BayesC)
    if (method==4) {
      for ( int isort = 0; isort < m; isort++) {
        int i = order[isort];
        sampleBetaCMt(i, nt, 
                    nmodels, models, cmodel, pi,
                    Ei, Bi,
                    ww, r, b, d, 
                    XXindices, XXvalues,
                    gen); 
      }
    }
    // Sample marker effects (BayesR)
    if (method==5) {
      for ( int isort = 0; isort < m; isort++) {
        int i = order[isort];
        sampleBetaR(i, nt, gamma,
                      dmarker, pimarker, pitrait,
                      E, B,
                      ww, r, b, d,
                      XXindices, XXvalues,
                      //LDindices, LDvalues,
                      gen);
      }
    }
    // Store values
    for (int t = 0; t < nt; t++) {
      for ( int i = 0; i < m; i++) {
        //if (d[t][i] == 1) {
        if (d[t][i] > 0) {
          if ((it > nburn) && (it % nthin == 0)) {
            dm[t][i] = dm[t][i] + 1.0;
            bm[t][i] = bm[t][i] + b[t][i];
          }
        }
      }
    }
    
    // Sample pi for Bayes C
    if(updatePi && method==4) {
      samplePi(cmodel, pi, gen);
      for (int k = 0; k<nmodels ; k++) {
        if(it>nburn) pis[k] = pis[k] + pi[k];
      }
    }
    
    // Sample pi for Bayes R
    if(updatePi && method==5) {
      samplePiMt(nt, pimarker, dmarker, gen);
      if(pimarker[0]<0.9) {
        pimarker[0]=0.9;
        pimarker[1]=1.0-pimarker[0];
      }
      samplePiR(nt, pitrait, d, gen);
      //Rcpp::Rcout << "pimarker[0]: " << pimarker[0] << std::endl;
      //Rcpp::Rcout << "pimarker[1]: " << pimarker[1] << std::endl;
      if(it>nburn) {
        for (int t = 0; t<nt ; t++) {
          for (size_t k = 0; k<gamma.size() ; k++) {
            pistrait[t][k] = pistrait[t][k] + pitrait[t][k];
          }
        }
        for (int k = 0; k<2 ; k++) {
          pismarker[k] = pismarker[k] + pimarker[k];
        }
      }
    }
    
    // Sample marker variance
    if(updateB && method==4) {
      sampleB(nt, m, nub, B, d, b, ssb_prior, gen);
      for (int t = 0; t < nt; t++) {
        vbs[t][it] = B(t,t);
      }
      for (int t1 = 0; t1 < nt; t1++) {
        for (int t2 = 0; t2 < nt; t2++) {
          if(it>nburn) cvbm[t1][t2] = cvbm[t1][t2] + B(t1,t2);
        }
      }
    }
    // Sample marker variance
    if(updateB && method==5) {
      sampleBR(nt, m, nub, gamma, B, d, b, ssb_prior, gen);
      for (int t = 0; t < nt; t++) {
        vbs[t][it] = B(t,t);
      }
      for (int t1 = 0; t1 < nt; t1++) {
        for (int t2 = 0; t2 < nt; t2++) {
          if(it>nburn) cvbm[t1][t2] = cvbm[t1][t2] + B(t1,t2);
        }
      }
    }
    //Update genetic variance
    computeG(nt, m, b, wy, r, n, G); 
    for (int t = 0; t < nt; t++) {
      vgs[t][it] = G(t,t);
    }
    for (int t1 = 0; t1 < nt; t1++  ) {
      for (int t2 = 0; t2 < nt; t2++) {
        if(it>nburn) cvgm[t1][t2] = cvgm[t1][t2] + G(t1,t2);
      } 
    } 
    
    // Sample residual variance 
    if(updateE) {
      sampleE(nt, m, nue, E, b, wy, r, sse_prior, yy, n, gen);
      for (int t = 0; t < nt; t++) {
        ves[t][it] = E(t,t);
      }
      for (int t1 = 0; t1 < nt; t1++) {
        for (int t2 = 0; t2 < nt; t2++) {
          if(it>nburn) cvem[t1][t2] = cvem[t1][t2] + E(t1,t2);
        } 
      } 
    }
    
    if ( (it > nburn) && (it % nthin == 0) ) {
      nsamples = nsamples + 1.0;
    }
    

  }
  // Summarize results
  std::vector<std::vector<std::vector<double>>> result;
  result.resize(20);
  
  result[0].resize(nt);
  result[1].resize(nt);
  result[2].resize(nt);
  result[3].resize(nt);
  result[4].resize(nt);
  result[5].resize(nt);
  result[6].resize(nt);
  result[7].resize(nt);
  result[8].resize(nt);
  result[9].resize(nt);
  result[10].resize(nt);
  result[11].resize(nt);
  result[12].resize(nt);
  result[13].resize(nt);
  result[14].resize(nt);
  result[15].resize(nt);
  result[16].resize(nt);
  result[17].resize(nt);
  result[18].resize(nt);
  result[19].resize(nt);
  
  for (int t=0; t < nt; t++) {
    result[0][t].resize(m);
    result[1][t].resize(m);
    result[2][t].resize(m);
    result[3][t].resize(m);
    result[4][t].resize(m);
    result[5][t].resize(m);
    result[6][t].resize(m);
    result[7][t].resize(nit+nburn);
    result[8][t].resize(nit+nburn);
    result[9][t].resize(nit+nburn);
    result[10][t].resize(nt);
    result[11][t].resize(nt);
    result[12][t].resize(nt);
    result[13][t].resize(nt);
    result[14][t].resize(nt);
    result[15][t].resize(nt);
    result[16][t].resize(nmodels);
    result[17][t].resize(nmodels);
    result[18][t].resize(4);
    result[19][t].resize(2);
  }
  
  for (int t=0; t < nt; t++) {
    for (int i=0; i < m; i++) {
      result[0][t][i] = bm[t][i]/nsamples;
      result[1][t][i] = dm[t][i]/nsamples;
      result[2][t][i] = wy[t][i];
      result[3][t][i] = r[t][i];
      result[4][t][i] = b[t][i];
      result[5][t][i] = d[t][i];
      result[6][t][i] = order[i];
    }
  }
  
  for (int t=0; t < nt; t++) {
    for (int i=0; i < nit+nburn; i++) {
      result[7][t][i] = vbs[t][i];
      result[8][t][i] = vgs[t][i];
      result[9][t][i] = ves[t][i];
    }
  }
  for (int t1=0; t1 < nt; t1++) {
    for (int t2=0; t2 < nt; t2++) {
      result[10][t1][t2] = cvbm[t1][t2]/nit;
      result[11][t1][t2] = cvgm[t1][t2]/nit;
      result[12][t1][t2] = cvem[t1][t2]/nit;
      result[13][t1][t2] = B(t1,t2);
      result[14][t1][t2] = G(t1,t2);
      result[15][t1][t2] = E(t1,t2);
    }
  }
  for (int t=0; t < nt; t++) {
    for (int i=0; i < nmodels; i++) {
      result[16][t][i] = pi[i];
      result[17][t][i] = pis[i]/nit;
    }
  }
  for (int t=0; t < nt; t++) {
    for (int i=0; i < 4; i++) {
      result[18][t][i] = pistrait[t][i]/nit;
    }
    for (int i=0; i < 2; i++) {
      result[19][t][i] = pismarker[i]/nit;
    }
  }
  return result;
}




// 
// // [[Rcpp::export]]
// std::vector<std::vector<std::vector<double>>>  mtsbayes(   std::vector<std::vector<double>> wy,
//                                                            std::vector<double> yy,
//                                                            std::vector<std::vector<double>> b,
//                                                            std::vector<std::vector<double>> LDvalues, 
//                                                            std::vector<std::vector<int>> LDindices, 
//                                                            arma::mat B,
//                                                            arma::mat E,
//                                                            std::vector<std::vector<double>> ssb_prior,
//                                                            std::vector<std::vector<double>> sse_prior,
//                                                            std::vector<std::vector<int>> models,
//                                                            std::vector<double> pi,
//                                                            double nub,
//                                                            double nue,
//                                                            bool updateB,
//                                                            bool updateE,
//                                                            bool updatePi,
//                                                            std::vector<int> n,
//                                                            int nit,
//                                                            int method) {
//   
//   // Define local variables
//   int nt = wy.size();
//   int m = wy[0].size();
//   int nmodels = models.size();
//   
//   double ssb, sse, dfb, u, logliksum, psum, detC, bxn, diff, cumprobc;
//   int mselect;
//   
//   std::vector<std::vector<double>> ww(nt, std::vector<double>(m, 0.0));
//   std::vector<std::vector<int>> d(nt, std::vector<int>(m, 0));
//   
//   std::vector<double> mu(nt), rhs(nt), conv(nt), dfe(nt);
//   std::vector<double> pmodel(nmodels), pcum(nmodels), loglik(nmodels), cmodel(nmodels);
//   std::vector<double> pis(nmodels);
//   
//   std::vector<std::vector<double>> bm(nt, std::vector<double>(m, 0.0));
//   std::vector<std::vector<double>> dm(nt, std::vector<double>(m, 0.0));
//   std::vector<std::vector<double>> ves(nt, std::vector<double>(nit, 0.0));
//   std::vector<std::vector<double>> vbs(nt, std::vector<double>(nit, 0.0));
//   std::vector<std::vector<double>> cvbm(nt, std::vector<double>(nt, 0.0));
//   std::vector<std::vector<double>> cvem(nt, std::vector<double>(nt, 0.0));
//   std::vector<std::vector<double>> mus(nt, std::vector<double>(nit, 0.0));
//   
//   std::vector<double> x2t(m);
//   std::vector<std::vector<double>> x2(nt, std::vector<double>(m, 0.0));
//   std::vector<std::vector<double>> r(nt, std::vector<double>(m, 0.0));
//   std::vector<int> order(m);
//   
//   
//   // Prior variance and degrees of freedom
//   // fix this dfe and n should be vectors of length nt
//   for (int t = 0; t < nt; t++) {
//     dfe[t] = n[t] + nue;
//   }
//   // Initialize variables
//   for (int i = 0; i < m; i++) {
//     for (int t = 0; t < nt; t++) {
//       // fix this dfe and n should be vectors of length nt
//       ww[t][i] = (double)n[t];
//       // fix this dfe and n should be vectors of length nt
//       //ww[t][i] = (double)n[t];;
//       //for (int j = 0; j < n; j++) {
//       //ww[t][i] = ww[t][i] + W[i][j]*W[i][j];
//       //}
//       r[t][i] = wy[t][i];
//       x2[t][i] = (wy[t][i]/ww[t][i])*(wy[t][i]/ww[t][i]);
//     }
//   }
//   
//   // Wy - W'Wb
//   for ( int i = 0; i < m; i++) {
//     for (int t = 0; t < nt; t++) {
//       if (b[t][i]!= 0.0) {
//         //      for (int j = 0; j < m; j++) {
//         bxn = double(n[t])*b[t][i];
//         for (size_t j = 0; j < LDindices[i].size(); j++) {
//           //r[t][LDindices[i][j]]=r[t][LDindices[i][j]] - LDvalues[i][j]*b[t][i];
//           r[t][LDindices[i][j]]=r[t][LDindices[i][j]] - LDvalues[i][j]*bxn;
//         }
//       }
//     }
//   }
//   
//   for (int i = 0; i < nmodels; i++) {
//     cmodel[i] = 1.0;
//     pis[i] = 0.0;
//   }
//   for ( int i = 0; i < m; i++) {
//     x2t[i] = 0.0;
//     for ( int t = 0; t < nt; t++) { 
//       x2t[i] = x2t[i] + x2[t][i];
//     }
//   }
//   
//   
//   // Establish order of markers as they are entered into the model
//   std::iota(order.begin(), order.end(), 0);
//   std::sort(  std::begin(order), 
//               std::end(order),
//               [&](int i1, int i2) { return x2t[i1] > x2t[i2]; } );
//   
//   
//   // Initialize (co)variance matrices
//   arma::mat C(nt,nt, fill::zeros);
//   arma::mat Bi = arma::inv(B);
//   arma::mat Ei = arma::inv(E);
//   
//   
//   // Start Gibbs sampler
//   std::random_device rd;
//   unsigned int seed;
//   seed = rd();
//   std::mt19937 gen(seed);
//   
//   for ( int it = 0; it < nit; it++) {
//     for ( int t = 0; t < nt; t++) {
//       conv[t] = 0.0;
//     }
//     
// 
//     // Sample marker effects (Mixed)
//     if (method==1) {
//     for ( int i = 0; i < m; i++) {
//         for ( int t = 0; t < nt; t++) {
//           rhs[t] = Ei(t,t)*r[t][i] + Ei(t,t)*ww[t][i]*b[t][i];
//         }
//         for ( int t = 0; t < nt; t++) { 
//           d[t][i] = 1;
//         }
//         arma::mat C(nt,nt, fill::zeros);
//         for ( int t1 = 0; t1 < nt; t1++) {
//           for ( int t2 = t1; t2 < nt; t2++) {
//             C(t1,t2) = Bi(t1,t2);
//             C(t2,t1) = Bi(t2,t1);
//           }
//           C(t1,t1) = C(t1,t1) + ww[t1][i]*Ei(t1,t1);       
//         }
//         arma::mat Ci = arma::inv(C);
//         arma::mat mub = mvrnormARMA(Ci);
//         for ( int t1 = 0; t1 < nt; t1++) {
//           for ( int t2 = 0; t2 < nt; t2++) {
//             mub(0,t1) = mub(0,t1) + Ci(t1,t2)*rhs[t2];
//           }
//         } 
//         
//         
//         for ( int t = 0; t < nt; t++) {
//           diff = (mub(0,t)-b[t][i])*double(n[t]);
//           for (size_t j = 0; j < LDindices[i].size(); j++) {
//             r[t][LDindices[i][j]]=r[t][LDindices[i][j]] - LDvalues[i][j]*diff;
//           }
//           conv[t] = conv[t] + diff*diff;
//           b[t][i] = mub(0,t);
//           dm[t][i] = dm[t][i] + 1.0;
//           bm[t][i] = bm[t][i] + b[t][i];     
//         }
//         
//       }
//     }
//     
//     
//     
//     // Sample marker effects (BayesC)
//     if (method==4) {
//       for ( int isort = 0; isort < m; isort++) {
//         int i = order[isort];
//         
//         for ( int t = 0; t < nt; t++) {
//           rhs[t] = Ei(t,t)*r[t][i] + Ei(t,t)*ww[t][i]*b[t][i];
//         }
//         for ( int k = 0; k < nmodels; k++) {
//           arma::mat C = Bi;
//           for ( int t1 = 0; t1 < nt; t1++) {
//             if(models[k][t1]==1) {
//               C(t1,t1) = C(t1,t1) + ww[t1][i]*Ei(t1,t1);
//             }
//           }
//           arma::mat Ci = arma::inv(C);
//           detC = arma::det(Ci);
//           loglik[k] = 0.5*std::log(detC) + std::log(pi[k]);
//           for ( int t1 = 0; t1 < nt; t1++) {
//             for ( int t2 = t1; t2 < nt; t2++) {
//               if(models[k][t1]==1 && models[k][t2]==1) {
//                 loglik[k] = loglik[k] + 0.5*rhs[t1]*rhs[t2]*Ci(t1,t2);
//                 if(t1!=t2) {
//                   loglik[k] = loglik[k] + 0.5*rhs[t2]*rhs[t1]*Ci(t2,t1);
//                 }
//               }
//             }
//           }
//         }
// 
// 
//         // Compute marker variance class probability
//         std::fill(pmodel.begin(), pmodel.end(), 0.0);
//         for (int k = 0; k<nmodels ; k++) {
//           logliksum = 0.0;
//           for (int l = 0; l<nmodels ; l++) {
//             logliksum += std::exp(loglik[l] - loglik[k]);
//           }
//           pmodel[k] = 1.0/logliksum;
//         }
//         
//         // // Compute marker variance class cumulative probability   
//         // std::fill(pcum.begin(), pcum.end(), 0.0);
//         // for (int k = 0; k<nmodels ; k++) {
//         //   psum = 0.0;
//         //   for (int l = 0; l<k+1 ; l++) {
//         //     psum += pmodel[l];
//         //   }
//         //   pcum[k] = psum;
//         // }
//         // 
//         // // Sample marker variance class indicator d  (note that lambda and d are 1-based)
//         // std::uniform_real_distribution<double> runif(0.0, 1.0);
//         // u = runif(gen);
//         // mselect = 0;
//         // for (int k = 0; k<nmodels ; k++) {
//         //   if (pcum[k]>u) {
//         //     mselect = k;
//         //     break;
//         //   }
//         // }
// 
//         // sample variance class indicator
//         std::uniform_real_distribution<double> runif(0.0, 1.0);
//         u = runif(gen);
//         cumprobc = 0.0;
//         for (int k = 0; k<nmodels ; k++) {
//           cumprobc += pmodel[k];
//           if(u < cumprobc){
//             mselect = k;
//             break;
//           }
//         }
//         
//         cmodel[mselect] = cmodel[mselect] + 1.0; 
//         for ( int t = 0; t < nt; t++) { 
//           d[t][i] = models[mselect][t];
//         }
//         arma::mat C = Bi;
//         for ( int t1 = 0; t1 < nt; t1++) {
//           if(models[mselect][t1]==1) {
//             C(t1,t1) = C(t1,t1) + ww[t1][i]*Ei(t1,t1);       
//           } 
//         }
//         arma::mat Ci = arma::inv(C);
//         arma::mat mub = mvrnormARMA(Ci);
//         for ( int t1 = 0; t1 < nt; t1++) {
//           mub(0,t1) = 0.0;
//           for ( int t2 = 0; t2 < nt; t2++) {
//             if(models[mselect][t2]==1) {
//               mub(0,t1) = mub(0,t1) + Ci(t1,t2)*rhs[t2];
//             }
//           }
//         } 
//         
//         for ( int t = 0; t < nt; t++) {
//           diff = (mub(0,t)-b[t][i])*double(n[t]);
//           for (size_t j = 0; j < LDindices[i].size(); j++) {
//             r[t][LDindices[i][j]]=r[t][LDindices[i][j]] - LDvalues[i][j]*diff;
//           }
//           conv[t] = conv[t] + diff*diff;
//           b[t][i] = mub(0,t);
//           if(d[t][i]==1) {
//             dm[t][i] = dm[t][i] + 1.0;
//             bm[t][i] = bm[t][i] + b[t][i];     
//           }
//         }
//         
//       }
//       // Sample pi for Bayes C
//       if(updatePi) {
//         for (int k = 0; k<nmodels ; k++) {
//           std::gamma_distribution<double> rgamma(cmodel[k],1.0);
//           double rg = rgamma(gen);
//           pi[k] = rg;
//         }
//         double psum = std::accumulate(pi.begin(), pi.end(), 0.0);
//         for (int k = 0; k<nmodels ; k++) {
//           pi[k] = pi[k]/psum;
//           pis[k] = pis[k] + pi[k];
//         }
//         std::fill(cmodel.begin(), cmodel.end(), 1.0);
//       }
//       // Sample pi for Bayes C
//       //if(updatePi) {
//       //  for (int k = 0; k<nmodels ; k++) {
//       //    std::gamma_distribution<double> rgamma(cmodel[k],1.0);
//       //    double rg = rgamma(gen);
//       //    pi[k] = rg/(double)m;
//       //    pis[k] = pis[k] + pi[k];
//       //  }
//       //  std::fill(cmodel.begin(), cmodel.end(), 1.0);
//       //}
//     }
//     
//     // Sample marker variance
//     if(updateB) {
//       arma::mat Sb(nt,nt, fill::zeros);
//       dfb = 0.0;
//       for (int t1 = 0; t1 < nt; t1++) {
//         for (int t2 = t1; t2 < nt; t2++) {
//           ssb = 0.0;
//           if (t1==t2) {
//             for (int i=0; i < m; i++) {
//               if( d[t1][i]==1 ) {
//                 ssb = ssb + b[t1][i]*b[t1][i];  
//                 dfb = dfb + 1.0;
//               }
//             }
//             Sb(t1,t1) = ssb + ssb_prior[t1][t1];
//           } 
//           if (t1!=t2) {
//             for (int i=0; i < m; i++) {
//               if( d[t1][i]==1 && d[t2][i]==1 ) {
//                 ssb = ssb + b[t1][i]*b[t2][i];  
//               }
//             }
//             Sb(t1,t2) = ssb + ssb_prior[t1][t2];
//             Sb(t2,t1) = ssb + ssb_prior[t2][t1];
//           } 
//         }
//       }
//       int dfSb = dfb/nt + nub;
//       arma::mat B = riwishart(dfSb, Sb);
//       for (int t = 0; t < nt; t++) {
//         vbs[t][it] = B(t,t);
//       }
//       for (int t1 = 0; t1 < nt; t1++) {
//         for (int t2 = 0; t2 < nt; t2++) {
//           cvbm[t1][t2] = cvbm[t1][t2] + B(t1,t2);
//         } 
//       } 
//       arma::mat Bi = arma::inv(B);
//     }
//     
//     // // Sample marker variance
//     // if(updateB) {
//     //   arma::mat Sb(nt,nt, fill::zeros);
//     //   dfb = 0.0;
//     //   for (int t1 = 0; t1 < nt; t1++) {
//     //     for (int t2 = t1; t2 < nt; t2++) {
//     //       ssb = 0.0;
//     //       if (t1==t2) {
//     //         for (int i=0; i < m; i++) {
//     //           if( (d[t1][i]==1) ) {
//     //             ssb = ssb + b[t1][i]*b[t2][i];  
//     //             dfb = dfb + 1.0;
//     //           }
//     //         }
//     //         Sb(t1,t2) = ssb + ssb_prior[t1][t2];
//     //       } 
//     //       if (t1!=t2) {
//     //         for (int i=0; i < m; i++) {
//     //           if( (d[t1][i]==1) && (d[t2][i]==1) ) {
//     //             ssb = ssb + b[t1][i]*b[t2][i];  
//     //           }
//     //         }
//     //         Sb(t1,t2) = ssb + ssb_prior[t1][t2];
//     //         Sb(t2,t1) = ssb + ssb_prior[t2][t1];
//     //       } 
//     //     }
//     //   }
//     //   int dfSb = dfb/nt + nub;
//     //   arma::mat B = riwishart(dfSb, Sb);
//     //   for (int t = 0; t < nt; t++) {
//     //     vbs[t][it] = B(t,t);
//     //   }
//     //   for (int t1 = 0; t1 < nt; t1++) {
//     //     for (int t2 = 0; t2 < nt; t2++) {
//     //       cvbm[t1][t2] = cvbm[t1][t2] + B(t1,t2)/(sqrt(B(t1,t1))*sqrt(B(t2,t2)));
//     //     } 
//     //   } 
//     //   arma::mat Bi = arma::inv(B);
//     // }
//     
//     // Sample residual variance
//     if(updateE) {
//       arma::mat Se(nt,nt, fill::zeros);
//       for ( int t1 = 0; t1 < nt; t1++) {
//         sse = 0.0;
//         for ( int i = 0; i < m; i++) {
//           sse = sse + b[t1][i] * (r[t1][i] + wy[t1][i]);
//         }
//         sse = yy[t1] - sse;
//       }
//       //int dfSe = dfe + nue;
//       int dfSe = dfe[0] + nue;
//       arma::mat E = riwishart(dfSe, Se);
//       arma::mat Ei = arma::inv(E);
//       for (int t = 0; t < nt; t++) {
//         ves[t][it] = E(t,t);
//       }
//       for (int t1 = 0; t1 < nt; t1++) {
//         for (int t2 = 0; t2 < nt; t2++) {
//           cvem[t1][t2] = cvem[t1][t2] + E(t1,t2)/(sqrt(E(t1,t1))*sqrt(E(t2,t2)));
//         } 
//       } 
//     }
// 
//         
//     
//   }
//   
//   // Summarize results
//   std::vector<std::vector<std::vector<double>>> result;
//   result.resize(15);
//   
//   result[0].resize(nt);
//   result[1].resize(nt);
//   result[2].resize(nt);
//   result[3].resize(nt);
//   result[4].resize(nt);
//   result[5].resize(nt);
//   result[6].resize(nt);
//   result[7].resize(nt);
//   result[8].resize(nt);
//   result[9].resize(nt);
//   result[10].resize(nt);
//   result[11].resize(nt);
//   result[12].resize(nt);
//   result[13].resize(nt);
//   result[14].resize(nt);
//   
//   for (int t=0; t < nt; t++) {
//     result[0][t].resize(m);
//     result[1][t].resize(m);
//     result[2][t].resize(nit);
//     result[3][t].resize(nit);
//     result[4][t].resize(nit);
//     result[5][t].resize(nt);
//     result[6][t].resize(nt);
//     result[7][t].resize(m);
//     result[8][t].resize(m);
//     result[9][t].resize(m);
//     result[10][t].resize(nt);
//     result[11][t].resize(nt);
//     result[12][t].resize(nmodels);
//     result[13][t].resize(nmodels);
//     result[14][t].resize(m);
//   }
//   
//   for (int t=0; t < nt; t++) {
//     for (int i=0; i < m; i++) {
//       result[0][t][i] = bm[t][i]/nit;
//       result[1][t][i] = dm[t][i]/nit;
//       result[9][t][i] = b[t][i];
//       result[14][t][i] = order[i];
//     }
//     for (int i=0; i < m; i++) {
//       result[7][t][i] = wy[t][i];
//       result[8][t][i] = r[t][i];
//     }
//   }
//   
//   for (int t=0; t < nt; t++) {
//     for (int i=0; i < nit; i++) {
//       result[2][t][i] = mus[t][i];
//       result[3][t][i] = vbs[t][i];
//       result[4][t][i] = ves[t][i];
//     }
//   }
//   for (int t1=0; t1 < nt; t1++) {
//     for (int t2=0; t2 < nt; t2++) {
//       result[5][t1][t2] = cvbm[t1][t2]/nit;
//       result[6][t1][t2] = cvem[t1][t2]/nit;
//       result[10][t1][t2] = B(t1,t2);
//       result[11][t1][t2] = E(t1,t2);
//     }
//   }
//   for (int t=0; t < nt; t++) {
//     for (int i=0; i < nmodels; i++) {
//       result[12][t][i] = pi[i];
//       result[13][t][i] = pis[i]/nit;
//     }
//   }  
//   return result;
// }


// // [[Rcpp::export]]
// std::vector<std::vector<std::vector<double>>>  mtblr(   std::vector<std::vector<double>> wy,
//                                                         std::vector<std::vector<double>> ww,
//                                                         std::vector<double> yy,
//                                                         std::vector<std::vector<double>> b,
//                                                         std::vector<std::vector<std::vector<double>>> XXvalues,
//                                                         std::vector<std::vector<int>> XXindices,
//                                                         arma::mat B,
//                                                         arma::mat E,
//                                                         std::vector<std::vector<double>> ssb_prior,
//                                                         std::vector<std::vector<double>> sse_prior,
//                                                         std::vector<std::vector<int>> models,
//                                                         std::vector<double> pi,
//                                                         double nub,
//                                                         double nue,
//                                                         bool updateB,
//                                                         bool updateE,
//                                                         bool updatePi,
//                                                         std::vector<int> n,
//                                                         int nit,
//                                                         int nburn,
//                                                         int nthin,
//                                                         int seed,
//                                                         int method) {
//   
//   // Define local variables
//   int nt = wy.size();
//   int m = wy[0].size();
//   int nmodels = models.size();
//   double nsamples=0.0;
//   
//   double ssb, sse, ssg, dfb, u, logliksum, detC, diff, cumprobc;
//   int mselect;
//   
//   std::vector<std::vector<int>> d(nt, std::vector<int>(m, 0));
//   
//   std::vector<double> mu(nt), rhs(nt), conv(nt);
//   std::vector<double> pmodel(nmodels), pcum(nmodels), loglik(nmodels), cmodel(nmodels);
//   std::vector<double> pis(nmodels);
//   
//   std::vector<std::vector<double>> bm(nt, std::vector<double>(m, 0.0));
//   std::vector<std::vector<double>> dm(nt, std::vector<double>(m, 0.0));
//   std::vector<std::vector<double>> ves(nt, std::vector<double>(nit+nburn, 0.0));
//   std::vector<std::vector<double>> vbs(nt, std::vector<double>(nit+nburn, 0.0));
//   std::vector<std::vector<double>> vgs(nt, std::vector<double>(nit+nburn, 0.0));
//   std::vector<std::vector<double>> cvbm(nt, std::vector<double>(nt, 0.0));
//   std::vector<std::vector<double>> cvem(nt, std::vector<double>(nt, 0.0));
//   std::vector<std::vector<double>> cvgm(nt, std::vector<double>(nt, 0.0));
//   std::vector<std::vector<double>> mus(nt, std::vector<double>(nit+nburn, 0.0));
//   
//   std::vector<double> x2t(m);
//   std::vector<std::vector<double>> x2(nt, std::vector<double>(m, 0.0));
//   std::vector<std::vector<double>> r(nt, std::vector<double>(m, 0.0));
//   std::vector<int> order(m);
//   
//   
//   // Initialize variables
//   for (int i = 0; i < m; i++) {
//     for (int t = 0; t < nt; t++) {
//       r[t][i] = wy[t][i];
//       x2[t][i] = (wy[t][i]/ww[t][i])*(wy[t][i]/ww[t][i]);
//     }
//   }
//   
//   // Wy - W'Wb
//   for ( int i = 0; i < m; i++) {
//     for (int t = 0; t < nt; t++) {
//       if (b[t][i]!= 0.0) {
//         for (size_t j = 0; j < XXindices[i].size(); j++) {
//           r[t][XXindices[i][j]]=r[t][XXindices[i][j]] - XXvalues[t][i][j]*b[t][i];
//         }
//       }
//     }
//   }
//   
//   // for ( int i = 0; i < m; i++) {
//   //   x2t[i] = 0.0;
//   //   for ( int t = 0; t < nt; t++) {
//   //     x2t[i] = x2t[i] + x2[t][i];
//   //   }
//   // }
//   
//   for ( int i = 0; i < m; i++) {
//     x2t[i] = 0.0;
//     for ( int t = 0; t < nt; t++) { 
//       if(x2[t][i]>x2t[i]) {
//         x2t[i] = x2[t][i];
//       }
//     }
//   }
//   
//   std::fill(cmodel.begin(), cmodel.end(), 1.0);
//   std::fill(pis.begin(), pis.end(), 0.0);
//   
//   
//   // Establish order of markers as they are entered into the model
//   std::iota(order.begin(), order.end(), 0);
//   std::sort(  std::begin(order),
//               std::end(order),
//               [&](int i1, int i2) { return x2t[i1] > x2t[i2]; } );
//   
//   
//   // Initialize (co)variance matrices
//   arma::mat C(nt,nt, fill::zeros);
//   arma::mat Bi = arma::inv(B);
//   arma::mat Ei = arma::inv(E);
//   arma::mat G(nt,nt, fill::zeros);
//   
//   arma::rowvec probs(nmodels, fill::zeros);
//   
//   // Start Gibbs sampler
//   std::random_device rd;
//   //unsigned int seed;
//   //seed = rd();
//   std::mt19937 gen(seed);
//   
//   for ( int it = 0; it < nit+nburn; it++) {
//     for ( int t = 0; t < nt; t++) {
//       conv[t] = 0.0;
//     }
//     
//     if ( (it > nburn) && (it % nthin == 0) ) {
//       nsamples = nsamples + 1.0;
//     }
//     
//     // Sample marker effects (BayesC)
//     if (method==4) {
//       for ( int isort = 0; isort < m; isort++) {
//         int i = order[isort];
//         
//         // compute rhs
//         for ( int t = 0; t < nt; t++) {
//           rhs[t] = Ei(t,t)*r[t][i] + Ei(t,t)*ww[t][i]*b[t][i];
//         }
//         
//         // Compute
//         for ( int k = 0; k < nmodels; k++) {
//           arma::mat C = Bi;
//           for ( int t1 = 0; t1 < nt; t1++) {
//             if(models[k][t1]==1) {
//               C(t1,t1) = C(t1,t1) + ww[t1][i]*Ei(t1,t1);
//             }
//           }
//           arma::mat Ci = arma::inv(C);
//           detC = arma::det(Ci);
//           loglik[k] = 0.5*std::log(detC) + std::log(pi[k]);
//           for ( int t1 = 0; t1 < nt; t1++) {
//             for ( int t2 = t1; t2 < nt; t2++) {
//               if(models[k][t1]==1 && models[k][t2]==1) {
//                 loglik[k] = loglik[k] + 0.5*rhs[t1]*rhs[t2]*Ci(t1,t2);
//                 if(t1!=t2) {
//                   loglik[k] = loglik[k] + 0.5*rhs[t2]*rhs[t1]*Ci(t2,t1);
//                 }
//               }
//             }
//           }
//         }
//         
//         // Compute marker variance class probability
//         std::fill(pmodel.begin(), pmodel.end(), 0.0);
//         for (int k = 0; k<nmodels ; k++) {
//           logliksum = 0.0;
//           for (int l = 0; l<nmodels ; l++) {
//             logliksum += std::exp(loglik[l] - loglik[k]);
//           }
//           pmodel[k] = 1.0/logliksum;
//         }
//         
//         // Sample variance class indicator
//         std::uniform_real_distribution<double> runif(0.0, 1.0);
//         u = runif(gen);
//         mselect=0;
//         cumprobc = 0.0;
//         
//         for (int k = 0; k<nmodels ; k++) {
//           cumprobc += pmodel[k];
//           if(u < cumprobc){
//             mselect = k;
//             break;
//           }
//         }
//         cmodel[mselect] = cmodel[mselect] + 1.0;
//         for ( int t = 0; t < nt; t++) {
//           d[t][i] = models[mselect][t];
//         }
//         
//         // Sample marker effect conditional on variance class indicator
//         arma::mat C = Bi;
//         for ( int t1 = 0; t1 < nt; t1++) {
//           if(models[mselect][t1]==1) {
//             C(t1,t1) = C(t1,t1) + ww[t1][i]*Ei(t1,t1);       
//           } 
//         }
//         arma::mat Ci = arma::inv(C);
//         arma::mat mub = mvrnormARMA(Ci);
//         for ( int t1 = 0; t1 < nt; t1++) {
//           if(models[mselect][t1]==0) {
//             mub(0,t1) = 0.0;
//           }
//           if(models[mselect][t1]==1) {
//             mub(0,t1) = mub(0,t1) + Ci(t1,t1)*rhs[t1];
//             for ( int t2 = 0; t2 < nt; t2++) {
//               if(t1!=t2 && models[mselect][t2]==1) {
//                 mub(0,t1) = mub(0,t1) + Ci(t1,t2)*rhs[t2];
//               }
//             }
//           }
//         }
//         
//         // Adjust residuals based on sample marker effects
//         for ( int t = 0; t < nt; t++) {
//           diff = (mub(0,t)-b[t][i]);
//           if(diff!=0.0) {
//             for (size_t j = 0; j < XXindices[i].size(); j++) {
//               r[t][XXindices[i][j]] = r[t][XXindices[i][j]] - XXvalues[t][i][j]*diff;
//             }
//           }
//           b[t][i] = mub(0,t);
//           if(d[t][i]==1) {
//             if ( (it > nburn) && (it % nthin == 0) ) {
//               dm[t][i] = dm[t][i] + 1.0;
//               bm[t][i] = bm[t][i] + b[t][i];
//             }
//           }
//         }
//       }
//     }
//     
//     // Sample pi for Bayes C
//     if(updatePi) {
//       updatePiAndCModel(cmodel, pi, gen);
//       // for (int k = 0; k<nmodels ; k++) {
//       //   std::gamma_distribution<double> rgamma(cmodel[k],1.0);
//       //   double rg = rgamma(gen);
//       //   pi[k] = rg;
//       // }
//       // double psum = std::accumulate(pi.begin(), pi.end(), 0.0);
//       for (int k = 0; k<nmodels ; k++) {
//         //pi[k] = pi[k]/psum;
//         if(it>nburn) pis[k] = pis[k] + pi[k];
//       }
//       std::fill(cmodel.begin(), cmodel.end(), 1.0);
//     }
//     
//     // Sample marker variance
//     if(updateB) {
//       arma::mat Sb(nt,nt, fill::zeros);
//       arma::mat corb(nt,nt, fill::zeros);
//       arma::mat dfB(nt,nt, fill::zeros);
//       for (int t1 = 0; t1 < nt; t1++) {
//         ssb = 0.0;
//         dfb = 0.0;
//         for (int i=0; i < m; i++) {
//           if( d[t1][i]==1 ) {
//             ssb = ssb + b[t1][i]*b[t1][i];
//             dfb = dfb + 1.0;
//           }
//         }
//         Sb(t1,t1) = ssb;
//         dfB(t1,t1) = dfb;
//         for (int t2 = t1; t2 < nt; t2++) {
//           ssb = 0.0;
//           dfb = 0.0;
//           if (t1!=t2) {
//             for (int i=0; i < m; i++) {
//               if( d[t1][i]==1 && d[t2][i]==1 ) {
//                 ssb = ssb + b[t1][i]*b[t2][i];
//                 dfb = dfb + 1.0;
//               }
//             }
//             dfB(t1,t2) = dfb;
//             dfB(t2,t1) = dfb;
//             Sb(t1,t2) = ssb;
//             Sb(t2,t1) = ssb;
//           }
//         }
//       }
//       std::vector<double> stdv(nt);
//       for (int t = 0; t < nt; t++) {
//         stdv[t] = sqrt(Sb(t,t));
//       }
//       for (int t1 = 0; t1 < nt; t1++) {
//         for (int t2 = 0; t2 < nt; t2++) {
//           if (t1==t2) {
//             corb(t1,t2) = 1.0;
//           }
//           if (t1!=t2 && stdv[t1]!=0.0 && stdv[t2]!=0.0 ) {
//             corb(t1,t2) = Sb(t1,t2)/(stdv[t1] * stdv[t2]);
//             corb(t2,t1) = corb(t1,t2);
//           }
//         }
//       }
//       
//       arma::mat B(nt,nt, fill::zeros);
//       for (int t = 0; t < nt; t++) {
//         std::chi_squared_distribution<double> rchisq(dfB(t,t) + nub);
//         double chi2 = rchisq(gen);
//         B(t,t) = (Sb(t,t) + nub*ssb_prior[t][t])/chi2;
//       }
//       for (int t1 = 0; t1 < nt; t1++) {
//         for (int t2 = 0; t2 < nt; t2++) {
//           if (t1!=t2) {
//             B(t1,t2) = corb(t1,t2)*sqrt(B(t1,t1))*sqrt(B(t2,t2));
//           }
//         }
//       }
//       for (int t1 = 0; t1 < nt; t1++) {
//         B(t1,t1) = B(t1,t1) + 0.0001;
//       }
//       for (int t = 0; t < nt; t++) {
//         vbs[t][it] = B(t,t);
//       }
//       for (int t1 = 0; t1 < nt; t1++) {
//         for (int t2 = 0; t2 < nt; t2++) {
//           if(it>nburn) cvbm[t1][t2] = cvbm[t1][t2] + B(t1,t2);
//         }
//       }
//       bool issym = B.is_symmetric();
//       if (!issym){
//         B = 0.5* (B + B.t());
//       }
//       arma::mat Bi(nt,nt, fill::zeros);
//       bool success;
//       success = inv_sympd(Bi,B,inv_opts::allow_approx);
//       if(!success) {
//         std::cerr << "Error: Condition is false." << std::endl;
//       }
//       //arma::mat Bi = arma::inv(B);
//     }
//     //Update genetic variance
//     arma::mat G(nt,nt, fill::zeros);
//     for ( int t1 = 0; t1 < nt; t1++) {
//       for ( int t2 = t1; t2 < nt; t2++) {
//         if (t1==t2) {
//           ssg = 0.0;
//           for ( int i = 0; i < m; i++) {
//             ssg = ssg + b[t1][i]*(wy[t1][i]-r[t1][i]);
//           }
//           G(t1,t2) = ssg/(std::sqrt((double)n[t1])*std::sqrt((double)n[t2]));
//         }
//         if (t1!=t2) {
//           ssg = 0.0;
//           for ( int i = 0; i < m; i++) {
//             ssg = ssg + b[t1][i]*(wy[t1][i]-r[t1][i]);
//             ssg = ssg + b[t2][i]*(wy[t2][i]-r[t2][i]);
//           }
//           ssg = ssg/2.0;
//           G(t1,t2) = ssg/(std::sqrt((double)n[t1])*std::sqrt((double)n[t2]));
//           G(t2,t1) = G(t1,t2);
//         }
//       }
//     }
//     for (int t = 0; t < nt; t++) {
//       vgs[t][it] = G(t,t);
//     }
//     for (int t1 = 0; t1 < nt; t1++) {
//       for (int t2 = 0; t2 < nt; t2++) {
//         if(it>nburn) cvgm[t1][t2] = cvgm[t1][t2] + G(t1,t2);
//       } 
//     } 
//     
//     // Sample residual variance 
//     // (current version assume uncorrelated residuals among traits traits)
//     if(updateE) {
//       arma::mat Se(nt,nt, fill::zeros);
//       for ( int t1 = 0; t1 < nt; t1++) {
//         for ( int t2 = t1; t2 < nt; t2++) {
//           if (t1==t2) {
//             sse = 0.0;
//             for ( int i = 0; i < m; i++) {
//               sse = sse + b[t1][i] * (r[t1][i] + wy[t1][i]);
//             }
//             sse = yy[t1] - sse;
//             Se(t1,t1) = sse + nue*sse_prior[t1][t1];
//           }
//         }
//       }
//       arma::mat E(nt,nt, fill::zeros);
//       for (int t = 0; t < nt; t++) {
//         std::chi_squared_distribution<double> rchisq(n[t] + nue);
//         double chi2 = rchisq(gen);
//         E(t,t) = (Se(t,t))/chi2;
//       }
//       //int dfSe = n[0] + nue;
//       //arma::mat E = riwishart(dfSe, Se);
//       arma::mat Ei = arma::inv(E);
//       for (int t = 0; t < nt; t++) {
//         ves[t][it] = E(t,t);
//       }
//       for (int t1 = 0; t1 < nt; t1++) {
//         for (int t2 = 0; t2 < nt; t2++) {
//           if(it>nburn) cvem[t1][t2] = cvem[t1][t2] + E(t1,t2);
//         } 
//       } 
//     }
//     
//     
//   }
//   // Summarize results
//   std::vector<std::vector<std::vector<double>>> result;
//   result.resize(18);
//   
//   result[0].resize(nt);
//   result[1].resize(nt);
//   result[2].resize(nt);
//   result[3].resize(nt);
//   result[4].resize(nt);
//   result[5].resize(nt);
//   result[6].resize(nt);
//   result[7].resize(nt);
//   result[8].resize(nt);
//   result[9].resize(nt);
//   result[10].resize(nt);
//   result[11].resize(nt);
//   result[12].resize(nt);
//   result[13].resize(nt);
//   result[14].resize(nt);
//   result[15].resize(nt);
//   result[16].resize(nt);
//   result[17].resize(nt);
//   
//   for (int t=0; t < nt; t++) {
//     result[0][t].resize(m);
//     result[1][t].resize(m);
//     result[2][t].resize(m);
//     result[3][t].resize(m);
//     result[4][t].resize(m);
//     result[5][t].resize(m);
//     result[6][t].resize(m);
//     result[7][t].resize(nit+nburn);
//     result[8][t].resize(nit+nburn);
//     result[9][t].resize(nit+nburn);
//     result[10][t].resize(nt);
//     result[11][t].resize(nt);
//     result[12][t].resize(nt);
//     result[13][t].resize(nt);
//     result[14][t].resize(nt);
//     result[15][t].resize(nt);
//     result[16][t].resize(nmodels);
//     result[17][t].resize(nmodels);
//   }
//   
//   for (int t=0; t < nt; t++) {
//     for (int i=0; i < m; i++) {
//       result[0][t][i] = bm[t][i]/nsamples;
//       result[1][t][i] = dm[t][i]/nsamples;
//       result[2][t][i] = wy[t][i];
//       result[3][t][i] = r[t][i];
//       result[4][t][i] = b[t][i];
//       result[5][t][i] = d[t][i];
//       result[6][t][i] = order[i];
//     }
//   }
//   
//   for (int t=0; t < nt; t++) {
//     for (int i=0; i < nit+nburn; i++) {
//       result[7][t][i] = vbs[t][i];
//       result[8][t][i] = vgs[t][i];
//       result[9][t][i] = ves[t][i];
//     }
//   }
//   for (int t1=0; t1 < nt; t1++) {
//     for (int t2=0; t2 < nt; t2++) {
//       result[10][t1][t2] = cvbm[t1][t2]/nit;
//       result[11][t1][t2] = cvgm[t1][t2]/nit;
//       result[12][t1][t2] = cvem[t1][t2]/nit;
//       result[13][t1][t2] = B(t1,t2);
//       result[14][t1][t2] = G(t1,t2);
//       result[15][t1][t2] = E(t1,t2);
//     }
//   }
//   for (int t=0; t < nt; t++) {
//     for (int i=0; i < nmodels; i++) {
//       result[16][t][i] = pi[i];
//       result[17][t][i] = pis[i]/nit;
//     }
//   }
//   return result;
// }
// 
