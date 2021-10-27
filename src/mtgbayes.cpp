// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <RcppDist.h>
using namespace Rcpp;
using namespace arma;


#include <vector>
#include <numeric>
#include <random>


// [[Rcpp::export]]
arma::mat mmult(arma::mat A, arma::mat B) {
  return A * B;
}




// [[Rcpp::export]]
arma::mat mvrnorm(arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(1, ncols);
  return (Y * arma::chol(sigma));
}

// arma::mat mvrnormArma(arma::vec mu, arma::mat sigma) {
//   int ncols = sigma.n_cols;
//   arma::mat Y = arma::randn(1, ncols);
//   return mu + Y * arma::chol(sigma);
// }


// // [[Rcpp::export]]
// arma::mat rinvwish(int v, arma::mat S){
//   RNGScope scope;
//   int p = S.n_rows;
//   arma::mat L = arma::chol(arma::inv_sympd(S), "lower");
//   arma::mat W(p, p,fill::zeros);
//   arma::mat A(p,p, fill::zeros);
//   for(int i = 0; i < p; i++){
//     int df = v - (i + 1) + 1; 
//     A(i,i) = std::sqrt(R::rchisq(df)); 
//   }
//   for(int row = 1; row < p; row++){
//     for(int col = 0; col < row; col++){
//       A(row, col) = R::rnorm(0,1);
//     }
//   }
//   //arma::mat LA_inv = arma::pinv(arma::trimatl(arma::trimatl(L) * arma::trimatl(A)), 0.001);
//   arma::mat LA_inv = arma::inv(arma::trimatl(arma::trimatl(L) * arma::trimatl(A)));
//   W = LA_inv.t() * LA_inv;
//   return(W);
// }


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
  
  //Rcout << "Number of traits: " << nt << "\n";
  //Rcout << "Number of individuals: " << m << "\n";
  //Rcout << "Number of markers: " << m << "\n";
  
  
  double ssb, sse, ssg, dfb, dfe, u, logliksum, psum, detC, diff;
  int mselect;
  
  std::vector<std::vector<double>> e(nt, std::vector<double>(n, 0.0));
  std::vector<std::vector<double>> g(nt, std::vector<double>(n, 0.0));
  std::vector<std::vector<double>> ww(nt, std::vector<double>(m, 0.0));
  std::vector<std::vector<int>> d(nt, std::vector<int>(m, 0));
  
  std::vector<double> mu(nt), rhs(nt), conv(nt);
  //std::vector<double> sse_prior(nt);
  //std::vector<std::vector<double>> ssb_prior(nt, std::vector<double>(nt, 0.0)); 
  
  std::vector<double> pmodel(nmodels), pcum(nmodels), loglik(nmodels), cmodel(nmodels);
  std::vector<double> pi_post_mean(nmodels);
  
  std::vector<std::vector<double>> b_post_mean(nt, std::vector<double>(m, 0.0));
  std::vector<std::vector<double>> d_post_mean(nt, std::vector<double>(m, 0.0));
  std::vector<std::vector<double>> ve_post(nt, std::vector<double>(nit, 0.0));
  std::vector<std::vector<double>> vb_post(nt, std::vector<double>(nit, 0.0));
  std::vector<std::vector<double>> covb_post_mean(nt, std::vector<double>(nt, 0.0));
  std::vector<std::vector<double>> cove_post_mean(nt, std::vector<double>(nt, 0.0));
  std::vector<std::vector<double>> mu_post(nt, std::vector<double>(nit, 0.0));
  
  std::vector<double> x2t(m);
  std::vector<std::vector<double>> x2(nt, std::vector<double>(m, 0.0));
  std::vector<std::vector<double>> wy(nt, std::vector<double>(m, 0.0));
  std::vector<std::vector<double>> r(nt, std::vector<double>(m, 0.0));
  std::vector<int> order(m);
  
  
  // Prior variance and degrees of freedom
  dfe = n;
  //for (int t1 = 0; t1 < nt; t1++) {
  //  for (int t2 = 0; t2 < nt; t2++) {
  //    ssb_prior[t1][t2] = (nub-2.0)/nub * B(t1,t2);
  //  }
  //  sse_prior[t1] = nue*E(t1,t1);
  //}

  
  // Initialize variables
  for (int i = 0; i < m; i++) {
    for (int t = 0; t < nt; t++) {
      ww[t][i] = 0.0;
      for (int j = 0; j < n; j++) {
        ww[t][i] = ww[t][i] + W[i][j]*W[i][j];
      }
    }
  }
  
  for (int i = 0; i < nmodels; i++) {
    cmodel[i] = 1.0;
    pi_post_mean[i] = 0.0;
  }
    
  // Mean adjust y
  for ( int t = 0; t < nt; t++) {
    mu[t] = std::accumulate(y[t].begin(), y[t].end(), 0.0)/n;
    for ( int i = 0; i < n; i++) {
      e[t][i] = y[t][i] - mu[t];
    }                
  }
  
  
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
  std::iota(order.begin(), order.end(), 0);
  std::sort(  std::begin(order), 
              std::end(order),
              [&](int i1, int i2) { return x2t[i1] > x2t[i2]; } );

    
  // Initialize (co)variance matrices
  arma::mat C(nt,nt, fill::zeros);

  arma::mat Bi = arma::inv(B);
  arma::mat Ei = arma::inv(E);
  

  // Start Gibbs sampler
  
  //Rcout << "  " << "\n";
  //Rcout << "Starting Gibbs sampler" << "\n";
  //Rcout << "  " << "\n";
  
  std::random_device rd;
  unsigned int seed;
  seed = rd();
  std::mt19937 gen(seed);
  
  for ( int it = 0; it < nit; it++) {
    for ( int t = 0; t < nt; t++) {
      conv[t] = 0.0;
    }

    //arma::mat Ei = arma::inv(E);
    //arma::mat Bi = arma::inv(B);
    
    // Sample marker effects (Mixed)
    if (method<2) {
      
      //for ( int i = 0; i < m; i++) {
      for ( int i0 = 0; i0 < m; i0++) {
        int i = order[i0];
        
        // Compute rhs                     
        for ( int t = 0; t < nt; t++) {
          rhs[t] = 0.0;
          for ( int j = 0; j < n; j++) {
            //rhs[t] = rhs[t] + Ei(t,t)*W[i][j]*e[t][j]; 
            rhs[t] = rhs[t] + W[i][j]*e[t][j]; 
          }
          rhs[t] = Ei(t,t)*rhs[t] + Ei(t,t)*ww[t][i]*b[t][i];
        }
        
        // Sample beta
        for ( int t = 0; t < nt; t++) { 
          d[t][i] = 1;
        }
        arma::mat C = Bi;
        for ( int t = 0; t < nt; t++) {
          C(t,t) = C(t,t) + ww[t][i]*Ei(t,t);       
        }
        arma::mat Ci = arma::inv(C);
        arma::mat mub = mvrnorm(Ci);

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
          d_post_mean[t][i] = d_post_mean[t][i] + 1.0;
          b_post_mean[t][i] = b_post_mean[t][i] + b[t][i];     
        }
        
      }
    }
    
    // Sample marker effects (BayesC)
    if (method==4) {
      for ( int i0 = 0; i0 < m; i0++) {
        int i = order[i0];
        // Compute rhs                     
        for ( int t = 0; t < nt; t++) {
          rhs[t] = 0.0;
          for ( int j = 0; j < n; j++) {
            rhs[t] = rhs[t] + Ei(t,t)*W[i][j]*e[t][j]; 
          }
          rhs[t] = rhs[t] + Ei(t,t)*ww[t][i]*b[t][i];
        }
        
        
        // Compute lhs
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
        
        // Compute variance class probability
        std::fill(pmodel.begin(), pmodel.end(), 0.0);
        for (int k = 0; k<nmodels ; k++) {
          logliksum = 0.0;
          for (int l = 0; l<nmodels ; l++) {
            logliksum += std::exp(loglik[l] - loglik[k]);
          }
          pmodel[k] = 1.0/logliksum;
        }
        
        // Compute variance class cumulative probability   
        std::fill(pcum.begin(), pcum.end(), 0.0);
        for (int k = 0; k<nmodels ; k++) {
          psum = 0.0;
          for (int l = 0; l<k+1 ; l++) {
            psum += pmodel[l];
          }
          pcum[k] = psum;
        }
        

        // Sample marker variance class indicator - lambda and d are 1-based
        std::uniform_real_distribution<double> runif(0.0, 1.0);
        u = runif(gen);
        mselect = 0;
        for (int k = 0; k<nmodels ; k++) {
          if (pcum[k]>u) {
            mselect = k;
            break;
          }
        }
        cmodel[mselect] = cmodel[mselect] + 1.0; 
        
        // Sample beta
        for ( int t = 0; t < nt; t++) { 
          d[t][i] = models[mselect][t];
        }
        arma::mat C = Bi;
        for ( int t = 0; t < nt; t++) {
          if(models[mselect][t]==1) {
            C(t,t) = C(t,t) + ww[t][i]*Ei(t,t);       
          } 
        }
        arma::mat Ci = arma::inv(C);
        arma::mat mub = mvrnorm(Ci);

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

        for ( int t = 0; t < nt; t++) {
          
          diff = mub(0,t)-b[t][i];
          for (int j = 0; j < n; j++) {
            e[t][j]=e[t][j] - W[i][j]*(diff);
          }
          conv[t] = conv[t] + diff*diff;
          b[t][i] = mub(0,t);
          if(d[t][i]==1) {
            d_post_mean[t][i] = d_post_mean[t][i] + 1.0;
            b_post_mean[t][i] = b_post_mean[t][i] + b[t][i];     
          }
        }
        
      }
      
      if(updatePi) {
      for (int k = 0; k<nmodels ; k++) {
        std::gamma_distribution<double> rgamma(cmodel[k],1.0);
        double rg = rgamma(gen);
        pi[k] = rg;
        //pi[k] = rg/(double)m;
      }
      double psum = std::accumulate(pi.begin(), pi.end(), 0.0);
      for (int k = 0; k<nmodels ; k++) {
        pi[k] = pi[k]/psum;
        pi_post_mean[k] = pi_post_mean[k] + pi[k];
      }
      std::fill(cmodel.begin(), cmodel.end(), 1.0);
      }
    }
    
    
    // Sample marker variance
    if(updateB) {
      //Rcout << "Sampling marker variance: " << it + 1 << "\n";
      arma::mat Sb(nt,nt, fill::zeros);
      dfb = 0.0;
      for (int t1 = 0; t1 < nt; t1++) {
        for (int t2 = t1; t2 < nt; t2++) {
          ssb = 0.0;
          if (t1==t2) {
            for (int i=0; i < m; i++) {
              if( (d[t1][i]==1) ) {
                ssb = ssb + b[t1][i]*b[t1][i];  
                dfb = dfb + 1.0;
              }
            }
            Sb(t1,t1) = ssb + ssb_prior[t1][t1];
          } 
          if (t1!=t2) {
            for (int i=0; i < m; i++) {
              if( (d[t1][i]==1) && (d[t2][i]==1) ) {
                ssb = ssb + b[t1][i]*b[t2][i];  
              }
            }
            Sb(t1,t2) = ssb + ssb_prior[t1][t2];
            Sb(t2,t1) = ssb + ssb_prior[t2][t1];
          } 
        }
      }
      
      int dfSb = dfb/nt + nub;
      arma::mat B = riwishart(dfSb, Sb);
      
      for (int t = 0; t < nt; t++) {
        vb_post[t][it] = B(t,t);
      }
      for (int t1 = 0; t1 < nt; t1++) {
        for (int t2 = 0; t2 < nt; t2++) {
          covb_post_mean[t1][t2] = covb_post_mean[t1][t2] + B(t1,t2);
        } 
      } 
      arma::mat Bi = arma::inv(B);
      //arma::mat Bi = arma::pinv(B);
    }


    // Sample residual variance
    if(updateE) {
      
      // version 1
      // arma::mat E(nt,nt, fill::zeros);
      // for ( int t = 0; t < nt; t++) {
      //   sse = 0.0;
      //   for ( int j = 0; j < n; j++) {
      //     sse = sse + e[t][j]*e[t][j];
      //   }
      //   std::chi_squared_distribution<double> rchisq_e(dfe);
      //   chi2 = rchisq_e(gen);
      //   ve_post[t][it] = (sse + sse_prior[t][t])/chi2;
      //   E(t,t) = (sse + sse_prior[t][t])/chi2;
      // }
      // arma::mat Ei = arma::inv(E);
      // end version 1
      
      // version 2
      arma::mat Se(nt,nt, fill::zeros);
      for ( int t1 = 0; t1 < nt; t1++) {
        for ( int t2 = t1; t2 < nt; t2++) {
          sse = 0.0;
          for ( int j = 0; j < n; j++) {
            sse = sse + e[t1][j]*e[t2][j];
          }
          Se(t1,t2) = sse + sse_prior[t1][t2];
          if (t1!=t2) {
            Se(t1,t2) = Se(t2,t1);
          }
          //std::chi_squared_distribution<double> rchisq_e(dfe);
          //chi2 = rchisq_e(gen);
          //ve_post[t][it] = (sse + sse_prior[t])/chi2;
          //E(t,t) = (sse + sse_prior[t][t])/chi2;
        }
      }
      //int dfSe = dfe + nue;
      int dfSe = dfe + nue;
      //arma::mat E = riwishart(dfSe, Se);
      arma::mat E = riwishart(dfSe, Se);
      arma::mat Ei = arma::inv(E);
      for (int t = 0; t < nt; t++) {
        ve_post[t][it] = E(t,t);
      }
      for (int t1 = 0; t1 < nt; t1++) {
        for (int t2 = 0; t2 < nt; t2++) {
          cove_post_mean[t1][t2] = cove_post_mean[t1][t2] + E(t1,t2);
        } 
      } 
      
      // end version 2
      
    }
    
    
    // Update mu and adjust residuals
    for ( int t = 0; t < nt; t++) {
      mu[t] = std::accumulate(e[t].begin(), e[t].end(), 0.0)/n;
      for ( int i = 0; i < n; i++) {
        e[t][i] = e[t][i] - mu[t];
      }                
      mu_post[t][it] = mu[t];
    }
    
    //Rcout << "Finished iteration: " << it + 1 << "\n";

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
      result[0][t][i] = b_post_mean[t][i]/nit;
      result[1][t][i] = d_post_mean[t][i]/nit;
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
      result[2][t][i] = mu_post[t][i];
      result[3][t][i] = vb_post[t][i];
      result[4][t][i] = ve_post[t][i];
    }
  }
  for (int t1=0; t1 < nt; t1++) {
    for (int t2=0; t2 < nt; t2++) {
      result[5][t1][t2] = covb_post_mean[t1][t2]/nit;
      result[6][t1][t2] = cove_post_mean[t1][t2]/nit;
    }
  }
  for (int t1=0; t1 < nt; t1++) {
    for (int t2=0; t2 < nt; t2++) {
      result[5][t1][t2] = covb_post_mean[t1][t2]/nit;
      result[6][t1][t2] = cove_post_mean[t1][t2]/nit;
      result[10][t1][t2] = B(t1,t2);
      //result[10][t2][t1] = B(t2,t1);
      result[11][t1][t2] = E(t1,t2);
      //result[11][t2][t1] = E(t2,t1);
    }
  }
  for (int t=0; t < nt; t++) {
    for (int i=0; i < nmodels; i++) {
      result[12][t][i] = pi[i];
      result[13][t][i] = pi_post_mean[i]/nit;
    }
  }  
  return result;

}




// [[Rcpp::export]]
std::vector<std::vector<std::vector<double>>>  mtsbayes(   std::vector<std::vector<double>> wy,
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
                                                                 int n,
                                                                 int nit,
                                                                 int method) {
  
  
  
  // Define local variables
  int nt = wy.size();
  //int n = y[0].size();
  int m = wy[0].size();
  int nmodels = models.size();
  
  //Rcout << "Number of traits: " << nt << "\n";
  //Rcout << "Number of individuals: " << m << "\n";
  //Rcout << "Number of markers: " << m << "\n";
  
  
  double ssb, sse, dfb, dfe, chi2, u, logliksum, psum, detC, diff;
  int mselect;
  
  //std::vector<std::vector<double>> e(nt, std::vector<double>(n, 0.0));
  std::vector<std::vector<double>> ww(nt, std::vector<double>(m, 0.0));
  std::vector<std::vector<int>> d(nt, std::vector<int>(m, 0));
  
  std::vector<double> mu(nt), rhs(nt), conv(nt);
  //std::vector<double> sse_prior(nt);
  //std::vector<std::vector<double>> ssb_prior(nt, std::vector<double>(nt, 0.0)); 
  
  std::vector<double> pmodel(nmodels), pcum(nmodels), loglik(nmodels), cmodel(nmodels);
  std::vector<double> pi_post_mean(nmodels);
  
  std::vector<std::vector<double>> b_post_mean(nt, std::vector<double>(m, 0.0));
  std::vector<std::vector<double>> d_post_mean(nt, std::vector<double>(m, 0.0));
  std::vector<std::vector<double>> ve_post(nt, std::vector<double>(nit, 0.0));
  std::vector<std::vector<double>> vb_post(nt, std::vector<double>(nit, 0.0));
  std::vector<std::vector<double>> covb_post_mean(nt, std::vector<double>(nt, 0.0));
  std::vector<std::vector<double>> cove_post_mean(nt, std::vector<double>(nt, 0.0));
  std::vector<std::vector<double>> mu_post(nt, std::vector<double>(nit, 0.0));
  
  std::vector<double> x2t(m);
  std::vector<std::vector<double>> x2(nt, std::vector<double>(m, 0.0));
  //std::vector<std::vector<double>> wy(nt, std::vector<double>(m, 0.0));
  std::vector<std::vector<double>> r(nt, std::vector<double>(m, 0.0));
  std::vector<int> order(m);
  
  
  // Prior variance and degrees of freedom
  // fix this dfe and n should be vectors of length nt
  dfe = n + nue;
  // fix this dfe and n should be vectors of length nt
  //for (int t1 = 0; t1 < nt; t1++) {
  //  for (int t2 = 0; t2 < nt; t2++) {
  //    ssb_prior[t1][t2] = (nub-2.0)/nub * B(t1,t2);
  //  }
  //  sse_prior[t1] = nue*E(t1,t1);
  //}
  
  
  // Initialize variables
  for (int i = 0; i < m; i++) {
    for (int t = 0; t < nt; t++) {
      // fix this dfe and n should be vectors of length nt
      ww[t][i] = (double)n;;
      // fix this dfe and n should be vectors of length nt
      //ww[t][i] = (double)n[t];;
      //for (int j = 0; j < n; j++) {
        //ww[t][i] = ww[t][i] + W[i][j]*W[i][j];
      //}
      r[t][i] = wy[t][i];
      x2[t][i] = (wy[t][i]/ww[t][i])*(wy[t][i]/ww[t][i]);
    }
  }
  
  // Wy - W'Wb
  for ( int i = 0; i < m; i++) {
    for (int t = 0; t < nt; t++) {
      if (b[t][i]!= 0.0) {
      //      for (int j = 0; j < m; j++) {
      for (size_t j = 0; j < LDindices[i].size(); j++) {
        //r[j]=r[j] - LD[i][j]*b[i];
        r[t][LDindices[i][j]]=r[t][LDindices[i][j]] - LDvalues[i][j]*b[t][i];
      }
      }
    }
  }
  
  
  for (int i = 0; i < nmodels; i++) {
    cmodel[i] = 1.0;
    pi_post_mean[i] = 0.0;
  }
  
  // // Mean adjust y
  // for ( int t = 0; t < nt; t++) {
  //   mu[t] = std::accumulate(y[t].begin(), y[t].end(), 0.0)/n;
  //   for ( int i = 0; i < n; i++) {
  //     e[t][i] = y[t][i] - mu[t];
  //   }                
  // }
  
  
  // for ( int t = 0; t < nt; t++) {
  //   for ( int i = 0; i < m; i++) {
  //     r[t][i] = 0.0;
  //     for (int j = 0; j < n; j++) {
  //       r[t][i] = r[t][i] + W[i][j]*e[t][j];
  //     }
  //     r[t][i] = ww[t][i]*b[t][i];
  //     wy[t][i] = r[t][i];
  //     x2[t][i] = (wy[t][i]/ww[t][i])*(wy[t][i]/ww[t][i]);
  //   }
  // }
  for ( int i = 0; i < m; i++) {
    x2t[i] = 0.0;
    for ( int t = 0; t < nt; t++) { 
      x2t[i] = x2t[i] + x2[t][i];
    }
  }
  std::iota(order.begin(), order.end(), 0);
  std::sort(  std::begin(order), 
              std::end(order),
              [&](int i1, int i2) { return x2t[i1] > x2t[i2]; } );
  
  
  // Initialize (co)variance matrices
  arma::mat C(nt,nt, fill::zeros);
  
  arma::mat Bi = arma::inv(B);
  arma::mat Ei = arma::inv(E);
  
  
  // Start Gibbs sampler
  
  //Rcout << "  " << "\n";
  //Rcout << "Starting Gibbs sampler" << "\n";
  //Rcout << "  " << "\n";
  
  std::random_device rd;
  unsigned int seed;
  seed = rd();
  std::mt19937 gen(seed);
  
  for ( int it = 0; it < nit; it++) {
    for ( int t = 0; t < nt; t++) {
      conv[t] = 0.0;
    }
    
    //arma::mat Ei = arma::inv(E);
    //arma::mat Bi = arma::inv(B);
    
    // Sample marker effects (Mixed)
    if (method==1) {
      
      
      for ( int i = 0; i < m; i++) {
        
        // Compute rhs                     
        for ( int t = 0; t < nt; t++) {
          //rhs[t] = 0.0;
          //for ( int j = 0; j < n; j++) {
          //  rhs[t] = rhs[t] + Ei(t,t)*W[i][j]*e[t][j]; 
          //}
          //rhs[t] = rhs[t] + Ei(t,t)*ww[t][i]*b[t][i];
          rhs[t] = Ei(t,t)*r[t][i] + Ei(t,t)*ww[t][i]*b[t][i];
        }
        

        // Sample beta
        for ( int t = 0; t < nt; t++) { 
          d[t][i] = 1;
        }
        arma::mat C(nt,nt, fill::zeros);
        for ( int t1 = 0; t1 < nt; t1++) {
          for ( int t2 = t1; t2 < nt; t2++) {
            C(t1,t2) = Bi(t1,t2);
            C(t2,t1) = Bi(t2,t1);
          }
          C(t1,t1) = C(t1,t1) + ww[t1][i]*Ei(t1,t1);       
        }
        arma::mat Ci = arma::inv(C);
        arma::mat mub = mvrnorm(Ci);
        
        
        for ( int t1 = 0; t1 < nt; t1++) {
          for ( int t2 = 0; t2 < nt; t2++) {
            mub(0,t1) = mub(0,t1) + Ci(t1,t2)*rhs[t2];
          }
        } 
        
        
        for ( int t = 0; t < nt; t++) {
          
          diff = mub(0,t)-b[t][i];
          //for (int j = 0; j < n; j++) {
          //  e[t][j]=e[t][j] - W[i][j]*(diff);
          //}
          for (size_t j = 0; j < LDindices[i].size(); j++) {
            r[t][LDindices[i][j]]=r[t][LDindices[i][j]] - LDvalues[i][j]*diff;
          }
          conv[t] = conv[t] + diff*diff;
          b[t][i] = mub(0,t);
          d_post_mean[t][i] = d_post_mean[t][i] + 1.0;
          b_post_mean[t][i] = b_post_mean[t][i] + b[t][i];     
        }
        
      }
    }
    
    // Sample marker effects (BayesC)
    if (method==4) {
      for ( int i0 = 0; i0 < m; i0++) {
        int i = order[i0];
        // Compute rhs                     
        for ( int t = 0; t < nt; t++) {
          //rhs[t] = 0.0;
          //for ( int j = 0; j < n; j++) {
          //  rhs[t] = rhs[t] + Ei(t,t)*W[i][j]*e[t][j]; 
          //}
          //rhs[t] = rhs[t] + Ei(t,t)*ww[t][i]*b[t][i];
          rhs[t] = Ei(t,t)*r[t][i] + Ei(t,t)*ww[t][i]*b[t][i];
        }
        
        
        // Compute lhs
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
        
        // Compute variance class probability
        std::fill(pmodel.begin(), pmodel.end(), 0.0);
        for (int k = 0; k<nmodels ; k++) {
          logliksum = 0.0;
          for (int l = 0; l<nmodels ; l++) {
            logliksum += std::exp(loglik[l] - loglik[k]);
          }
          pmodel[k] = 1.0/logliksum;
        }
        
        // Compute variance class cumulative probability   
        std::fill(pcum.begin(), pcum.end(), 0.0);
        for (int k = 0; k<nmodels ; k++) {
          psum = 0.0;
          for (int l = 0; l<k+1 ; l++) {
            psum += pmodel[l];
          }
          pcum[k] = psum;
        }
        
        
        // Sample marker variance class indicator - lambda and d are 1-based
        std::uniform_real_distribution<double> runif(0.0, 1.0);
        u = runif(gen);
        mselect = 0;
        for (int k = 0; k<nmodels ; k++) {
          if (pcum[k]>u) {
            mselect = k;
            break;
          }
        }
        cmodel[mselect] = cmodel[mselect] + 1.0; 
        
        // Sample beta
        for ( int t = 0; t < nt; t++) { 
          d[t][i] = models[mselect][t];
        }
        arma::mat C = Bi;
        for ( int t1 = 0; t1 < nt; t1++) {
          if(models[mselect][t1]==1) {
            C(t1,t1) = C(t1,t1) + ww[t1][i]*Ei(t1,t1);       
          } 
        }
        arma::mat Ci = arma::inv(C);
        arma::mat mub = mvrnorm(Ci);
        
        
        for ( int t1 = 0; t1 < nt; t1++) {
          mub(0,t1) = 0.0;
          for ( int t2 = 0; t2 < nt; t2++) {
            if(models[mselect][t2]==1) {
              mub(0,t1) = mub(0,t1) + Ci(t1,t2)*rhs[t2];
            }
          }
        } 
        
        for ( int t = 0; t < nt; t++) {
          
          diff = mub(0,t)-b[t][i];
          //for (int j = 0; j < n; j++) {
          //  e[t][j]=e[t][j] - W[i][j]*(diff);
          //}
          for (size_t j = 0; j < LDindices[i].size(); j++) {
            r[t][LDindices[i][j]]=r[t][LDindices[i][j]] - LDvalues[i][j]*diff;
          }
          conv[t] = conv[t] + diff*diff;
          b[t][i] = mub(0,t);
          if(d[t][i]==1) {
            d_post_mean[t][i] = d_post_mean[t][i] + 1.0;
            b_post_mean[t][i] = b_post_mean[t][i] + b[t][i];     
          }
        }
        
      }
      
      if(updatePi) {
        for (int k = 0; k<nmodels ; k++) {
          std::gamma_distribution<double> rgamma(cmodel[k],1.0);
          double rg = rgamma(gen);
          pi[k] = rg/(double)m;
          pi_post_mean[k] = pi_post_mean[k] + pi[k];
        }
        std::fill(cmodel.begin(), cmodel.end(), 1.0);
      }
    }
    
    
    // Sample marker variance
    if(updateB) {
      //Rcout << "Sampling marker variance: " << it + 1 << "\n";
      arma::mat Sb(nt,nt, fill::zeros);
      dfb = 0.0;
      for (int t1 = 0; t1 < nt; t1++) {
        for (int t2 = t1; t2 < nt; t2++) {
          ssb = 0.0;
          if (t1==t2) {
            for (int i=0; i < m; i++) {
              if( (d[t1][i]==1) ) {
                ssb = ssb + b[t1][i]*b[t2][i];  
                dfb = dfb + 1.0;
              }
            }
            Sb(t1,t2) = ssb + ssb_prior[t1][t2];
          } 
          if (t1!=t2) {
            for (int i=0; i < m; i++) {
              if( (d[t1][i]==1) && (d[t2][i]==1) ) {
                ssb = ssb + b[t1][i]*b[t2][i];  
              }
            }
            Sb(t1,t2) = ssb + ssb_prior[t1][t2];
            Sb(t2,t1) = ssb + ssb_prior[t2][t1];
          } 
        }
      }
      int dfSb = dfb/nt + nub;
      //arma::mat B = rinvwish(dfSb, Sb);
      //arma::mat B = riwish(dfSb, Sb);
      //arma::mat B = riwishart(dfSb, Sb);
      arma::mat B = riwishart(dfSb, Sb);
      
      for (int t = 0; t < nt; t++) {
        vb_post[t][it] = B(t,t);
      }
      for (int t1 = 0; t1 < nt; t1++) {
        for (int t2 = 0; t2 < nt; t2++) {
          covb_post_mean[t1][t2] = covb_post_mean[t1][t2] + B(t1,t2)/(sqrt(B(t1,t1))*sqrt(B(t2,t2)));
        } 
      } 
      arma::mat Bi = arma::inv(B);
      //arma::mat Bi = arma::pinv(B);
    }
    
    // Sample residual variance
    if(updateE) {
      
      // version 1
      // arma::mat E(nt,nt, fill::zeros);
      // for ( int t = 0; t < nt; t++) {
      //   sse = 0.0;
      //   for ( int j = 0; j < n; j++) {
      //     sse = sse + e[t][j]*e[t][j];
      //   }
      //   std::chi_squared_distribution<double> rchisq_e(dfe);
      //   chi2 = rchisq_e(gen);
      //   ve_post[t][it] = (sse + sse_prior[t][t])/chi2;
      //   E(t,t) = (sse + sse_prior[t][t])/chi2;
      // }
      // arma::mat Ei = arma::inv(E);
      // end version 1
      
      // version 2
      arma::mat Se(nt,nt, fill::zeros);
      for ( int t1 = 0; t1 < nt; t1++) {
        sse = 0.0;
        for ( int i = 0; i < m; i++) {
          sse = sse + b[t1][i] * (r[t1][i] + wy[t1][i]);
        }
        sse = yy[t1] - sse;
        // for ( int t2 = t1; t2 < nt; t2++) {
        //   sse = 0.0;
        //   for ( int j = 0; j < n; j++) {
        //     sse = sse + e[t1][j]*e[t2][j];
        //   }
        //   Se(t1,t2) = sse + sse_prior[t1][t2];
        //   if (t1!=t2) {
        //     Se(t1,t2) = Se(t2,t1);
        //   }
        //   //std::chi_squared_distribution<double> rchisq_e(dfe);
        //   //chi2 = rchisq_e(gen);
        //   //ve_post[t][it] = (sse + sse_prior[t])/chi2;
        //   //E(t,t) = (sse + sse_prior[t][t])/chi2;
        // }
      }
      //int dfSe = dfe + nue;
      int dfSe = dfe + nue;
      //arma::mat E = riwishart(dfSe, Se);
      arma::mat E = riwishart(dfSe, Se);
      arma::mat Ei = arma::inv(E);
      for (int t = 0; t < nt; t++) {
        ve_post[t][it] = E(t,t);
      }
      for (int t1 = 0; t1 < nt; t1++) {
        for (int t2 = 0; t2 < nt; t2++) {
          cove_post_mean[t1][t2] = cove_post_mean[t1][t2] + E(t1,t2)/(sqrt(E(t1,t1))*sqrt(E(t2,t2)));
        } 
      } 
      
      // end version 2
      
    }
    
    
    // Update mu and adjust residuals
    //for ( int t = 0; t < nt; t++) {
    //  mu[t] = std::accumulate(e[t].begin(), e[t].end(), 0.0)/n;
    //  for ( int i = 0; i < n; i++) {
    //    e[t][i] = e[t][i] - mu[t];
    //  }                
    //  mu_post[t][it] = mu[t];
    //}
    
    //Rcout << "Finished iteration: " << it + 1 << "\n";
    
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
    result[7][t].resize(m);
    result[8][t].resize(m);
    result[9][t].resize(m);
    result[10][t].resize(nt);
    result[11][t].resize(nt);
    result[12][t].resize(nmodels);
    result[13][t].resize(nmodels);
    result[14][t].resize(m);
  }
  
  for (int t=0; t < nt; t++) {
    for (int i=0; i < m; i++) {
      result[0][t][i] = b_post_mean[t][i]/nit;
      result[1][t][i] = d_post_mean[t][i]/nit;
      result[9][t][i] = b[t][i];
      result[14][t][i] = order[i];
    }
    for (int i=0; i < m; i++) {
      result[7][t][i] = wy[t][i];
      result[8][t][i] = r[t][i];
    }
  }
  
  for (int t=0; t < nt; t++) {
    for (int i=0; i < nit; i++) {
      result[2][t][i] = mu_post[t][i];
      result[3][t][i] = vb_post[t][i];
      result[4][t][i] = ve_post[t][i];
    }
  }
  for (int t1=0; t1 < nt; t1++) {
    for (int t2=0; t2 < nt; t2++) {
      result[5][t1][t2] = covb_post_mean[t1][t2]/nit;
      result[6][t1][t2] = cove_post_mean[t1][t2]/nit;
      result[10][t1][t2] = B(t1,t2);
      //result[10][t2][t1] = B(t2,t1);
      result[11][t1][t2] = E(t1,t2);
      //result[11][t2][t1] = E(t2,t1);
    }
  }
  for (int t=0; t < nt; t++) {
    for (int i=0; i < nmodels; i++) {
      result[12][t][i] = pi[i];
      result[13][t][i] = pi_post_mean[i]/nit;
    }
  }  
  return result;
  
}
