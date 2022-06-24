#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
std::vector<std::vector<double>> solvebed( const char* file, 
                             			int n, 
                             			std::vector<int> cls, 
                             			int nit, 
                             			std::vector<double> af, 
                             			std::vector<double> b, 
                             			std::vector<double> lambda, 
                             			std::vector<double> y) {
  
  
  FILE *file_stream = fopen( file, "rb" );
  
  int m = cls.size();
  
  size_t nbytes_read;
  size_t nbytes = ( n + 3 ) / 4;
  
  unsigned char *buffer = (unsigned char *) malloc( nbytes );
  unsigned char buf_k; 
  
  std::vector<std::vector<double>> g(m, std::vector<double>(n, 0.0));
  std::vector<double> map(4);
  std::vector<double> dxx(m);
  for ( int i = 0; i < m; i++) {
    dxx[i] = 0.0;
  }                
  
  
  ////////////////////////////////////////////////////////////////////////
  //  00 01 10 11         bit level  corresponds to
  //  0  1  2  3          xij level  corresponds to
  //  2  NA  1  0         number of copies of first allele in bim file
  ////////////////////////////////////////////////////////////////////////
  
  //std::cout << "  " << "\n";
  //std::cout << "Reading genotypes" << "\n";
  //std::cout << "  " << "\n";
  
  for (int i = 0; i < m; i++) {
    // cls[i] is 1-based
    long int offset = (cls[i]-1)*nbytes + 3;
    fseek( file_stream, offset, SEEK_SET );
    nbytes_read = fread( buffer, sizeof(unsigned char), nbytes, file_stream );
    //if (nbytes_read != nbytes) {
    //  std::cout << "Error reading data: nbytes_read != nbytes" << "\n";
    //}
    int j = 0; 
    map[0] = 2.0 - 2.0*af[i];
    map[1] = 0.0;
    map[2] = 1.0 - 2.0*af[i];
    map[3] = -2.0*af[i];
    for (size_t k = 0; k < nbytes; k++) {
      buf_k = buffer[k];
      for (int pos = 0; pos < 4; pos++, j++) {
        if (j < n) {
          g[i][j] = map[buf_k & 3];
          dxx[i] = dxx[i] + g[i][j]*g[i][j];
          buf_k = buf_k >> 2;
        } 
      }
    }
  }
  free( buffer );
  fclose( file_stream );
  
  
  //std::cout << "  " << "\n";
  //std::cout << "Starting solver" << "\n";
  //std::cout << "  " << "\n";
  
  std::vector<double> e(n);
  for ( int i = 0; i < n; i++) {
    e[i] = y[i];
  }                
  
  double rhs, lhs, bnew, conv, diff;
  
  for ( int it = 0; it < nit; it++) {
    conv = 0.0;
    for ( int i = 0; i < m; i++) {
      lhs = dxx[i] + lambda[i];
      rhs = 0.0;
      for ( int j = 0; j < n; j++) {
        rhs = rhs + g[i][j]*e[j]; 
      }
      rhs = rhs + dxx[i]*b[i];
      bnew = rhs/lhs;
      diff = bnew-b[i];
      for (int j = 0; j < n; j++) {
        e[j]=e[j] - g[i][j]*(diff);
      }
      conv = conv + diff*diff;
      b[i] = bnew;
    }
    //std::cout << "Finished iteration: " << it + 1 << "\n";
    //std::cout << "Convergence: " << conv << "\n";
    
  }
  
  fclose( file_stream );
  
  // Summarize results
  std::vector<std::vector<double>> result(2);
  result[0].resize(m);
  result[1].resize(n);
  
  for (int i=0; i < m; i++) {
    result[0][i] = b[i];
  }
  for (int i=0; i < n; i++) {
    result[1][i] = e[i];
  }
  return result;

}

// [[Rcpp::export]]
std::vector<std::vector<std::vector<double>>> mtsolvebed( const char* file, 
                                            int n, 
                                            std::vector<int> cls, 
                                            int nit, 
                                            std::vector<double> af, 
                                            std::vector<std::vector<double>> b, 
                                            std::vector<std::vector<double>> lambda, 
                                            std::vector<std::vector<double>> y) {
  
  FILE *file_stream = fopen( file, "rb" );
  
  
  int m = cls.size();
  
  size_t nbytes_read;
  size_t nbytes = ( n + 3 ) / 4;
  
  unsigned char *buffer = (unsigned char *) malloc( nbytes );
  unsigned char buf_k; 
  
  std::vector<std::vector<double>> g(m, std::vector<double>(n, 0.0));
  std::vector<double> map(4);
  std::vector<double> dxx(m);
  for ( int i = 0; i < m; i++) {
    dxx[i] = 0.0;
  }                
  
  
  ////////////////////////////////////////////////////////////////////////
  //  00 01 10 11         bit level  corresponds to
  //  0  1  2  3          xij level  corresponds to
  //  2  NA  1  0         number of copies of first allele in bim file
  ////////////////////////////////////////////////////////////////////////
  
  //std::cout << "  " << "\n";
  //std::cout << "Reading genotypes" << "\n";
  //std::cout << "  " << "\n";
  
  for (int i = 0; i < m; i++) {
    // cls[i] is 1-based
    long int offset = (cls[i]-1)*nbytes + 3;
    fseek( file_stream, offset, SEEK_SET );
    nbytes_read = fread( buffer, sizeof(unsigned char), nbytes, file_stream );
    //if (nbytes_read != nbytes) {
    //  std::cout << "Error reading data: nbytes_read != nbytes" << "\n";
    //}
    int j = 0; 
    map[0] = 2.0 - 2.0*af[i];
    map[1] = 0.0;
    map[2] = 1.0 - 2.0*af[i];
    map[3] = -2.0*af[i];
    for (size_t k = 0; k < nbytes; k++) {
      buf_k = buffer[k];
      for (int pos = 0; pos < 4; pos++, j++) {
        if (j < n) {
          g[i][j] = map[buf_k & 3];
          dxx[i] = dxx[i] + g[i][j]*g[i][j];
          buf_k = buf_k >> 2;
        }
      }
    }
  }
  free( buffer );
  fclose( file_stream );
  
  
  //std::cout << "  " << "\n";
  //std::cout << "Starting solver" << "\n";
  //std::cout << "  " << "\n";
  
  int nt = y.size();
  std::vector<std::vector<double>> e(nt, std::vector<double>(n, 0.0));
  
  for ( int t = 0; t < nt; t++) {
    for ( int i = 0; i < n; i++) {
      e[t][i] = y[t][i];
    }                
  }
  
  std::vector<double> rhs(nt), lhs(nt), bnew(nt), conv(nt), diff(nt);
  
  for ( int it = 0; it < nit; it++) {
    for ( int t = 0; t < nt; t++) {
      conv[t] = 0.0;
    }
    for ( int i = 0; i < m; i++) {
      for ( int t = 0; t < nt; t++) {
        lhs[t] = dxx[i] + lambda[t][i];
        rhs[t] = 0.0;
        for ( int j = 0; j < n; j++) {
          rhs[t] = rhs[t] + g[i][j]*e[t][j]; 
        }
        rhs[t] = rhs[t] + dxx[i]*b[t][i];
        bnew[t] = rhs[t]/lhs[t];
        diff[t] = bnew[t]-b[t][i];
        for ( int j = 0; j < n; j++) {
          e[t][j]=e[t][j] - g[i][j]*(diff[t]);
        }
        conv[t] = conv[t] + diff[t]*diff[t];
        b[t][i] = bnew[t];
      }
    }
    //std::cout << "Finished iteration: " << it + 1 << "\n";
    //for ( int t = 0; t < nt; t++) {
    //  std::cout << "Convergence: " << conv[t] << "\n";
    //}
  }
  
  // Summarize results
  std::vector<std::vector<std::vector<double>>> result;
  result.resize(2);
  
  result[0].resize(nt);
  result[1].resize(nt);
  
  for (int t=0; t < nt; t++) {
    result[0][t].resize(m);
    result[1][t].resize(n);
  }
  
  for (int t=0; t < nt; t++) {
    for (int i=0; i < m; i++) {
      result[0][t][i] = b[t][i];
    }
  }
  
  for (int t=0; t < nt; t++) {
    for (int i=0; i < n; i++) {
      result[1][t][i] = e[t][i];
    }
  }
  
  return result;

}
