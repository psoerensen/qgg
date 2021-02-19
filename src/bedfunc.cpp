#include <Rcpp.h>
#include <cerrno>
using namespace Rcpp;


// [[Rcpp::export]]
IntegerMatrix   readbed( const char* file, 
                         int n, 
                         std::vector<int> cls) {
  
  FILE *file_stream = fopen ( file, "rb" );
  
  int m = cls.size();
  
  size_t nbytes = ( n + 3 ) / 4;
  size_t nbytes_read;
  
  unsigned char *buffer = (unsigned char *) malloc( nbytes );
  unsigned char buf_k; 
  
  //std::vector<int> map(4);
  //map[0] = 2;
  //map[1] = NA_INTEGER;
  //map[2] = 1;
  //map[3] = 0;
  
  int map[4];
  map[0] = 2;
  map[1] = NA_INTEGER;
  map[2] = 1;
  map[3] = 0;
  
  IntegerMatrix X(n, m);
  
  
  ////////////////////////////////////////////////////////////////////////
  //  00 01 10 11         bit level  corresponds to
  //  0  1  2  3          xij level  corresponds to
  //  2  NA  1  0         number of copies of first allele in bim file
  ////////////////////////////////////////////////////////////////////////
  
  for ( int i = 0; i < m; i++) {
    // cls[i] is 1-based
    long int offset = (cls[i]-1)*nbytes + 3;
    fseek( file_stream, offset, SEEK_SET );
    nbytes_read = fread( buffer, sizeof(unsigned char), nbytes, file_stream );
    int j = 0; 
    for ( size_t k = 0; k < nbytes; k++) {
      buf_k = buffer[k];
      for ( int pos = 0; pos < 4; pos++, j++) {
        if (j < n) {
          X(j,i) = map[buf_k & 3];
          buf_k = buf_k >> 2;
        } 
      }
    }
  }
  
  free( buffer );
  
  return X;
}


// [[Rcpp::export]]
NumericMatrix readW( const char* file, 
                        int n, 
                        std::vector<int> cls, 
                        std::vector<float> af) {
  
  
  FILE *file_stream = fopen( file, "rb" );
  
  int m = cls.size();
  
  size_t nbytes_read;
  size_t nbytes = ( n + 3 ) / 4;
  
  unsigned char *buffer = (unsigned char *) malloc( nbytes );
  unsigned char buf_k; 
  
  //std::vector<std::vector<float>> g(m, std::vector<float>(n, 0.0));
  NumericMatrix W(n, m);
  
  std::vector<float> map(4);
  
  
  ////////////////////////////////////////////////////////////////////////
  //  00 01 10 11         bit level  corresponds to
  //  0  1  2  3          xij level  corresponds to
  //  2  NA  1  0         number of copies of first allele in bim file
  ////////////////////////////////////////////////////////////////////////
  
  std::cout << "  " << "\n";
  std::cout << "Reading genotypes" << "\n";
  std::cout << "  " << "\n";
  
  for (int i = 0; i < m; i++) {
    // cls[i] is 1-based
    long int offset = (cls[i]-1)*nbytes + 3;
    fseek( file_stream, offset, SEEK_SET );
    nbytes_read = fread( buffer, sizeof(unsigned char), nbytes, file_stream );
    int j = 0; 
    map[0] = 2.0 - 2.0*af[i];
    map[1] = 0.0;
    map[2] = 1.0 - 2.0*af[i];
    map[3] = -2.0*af[i];
    for (size_t k = 0; k < nbytes; k++) {
      buf_k = buffer[k];
      for (int pos = 0; pos < 4; pos++, j++) {
        if (j < n) {
          //g[i][j] = map[buf_k & 3];
          W(j,i) = map[buf_k & 3];
          buf_k = buf_k >> 2;
        } 
      }
    }
  }
  free( buffer );
  return W;
}


// [[Rcpp::export]]
std::vector<std::vector<float>> getWlist( const char* file, 
                                          int n, 
                                          std::vector<int> cls, 
                                          std::vector<float> af) {
  
  
  FILE *file_stream = fopen( file, "rb" );
  
  int m = cls.size();
  
  size_t nbytes_read;
  size_t nbytes = ( n + 3 ) / 4;
  
  unsigned char *buffer = (unsigned char *) malloc( nbytes );
  unsigned char buf_k; 
  
  std::vector<std::vector<float>> g(m, std::vector<float>(n, 0.0));
  std::vector<float> map(4);
  
  
  ////////////////////////////////////////////////////////////////////////
  //  00 01 10 11         bit level  corresponds to
  //  0  1  2  3          xij level  corresponds to
  //  2  NA  1  0         number of copies of first allele in bim file
  ////////////////////////////////////////////////////////////////////////
  
  std::cout << "  " << "\n";
  std::cout << "Reading genotypes" << "\n";
  std::cout << "  " << "\n";
  
  for (int i = 0; i < m; i++) {
    // cls[i] is 1-based
    long int offset = (cls[i]-1)*nbytes + 3;
    fseek( file_stream, offset, SEEK_SET );
    nbytes_read = fread( buffer, sizeof(unsigned char), nbytes, file_stream );
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
          buf_k = buf_k >> 2;
        } 
      }
    }
  }
  free( buffer );
  return g;
}



// [[Rcpp::export]]
IntegerMatrix freqbed( const char* file, 
                       int n, 
                       std::vector<int> cls) {
  
  FILE *file_stream = fopen( file, "rb" );
  
  int m = cls.size();
  
  size_t nbytes_read;
  size_t nbytes = ( n + 3 ) / 4;
  
  unsigned char buf_k, xij;
  unsigned char *buffer = (unsigned char *) malloc( nbytes );
  
  std::vector<int> map(4);
  map[0] = 0;
  map[1] = 1;
  map[2] = 2;
  map[3] = 3;
  
  ////////////////////////////////////////////////////////////////////////
  //  00 01 10 11         bit level  corresponds to
  //  0  1  2  3          xij level  corresponds to
  //  2  NA  1  0         number of copies of first allele in bim file
  ////////////////////////////////////////////////////////////////////////
  
  IntegerMatrix X(4, m);
  
  for (int i = 0; i < m; i++) {
    // cls[i] is 1-based
    long int offset = (cls[i]-1)*nbytes + 3;
    fseek( file_stream, offset, SEEK_SET );
    nbytes_read = fread( buffer, sizeof(unsigned char), nbytes, file_stream );
    int j = 0;
    for (size_t k = 0; k < nbytes; k++) {
      buf_k = buffer[k];
      for (int pos = 0; pos < 4; pos++, j++) {
        if (j < n) {
          X(map[buf_k & 3], i) +=  1;
          buf_k = buf_k >> 2;
        } 
      }
    }
  }
  free( buffer );
  
  return X;
}



// [[Rcpp::export]]
std::vector<std::vector<std::vector<float>>> summarybed( const char* file, 
                                                         int n, 
                                                         std::vector<int> cls,
                                                         std::vector<float> af, 
                                                         std::vector<std::vector<float>> y) {
  
  FILE *file_stream = fopen( file, "rb" );
  
  size_t nbytes = ( n + 3 ) / 4;
  size_t nbytes_read;
  
  unsigned char buf_k; 
  unsigned char *buffer = (unsigned char *) malloc( nbytes );
  
  int m = cls.size();
  int nt = y.size();
  
  std::vector<std::vector<float>> xy(nt, std::vector<float>(m, 0.0));
  std::vector<std::vector<float>> xx(nt, std::vector<float>(m, 0.0));
  std::vector<float> map(4);
  std::vector<float> x(n);
  
  
  ////////////////////////////////////////////////////////////////////////
  //  00 01 10 11         bit level  corresponds to
  //  0  1  2  3          xij level  corresponds to
  //  2  NA  1  0         number of copies of first allele in bim file
  ////////////////////////////////////////////////////////////////////////
  
  
  for (int i = 0; i < m; i++) {
    // cls[i] is 1-based
    long int offset = (cls[i]-1)*nbytes + 3;
    fseek( file_stream, offset, SEEK_SET );
    nbytes_read = fread( buffer, sizeof(unsigned char), nbytes, file_stream );
    int j = 0; 
    map[0] = 2.0 - 2.0*af[i];
    map[1] = 0.0;
    map[2] = 1.0 - 2.0*af[i];
    map[3] = -2.0*af[i];
    for (size_t k = 0; k < nbytes; k++) {
      buf_k = buffer[k];
      for (int pos = 0; pos < 4; pos++, j++) {
        if (j < n) {
          x[j] = map[buf_k & 3];
          buf_k = buf_k >> 2;
        }
      }
    }
    for (int t = 0; t < nt; t++) {
      xx[t][i] = std::inner_product(std::begin(x), std::end(x), std::begin(x), 0.0);
      xy[t][i] = std::inner_product(std::begin(x), std::end(x), std::begin(y[t]), 0.0);
    }
    
  }
  free( buffer );
  
  std::vector<std::vector<std::vector<float>>> result;
  result.resize(2);
  result[0].resize(nt);
  result[1].resize(nt);
  
  for (int t = 0; t < nt; t++) {
    result[0][t].resize(m);
    result[1][t].resize(m);
    for (int i = 0; i < m; i++) {
      result[0][t][i] = xx[t][i];        
      result[1][t][i] = xy[t][i];        
    }
  }                
  
  return result;
  
}
