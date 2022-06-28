#include <Rcpp.h>
#include <cerrno>
using namespace Rcpp;

// [[Rcpp::export]]
std::vector<double> grsbed( const char* file, 
                           int n, 
                           std::vector<int> cls, 
                           std::vector<double> af, 
                           std::vector<double> b) {
  
  FILE *file_stream = fopen( file, "rb" );
  
  size_t nbytes_read;
  size_t nbytes = ( n + 3 ) / 4;
  
  int m = cls.size();
  
  unsigned char buf_k; 
  unsigned char *buffer = (unsigned char *) malloc( nbytes );
  
  std::vector<double> map(4);
  std::vector<double> grs(n);
  for ( int i = 0; i < n; i++) {
    grs[i] = 0.0;
  }                
  
  
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
    if (nbytes_read != nbytes) {
      Rcerr << "Error reading data: nbytes_read != nbytes" << "\n";
    }
    int j = 0;
    map[0] = b[i]*2.0;
    map[1] = 2.0*af[i]*b[i];
    map[2] = b[i];
    map[3] = 0.0;
    for (size_t k = 0; k < nbytes; k++) {
      buf_k = buffer[k];
      for (int pos = 0; pos < 4; pos++, j++) {
        if (j < n) {
          grs[j] = grs[j] + map[buf_k & 3];
          buf_k = buf_k >> 2;
        } 
      }
    }
  }
  free( buffer );
  fclose( file_stream );
  
  return grs;
}


// [[Rcpp::export]]
std::vector<std::vector<double>> mtgrsbed( const char* file, 
                                          int n, 
                                          std::vector<int> cls, 
                                          std::vector<double> af,
                                          bool scale,
                                          std::vector<std::vector<double>> b) {
  
  
  FILE *file_stream = fopen( file, "rb" );
  
  int nt = b.size();
  int m = cls.size();
  
  size_t nbytes_read;
  size_t nbytes = ( n + 3 ) / 4;
  
  unsigned char buf_k; 
  std::vector<double> x(n); 
  unsigned char *buffer = (unsigned char *) malloc( nbytes );
  
  std::vector<double> map(4);
  std::vector<std::vector<double>> grs(nt, std::vector<double>(n, 0.0));

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
    if (nbytes_read != nbytes) {
      Rcerr << "Error reading data: nbytes_read != nbytes" << "\n";
    }
    int j = 0; 
    if(scale) {
      map[0] = (2.0 - 2.0*af[i])/sqrtf(2.0*af[i]*(1.0-af[i]));
      map[1] = 0.0;
      map[2] = (1.0 - 2.0*af[i])/sqrtf(2.0*af[i]*(1.0-af[i]));
      map[3] = (-2.0*af[i])/sqrtf(2.0*af[i]*(1.0-af[i]));
    } else {
      map[0] = 2.0;
      map[1] = - 2.0*af[i];
      map[2] = 1.0 ;
      map[3] = 0.0;
    }
    for (size_t k = 0; k < nbytes; k++) {
      buf_k = buffer[k];
      for (int pos = 0; pos < 4; pos++, j++) {
        if (j < n) {
          x[j] = map[buf_k & 3];
          buf_k = buf_k >> 2;
        } 
      }
    }
    for ( int t = 0; t < nt; t++) {
      for ( int j = 0; j < n; j++) {
        grs[t][j] = grs[t][j] + b[t][i]*x[j];
      }
    }
  }
  free( buffer );
  fclose( file_stream );
  
  return grs;
}
