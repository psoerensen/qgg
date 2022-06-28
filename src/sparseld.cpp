//[[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
std::vector<int> pruneld( const char* file,
                          int ldsize,
                          std::vector<int> cls,
                          std::vector<float> p,
                          float threshold,
                          float r2) {
  
  FILE *file_stream = fopen( file, "rb" );
  
  int m = cls.size();
  std::vector<int> mask1(m), mask2(m);
  for ( int i = 0; i < m; i++) { 
    mask1[i] = 0;
    mask2[i] = 0;
  }
  
  
  size_t nbytes_read;
  size_t nbytes=ldsize*2+1; 
  
  float *buffer = (float *) malloc( nbytes*4 );
  
  for ( int i = 0; i < m; i++) {
    // cls[i] is 1-based
    int i0 = cls[i] - 1;
    if(mask1[i0]==0 && mask2[i0]==0 && p[i0]<threshold) {
      mask2[i0]=1;
      long int offset = i0*nbytes*4;
      fseek( file_stream, offset, SEEK_SET );
      nbytes_read = fread( buffer, sizeof(float), nbytes, file_stream );
      if (nbytes_read != nbytes) {
        Rcerr << "Error reading data: nbytes_read != nbytes" << "\n";
      }
      int k0=0;
      for ( size_t j = 0; j < nbytes; j++) {
        int k1 = i0 - ldsize + k0;
        if (k1>=0 && k1<m && k1!=i0) {
          float r2ij=buffer[j]*buffer[j];
          if (r2ij>=r2) mask1[k1] = 1; 
        }
        k0 = k0 + 1;
      }
    }
  }
  free( buffer );
  fclose( file_stream );

  
  return mask2;
}

// [[Rcpp::export]]
std::vector<std::vector<std::vector<int>>> pruneldmat( const char* file,
                                                       int ldsize,
                                                       std::vector<std::vector<float>> p,
                                                       std::vector<float> threshold,
                                                       float r2) {
  
  FILE *file_stream = fopen( file, "rb" );
  
  int np = p.size();
  int m = p[0].size();
  
  std::vector<int> cls(m);
  
  
  int nthold = threshold.size();
  std::vector<int> mask1(m), mask2(m);
  std::vector<std::vector<float>> ld(m, std::vector<float>(ldsize*2+1, 0.0));
  
  std::vector<std::vector<std::vector<int>>> mask(np, std::vector<std::vector<int>>(nthold, std::vector<int>(m)));
  
  
  size_t nbytes_read;
  size_t nbytes=ldsize*2+1; 
  
  float *buffer = (float *) malloc( nbytes*4 );
  for ( int i = 0; i < m; i++) {
    nbytes_read = fread( buffer, sizeof(float), nbytes, file_stream );
    if (nbytes_read != nbytes) {
      Rcerr << "Error reading data: nbytes_read != nbytes" << "\n";
    }
    for ( size_t j = 0; j < nbytes; j++) {
      ld[i][j] = buffer[j];
    }
  }
  free( buffer );
  fclose( file_stream );
  
  for ( int t1 = 0; t1 < np; t1++) {
    
    std::iota(cls.begin(), cls.end(), 0);
    std::sort(  std::begin(cls), 
                std::end(cls),
                [&](int i1, int i2) { return p[t1][i1] < p[t1][i2]; } );
    
    //std::cout << "Pruning column: " << t1+1 << "\n";
    
    
    for ( int t2 = 0; t2 < nthold; t2++) {
      
      for ( int i = 0; i < m; i++) { 
        mask1[i] = 0;
        mask2[i] = 0;
      }
      
      for ( int i = 0; i < m; i++) {
        int j = cls[i] ;
        if(mask1[j]==0 && mask2[j]==0 && p[t1][j]<threshold[t2]) {
          mask2[j]=1;
          int k0=0;
          for ( size_t k = 0; k < nbytes; k++) {
            int k1 = j - ldsize + k0;
            if (k1>=0 && k1<m && k1!=j) {
              float r2ij=ld[j][k]*ld[j][k];
              if (r2ij>=r2) mask1[k1] = 1; 
            }
            k0 = k0 + 1;
          }
        }
      }
      for ( int i = 0; i < m; i++) { 
        mask[t1][t2][i] = mask2[i];
      }
      //std::cout << "Finished threshold: " << threshold[t2] << "\n";
      
    }
  }
  return mask;
}
