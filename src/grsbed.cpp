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

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(openmp)]]

#include <Rcpp.h>
#include <cstdio>
#include <vector>
#include <cmath>
#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

// [[Rcpp::export]]
std::vector<std::vector<double>> mtgrsbed_omp(
    const char* file,
    int n,
    const std::vector<int>& cls,                     // 1-based SNP indices in BED
    const std::vector<double>& af,
    bool scale,
    const std::vector<std::vector<double>>& b,       // b[trait][snp]
    int nthreads = 1
) {
  FILE* file_stream = std::fopen(file, "rb");
  if (!file_stream) {
    stop("Could not open BED file.");
  }
  
  const int nt = static_cast<int>(b.size());
  if (nt == 0) {
    std::fclose(file_stream);
    return {};
  }

  const int m = static_cast<int>(cls.size());
  if ((int)af.size() != m) {
    std::fclose(file_stream);
    stop("af.size() != cls.size()");
  }

  for (int t = 0; t < nt; ++t) {
    if ((int)b[t].size() != m) {
      std::fclose(file_stream);
      stop("All rows of b must have length m = cls.size()");
    }
  }

  const size_t nbytes = (n + 3) / 4;

  // BED uses 2 bits/genotype:
  // 00 01 10 11  ->  2, NA, 1, 0 copies of first allele in .bim
  unsigned char* buffer = (unsigned char*) std::malloc(nbytes);
  if (!buffer) {
    std::fclose(file_stream);
    stop("Failed to allocate BED buffer.");
  }

  // Decoded genotype values for current SNP
  std::vector<double> x(n);

  // Repack b into SNP-major layout: b_snp[i * nt + t]
  // so all 500 trait weights for one SNP are contiguous.
  std::vector<double> b_snp((size_t)m * nt);
  for (int t = 0; t < nt; ++t) {
    for (int i = 0; i < m; ++i) {
      b_snp[(size_t)i * nt + t] = b[t][i];
    }
  }

  // Store scores individual-major: grs_flat[j * nt + t]
  // so all 500 trait scores for one individual are contiguous.
  std::vector<double> grs_flat((size_t)n * nt, 0.0);

#ifdef _OPENMP
  if (nthreads > 0) omp_set_num_threads(nthreads);
#endif

#pragma omp parallel
{
  std::vector<double> map(4);
  
  for (int i = 0; i < m; ++i) {
    
    // One thread reads + decodes current SNP into x[]
#pragma omp single
{
  long long offset = (long long)(cls[i] - 1) * (long long)nbytes + 3LL;
  
  if (std::fseek(file_stream, offset, SEEK_SET) != 0) {
    std::free(buffer);
    std::fclose(file_stream);
    stop("fseek failed.");
  }
  
  size_t nbytes_read = std::fread(buffer, sizeof(unsigned char), nbytes, file_stream);
  if (nbytes_read != nbytes) {
    std::free(buffer);
    std::fclose(file_stream);
    stop("Error reading BED data: nbytes_read != nbytes");
  }
  
  const double p = af[i];
  if (scale) {
    double denom = std::sqrt(2.0 * p * (1.0 - p));
    // guard against monomorphic/nearly monomorphic markers
    if (denom <= 0.0 || !std::isfinite(denom)) {
      map[0] = 0.0;
      map[1] = 0.0;
      map[2] = 0.0;
      map[3] = 0.0;
    } else {
      map[0] = (2.0 - 2.0 * p) / denom;
      map[1] = 0.0; // missing
      map[2] = (1.0 - 2.0 * p) / denom;
      map[3] = (-2.0 * p) / denom;
    }
  } else {
    map[0] = 2.0;
    map[1] = -2.0 * p; // missing imputed to mean
    map[2] = 1.0;
    map[3] = 0.0;
  }
  
  int j = 0;
  for (size_t k = 0; k < nbytes; ++k) {
    unsigned char buf_k = buffer[k];
    for (int pos = 0; pos < 4 && j < n; ++pos, ++j) {
      x[j] = map[buf_k & 3];
      buf_k >>= 2;
    }
  }
} // implicit barrier here

const double* bi = &b_snp[(size_t)i * nt];

// Parallelize across individuals.
// Each thread updates a disjoint chunk of grs_flat.
#pragma omp for schedule(static)
for (int j = 0; j < n; ++j) {
  const double xj = x[j];
  double* gj = &grs_flat[(size_t)j * nt];
  
  // update all traits for this individual
  for (int t = 0; t < nt; ++t) {
    gj[t] += bi[t] * xj;
  }
}
// implicit barrier here before next SNP
  }
}

std::free(buffer);
std::fclose(file_stream);

// Convert back to R-friendly trait-major nested vectors: grs[t][j]
std::vector<std::vector<double>> grs(nt, std::vector<double>(n));
for (int j = 0; j < n; ++j) {
  const double* gj = &grs_flat[(size_t)j * nt];
  for (int t = 0; t < nt; ++t) {
    grs[t][j] = gj[t];
  }
}

return grs;
}
