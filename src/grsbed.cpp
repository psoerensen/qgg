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

// // [[Rcpp::plugins(cpp11)]]
// // [[Rcpp::plugins(openmp)]]
// 
// #include <Rcpp.h>
// #include <cstdio>
// #include <vector>
// #include <cmath>
// #include <algorithm>
// 
// #ifdef _OPENMP
// #include <omp.h>
// #endif
// 
// using namespace Rcpp;
// 
// // [[Rcpp::export]]
// std::vector<std::vector<double>> mtgrsbed_omp(
//     const char* file,
//     int n,
//     const std::vector<int>& cls,                     // 1-based SNP indices in BED
//     const std::vector<double>& af,
//     bool scale,
//     const std::vector<std::vector<double>>& b,       // b[trait][snp]
//     int nthreads = 1
// ) {
//   FILE* file_stream = std::fopen(file, "rb");
//   if (!file_stream) {
//     stop("Could not open BED file.");
//   }
//   
//   const int nt = static_cast<int>(b.size());
//   if (nt == 0) {
//     std::fclose(file_stream);
//     return {};
//   }
// 
//   const int m = static_cast<int>(cls.size());
//   if ((int)af.size() != m) {
//     std::fclose(file_stream);
//     stop("af.size() != cls.size()");
//   }
// 
//   for (int t = 0; t < nt; ++t) {
//     if ((int)b[t].size() != m) {
//       std::fclose(file_stream);
//       stop("All rows of b must have length m = cls.size()");
//     }
//   }
// 
//   const size_t nbytes = (n + 3) / 4;
// 
//   // BED uses 2 bits/genotype:
//   // 00 01 10 11  ->  2, NA, 1, 0 copies of first allele in .bim
//   unsigned char* buffer = (unsigned char*) std::malloc(nbytes);
//   if (!buffer) {
//     std::fclose(file_stream);
//     stop("Failed to allocate BED buffer.");
//   }
// 
//   // Decoded genotype values for current SNP
//   std::vector<double> x(n);
// 
//   // Repack b into SNP-major layout: b_snp[i * nt + t]
//   // so all 500 trait weights for one SNP are contiguous.
//   std::vector<double> b_snp((size_t)m * nt);
//   for (int t = 0; t < nt; ++t) {
//     for (int i = 0; i < m; ++i) {
//       b_snp[(size_t)i * nt + t] = b[t][i];
//     }
//   }
// 
//   // Store scores individual-major: grs_flat[j * nt + t]
//   // so all 500 trait scores for one individual are contiguous.
//   std::vector<double> grs_flat((size_t)n * nt, 0.0);
// 
// #ifdef _OPENMP
//   if (nthreads > 0) omp_set_num_threads(nthreads);
// #endif
// 
// #pragma omp parallel
// {
//   std::vector<double> map(4);
//   
//   for (int i = 0; i < m; ++i) {
//     
//     // One thread reads + decodes current SNP into x[]
// #pragma omp single
// {
//   long long offset = (long long)(cls[i] - 1) * (long long)nbytes + 3LL;
//   
//   if (std::fseek(file_stream, offset, SEEK_SET) != 0) {
//     std::free(buffer);
//     std::fclose(file_stream);
//     stop("fseek failed.");
//   }
//   
//   size_t nbytes_read = std::fread(buffer, sizeof(unsigned char), nbytes, file_stream);
//   if (nbytes_read != nbytes) {
//     std::free(buffer);
//     std::fclose(file_stream);
//     stop("Error reading BED data: nbytes_read != nbytes");
//   }
//   
//   const double p = af[i];
//   if (scale) {
//     double denom = std::sqrt(2.0 * p * (1.0 - p));
//     // guard against monomorphic/nearly monomorphic markers
//     if (denom <= 0.0 || !std::isfinite(denom)) {
//       map[0] = 0.0;
//       map[1] = 0.0;
//       map[2] = 0.0;
//       map[3] = 0.0;
//     } else {
//       map[0] = (2.0 - 2.0 * p) / denom;
//       map[1] = 0.0; // missing
//       map[2] = (1.0 - 2.0 * p) / denom;
//       map[3] = (-2.0 * p) / denom;
//     }
//   } else {
//     map[0] = 2.0;
//     map[1] = -2.0 * p; // missing imputed to mean
//     map[2] = 1.0;
//     map[3] = 0.0;
//   }
//   
//   int j = 0;
//   for (size_t k = 0; k < nbytes; ++k) {
//     unsigned char buf_k = buffer[k];
//     for (int pos = 0; pos < 4 && j < n; ++pos, ++j) {
//       x[j] = map[buf_k & 3];
//       buf_k >>= 2;
//     }
//   }
// } // implicit barrier here
// 
// const double* bi = &b_snp[(size_t)i * nt];
// 
// // Parallelize across individuals.
// // Each thread updates a disjoint chunk of grs_flat.
// #pragma omp for schedule(static)
// for (int j = 0; j < n; ++j) {
//   const double xj = x[j];
//   double* gj = &grs_flat[(size_t)j * nt];
//   
//   // update all traits for this individual
//   for (int t = 0; t < nt; ++t) {
//     gj[t] += bi[t] * xj;
//   }
// }
// // implicit barrier here before next SNP
//   }
// }
// 
// std::free(buffer);
// std::fclose(file_stream);
// 
// // Convert back to R-friendly trait-major nested vectors: grs[t][j]
// std::vector<std::vector<double>> grs(nt, std::vector<double>(n));
// for (int j = 0; j < n; ++j) {
//   const double* gj = &grs_flat[(size_t)j * nt];
//   for (int t = 0; t < nt; ++t) {
//     grs[t][j] = gj[t];
//   }
// }
// 
// return grs;
// }


// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(openmp)]]

#include <Rcpp.h>
#include <cstdio>
#include <vector>
#include <cmath>
#include <algorithm>
#include <stdexcept>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

// -----------------------------------------------------------------------------
// Pure C++ core
// b_snp is SNP-major: b_snp[i * nt + t]
// grs_flat is individual-major: grs_flat[j * nt + t]
// -----------------------------------------------------------------------------
// void mtgrsbed_core(
//     const char* file,
//     int n,
//     const int* cls,
//     const double* af,
//     int m,
//     int nt,
//     bool scale,
//     const double* b_snp,
//     double* grs_flat,
//     int nthreads
// ) {
//   FILE* file_stream = std::fopen(file, "rb");
//   if (!file_stream) {
//     throw std::runtime_error("Could not open BED file.");
//   }
//   
//   const size_t nbytes = (n + 3) / 4;
//   
//   unsigned char* buffer = (unsigned char*) std::malloc(nbytes);
//   if (!buffer) {
//     std::fclose(file_stream);
//     throw std::runtime_error("Failed to allocate BED buffer.");
//   }
//   
//   std::vector<double> x(n);
//   
// #ifdef _OPENMP
//   if (nthreads > 0) omp_set_num_threads(nthreads);
// #endif
//   
// #pragma omp parallel
// {
//   std::vector<double> map(4);
//   
//   for (int i = 0; i < m; ++i) {
//     
// #pragma omp single
// {
//   long long offset = (long long)(cls[i] - 1) * (long long)nbytes + 3LL;
//   
//   if (std::fseek(file_stream, offset, SEEK_SET) != 0) {
//     std::free(buffer);
//     std::fclose(file_stream);
//     throw std::runtime_error("fseek failed.");
//   }
//   
//   size_t nbytes_read = std::fread(buffer, sizeof(unsigned char), nbytes, file_stream);
//   if (nbytes_read != nbytes) {
//     std::free(buffer);
//     std::fclose(file_stream);
//     throw std::runtime_error("Error reading BED data: nbytes_read != nbytes.");
//   }
//   
//   const double p = af[i];
//   if (scale) {
//     const double denom = std::sqrt(2.0 * p * (1.0 - p));
//     if (denom <= 0.0 || !std::isfinite(denom)) {
//       map[0] = 0.0;
//       map[1] = 0.0;
//       map[2] = 0.0;
//       map[3] = 0.0;
//     } else {
//       map[0] = (2.0 - 2.0 * p) / denom;
//       map[1] = 0.0;
//       map[2] = (1.0 - 2.0 * p) / denom;
//       map[3] = (-2.0 * p) / denom;
//     }
//   } else {
//     map[0] = 2.0;
//     map[1] = -2.0 * p;
//     map[2] = 1.0;
//     map[3] = 0.0;
//   }
//   
//   int j = 0;
//   for (size_t k = 0; k < nbytes; ++k) {
//     unsigned char buf_k = buffer[k];
//     for (int pos = 0; pos < 4 && j < n; ++pos, ++j) {
//       x[j] = map[buf_k & 3];
//       buf_k >>= 2;
//     }
//   }
// }
// 
// const double* bi = &b_snp[(size_t)i * nt];
// 
// #pragma omp for schedule(static)
// for (int j = 0; j < n; ++j) {
//   const double xj = x[j];
//   double* gj = &grs_flat[(size_t)j * nt];
//   
// #pragma omp simd
//   for (int t = 0; t < nt; ++t) {
//     gj[t] += bi[t] * xj;
//   }
// }
//   }
// }
// 
// std::free(buffer);
// std::fclose(file_stream);
// }

void mtgrsbed_core(
    const char* file,
    int n,
    const int* cls,
    const double* af,
    int m,
    int nt,
    bool scale,
    const double* b_snp,
    double* grs_flat,
    int nthreads,
    int MG = 64,
    int JB = 2048,
    int TB = 64
) {
  FILE* file_stream = std::fopen(file, "rb");
  if (!file_stream) {
    throw std::runtime_error("Could not open BED file.");
  }
  
  const size_t nbytes = (n + 3) / 4;
  
  // Raw BED bytes for a block of MG markers
  unsigned char* block_buffer =
    (unsigned char*) std::malloc((size_t)MG * nbytes);
  if (!block_buffer) {
    std::fclose(file_stream);
    throw std::runtime_error("Failed to allocate BED block buffer.");
  }
  
  // Shared maps for current marker block
  std::vector<double> map0(MG), map1(MG), map2(MG), map3(MG);
  
#ifdef _OPENMP
  if (nthreads > 0) omp_set_num_threads(nthreads);
#endif
  
#pragma omp parallel
{
  for (int i0 = 0; i0 < m; i0 += MG) {
    const int imax = std::min(i0 + MG, m);
    const int mlen = imax - i0;
    
#pragma omp single
{
  for (int ii = 0; ii < mlen; ++ii) {
    const int i = i0 + ii;
    
    const long long offset =
      (long long)(cls[i] - 1) * (long long)nbytes + 3LL;
    
    if (std::fseek(file_stream, offset, SEEK_SET) != 0) {
      std::free(block_buffer);
      std::fclose(file_stream);
      throw std::runtime_error("fseek failed.");
    }
    
    unsigned char* buf_i = block_buffer + (size_t)ii * nbytes;
    
    const size_t nbytes_read =
      std::fread(buf_i, sizeof(unsigned char), nbytes, file_stream);
    
    if (nbytes_read != nbytes) {
      std::free(block_buffer);
      std::fclose(file_stream);
      throw std::runtime_error("Error reading BED data: nbytes_read != nbytes.");
    }
    
    const double p = af[i];
    if (scale) {
      const double denom = std::sqrt(2.0 * p * (1.0 - p));
      if (denom <= 0.0 || !std::isfinite(denom)) {
        map0[ii] = 0.0;
        map1[ii] = 0.0;
        map2[ii] = 0.0;
        map3[ii] = 0.0;
      } else {
        map0[ii] = (2.0 - 2.0 * p) / denom;
        map1[ii] = 0.0;                    // missing
        map2[ii] = (1.0 - 2.0 * p) / denom;
        map3[ii] = (-2.0 * p) / denom;
      }
    } else {
      map0[ii] = 2.0;
      map1[ii] = -2.0 * p;                // mean-imputed centered coding
      map2[ii] = 1.0;
      map3[ii] = 0.0;
    }
  }
} // implicit barrier here

#pragma omp for schedule(static)
for (int j0 = 0; j0 < n; j0 += JB) {
  const int jmax = std::min(j0 + JB, n);
  
  const int byte0 = j0 >> 2;
  const int byte1 = (jmax + 3) >> 2;
  
  for (int t0 = 0; t0 < nt; t0 += TB) {
    const int tmax = std::min(t0 + TB, nt);
    const int tlen = tmax - t0;
    
    // Loop over markers in current block
    for (int ii = 0; ii < mlen; ++ii) {
      const int i = i0 + ii;
      const unsigned char* buf_i = block_buffer + (size_t)ii * nbytes;
      const double* bi = &b_snp[(size_t)i * nt + t0];
      
      const double m0 = map0[ii];
      const double m1 = map1[ii];
      const double m2 = map2[ii];
      const double m3 = map3[ii];
      
      for (int kb = byte0; kb < byte1; ++kb) {
        unsigned char buf_k = buf_i[kb];
        const int jbase = kb << 2;
        
        for (int pos = 0; pos < 4; ++pos) {
          const int j = jbase + pos;
          if (j < j0 || j >= jmax || j >= n) {
            buf_k >>= 2;
            continue;
          }
          
          double xj;
          switch (buf_k & 3u) {
          case 0u: xj = m0; break;   // BED 00 -> 2 copies
          case 1u: xj = m1; break;   // BED 01 -> missing
          case 2u: xj = m2; break;   // BED 10 -> 1 copy
          default: xj = m3; break;   // BED 11 -> 0 copies
          }
          buf_k >>= 2;
          
          double* gj = &grs_flat[(size_t)j * nt + t0];
          
#pragma omp simd
          for (int t = 0; t < tlen; ++t) {
            gj[t] += bi[t] * xj;
          }
        }
      }
    }
  }
}
  }
}

std::free(block_buffer);
std::fclose(file_stream);
}

void mtgrsbed_core_fast01(
    const char* file,
    int n,
    const int* cls,
    int m,
    int nt,
    const double* b_snp,
    double* grs_flat,
    int nthreads
) {
  FILE* file_stream = std::fopen(file, "rb");
  if (!file_stream) {
    throw std::runtime_error("Could not open BED file.");
  }
  
  const size_t nbytes = (n + 3) / 4;
  unsigned char* buffer = (unsigned char*) std::malloc(nbytes);
  if (!buffer) {
    std::fclose(file_stream);
    throw std::runtime_error("Failed to allocate BED buffer.");
  }
  
  const int JB = 2048;
  const int TB = 64;
  
#ifdef _OPENMP
  if (nthreads > 0) omp_set_num_threads(nthreads);
#endif
  
#pragma omp parallel
{
  for (int i = 0; i < m; ++i) {
    
#pragma omp single
{
  const long long offset =
    (long long)(cls[i] - 1) * (long long)nbytes + 3LL;
  
  if (std::fseek(file_stream, offset, SEEK_SET) != 0) {
    std::free(buffer);
    std::fclose(file_stream);
    throw std::runtime_error("fseek failed.");
  }
  
  const size_t nbytes_read =
    std::fread(buffer, sizeof(unsigned char), nbytes, file_stream);
  
  if (nbytes_read != nbytes) {
    std::free(buffer);
    std::fclose(file_stream);
    throw std::runtime_error("Error reading BED data: nbytes_read != nbytes.");
  }
}

const double* bi = &b_snp[(size_t)i * nt];

#pragma omp for schedule(static)
for (int j0 = 0; j0 < n; j0 += JB) {
  const int jmax = std::min(j0 + JB, n);
  
  const int byte0 = j0 >> 2;
  const int byte1 = (jmax + 3) >> 2;
  
  for (int t0 = 0; t0 < nt; t0 += TB) {
    const int tmax = std::min(t0 + TB, nt);
    const int tlen = tmax - t0;
    const double* bij = bi + t0;
    
    for (int kb = byte0; kb < byte1; ++kb) {
      unsigned char buf_k = buffer[kb];
      const int jbase = kb << 2;
      
      for (int pos = 0; pos < 4; ++pos) {
        const int j = jbase + pos;
        if (j < j0 || j >= jmax || j >= n) {
          buf_k >>= 2;
          continue;
        }
        
        const unsigned code = (unsigned)(buf_k & 3u);
        buf_k >>= 2;
        
        // BED:
        // 00 -> 2
        // 01 -> missing
        // 10 -> 1
        // 11 -> 0
        
        if (code == 3u) {
          continue; // genotype 0, skip
        }
        
        double* gj = &grs_flat[(size_t)j * nt + t0];
        
        if (code == 2u) {
          // genotype 1
#pragma omp simd
          for (int t = 0; t < tlen; ++t) {
            gj[t] += bij[t];
          }
        } else if (code == 0u) {
          // genotype 2
#pragma omp simd
          for (int t = 0; t < tlen; ++t) {
            gj[t] += 2.0 * bij[t];
          }
        } else {
          // missing (01)
          // if you want mean imputation on raw scale, add 2p * b here.
          // for the transformed-weight trick, it's often simplest to leave
          // this as zero and handle consistency separately if needed.
        }
      }
    }
  }
}
  }
}

std::free(buffer);
std::fclose(file_stream);
}


// -----------------------------------------------------------------------------
// Thin R wrapper
// Accepts b as list of trait-vectors: b[[t]][i]
// Returns list of trait-vectors: grs[[t]][j]
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
std::vector<std::vector<double>> mtgrsbed_omp(
    std::string file,
    int n,
    const std::vector<int>& cls,
    const std::vector<double>& af,
    bool scale,
    const std::vector<std::vector<double>>& b,
    int nthreads = 1
) {
  const int nt = static_cast<int>(b.size());
  if (nt == 0) {
    return {};
  }
  
  const int m = static_cast<int>(cls.size());
  if ((int)af.size() != m) {
    stop("af.size() != cls.size()");
  }
  
  for (int t = 0; t < nt; ++t) {
    if ((int)b[t].size() != m) {
      stop("All rows of b must have length m = cls.size()");
    }
  }
  
  // Repack to SNP-major contiguous layout
  std::vector<double> b_snp((size_t)m * nt);
  for (int t = 0; t < nt; ++t) {
    for (int i = 0; i < m; ++i) {
      b_snp[(size_t)i * nt + t] = b[t][i];
    }
  }
  
  // Output as individual-major flat array
  std::vector<double> grs_flat((size_t)n * nt, 0.0);
  
  try {
    mtgrsbed_core(
      file.c_str(),
      n,
      cls.data(),
      af.data(),
      m,
      nt,
      scale,
      b_snp.data(),
      grs_flat.data(),
      nthreads
    );
  } catch (const std::exception& e) {
    stop(e.what());
  }
  
  // Convert back to current R-facing format: list of trait vectors
  std::vector<std::vector<double>> grs(nt, std::vector<double>(n));
  for (int j = 0; j < n; ++j) {
    const double* gj = &grs_flat[(size_t)j * nt];
    for (int t = 0; t < nt; ++t) {
      grs[t][j] = gj[t];
    }
  }
  
  return grs;
}


// [[Rcpp::export]]
Rcpp::NumericMatrix mtgrsbed_matrix(
    std::string file,
    int n,
    const std::vector<int>& cls,
    const std::vector<double>& af,
    bool scale,
    Rcpp::NumericMatrix S,   // m x nt (SNP x trait)
    int nthreads = 1,
    int MG = 64,             // marker block size
    int JB = 2048,           // individual block size
    int TB = 32              // trait block size
) {
  const int m  = S.nrow();   // SNPs
  const int nt = S.ncol();   // traits
  
  if ((int)cls.size() != m) {
    Rcpp::stop("cls.size() must match number of rows in S");
  }
  
  if ((int)af.size() != m) {
    Rcpp::stop("af.size() must match number of rows in S");
  }
  
  // Pointer to R matrix (column-major: S[i + t*m])
  const double* s = REAL(S);
  
  // Repack to SNP-major layout: b_snp[i*nt + t]
  std::vector<double> b_snp((size_t)m * nt);
  
  for (int t = 0; t < nt; ++t) {
    for (int i = 0; i < m; ++i) {
      b_snp[(size_t)i * nt + t] = s[i + (size_t)t * m];
    }
  }
  
  // Allocate output (individual-major flat)
  std::vector<double> grs_flat((size_t)n * nt, 0.0);
  
  try {
    mtgrsbed_core(
      file.c_str(),
      n,
      cls.data(),
      af.data(),
      m,
      nt,
      scale,
      b_snp.data(),
      grs_flat.data(),
      nthreads,
      MG,
      JB,
      TB
    );
  } catch (const std::exception& e) {
    Rcpp::stop(e.what());
  }
  
  // Convert to R matrix (n x nt)
  Rcpp::NumericMatrix out(n, nt);
  
  for (int t = 0; t < nt; ++t) {
    for (int j = 0; j < n; ++j) {
      out(j, t) = grs_flat[(size_t)j * nt + t];
    }
  }
  
  return out;
}


// [[Rcpp::export]]
Rcpp::NumericMatrix mtgrsbed_matrix_fast01(
    std::string file,
    int n,
    const std::vector<int>& cls,
    const std::vector<double>& af,
    bool scale,
    Rcpp::NumericMatrix S,
    int nthreads = 1
) {
  if (scale) {
    Rcpp::stop("mtgrsbed_matrix_fast01 requires scale = FALSE");
  }
  
  const int m  = S.nrow();
  const int nt = S.ncol();
  
  if ((int)cls.size() != m) {
    Rcpp::stop("cls.size() must match number of rows in S");
  }
  
  if ((int)af.size() != m) {
    Rcpp::stop("af.size() must match number of rows in S");
  }
  
  // pointer to matrix
  const double* s = REAL(S);
  
  // convert to SNP-major layout
  std::vector<double> b_snp((size_t)m * nt);
  
  for (int t = 0; t < nt; ++t) {
    for (int i = 0; i < m; ++i) {
      b_snp[(size_t)i * nt + t] = s[i + (size_t)t * m];
    }
  }
  
  // output buffer
  std::vector<double> grs_flat((size_t)n * nt, 0.0);
  
  try {
    mtgrsbed_core_fast01(
      file.c_str(),
      n,
      cls.data(),
      m,
      nt,
      b_snp.data(),
      grs_flat.data(),
      nthreads
    );
  } catch (const std::exception& e) {
    Rcpp::stop(e.what());
  }
  
  // return matrix
  Rcpp::NumericMatrix out(n, nt);
  
  for (int t = 0; t < nt; ++t) {
    for (int j = 0; j < n; ++j) {
      out(j, t) = grs_flat[(size_t)j * nt + t];
    }
  }
  
  return out;
}


  
  