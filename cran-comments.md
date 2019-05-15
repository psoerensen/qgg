###########################################################################
## Test Platform (using devtools::rhub_check)
###########################################################################

* Windows Server 2008 R2 SP1, R-devel, 32/64 bit
* Fedora Linux, R-devel, clang, gfortran
* Debian Linux, R-devel, GCC ASAN/UBSAN

###########################################################################
## R CMD check results
###########################################################################

Notes 1 relates to new submission which is correct.

Notes 2 (which is a warning on the Windows server) relates to the use of 
Fortran I/O. 
We use Fortran stream based I/O to read/write genotype files. 
Currently a large number of really big cohort studies are becoming availabe
and therefore these files can be huge (>1000Gb). 

The qgg package provides an infrastructure for efficient processing of 
large-scale genetic and phenotypic data including core functions for:
* fitting linear mixed models
* construction of genomic relationship matrices
* estimating genetic parameters (heritability and correlation)
* genomic prediction
* single marker association analysis
* gene set enrichment analysis

qgg handles large-scale data by taking advantage of:
* multi-core processing using openMP
* multithreaded matrix operations implemented in BLAS libraries (e.g. OpenBLAS, ATLAS or MKL)
* fast and memory-efficient batch processing of genotype data stored in binary files (e.g. PLINK bedfiles)     


There were no ERRORs or WARNINGs or NOTES on Debian Linux 
##########################################################


There were no ERRORs or WARNINGs and 2 NOTES on Fedora Linux 
#############################################################

NOTES:

* checking CRAN incoming feasibility ... NOTE
Maintainer: ‘Peter Soerensen <pso@mbg.au.dk>’

New submission


* checking compiled code ... NOTE
File ‘qgg/libs/qgg.so’:
  Found ‘_gfortran_st_close’, possibly from ‘close’ (Fortran)
    Objects: ‘bedfuncs.o’, ‘bigreml.o’
  Found ‘_gfortran_st_open’, possibly from ‘open’ (Fortran)
    Objects: ‘bedfuncs.o’, ‘bigreml.o’
  Found ‘_gfortran_st_read’, possibly from ‘read’ (Fortran)
    Objects: ‘bedfuncs.o’, ‘bigreml.o’
  Found ‘_gfortran_st_write’, possibly from ‘write’ (Fortran), ‘print’
    (Fortran)
    Object: ‘bedfuncs.o’

Compiled code should not call entry points which might terminate R nor
write to stdout/stderr instead of to the console, nor use Fortran I/O
nor system RNGs.

See ‘Writing portable packages’ in the ‘Writing R Extensions’ manual.


There were no ERRORs, 1 WARNINGS, and 1 NOTES Windows Server 
#############################################################

WARNINGS: 

* checking compiled code ... WARNING
File 'qgg/libs/i386/qgg.dll':
  Found '_gfortran_st_close', possibly from 'close' (Fortran)
    Objects: 'bedfuncs.o', 'bigreml.o'
File 'qgg/libs/x64/qgg.dll':
  Found '_gfortran_st_close', possibly from 'close' (Fortran)
  Found '_gfortran_st_open', possibly from 'open' (Fortran)
  Found '_gfortran_st_open', possibly from 'open' (Fortran)
    Objects: 'bedfuncs.o', 'bigreml.o'
    Objects: 'bedfuncs.o', 'bigreml.o'
    Objects: 'bedfuncs.o', 'bigreml.o'

Compiled code should not call entry points which might terminate R nor
write to stdout/stderr instead of to the console, nor use Fortran I/O
nor system RNGs.

See 'Writing portable packages' in the 'Writing R Extensions' manual.

NOTES: 

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Peter Soerensen <pso@mbg.au.dk>'
New submission

