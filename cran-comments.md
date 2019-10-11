## Resubmission
This is a resubmission where I have:

* Corrected Date Field.

## Resubmission
This is a resubmission where I have:

* Corrected description text.
* Removed print statements in functions.
* Ensured that functions do not write in users home filespace by default.


## Resubmission
This is a resubmission where I have:

* Removed capital letters in description text.
* Added more details about the package functionality in the DESCRIPTION file including a reference.
* Unwrapped test examples and ensured they can be executed in less than 5 sec.
* Ensured that examples only use 1 core.
* Ensured that functions can use files not in working directory.
* Included a small data set in inst/extdata.

There was a note about:

  Possibly mis-spelled words in DESCRIPTION:
    al (14:14)
    et (14:11)
    heritability (13:281)
    qgg (13:18)
    Rohde (14:5)

I believe these are not mis-spelled.


## Resubmission
This is a resubmission. 

There were two notes. 

Note 1 relates to the package is a new submission which is correct

Note 2 was related to 'cran-comments.md' not being ignored. 

This has now been fixed by adding '^cran-comments\.md$' to the .Rbuildignore file.


###########################################################################
## Test Platform (using devtools::rhub_check)
###########################################################################

* Windows Server 2008 R2 SP1, R-devel, 32/64 bit
* Ubuntu Linux 16.04 LTS, R-release, GCC
* Fedora Linux, R-devel, clang, gfortran
* Debian Linux, R-devel, GCC ASAN/UBSAN

   https://builder.r-hub.io/status/qgg_1.0.0.tar.gz-6892dda9117e451c8040d51af6f1eceb
   https://builder.r-hub.io/status/qgg_1.0.0.tar.gz-a264ecf4dfa94a4d9457a479c08dfe4d
   https://builder.r-hub.io/status/qgg_1.0.0.tar.gz-8d880f8d18da448196a7d486e8a9759a
   https://builder.r-hub.io/status/qgg_1.0.0.tar.gz-8f009b45da004cd6abf269c662518f10

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

In addition to the examples included in the package submission we have prepared an extensive set of examples
avaliable on our github repository:

http://psoerensen.github.io/qgg/index.html
http://psoerensen.github.io/qgg/articles/qgg.html


###########################################################################
## R CMD check results
###########################################################################

On Windows and Ubuntu linux the following note appeared:

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Peter Soerensen <pso@mbg.au.dk>'

New submission


On Windows only the following note appeared:

* checking for non-standard things in the check directory ... NOTE
Found the following files/directories:
  'examples_i386' 'examples_x64' 'qgg-Ex_i386.Rout' 'qgg-Ex_x64.Rout'



Notes 1 relates to new submission which is correct.

Notes 2 only appeared on Windows and seems to be related to the r-hub builder.


