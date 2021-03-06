## Submission
This is a new submission where:

* Fixed issue with module name multiple defined detected on solaris 
* Check for CRAN was OK except for the window build related to data.table not available but assume this is a temporary error.
* Check on solaris was OK.

Results for check_on_solaris:
https://builder.r-hub.io/status/qgg_1.0.4.tar.gz-032fd7d9e906414cb679bc79d54e66f7

Results for check_for_cran():
   https://builder.r-hub.io/status/qgg_1.0.4.tar.gz-23b75abdfecd459ea0141e9e58c65887
   https://builder.r-hub.io/status/qgg_1.0.4.tar.gz-6a5bbd67917c468d802420182af4ce7d
   https://builder.r-hub.io/status/qgg_1.0.4.tar.gz-187feeb12c2d45228f8b0c35317053e0
   https://builder.r-hub.io/status/qgg_1.0.4.tar.gz-9be5d8e78be94214a5d8bc00ce476cf4

## Submission
This is a new submission where:

* Fixed issue with LTO (in fortran 2 c interface) detected on solaris 
* Added a few user requested features.
* Check for CRAN was OK except for the window build related to data.table not available but assume this is a temporary error.

 
Results for check_for_cran():
   https://builder.r-hub.io/status/qgg_1.0.4.tar.gz-5c9f861803a94d4ab7bf4ed23efc2ec3
   https://builder.r-hub.io/status/qgg_1.0.4.tar.gz-c3230d7bb1064c2b809cb31f08e3a067
   https://builder.r-hub.io/status/qgg_1.0.4.tar.gz-b7914ccce695408a9dc01e3ec7682543
   https://builder.r-hub.io/status/qgg_1.0.4.tar.gz-a31493d3a3f8461a90274493d12c78f9


## Submission
This is a new submission where:

* Dependencies on openmp is removed as required by CRAN.
* There was a note about misspelling, but checked that spelling is correct.
* There was a note that package was archived on CRAN due to the issues with openmp not being corrected in time.

 
Results for check_for_cran():

   https://builder.r-hub.io/status/qgg_1.0.3.tar.gz-bfc21379cb454b4891678639fe1de18c
   https://builder.r-hub.io/status/qgg_1.0.3.tar.gz-44c5b58cf6b145baac2708f2ccef234d
   https://builder.r-hub.io/status/qgg_1.0.3.tar.gz-6040ed0410ac47f2824f272a5caa67e6
   https://builder.r-hub.io/status/qgg_1.0.3.tar.gz-b4c0e33aba66401fac860801337ae75e

## Submission
This is a new submission where I have fixed the following issues encountered on CRAN checks:

* Fixed additional issues related to build errors on fedora and solaris.
* There was a note about misspelling, but checked that spelling is correct.
* There was a note that package was archived on CRAN due to "reported installation issued was ignored in update".
* Clearly last updates did not fix the installation issues. 
* New updates appear to have fixed previous installation issued which can be verified in the check links provided below.  
 
Results for check_for_cran():
   https://builder.r-hub.io/status/qgg_1.0.2.tar.gz-1dcb3475744a4294b8114d6aa8936b98
   https://builder.r-hub.io/status/qgg_1.0.2.tar.gz-057edb786f4a4d6dbb9f449d9492dc9d
   https://builder.r-hub.io/status/qgg_1.0.2.tar.gz-37b08a62e5474ce6be940589d6f0ca1b
   https://builder.r-hub.io/status/qgg_1.0.2.tar.gz-5fe973bab2564cf39172fbded9fc5cc7

Results for check_on_solaris():
   https://builder.r-hub.io/status/qgg_1.0.2.tar.gz-6759129787ba448880e6aa81de6044dc

Results for check_on_fedora():
   https://builder.r-hub.io/status/qgg_1.0.2.tar.gz-36362d7179774f43b766b48b5f2eb752


## Submission
This is a new submission where I have fixed the following issues encountered on CRAN checks:

* Fixed issues related to using fortran modules when building in parallel.
* Removed unused variables in fortran code.
* Fixed issue with a type mismatch in fortran code (LTO).
* The revised version passes check_for_cran() checks on all platforms.
* There was a note about misspelling, but checked that spelling is correct. 


## Resubmission
This is a resubmission where I have:

* Changed maintainer email address.
* Fixed minor bugs related to multi-core processing of bed files.
* Removed knitr/rmarkdown from suggest.


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


