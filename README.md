
<!-- README.md is generated from README.Rmd. Please edit that file -->

### An R package for Quantitative Genetic and Genomic analyses

The **qgg** package was developed based on the hypothesis that certain
regions on the genome, so-called *genomic features*, may be enriched for
causal variants affecting the trait. Several genomic feature classes can
be formed based on previous studies and different sources of information
such as genes, chromosomes or biological pathways.

**qgg** provides an infrastructure for efficient processing of
large-scale genetic and phenotypic data including core functions for:

- fitting linear mixed models  
- construction of genomic relationship matrices  
- estimating genetic parameters (heritability and correlation)  
- genomic prediction  
- single marker association analysis  
- gene set enrichment analysis

**qgg** handles large-scale data by taking advantage of:

- multi-core processing using [openMP](https://www.openmp.org/)  
- multithreaded matrix operations implemented in BLAS libraries
  (e.g. [OpenBLAS](https://www.openblas.net/),
  [ATLAS](https://math-atlas.sourceforge.net/) or
  [MKL](https://en.wikipedia.org/wiki/Math_Kernel_Library))  
- fast and memory-efficient batch processing of genotype data stored in
  binary files (e.g. [PLINK](https://www.cog-genomics.org/plink2)
  bedfiles)

The **qgg** package provides a range of genomic feature modeling
approaches, including genomic feature best linear unbiased prediction
(GFBLUP) models, implemented using likelihood or Bayesian methods.
Multiple features and multiple traits can be included in these models
and different genetic models (e.g. additive, dominance, gene by gene and
gene by environment interactions) can be used. Further extensions
include a weighted GFBLUP model using differential weighting of the
individual genetic marker relationships. Marker set tests, which are
computationally very fast, can be performed. These marker set tests
allow the rapid analyses of different layers of genomic feature classes
to discover genomic features potentially enriched for causal variants.
Marker set tests can thus facilitate more accurate prediction models.

### Install

You can install qgg from CRAN with:

``` r
install.packages("qgg")
```

The most recent version of `qgg` can be obtained from github:

``` r
library(devtools)
devtools::install_github("psoerensen/qgg")
```

### Tutorials

Below is a set of tutorials used for the qgg package:

This tutorial provides a brief introduction to R package qgg using small
simulated data examples.  
[Practicals_brief_introduction](https://psoerensen.github.io/qgtutorials/Quick-tutorials-for-qgg-package.pdf)

This tutorial provides an introduction to R package qgg using 1000G
data.  
[Practicals_1000G_tutorials](https://psoerensen.github.io/qgtutorials/1000G-tutorials-for-qgg-package.pdf)

This tutorial provide a simple introduction to polygenic risk scoring
(PRS) of complex traits and diseases using simulated data. The practical
will be a mix of theoretical and practical exercises in R that are used
for illustrating/applying the theory presented in the corresponding
lecture notes on polygenic risk scoring.  
[Practicals_human_example](https://psoerensen.github.io/qgtutorials/Practicals_human_example.pdf)

In this tutorial we will be analysing quantitative traits observed in a
mice population. The mouse data consist of phenotypes for traits related
to growth and obesity (e.g. body weight, glucose levels in blood),
pedigree information, and genetic marker data.  
[Practicals_mouse_example](https://psoerensen.github.io/qgtutorials/Practicals_mouse_example.pdf)

### Notes

Below is a set of notes for the quantitative genetic theory, statistical
models and methods implemented in the qgg package:

[Quantitative Genetics
Theory](https://psoerensen.github.io/qgnotes/Quantitative-Genetics-Theory.pdf)

[Estimation of Genetic
Predisposition](https://psoerensen.github.io/qgnotes/Estimation-of-Genetic-Predisposition.pdf)

[Estimation of Genetic
Parameters](https://psoerensen.github.io/qgnotes/Estimation-of-Genetic-Parameters.pdf)

[Linear Mixed Models](https://psoerensen.github.io/qgnotes/LMM.pdf)

[Best Linear Unbiased Prediction
Models](https://psoerensen.github.io/qgnotes/BLUP.pdf)

[REstricted Maximum Likelihood
Methods](https://psoerensen.github.io/qgnotes/REML.pdf)

[Gene Set Enrichment
Analysis](https://psoerensen.github.io/qgnotes/GSEA.pdf)

[Bayesian Linear Regression
Models](https://psoerensen.github.io/qgnotes/BLR.pdf)

#### References

1.  Edwards SM, Thomsen B, Madsen P, Sørensen P. 2015. Partitioning of
    genomic variance reveals biological pathways associated with udder
    health and milk production traits in dairy cattle. *Genet Sel Evol*
    47:60. <doi:10.1186/s12711-015-0132-6>  
2.  Edwards SM, Sørensen IF, Sarup P, Mackay TFC, Sørensen P. 2016.
    Genomic prediction for quantitative traits is improved by mapping
    variants to gene ontology categories in *Drosophila melanogaster*.
    *Genetics* 203:1871–1883. <doi:10.1534/genetics.116.187161>
3.  Ehsani A, Janss L, Pomp D, Sørensen P. 2015. Decomposing genomic
    variance using information from GWA, GWE and eQTL analysis. *Anim
    Genet* 47:165–173. <doi:10.1111/age.12396>
4.  Fang L, Sahana G, Ma P, Su G, Yu Y, Zhang S, Lund MS,
    Sørensen P. 2017. Exploring the genetic architecture and improving
    genomic prediction accuracy for mastitis and milk production traits
    in dairy cattle by mapping variants to hepatic transcriptomic
    regions responsive to intra-mammary infection. *Genet Sel Evol*
    49:1–18. <doi:10.1186/s12711-017-0319-0>
5.  Fang L, Sahana G, Su G, Yu Y, Zhang S, Lund MS, Sørensen P. 2017.
    Integrating sequence-based GWAS and RNA-seq provides novel insights
    into the genetic basis of mastitis and milk production in dairy
    cattle. *Sci Rep* 7:45560. <doi:10.1038/srep45560>
6.  Fang L, Sørensen P, Sahana G, Panitz F, Su G, Zhang S, Yu Y, Li B,
    Ma L, Liu G, Lund MS, Thomsen B. 2018. MicroRNA-guided
    prioritization of genome-wide association signals reveals the
    importance of microRNA-target gene networks for complex traits in
    cattle. *Sci Rep* 8:1–14. <doi:10.1038/s41598-018-27729-y>
7.  Ørsted M, Rohde PD, Hoffmann AA, Sørensen P, Kristensen TN. 2017.
    Environmental variation partitioned into separate heritable
    components. *Evolution* (N Y) 72:136–152. <doi:10.1111/evo.13391>
8.  Ørsted M, Hoffmann AA, Rohde PD, Sørensen P, Kristensen TN. 2018.
    Strong impact of thermal environment on the quantitative genetic
    basis of a key stress tolerance trait. *Heredity* (Edinb).
    <doi:10.1038/s41437-018-0117-7>
9.  Rohde PD, Krag K, Loeschcke V, Overgaard J, Sørensen P, Kristensen
    TN. 2016. A quantitative genomic approach for analysis of fitness
    and stress related traits in a *Drosophila melanogaster model
    population*. *Int J Genomics* 2016:1–11.
10. Rohde PD, Demontis D, Cuyabano BCD, The GEMS Group, Børglum AD,
    Sørensen P. 2016. Covariance Association Test (CVAT) identify
    genetic markers associated with schizophrenia in functionally
    associated biological processes. *Genetics* 203:1901–1913.
    <doi:10.1534/genetics.116.189498>
11. Rohde PD, Gaertner B, Ward K, Sørensen P, Mackay TFC. 2017. Genomic
    analysis of genotype-by-social environment interaction for
    *Drosophila melanogaster*. *Genetics* 206:1969–1984.
    <doi:10.1534/genetics.117.200642/-/DC1.1>
12. Rohde PD, Østergaard S, Kristensen TN, Sørensen P, Loeschcke V,
    Mackay TFC, Sarup P. 2018. Functional validation of candidate genes
    detected by genomic feature models. *G3 Genes, Genomes, Genet*
    8:1659–1668. <doi:10.1534/g3.118.200082>
13. Sarup P, Jensen J, Ostersen T, Henryon M, Sørensen P. 2016.
    Increased prediction accuracy using a genomic feature model
    including prior information on quantitative trait locus regions in
    purebred Danish Duroc pigs. *BMC Genet* 17:11.
    <doi:10.1186/s12863-015-0322-9>
14. Sørensen P, de los Campos G, Morgante F, Mackay TFC,
    Sorensen D. 2015. Genetic control of environmental variation of two
    quantitative traits of *Drosophila melanogaster* revealed by
    whole-genome sequencing. *Genetics* 201:487–497.
    <doi:10.1534/genetics.115.180273>
15. Sørensen IF, Edwards SM, Rohde PD, Sørensen P. 2017. Multiple trait
    covariance association test identifies gene ontology categories
    associated with chill coma recovery time in *Drosophila
    melanogaster*. *Sci Rep* 7:2413. <doi:10.1038/s41598-017-02281-3>
