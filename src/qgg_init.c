#include <R_ext/RS.h>
#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>


/* tools::package_native_routine_registration_skeleton('C://Users//au223366//Documents//GitHub//qgg',,,FALSE)
*/


/* .Call calls */
extern SEXP _qgg_bayes(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _qgg_cp(SEXP);
extern SEXP _qgg_freqbed(SEXP, SEXP, SEXP, SEXP);
extern SEXP _qgg_readW(SEXP, SEXP, SEXP, SEXP);
extern SEXP _qgg_getWlist(SEXP, SEXP, SEXP, SEXP);
extern SEXP _qgg_grsbed(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _qgg_mmult(SEXP, SEXP);
extern SEXP _qgg_mtbayes(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _qgg_mtsbayes(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _qgg_mtgrsbed(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _qgg_mtsolvebed(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _qgg_mvrnorm(SEXP);
extern SEXP _qgg_pruneld(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _qgg_pruneldmat(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _qgg_psets(SEXP, SEXP, SEXP, SEXP);
extern SEXP _qgg_rcpparma_bothproducts(SEXP);
extern SEXP _qgg_rcpparma_hello_world();
extern SEXP _qgg_rcpparma_innerproduct(SEXP);
extern SEXP _qgg_rcpparma_outerproduct(SEXP);
extern SEXP _qgg_readG(SEXP, SEXP, SEXP);
extern SEXP _qgg_riwishart(SEXP, SEXP);
extern SEXP _qgg_rwishart(SEXP, SEXP);
extern SEXP _qgg_sbayes(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _qgg_sbayes_spa(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _qgg_solvebed(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _qgg_summarybed(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);


static const R_CallMethodDef CallEntries[] = {
    {"_qgg_bayes",                 (DL_FUNC) &_qgg_bayes,                 17},
    {"_qgg_cp",                    (DL_FUNC) &_qgg_cp,                     1},
    {"_qgg_freqbed",               (DL_FUNC) &_qgg_freqbed,                4},
    {"_qgg_readW",                 (DL_FUNC) &_qgg_readW,                  4},
    {"_qgg_getWlist",              (DL_FUNC) &_qgg_getWlist,               4},
    {"_qgg_grsbed",                (DL_FUNC) &_qgg_grsbed,                 5},
    {"_qgg_mmult",                 (DL_FUNC) &_qgg_mmult,                  2},
    {"_qgg_mtbayes",               (DL_FUNC) &_qgg_mtbayes,               16},
    {"_qgg_mtsbayes",               (DL_FUNC) &_qgg_mtsbayes,               19},
    {"_qgg_mtgrsbed",              (DL_FUNC) &_qgg_mtgrsbed,               6},
    {"_qgg_mtsolvebed",            (DL_FUNC) &_qgg_mtsolvebed,             8},
    {"_qgg_mvrnorm",               (DL_FUNC) &_qgg_mvrnorm,                1},
    {"_qgg_pruneld",               (DL_FUNC) &_qgg_pruneld,                6},
    {"_qgg_pruneldmat",            (DL_FUNC) &_qgg_pruneldmat,             5},
    {"_qgg_psets",                 (DL_FUNC) &_qgg_psets,                  4},
    {"_qgg_rcpparma_bothproducts", (DL_FUNC) &_qgg_rcpparma_bothproducts,  1},
    {"_qgg_rcpparma_hello_world",  (DL_FUNC) &_qgg_rcpparma_hello_world,   0},
    {"_qgg_rcpparma_innerproduct", (DL_FUNC) &_qgg_rcpparma_innerproduct,  1},
    {"_qgg_rcpparma_outerproduct", (DL_FUNC) &_qgg_rcpparma_outerproduct,  1},
    {"_qgg_readG",               (DL_FUNC) &_qgg_readG,                3},
    {"_qgg_riwishart",             (DL_FUNC) &_qgg_riwishart,              2},
    {"_qgg_rwishart",              (DL_FUNC) &_qgg_rwishart,               2},
    {"_qgg_sbayes",                (DL_FUNC) &_qgg_sbayes,                19},
    {"_qgg_sbayes_spa",                (DL_FUNC) &_qgg_sbayes_spa,                21},
    {"_qgg_solvebed",              (DL_FUNC) &_qgg_solvebed,               8},
    {"_qgg_summarybed",            (DL_FUNC) &_qgg_summarybed,             6},
    {NULL, NULL, 0}
};


/* .Fortran calls */
extern void F77_NAME(eiggrm)(int *n, double *GRM, double *evals, int *ncores);
extern void F77_NAME(grmbed)(int *n, int *nr, int *rws, int *nc, int *cls1, int *cls2, int *scale, int *nbytes, int *fnRAWCHAR, int *nchars, int *msize, int *ncores, int *fnGCHAR, int *ncharsg, int *gmodel);
extern void F77_NAME(reml)(int *n, int *nf, int *nr, double *tol, int *maxit, int *ncores, int *ngr, int *indx, double *y, double *X, double *theta, double *ai, double *b, double *varb, double *u, double *Vy, double *Py, double *llik, double *trPG, double *trVG, int *ncharsg, int *fnGCHAR);

static const R_FortranMethodDef FortranEntries[] = {
    {"eiggrm",     (DL_FUNC) &F77_NAME(eiggrm),      4},
    {"grmbed",     (DL_FUNC) &F77_NAME(grmbed),     15},
    {"reml",       (DL_FUNC) &F77_NAME(reml),       22},
    {NULL, NULL, 0}
};

void R_init_qgg(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}