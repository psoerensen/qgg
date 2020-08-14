#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Fortran calls */
extern void F77_NAME(bed2raw)(int *m, int *cls, int *nbytes, int *append, int *fnBEDCHAR, int *fnRAWCHAR, int *ncharbed, int *ncharraw);
extern void F77_NAME(eiggrm)(int *n, double *GRM, double *evals, int *ncores);
extern void F77_NAME(grmbed)(int *n, int *nr, int *rws, int *nc, int *cls1, int *cls2, int *scale, int *nbytes, int *fnRAWCHAR, int *nchars, int *msize, int *ncores, int *fnGCHAR, int *ncharsg, int *gmodel);
extern void F77_NAME(mpgrs)(int *n, int *nr, int *rws, int *nc, int *cls, int *nbytes, int *fnRAWCHAR, int *nchars, int *nprs, double *s, double *prs, double *af, int *impute, int *direction);
extern void F77_NAME(psets)(int *m, double *stat, int *nsets, double *setstat, int *msets, int *p, int *np, int *ncores);
extern void F77_NAME(readbed)(int *n, int *nr, int *rws, int *nc, int *cls, int *impute, int *scale, int *direction, double *W, int *nbytes, int *nchars, int *fnRAWCHAR);
extern void F77_NAME(reml)(int *n, int *nf, int *nr, double *tol, int *maxit, int *ncores, int *ngr, int *indx, double *y, double *X, double *theta, double *ai, double *b, double *varb, double *u, double *Vy, double *Py, double *llik, double *trPG, double *trVG, int *ncharsg, int *fnGCHAR);
extern void F77_NAME(solvebed)(int *n, int *nr, int *rws, int *nc, int *cls, int *scale, int *nbytes, int *fnRAWCHAR, int *nchars, int *ncores, int *nit, double *lambda, double *tol, double *y, double *g, double *e, double *s, double *mean, double *sd);
extern void F77_NAME(summarybed)(int *n, int *nr, int *rws, int *nc, int *cls, double *af, double *nmiss, double *n0, double *n1, double *n2, int *nbytes, int *fnRAWCHAR, int *nchars, int *ncores);
extern void F77_NAME(freqbed)(int *n, int *nr, int *rws, int *nc, int *cls, double *af, double *nmiss, double *n0, double *n1, double *n2, int *nbytes, int *fnRAWCHAR, int *nchars);




static const R_FortranMethodDef FortranEntries[] = {
    {"bed2raw",    (DL_FUNC) &F77_NAME(bed2raw),     8},
    {"eiggrm",     (DL_FUNC) &F77_NAME(eiggrm),      4},
    {"grmbed",     (DL_FUNC) &F77_NAME(grmbed),     15},
    {"mpgrs",      (DL_FUNC) &F77_NAME(mpgrs),      14},
    {"psets",      (DL_FUNC) &F77_NAME(psets),       8},
    {"readbed",    (DL_FUNC) &F77_NAME(readbed),    12},
    {"reml",       (DL_FUNC) &F77_NAME(reml),       22},
    {"solvebed",   (DL_FUNC) &F77_NAME(solvebed),   19},
    {"summarybed", (DL_FUNC) &F77_NAME(summarybed), 14},
    {"freqbed", (DL_FUNC) &F77_NAME(summarybed), 13},
    {NULL, NULL, 0}
};

void R_init_qgg(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}