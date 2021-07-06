// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// readG
IntegerMatrix readG(const char* file, int n, std::vector<int> cls);
RcppExport SEXP _qgg_readG(SEXP fileSEXP, SEXP nSEXP, SEXP clsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const char* >::type file(fileSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type cls(clsSEXP);
    rcpp_result_gen = Rcpp::wrap(readG(file, n, cls));
    return rcpp_result_gen;
END_RCPP
}
// readW
NumericMatrix readW(const char* file, int n, std::vector<int> cls, std::vector<float> af);
RcppExport SEXP _qgg_readW(SEXP fileSEXP, SEXP nSEXP, SEXP clsSEXP, SEXP afSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const char* >::type file(fileSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type cls(clsSEXP);
    Rcpp::traits::input_parameter< std::vector<float> >::type af(afSEXP);
    rcpp_result_gen = Rcpp::wrap(readW(file, n, cls, af));
    return rcpp_result_gen;
END_RCPP
}
// getWlist
std::vector<std::vector<float>> getWlist(const char* file, int n, std::vector<int> cls, std::vector<float> af);
RcppExport SEXP _qgg_getWlist(SEXP fileSEXP, SEXP nSEXP, SEXP clsSEXP, SEXP afSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const char* >::type file(fileSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type cls(clsSEXP);
    Rcpp::traits::input_parameter< std::vector<float> >::type af(afSEXP);
    rcpp_result_gen = Rcpp::wrap(getWlist(file, n, cls, af));
    return rcpp_result_gen;
END_RCPP
}
// freqbed
IntegerMatrix freqbed(const char* file, int n, std::vector<int> mask, std::vector<int> cls);
RcppExport SEXP _qgg_freqbed(SEXP fileSEXP, SEXP nSEXP, SEXP maskSEXP, SEXP clsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const char* >::type file(fileSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type mask(maskSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type cls(clsSEXP);
    rcpp_result_gen = Rcpp::wrap(freqbed(file, n, mask, cls));
    return rcpp_result_gen;
END_RCPP
}
// summarybed
std::vector<std::vector<std::vector<float>>> summarybed(const char* file, int n, std::vector<int> cls, std::vector<float> af, std::vector<std::vector<float>> weights, std::vector<std::vector<float>> y);
RcppExport SEXP _qgg_summarybed(SEXP fileSEXP, SEXP nSEXP, SEXP clsSEXP, SEXP afSEXP, SEXP weightsSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const char* >::type file(fileSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type cls(clsSEXP);
    Rcpp::traits::input_parameter< std::vector<float> >::type af(afSEXP);
    Rcpp::traits::input_parameter< std::vector<std::vector<float>> >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< std::vector<std::vector<float>> >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(summarybed(file, n, cls, af, weights, y));
    return rcpp_result_gen;
END_RCPP
}
// bayes
std::vector<std::vector<double>> bayes(std::vector<double> y, std::vector<std::vector<double>> W, std::vector<double> b, std::vector<double> lambda, double pi, double vara, double varb, double vare, double ssb_prior, double sse_prior, double nub, double nue, bool updateB, bool updateE, bool updatePi, int nit, int method);
RcppExport SEXP _qgg_bayes(SEXP ySEXP, SEXP WSEXP, SEXP bSEXP, SEXP lambdaSEXP, SEXP piSEXP, SEXP varaSEXP, SEXP varbSEXP, SEXP vareSEXP, SEXP ssb_priorSEXP, SEXP sse_priorSEXP, SEXP nubSEXP, SEXP nueSEXP, SEXP updateBSEXP, SEXP updateESEXP, SEXP updatePiSEXP, SEXP nitSEXP, SEXP methodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<double> >::type y(ySEXP);
    Rcpp::traits::input_parameter< std::vector<std::vector<double>> >::type W(WSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type b(bSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type pi(piSEXP);
    Rcpp::traits::input_parameter< double >::type vara(varaSEXP);
    Rcpp::traits::input_parameter< double >::type varb(varbSEXP);
    Rcpp::traits::input_parameter< double >::type vare(vareSEXP);
    Rcpp::traits::input_parameter< double >::type ssb_prior(ssb_priorSEXP);
    Rcpp::traits::input_parameter< double >::type sse_prior(sse_priorSEXP);
    Rcpp::traits::input_parameter< double >::type nub(nubSEXP);
    Rcpp::traits::input_parameter< double >::type nue(nueSEXP);
    Rcpp::traits::input_parameter< bool >::type updateB(updateBSEXP);
    Rcpp::traits::input_parameter< bool >::type updateE(updateESEXP);
    Rcpp::traits::input_parameter< bool >::type updatePi(updatePiSEXP);
    Rcpp::traits::input_parameter< int >::type nit(nitSEXP);
    Rcpp::traits::input_parameter< int >::type method(methodSEXP);
    rcpp_result_gen = Rcpp::wrap(bayes(y, W, b, lambda, pi, vara, varb, vare, ssb_prior, sse_prior, nub, nue, updateB, updateE, updatePi, nit, method));
    return rcpp_result_gen;
END_RCPP
}
// fbayes
std::vector<std::vector<double>> fbayes(std::vector<double> y, std::vector<std::vector<double>> W, std::vector<std::vector<double>> LD, std::vector<double> b, std::vector<double> lambda, double pi, double vara, double varb, double vare, double nub, double nue, bool updateB, bool updateE, bool updatePi, int nit, int method);
RcppExport SEXP _qgg_fbayes(SEXP ySEXP, SEXP WSEXP, SEXP LDSEXP, SEXP bSEXP, SEXP lambdaSEXP, SEXP piSEXP, SEXP varaSEXP, SEXP varbSEXP, SEXP vareSEXP, SEXP nubSEXP, SEXP nueSEXP, SEXP updateBSEXP, SEXP updateESEXP, SEXP updatePiSEXP, SEXP nitSEXP, SEXP methodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<double> >::type y(ySEXP);
    Rcpp::traits::input_parameter< std::vector<std::vector<double>> >::type W(WSEXP);
    Rcpp::traits::input_parameter< std::vector<std::vector<double>> >::type LD(LDSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type b(bSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type pi(piSEXP);
    Rcpp::traits::input_parameter< double >::type vara(varaSEXP);
    Rcpp::traits::input_parameter< double >::type varb(varbSEXP);
    Rcpp::traits::input_parameter< double >::type vare(vareSEXP);
    Rcpp::traits::input_parameter< double >::type nub(nubSEXP);
    Rcpp::traits::input_parameter< double >::type nue(nueSEXP);
    Rcpp::traits::input_parameter< bool >::type updateB(updateBSEXP);
    Rcpp::traits::input_parameter< bool >::type updateE(updateESEXP);
    Rcpp::traits::input_parameter< bool >::type updatePi(updatePiSEXP);
    Rcpp::traits::input_parameter< int >::type nit(nitSEXP);
    Rcpp::traits::input_parameter< int >::type method(methodSEXP);
    rcpp_result_gen = Rcpp::wrap(fbayes(y, W, LD, b, lambda, pi, vara, varb, vare, nub, nue, updateB, updateE, updatePi, nit, method));
    return rcpp_result_gen;
END_RCPP
}
// sbayes
std::vector<std::vector<double>> sbayes(std::vector<double> wy, std::vector<std::vector<double>> LD, std::vector<double> b, std::vector<double> lambda, double yy, double pi, double vara, double varb, double vare, double ssb_prior, double sse_prior, double nub, double nue, bool updateB, bool updateE, bool updatePi, int n, int nit, int method);
RcppExport SEXP _qgg_sbayes(SEXP wySEXP, SEXP LDSEXP, SEXP bSEXP, SEXP lambdaSEXP, SEXP yySEXP, SEXP piSEXP, SEXP varaSEXP, SEXP varbSEXP, SEXP vareSEXP, SEXP ssb_priorSEXP, SEXP sse_priorSEXP, SEXP nubSEXP, SEXP nueSEXP, SEXP updateBSEXP, SEXP updateESEXP, SEXP updatePiSEXP, SEXP nSEXP, SEXP nitSEXP, SEXP methodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<double> >::type wy(wySEXP);
    Rcpp::traits::input_parameter< std::vector<std::vector<double>> >::type LD(LDSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type b(bSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type yy(yySEXP);
    Rcpp::traits::input_parameter< double >::type pi(piSEXP);
    Rcpp::traits::input_parameter< double >::type vara(varaSEXP);
    Rcpp::traits::input_parameter< double >::type varb(varbSEXP);
    Rcpp::traits::input_parameter< double >::type vare(vareSEXP);
    Rcpp::traits::input_parameter< double >::type ssb_prior(ssb_priorSEXP);
    Rcpp::traits::input_parameter< double >::type sse_prior(sse_priorSEXP);
    Rcpp::traits::input_parameter< double >::type nub(nubSEXP);
    Rcpp::traits::input_parameter< double >::type nue(nueSEXP);
    Rcpp::traits::input_parameter< bool >::type updateB(updateBSEXP);
    Rcpp::traits::input_parameter< bool >::type updateE(updateESEXP);
    Rcpp::traits::input_parameter< bool >::type updatePi(updatePiSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type nit(nitSEXP);
    Rcpp::traits::input_parameter< int >::type method(methodSEXP);
    rcpp_result_gen = Rcpp::wrap(sbayes(wy, LD, b, lambda, yy, pi, vara, varb, vare, ssb_prior, sse_prior, nub, nue, updateB, updateE, updatePi, n, nit, method));
    return rcpp_result_gen;
END_RCPP
}
// grsbed
std::vector<float> grsbed(const char* file, int n, std::vector<int> cls, std::vector<float> af, std::vector<float> b);
RcppExport SEXP _qgg_grsbed(SEXP fileSEXP, SEXP nSEXP, SEXP clsSEXP, SEXP afSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const char* >::type file(fileSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type cls(clsSEXP);
    Rcpp::traits::input_parameter< std::vector<float> >::type af(afSEXP);
    Rcpp::traits::input_parameter< std::vector<float> >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(grsbed(file, n, cls, af, b));
    return rcpp_result_gen;
END_RCPP
}
// mtgrsbed
std::vector<std::vector<float>> mtgrsbed(const char* file, int n, std::vector<int> cls, std::vector<float> af, bool scale, std::vector<std::vector<float>> b);
RcppExport SEXP _qgg_mtgrsbed(SEXP fileSEXP, SEXP nSEXP, SEXP clsSEXP, SEXP afSEXP, SEXP scaleSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const char* >::type file(fileSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type cls(clsSEXP);
    Rcpp::traits::input_parameter< std::vector<float> >::type af(afSEXP);
    Rcpp::traits::input_parameter< bool >::type scale(scaleSEXP);
    Rcpp::traits::input_parameter< std::vector<std::vector<float>> >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(mtgrsbed(file, n, cls, af, scale, b));
    return rcpp_result_gen;
END_RCPP
}
// mmult
arma::mat mmult(arma::mat A, arma::mat B);
RcppExport SEXP _qgg_mmult(SEXP ASEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::mat >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(mmult(A, B));
    return rcpp_result_gen;
END_RCPP
}
// mvrnorm
arma::mat mvrnorm(arma::mat sigma);
RcppExport SEXP _qgg_mvrnorm(SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type sigma(sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(mvrnorm(sigma));
    return rcpp_result_gen;
END_RCPP
}
// rwishart
arma::mat rwishart(unsigned int df, const arma::mat& S);
RcppExport SEXP _qgg_rwishart(SEXP dfSEXP, SEXP SSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type df(dfSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type S(SSEXP);
    rcpp_result_gen = Rcpp::wrap(rwishart(df, S));
    return rcpp_result_gen;
END_RCPP
}
// riwishart
arma::mat riwishart(unsigned int df, const arma::mat& S);
RcppExport SEXP _qgg_riwishart(SEXP dfSEXP, SEXP SSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type df(dfSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type S(SSEXP);
    rcpp_result_gen = Rcpp::wrap(riwishart(df, S));
    return rcpp_result_gen;
END_RCPP
}
// mtbayes
std::vector<std::vector<std::vector<double>>> mtbayes(std::vector<std::vector<double>> y, std::vector<std::vector<double>> W, std::vector<std::vector<double>> b, arma::mat B, arma::mat E, std::vector<std::vector<double>> ssb_prior, std::vector<std::vector<double>> sse_prior, std::vector<std::vector<int>> models, std::vector<double> pi, double nub, double nue, bool updateB, bool updateE, bool updatePi, int nit, int method);
RcppExport SEXP _qgg_mtbayes(SEXP ySEXP, SEXP WSEXP, SEXP bSEXP, SEXP BSEXP, SEXP ESEXP, SEXP ssb_priorSEXP, SEXP sse_priorSEXP, SEXP modelsSEXP, SEXP piSEXP, SEXP nubSEXP, SEXP nueSEXP, SEXP updateBSEXP, SEXP updateESEXP, SEXP updatePiSEXP, SEXP nitSEXP, SEXP methodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<std::vector<double>> >::type y(ySEXP);
    Rcpp::traits::input_parameter< std::vector<std::vector<double>> >::type W(WSEXP);
    Rcpp::traits::input_parameter< std::vector<std::vector<double>> >::type b(bSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type B(BSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type E(ESEXP);
    Rcpp::traits::input_parameter< std::vector<std::vector<double>> >::type ssb_prior(ssb_priorSEXP);
    Rcpp::traits::input_parameter< std::vector<std::vector<double>> >::type sse_prior(sse_priorSEXP);
    Rcpp::traits::input_parameter< std::vector<std::vector<int>> >::type models(modelsSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type pi(piSEXP);
    Rcpp::traits::input_parameter< double >::type nub(nubSEXP);
    Rcpp::traits::input_parameter< double >::type nue(nueSEXP);
    Rcpp::traits::input_parameter< bool >::type updateB(updateBSEXP);
    Rcpp::traits::input_parameter< bool >::type updateE(updateESEXP);
    Rcpp::traits::input_parameter< bool >::type updatePi(updatePiSEXP);
    Rcpp::traits::input_parameter< int >::type nit(nitSEXP);
    Rcpp::traits::input_parameter< int >::type method(methodSEXP);
    rcpp_result_gen = Rcpp::wrap(mtbayes(y, W, b, B, E, ssb_prior, sse_prior, models, pi, nub, nue, updateB, updateE, updatePi, nit, method));
    return rcpp_result_gen;
END_RCPP
}
// solvebed
std::vector<std::vector<double>> solvebed(const char* file, int n, std::vector<int> cls, int nit, std::vector<float> af, std::vector<float> b, std::vector<float> lambda, std::vector<float> y);
RcppExport SEXP _qgg_solvebed(SEXP fileSEXP, SEXP nSEXP, SEXP clsSEXP, SEXP nitSEXP, SEXP afSEXP, SEXP bSEXP, SEXP lambdaSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const char* >::type file(fileSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type cls(clsSEXP);
    Rcpp::traits::input_parameter< int >::type nit(nitSEXP);
    Rcpp::traits::input_parameter< std::vector<float> >::type af(afSEXP);
    Rcpp::traits::input_parameter< std::vector<float> >::type b(bSEXP);
    Rcpp::traits::input_parameter< std::vector<float> >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< std::vector<float> >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(solvebed(file, n, cls, nit, af, b, lambda, y));
    return rcpp_result_gen;
END_RCPP
}
// mtsolvebed
std::vector<std::vector<std::vector<double>>> mtsolvebed(const char* file, int n, std::vector<int> cls, int nit, std::vector<float> af, std::vector<std::vector<float>> b, std::vector<std::vector<float>> lambda, std::vector<std::vector<float>> y);
RcppExport SEXP _qgg_mtsolvebed(SEXP fileSEXP, SEXP nSEXP, SEXP clsSEXP, SEXP nitSEXP, SEXP afSEXP, SEXP bSEXP, SEXP lambdaSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const char* >::type file(fileSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type cls(clsSEXP);
    Rcpp::traits::input_parameter< int >::type nit(nitSEXP);
    Rcpp::traits::input_parameter< std::vector<float> >::type af(afSEXP);
    Rcpp::traits::input_parameter< std::vector<std::vector<float>> >::type b(bSEXP);
    Rcpp::traits::input_parameter< std::vector<std::vector<float>> >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< std::vector<std::vector<float>> >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(mtsolvebed(file, n, cls, nit, af, b, lambda, y));
    return rcpp_result_gen;
END_RCPP
}
// pruneld
std::vector<int> pruneld(const char* file, int ldsize, std::vector<int> cls, std::vector<float> p, float threshold, float r2);
RcppExport SEXP _qgg_pruneld(SEXP fileSEXP, SEXP ldsizeSEXP, SEXP clsSEXP, SEXP pSEXP, SEXP thresholdSEXP, SEXP r2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const char* >::type file(fileSEXP);
    Rcpp::traits::input_parameter< int >::type ldsize(ldsizeSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type cls(clsSEXP);
    Rcpp::traits::input_parameter< std::vector<float> >::type p(pSEXP);
    Rcpp::traits::input_parameter< float >::type threshold(thresholdSEXP);
    Rcpp::traits::input_parameter< float >::type r2(r2SEXP);
    rcpp_result_gen = Rcpp::wrap(pruneld(file, ldsize, cls, p, threshold, r2));
    return rcpp_result_gen;
END_RCPP
}
// pruneldmat
std::vector<std::vector<std::vector<int>>> pruneldmat(const char* file, int ldsize, std::vector<std::vector<float>> p, std::vector<float> threshold, float r2);
RcppExport SEXP _qgg_pruneldmat(SEXP fileSEXP, SEXP ldsizeSEXP, SEXP pSEXP, SEXP thresholdSEXP, SEXP r2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const char* >::type file(fileSEXP);
    Rcpp::traits::input_parameter< int >::type ldsize(ldsizeSEXP);
    Rcpp::traits::input_parameter< std::vector<std::vector<float>> >::type p(pSEXP);
    Rcpp::traits::input_parameter< std::vector<float> >::type threshold(thresholdSEXP);
    Rcpp::traits::input_parameter< float >::type r2(r2SEXP);
    rcpp_result_gen = Rcpp::wrap(pruneldmat(file, ldsize, p, threshold, r2));
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_hello_world
arma::mat rcpparma_hello_world();
RcppExport SEXP _qgg_rcpparma_hello_world() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rcpparma_hello_world());
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_outerproduct
arma::mat rcpparma_outerproduct(const arma::colvec& x);
RcppExport SEXP _qgg_rcpparma_outerproduct(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_outerproduct(x));
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_innerproduct
double rcpparma_innerproduct(const arma::colvec& x);
RcppExport SEXP _qgg_rcpparma_innerproduct(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_innerproduct(x));
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_bothproducts
Rcpp::List rcpparma_bothproducts(const arma::colvec& x);
RcppExport SEXP _qgg_rcpparma_bothproducts(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_bothproducts(x));
    return rcpp_result_gen;
END_RCPP
}
// cp
arma::mat cp(arma::mat& W);
RcppExport SEXP _qgg_cp(SEXP WSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type W(WSEXP);
    rcpp_result_gen = Rcpp::wrap(cp(W));
    return rcpp_result_gen;
END_RCPP
}
// psets
std::vector<int> psets(std::vector<int> msets, std::vector<float> setstat, std::vector<float> stat, int np);
RcppExport SEXP _qgg_psets(SEXP msetsSEXP, SEXP setstatSEXP, SEXP statSEXP, SEXP npSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<int> >::type msets(msetsSEXP);
    Rcpp::traits::input_parameter< std::vector<float> >::type setstat(setstatSEXP);
    Rcpp::traits::input_parameter< std::vector<float> >::type stat(statSEXP);
    Rcpp::traits::input_parameter< int >::type np(npSEXP);
    rcpp_result_gen = Rcpp::wrap(psets(msets, setstat, stat, np));
    return rcpp_result_gen;
END_RCPP
}
