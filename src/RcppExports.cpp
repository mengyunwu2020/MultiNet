// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// sqrt_regression
Eigen::ArrayXd sqrt_regression(Eigen::ArrayXd ww, Eigen::MatrixXd X1, Eigen::ArrayXd Y, NumericVector groupLen, NumericVector rangeGroupInd, NumericVector lambda2, double lambda1, int nlambda, int d, int kG, int num_relaxation_round, int max_iter, int innerIter, int standard);
RcppExport SEXP _MultiNet_sqrt_regression(SEXP wwSEXP, SEXP X1SEXP, SEXP YSEXP, SEXP groupLenSEXP, SEXP rangeGroupIndSEXP, SEXP lambda2SEXP, SEXP lambda1SEXP, SEXP nlambdaSEXP, SEXP dSEXP, SEXP kGSEXP, SEXP num_relaxation_roundSEXP, SEXP max_iterSEXP, SEXP innerIterSEXP, SEXP standardSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::ArrayXd >::type ww(wwSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X1(X1SEXP);
    Rcpp::traits::input_parameter< Eigen::ArrayXd >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type groupLen(groupLenSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type rangeGroupInd(rangeGroupIndSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lambda2(lambda2SEXP);
    Rcpp::traits::input_parameter< double >::type lambda1(lambda1SEXP);
    Rcpp::traits::input_parameter< int >::type nlambda(nlambdaSEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    Rcpp::traits::input_parameter< int >::type kG(kGSEXP);
    Rcpp::traits::input_parameter< int >::type num_relaxation_round(num_relaxation_roundSEXP);
    Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< int >::type innerIter(innerIterSEXP);
    Rcpp::traits::input_parameter< int >::type standard(standardSEXP);
    rcpp_result_gen = Rcpp::wrap(sqrt_regression(ww, X1, Y, groupLen, rangeGroupInd, lambda2, lambda1, nlambda, d, kG, num_relaxation_round, max_iter, innerIter, standard));
    return rcpp_result_gen;
END_RCPP
}
// nodewise
List nodewise(Eigen::Map<Eigen::MatrixXd> data, int K, int nrow, int p, int num_relaxation_round, Eigen::Map<Eigen::MatrixXd> residual, int max_iter, int innerIter, Eigen::Map<Eigen::MatrixXd> tmp, NumericVector index, NumericVector lambda1, NumericVector lambda2, int nlam, NumericVector groupLen, NumericVector rangeGroupInd, Eigen::Map<Eigen::MatrixXd> beta, Function f, double tau, int standard, Eigen::Map<Eigen::MatrixXd> cc);
RcppExport SEXP _MultiNet_nodewise(SEXP dataSEXP, SEXP KSEXP, SEXP nrowSEXP, SEXP pSEXP, SEXP num_relaxation_roundSEXP, SEXP residualSEXP, SEXP max_iterSEXP, SEXP innerIterSEXP, SEXP tmpSEXP, SEXP indexSEXP, SEXP lambda1SEXP, SEXP lambda2SEXP, SEXP nlamSEXP, SEXP groupLenSEXP, SEXP rangeGroupIndSEXP, SEXP betaSEXP, SEXP fSEXP, SEXP tauSEXP, SEXP standardSEXP, SEXP ccSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type data(dataSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type nrow(nrowSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type num_relaxation_round(num_relaxation_roundSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type residual(residualSEXP);
    Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< int >::type innerIter(innerIterSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type tmp(tmpSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type index(indexSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lambda1(lambda1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lambda2(lambda2SEXP);
    Rcpp::traits::input_parameter< int >::type nlam(nlamSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type groupLen(groupLenSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type rangeGroupInd(rangeGroupIndSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< Function >::type f(fSEXP);
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< int >::type standard(standardSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type cc(ccSEXP);
    rcpp_result_gen = Rcpp::wrap(nodewise(data, K, nrow, p, num_relaxation_round, residual, max_iter, innerIter, tmp, index, lambda1, lambda2, nlam, groupLen, rangeGroupInd, beta, f, tau, standard, cc));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_MultiNet_sqrt_regression", (DL_FUNC) &_MultiNet_sqrt_regression, 14},
    {"_MultiNet_nodewise", (DL_FUNC) &_MultiNet_nodewise, 20},
    {NULL, NULL, 0}
};

RcppExport void R_init_MultiNet(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
