// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// bjimpute
Rcpp::NumericVector bjimpute(SEXP y, SEXP cen, SEXP x, SEXP inibeta);
RcppExport SEXP _ELMSurv_bjimpute(SEXP ySEXP, SEXP cenSEXP, SEXP xSEXP, SEXP inibetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type y(ySEXP);
    Rcpp::traits::input_parameter< SEXP >::type cen(cenSEXP);
    Rcpp::traits::input_parameter< SEXP >::type x(xSEXP);
    Rcpp::traits::input_parameter< SEXP >::type inibeta(inibetaSEXP);
    rcpp_result_gen = Rcpp::wrap(bjimpute(y, cen, x, inibeta));
    return rcpp_result_gen;
END_RCPP
}
// kernmat
Rcpp::NumericMatrix kernmat(SEXP xtrain, SEXP kernel_type, SEXP kernel_para, SEXP xtest);
RcppExport SEXP _ELMSurv_kernmat(SEXP xtrainSEXP, SEXP kernel_typeSEXP, SEXP kernel_paraSEXP, SEXP xtestSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type xtrain(xtrainSEXP);
    Rcpp::traits::input_parameter< SEXP >::type kernel_type(kernel_typeSEXP);
    Rcpp::traits::input_parameter< SEXP >::type kernel_para(kernel_paraSEXP);
    Rcpp::traits::input_parameter< SEXP >::type xtest(xtestSEXP);
    rcpp_result_gen = Rcpp::wrap(kernmat(xtrain, kernel_type, kernel_para, xtest));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_ELMSurv_bjimpute", (DL_FUNC) &_ELMSurv_bjimpute, 4},
    {"_ELMSurv_kernmat", (DL_FUNC) &_ELMSurv_kernmat, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_ELMSurv(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
