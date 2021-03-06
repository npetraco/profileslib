// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// computeEnvelope_lb
NumericMatrix computeEnvelope_lb(NumericVector array, int constraint);
RcppExport SEXP _profileslib_computeEnvelope_lb(SEXP arraySEXP, SEXP constraintSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type array(arraySEXP);
    Rcpp::traits::input_parameter< int >::type constraint(constraintSEXP);
    rcpp_result_gen = Rcpp::wrap(computeEnvelope_lb(array, constraint));
    return rcpp_result_gen;
END_RCPP
}
// LB_Keogh_lb
double LB_Keogh_lb(NumericVector candidateArray, NumericVector upperBound, NumericVector lowerBound);
RcppExport SEXP _profileslib_LB_Keogh_lb(SEXP candidateArraySEXP, SEXP upperBoundSEXP, SEXP lowerBoundSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type candidateArray(candidateArraySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type upperBound(upperBoundSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lowerBound(lowerBoundSEXP);
    rcpp_result_gen = Rcpp::wrap(LB_Keogh_lb(candidateArray, upperBound, lowerBound));
    return rcpp_result_gen;
END_RCPP
}
// LB_Improved
double LB_Improved(NumericVector candidateArray, NumericVector queryArray, NumericVector ubQuery, NumericVector lbQuery, int windowSize);
RcppExport SEXP _profileslib_LB_Improved(SEXP candidateArraySEXP, SEXP queryArraySEXP, SEXP ubQuerySEXP, SEXP lbQuerySEXP, SEXP windowSizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type candidateArray(candidateArraySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type queryArray(queryArraySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ubQuery(ubQuerySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lbQuery(lbQuerySEXP);
    Rcpp::traits::input_parameter< int >::type windowSize(windowSizeSEXP);
    rcpp_result_gen = Rcpp::wrap(LB_Improved(candidateArray, queryArray, ubQuery, lbQuery, windowSize));
    return rcpp_result_gen;
END_RCPP
}
// computeEnvelope
NumericMatrix computeEnvelope(NumericVector array, int constraint);
RcppExport SEXP _profileslib_computeEnvelope(SEXP arraySEXP, SEXP constraintSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type array(arraySEXP);
    Rcpp::traits::input_parameter< int >::type constraint(constraintSEXP);
    rcpp_result_gen = Rcpp::wrap(computeEnvelope(array, constraint));
    return rcpp_result_gen;
END_RCPP
}
// LB_Keogh2
double LB_Keogh2(const NumericVector& candidateArray, const NumericVector& upperBound, const NumericVector& lowerBound);
RcppExport SEXP _profileslib_LB_Keogh2(SEXP candidateArraySEXP, SEXP upperBoundSEXP, SEXP lowerBoundSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type candidateArray(candidateArraySEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type upperBound(upperBoundSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type lowerBound(lowerBoundSEXP);
    rcpp_result_gen = Rcpp::wrap(LB_Keogh2(candidateArray, upperBound, lowerBound));
    return rcpp_result_gen;
END_RCPP
}
// LB_Keogh
double LB_Keogh(NumericVector candidateArray, NumericVector upperBound, NumericVector lowerBound);
RcppExport SEXP _profileslib_LB_Keogh(SEXP candidateArraySEXP, SEXP upperBoundSEXP, SEXP lowerBoundSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type candidateArray(candidateArraySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type upperBound(upperBoundSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lowerBound(lowerBoundSEXP);
    rcpp_result_gen = Rcpp::wrap(LB_Keogh(candidateArray, upperBound, lowerBound));
    return rcpp_result_gen;
END_RCPP
}
// readSurFile
List readSurFile(std::string path);
RcppExport SEXP _profileslib_readSurFile(SEXP pathSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type path(pathSEXP);
    rcpp_result_gen = Rcpp::wrap(readSurFile(path));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_profileslib_computeEnvelope_lb", (DL_FUNC) &_profileslib_computeEnvelope_lb, 2},
    {"_profileslib_LB_Keogh_lb", (DL_FUNC) &_profileslib_LB_Keogh_lb, 3},
    {"_profileslib_LB_Improved", (DL_FUNC) &_profileslib_LB_Improved, 5},
    {"_profileslib_computeEnvelope", (DL_FUNC) &_profileslib_computeEnvelope, 2},
    {"_profileslib_LB_Keogh2", (DL_FUNC) &_profileslib_LB_Keogh2, 3},
    {"_profileslib_LB_Keogh", (DL_FUNC) &_profileslib_LB_Keogh, 3},
    {"_profileslib_readSurFile", (DL_FUNC) &_profileslib_readSurFile, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_profileslib(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
