// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// anova_test
int anova_test(IntegerVector treatment, NumericVector response, int K, double alpha);
RcppExport SEXP _restricted_anova_test(SEXP treatmentSEXP, SEXP responseSEXP, SEXP KSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type treatment(treatmentSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type response(responseSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(anova_test(treatment, response, K, alpha));
    return rcpp_result_gen;
END_RCPP
}
// anova_test_per_subject
List anova_test_per_subject(List treatment_list, List response_list, int K, double alpha, bool time_drift);
RcppExport SEXP _restricted_anova_test_per_subject(SEXP treatment_listSEXP, SEXP response_listSEXP, SEXP KSEXP, SEXP alphaSEXP, SEXP time_driftSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type treatment_list(treatment_listSEXP);
    Rcpp::traits::input_parameter< List >::type response_list(response_listSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< bool >::type time_drift(time_driftSEXP);
    rcpp_result_gen = Rcpp::wrap(anova_test_per_subject(treatment_list, response_list, K, alpha, time_drift));
    return rcpp_result_gen;
END_RCPP
}
// restricted_one_simulation
List restricted_one_simulation(NumericVector w, int nsbj, std::string procedure, double parameter);
RcppExport SEXP _restricted_restricted_one_simulation(SEXP wSEXP, SEXP nsbjSEXP, SEXP procedureSEXP, SEXP parameterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type w(wSEXP);
    Rcpp::traits::input_parameter< int >::type nsbj(nsbjSEXP);
    Rcpp::traits::input_parameter< std::string >::type procedure(procedureSEXP);
    Rcpp::traits::input_parameter< double >::type parameter(parameterSEXP);
    rcpp_result_gen = Rcpp::wrap(restricted_one_simulation(w, nsbj, procedure, parameter));
    return rcpp_result_gen;
END_RCPP
}
// restricted_multiple_simulations
List restricted_multiple_simulations(NumericVector w, int nsbj, std::string procedure, double parameter, int nsim);
RcppExport SEXP _restricted_restricted_multiple_simulations(SEXP wSEXP, SEXP nsbjSEXP, SEXP procedureSEXP, SEXP parameterSEXP, SEXP nsimSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type w(wSEXP);
    Rcpp::traits::input_parameter< int >::type nsbj(nsbjSEXP);
    Rcpp::traits::input_parameter< std::string >::type procedure(procedureSEXP);
    Rcpp::traits::input_parameter< double >::type parameter(parameterSEXP);
    Rcpp::traits::input_parameter< int >::type nsim(nsimSEXP);
    rcpp_result_gen = Rcpp::wrap(restricted_multiple_simulations(w, nsbj, procedure, parameter, nsim));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_restricted_anova_test", (DL_FUNC) &_restricted_anova_test, 4},
    {"_restricted_anova_test_per_subject", (DL_FUNC) &_restricted_anova_test_per_subject, 5},
    {"_restricted_restricted_one_simulation", (DL_FUNC) &_restricted_restricted_one_simulation, 4},
    {"_restricted_restricted_multiple_simulations", (DL_FUNC) &_restricted_restricted_multiple_simulations, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_restricted(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
