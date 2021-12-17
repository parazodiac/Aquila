// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// Oxidized_MoransI
bool Oxidized_MoransI(SEXP weight_folder_path, SEXP values_folder_path, SEXP out_file_path, SEXP num_threads);
RcppExport SEXP _aquila_Oxidized_MoransI(SEXP weight_folder_pathSEXP, SEXP values_folder_pathSEXP, SEXP out_file_pathSEXP, SEXP num_threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type weight_folder_path(weight_folder_pathSEXP);
    Rcpp::traits::input_parameter< SEXP >::type values_folder_path(values_folder_pathSEXP);
    Rcpp::traits::input_parameter< SEXP >::type out_file_path(out_file_pathSEXP);
    Rcpp::traits::input_parameter< SEXP >::type num_threads(num_threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(Oxidized_MoransI(weight_folder_path, values_folder_path, out_file_path, num_threads));
    return rcpp_result_gen;
END_RCPP
}
// Oxidized_GearysC
bool Oxidized_GearysC(SEXP weight_folder_path, SEXP values_folder_path, SEXP out_file_path, SEXP num_threads);
RcppExport SEXP _aquila_Oxidized_GearysC(SEXP weight_folder_pathSEXP, SEXP values_folder_pathSEXP, SEXP out_file_pathSEXP, SEXP num_threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type weight_folder_path(weight_folder_pathSEXP);
    Rcpp::traits::input_parameter< SEXP >::type values_folder_path(values_folder_pathSEXP);
    Rcpp::traits::input_parameter< SEXP >::type out_file_path(out_file_pathSEXP);
    Rcpp::traits::input_parameter< SEXP >::type num_threads(num_threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(Oxidized_GearysC(weight_folder_path, values_folder_path, out_file_path, num_threads));
    return rcpp_result_gen;
END_RCPP
}
// Oxidized_NNHelper
bool Oxidized_NNHelper(SEXP mat_folder_path, SEXP out_file_path, SEXP num_threads);
RcppExport SEXP _aquila_Oxidized_NNHelper(SEXP mat_folder_pathSEXP, SEXP out_file_pathSEXP, SEXP num_threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type mat_folder_path(mat_folder_pathSEXP);
    Rcpp::traits::input_parameter< SEXP >::type out_file_path(out_file_pathSEXP);
    Rcpp::traits::input_parameter< SEXP >::type num_threads(num_threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(Oxidized_NNHelper(mat_folder_path, out_file_path, num_threads));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_aquila_Oxidized_MoransI", (DL_FUNC) &_aquila_Oxidized_MoransI, 4},
    {"_aquila_Oxidized_GearysC", (DL_FUNC) &_aquila_Oxidized_GearysC, 4},
    {"_aquila_Oxidized_NNHelper", (DL_FUNC) &_aquila_Oxidized_NNHelper, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_aquila(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
