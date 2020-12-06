#define R_NO_REMAP
#define STRICT_R_HEADERS
#include <Rinternals.h>
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

extern "C" {
	char * string_from_rust();
}

// [[Rcpp::export]]
SEXP hello_wrapper(){
  return Rf_ScalarString(Rf_mkCharCE(string_from_rust(), CE_UTF8));
}

