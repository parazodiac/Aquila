#define R_NO_REMAP
#define STRICT_R_HEADERS
#include <Rinternals.h>
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

extern "C" {
	bool morans_i(const char* weight_folder_path, const char* values_folder_path, const char* out_file_path, const char* num_threads);
	bool gearys_c(const char* weight_folder_path, const char* values_folder_path, const char* out_file_path, const char* num_threads);
	bool nnhelper_rust(const char* mat_folder_path, const char* out_file_path, const char* num_threads);
}

// [[Rcpp::export]]
bool Oxidized_MoransI(SEXP weight_folder_path, SEXP values_folder_path, SEXP out_file_path, SEXP num_threads){
  return morans_i(CHAR(STRING_ELT(weight_folder_path, 0)), 
                  CHAR(STRING_ELT(values_folder_path, 0)),
                  CHAR(STRING_ELT(out_file_path, 0)),
                  CHAR(STRING_ELT(num_threads, 0))
  );
}

// [[Rcpp::export]]
bool Oxidized_GearysC(SEXP weight_folder_path, SEXP values_folder_path, SEXP out_file_path, SEXP num_threads){
  return gearys_c(CHAR(STRING_ELT(weight_folder_path, 0)), 
                  CHAR(STRING_ELT(values_folder_path, 0)),
                  CHAR(STRING_ELT(out_file_path, 0)),
                  CHAR(STRING_ELT(num_threads, 0))
  );
}

// [[Rcpp::export]]
bool Oxidized_NNHelper(SEXP mat_folder_path, SEXP out_file_path, SEXP num_threads){
  return nnhelper_rust(CHAR(STRING_ELT(mat_folder_path, 0)), 
                       CHAR(STRING_ELT(out_file_path, 0)),
                       CHAR(STRING_ELT(num_threads, 0))
  );
}
