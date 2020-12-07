#define R_NO_REMAP
#define STRICT_R_HEADERS
#include <Rinternals.h>
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

extern "C" {
	bool morans_i(const char* weight_folder_path, const char* values_folder_path, const char* out_file_path);
	bool geary_c(const char* weight_folder_path, const char* values_folder_path, const char* out_file_path);
}

// [[Rcpp::export]]
bool Oxidized_MoransI(SEXP weight_folder_path, SEXP values_folder_path, SEXP out_file_path){
  return morans_i(CHAR(STRING_ELT(weight_folder_path, 0)), 
          CHAR(STRING_ELT(values_folder_path, 0)),
          CHAR(STRING_ELT(out_file_path, 0))
  );
}


// [[Rcpp::export]]
bool Oxidized_GearyC(SEXP weight_folder_path, SEXP values_folder_path, SEXP out_file_path){
  return morans_i(CHAR(STRING_ELT(weight_folder_path, 0)), 
                  CHAR(STRING_ELT(values_folder_path, 0)),
                  CHAR(STRING_ELT(out_file_path, 0))
  );
}