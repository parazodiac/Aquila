#' @useDynLib aquila
#' @export
MoransI <- function(values, weights, temp_dir) {
  values_path <- file.path(temp_dir, "values")
  dir.create(values_path, showWarnings = T)
  
  weights_path <- file.path(temp_dir, "weights")
  dir.create(weights_path, showWarnings = T)
  
  print(paste0("Writing values at: ", values_path))
  values <- Matrix(values, sparse = TRUE)
  WriteCRMatrix(values_path, values)
  print("Done writing values")
  
  print(paste0("Writing weights at: ", weights_path))
  weights <- Matrix(weights, sparse = TRUE)
  WriteCRMatrix(weights_path, weights)
  print("Done writing weights")
  
  print("Calculating MoransI")
  out_path <- file.path(temp_dir, "moransi.is")
  Oxidized_MoransI(weights_path, values_path, out_path)
  
  unlink(values_path, recursive=TRUE)
  unlink(weights_path, recursive=TRUE)
  
  print("Reading back stats")
  read.table(out_path)
}

WriteCRMatrix <- function(base.path, mat) {
  write.table(rownames(mat), paste0(base.path, "/genes.tsv"), 
              quote = F, row.names = F, 
              col.names = F)
  write.table(colnames(mat), paste0(base.path, "/barcodes.tsv"), 
              quote = F, row.names = F, 
              col.names = F)
  writeMM(mat, file = paste0(base.path, "/matrix.mtx"))
}
