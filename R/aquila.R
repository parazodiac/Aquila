WriteIndus <- function(base.path, mat) {
  write.table(rownames(mat), paste0(base.path, "/genes.tsv"), 
              quote = F, row.names = F, 
              col.names = F)
  write.table(colnames(mat), paste0(base.path, "/barcodes.tsv"), 
              quote = F, row.names = F, 
              col.names = F)
  writeMM(mat, file = paste0(base.path, "/matrix.mtx"))
}


MoransI <- function(values, weights, temp_dir) {
  print("Writing values")
  values_path <- paste0(temp_dir, "/", "values")
  WriteIndus(values_path, values)
  
  print("Writing weights")
  weights_path <- paste0(temp_dir, "/", "weights")
  WriteIndus(weights_path, weights)
  
  print("Calculating MoransI")
  out_path <- paste0(temp_dir, "/", "moransi.is")
  Oxidized_MoransI(weights_path, values_path, out_path)
  
  print("Reading back stats")
  read.table(out_path)
}