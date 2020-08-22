#' Create Signac Object
#'
#' This function loads the mRNA and peaks file as a Signac object.
#'
#' @param mrna.h5.file Path to the input mrna h5 file
#' @param peaks.fragments.file Path to the input peaks file
#' @return A Signac object
#' @export
GetSignacObject <- function(mrna.h5.file, peaks.fragments.file) {
  annotations <- GetGRangesFromEnsDb(EnsDb.Hsapiens.v86)
  seqlevelsStyle(annotations) <- "UCSC"
  genome(annotations) <- "hg38"
  annotations

  # create Seurat object from gene counts
  counts <- Read10X_h5(mrna.h5.file)
  pbmc <- CreateSeuratObject(counts = counts$`Gene Expression`,
                             project = "coassay",
                             assay = "RNA")

  # create chromatin assay from the peak counts
  pbmc[["ATAC"]] <- CreateChromatinAssay(counts = counts$Peaks,
                                         min.cells = 5,
                                         sep = c(":", "-"),
                                         genome = "hg38",
                                         fragments = peaks.fragments.file,
                                         annotation = annotations)
  pbmc
}
