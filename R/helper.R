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
  object <- CreateSeuratObject(counts = counts$`Gene Expression`,
                             project = "coassay",
                             assay = "RNA")

  # create chromatin assay from the peak counts
  object[["ATAC"]] <- CreateChromatinAssay(counts = counts$Peaks,
                                         min.cells = 5,
                                         sep = c(":", "-"),
                                         genome = "hg38",
                                         fragments = peaks.fragments.file,
                                         annotation = annotations)
  object
}


#' Process Signac Object for ATAC assay
#'
#' @param object Path to the signac object
#' @return Processed Signac object
#' @export
ProcessATAC <- function(object) {
  DefaultAssay(object) <- "ATAC"
  object <- NucleosomeSignal(object)
  object <- TSSEnrichment(object)

  object$blacklist <- FractionCountsInRegion(object, regions = blacklist_hg38_unified)
  object <- FindTopFeatures(object, min.cutoff = 10)
  object <- RunTFIDF(object)
  object <- RunSVD(object)
  object <- RunUMAP(object,
                    reduction = "lsi",
                    dims = 2:40,
                    reduction.name = "umap.atac")
  object
}

#' Process Signac Object for RNA assay
#'
#' @param object Path to the signac object
#' @return Processed Signac object
#' @export
ProcessRNA <- function(object) {
  # perform PCA and generate UMAP for RNA data
  DefaultAssay(object) <- "RNA"

  object <- SCTransform(object, ncells = 5000)
  object <- RunPCA(object)
  object <- RunUMAP(object, dims = 1:40, reduction.name = "umap.rna")
  object
}


#' Cluster Signac Object for ATAC assay
#'
#' @param object Path to the signac object
#' @return Cluster Signac object
#' @export
ClusterATAC <- function(object) {
  object <- FindNeighbors(object = object, reduction = 'lsi', dims = 2:30)
  object <- FindClusters(object = object, verbose = FALSE, algorithm = 3)
}


#' Transfer PBMC RNA labels to the data
#'
#' @param object Path to the signac object
#' @return RNA annotated Signac object
#' @export
TransferRnaPBMC <- function(object) {
  # preprocess RNA annotated data
  pbmc_rna <- readRDS("~/annotations/pbmc_10k_v3_annotated.rds")
  pbmc_rna <- SCTransform(pbmc_rna, ncells = 5000)

  anchors <- FindTransferAnchors(reference = pbmc_rna, query = object,
                                 normalization.method = "SCT")
  labels <- TransferData(anchorset = anchors, refdata = pbmc_rna$annotated,
                         weight.reduction = object[["pca"]])
  object <- AddMetaData(object, metadata = labels)
  object
}

AssayFromCounts <- function(object, counts){
  CreateChromatinAssay(
    counts = counts,
    genome = "hg38",
    fragments = Fragments(object),
    annotation = Annotation(object)
  )
}

NewPeaks <- function(macs.peak.file) {
  df <- read.table(
    file = macs.peak.file,
    col.names = c("chr", "start", "end", "name", "score", "strand", "one", "two", "three", "four")
  )

  gr <- makeGRangesFromDataFrame(df = df)
  gr
}

NewPeaksSeparated <- function(macs.peak.folder) {
  peakfiles <- list.files(path = macs.peak.folder,
                          pattern = "*.narrowPeak", full.names = TRUE)
  peaks.sep <- sapply(X = peakfiles, FUN = function(x){
    df <- read.table(
      file = x,
      col.names = c("chr", "start", "end", "name", "score", "dot", "v7", "v8", "v9", "v10")
    )
    gr <- makeGRangesFromDataFrame(df = df)
    return(gr)
  })

  unified <- Reduce(f = c, x = peaks.sep)
  unified <- reduce(x = unified)
  unified <- keepStandardChromosomes(unified, pruning.mode = "coarse")
  unified
}
