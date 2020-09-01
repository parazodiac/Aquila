#' @useDynLib Circus
NULL

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
ProcessATAC <- function(object, assay.name) {
  DefaultAssay(object) <- assay.name
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

#' Create new Chromatin Assay from new peaks
#'
#' @param object Path to the signac object
#' @param peaks.file Path to the new peaks file
#' @param assay.name Name of the new assay
#' @return new chromain assay Signac object
#' @export
AddAssayFromPeaks <- function(object, peaks.file, assay.name) {
  gr <- NewPeaks(peaks.file)
  counts <- FeatureMatrix(
    fragments = Fragments(object),
    features = gr,
    cells = colnames(object)
  )
  object[[assay.name]] <-AssayFromCounts(pbmc, counts)
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

#' @export
GetCorrelationGamma <- function(peak.counts, gene.counts, pbmc, ctype) {
  pranges <- StringToGRanges(rownames(peak.counts), sep = c("-", "-"))
  pranges <- Extend(x = pranges, upstream = 100000, downstream = 100000)

  annotation <- pbmc@assays$macs.atac@annotation
  overlaps <- data.frame(findOverlaps(pranges, annotation))
  overlaps$queryHits <- rownames(peak.counts)[overlaps$queryHits]
  overlaps$subjectHits <- annotation[overlaps$subjectHits]$gene_name
  overlaps <- unique(overlaps)
  overlaps <- overlaps[overlaps$subjectHits %in% rownames(gene.counts), ]

  if (ctype == "ident") {
    return (overlaps)
  }

  num.rows = nrow(overlaps)
  corrs <- numeric(num.rows)
  for (i in 1:num.rows) {
    if (i %% 1000 == 0) {
      cat("\r", i)
    }

    pname <- overlaps[i, "queryHits"]
    gname <- overlaps[i, "subjectHits"]
    arow <- peak.counts[pname, ]
    rrow <- gene.counts[gname, ]
    if (sum(arow) > 10 && sum(rrow) > 10) {
      if (ctype == "dot") {
        corrs[i] <- sum(arow * rrow)
      } else if (ctype == "dot_depth") {
        corrs[i] <- sum((arow / sum(arow)) * (rrow / sum(rrow)))
      } else if (ctype == "ks") {
        corrs[i] <- ks.test(arow, rrow)[[1]]
      } else if (ctype == "spearman") {
        corrs[i] <- cor(arow, rrow, method = "spearman")
      } else {
        print("ERROR")
        return (list())
      }
    }
  }

  overlaps["corr"] <- corrs / sum(abs(corrs))
  overlaps <- overlaps[overlaps$corr != 0, ]
  overlaps
}

#' @export
PrincipleTime <- function(cell.embeddings) {
  rna.princurve <- principal_curve(cell.embeddings)
  rna.pseudo.time <- rna.princurve$lambda
  rna.pseudo.time <- rna.pseudo.time / max(rna.pseudo.time)
  rna.pseudo.time
}

#' @export
SubsetGamma <- function(object, keep.cells, atac.counts, rna.counts, ctype) {
  pseudotime <- PrincipleTime(object@reductions$umap.rna@cell.embeddings[keep.cells, ])
  cell.order <- names(pseudotime[order(pseudotime)])

  d.atac <- as.matrix(atac.counts[, cell.order])
  d.rna <- as.matrix(rna.counts[, cell.order])
  GetCorrelationGamma(d.atac, d.rna, object, ctype)
}

#' @export
SubsetCounts <- function(chr.name, object) {
  annotation <- pbmc@assays$macs.atac@annotation
  keep.gnames <- unique(annotation[annotation@seqnames == chr.name]$gene_name)
  keep.gnames <- intersect(keep.gnames, rownames(pbmc@assays$RNA@counts))
  keep.peaks <- grep(paste0(chr.name, "-"), rownames(pbmc@assays$macs.atac@counts))

  atac.counts <- pbmc@assays$macs.atac@counts[keep.peaks, ]
  rna.counts <- pbmc@assays$RNA@counts[keep.gnames, ]
  list(atac.counts, rna.counts)
}

#' @export
PlotATACnRNA <- function(atac, rna, pname, gname, cells) {
  df <- data.frame(cbind(atac[pname, ], rna[gname, ]))
  pname <- gsub("-", "_", pname)
  colnames(df) <- c(pname, gname)
  df["num"] <- seq(dim(df)[1])
  df <- df[cells, ]

  p1 <- ggplot(data=df, aes_string(x="num", y=pname)) + geom_point()  + geom_smooth(method="gam")
  p2 <- ggplot(data=df, aes_string(x="num", y=gname)) + geom_point() + geom_smooth(method="gam")
  p1 / p2
}

#' @export
GetConnectedRegion <- function(overlaps, gname){
  ghits <- c(gname)
  num.ghits <- 1
  gchange <- TRUE

  num.qhits <- 0
  qchange <- TRUE

  while (gchange || qchange) {
    qhits <- unique(overlaps[overlaps$subjectHits %in% ghits, ]$queryHits)
    ghits <- unique(overlaps[overlaps$queryHits %in% qhits, ]$subjectHits)
    new.num.qhits <- length(qhits)
    new.num.ghits <- length(ghits)

    if (new.num.ghits != num.ghits) {
      num.ghits <- new.num.ghits
      gchange = TRUE
    } else {
      gchange = FALSE
    }

    if (new.num.qhits != num.qhits) {
      num.qhits <- new.num.qhits
      qchange = TRUE
    } else {
      qchange = FALSE
    }
  }

  return(ghits)
}

#' @export
GetCellTypeGamma <- function(atac.counts, rna.counts, keep.cells, gname = "TCL1A") {
  alpha <- as.matrix(atac.counts[, keep.cells])
  beta <- as.matrix(rna.counts[, keep.cells])
  overlaps <- GetCorrelationGamma(alpha, beta, pbmc, "ident")
  GibbsGamma(overlaps, gname, alpha, beta)
}

#' @export
GibbsGamma <- function(overlaps, gname, alpha, beta, max.it = 1000000) {
  ghits <- GetConnectedRegion(overlaps, gname)
  phits <- unique(overlaps[overlaps$subjectHits %in% ghits, ]$queryHits)

  alpha <- alpha[phits, ]
  beta <- beta[ghits, ]

  phits <- c(phits, "pnoise")
  ghits <- c(ghits, "gnoise")
  colnames(overlaps) <- c("peaks", "genes")

  gamma <- matrix(0, nrow = length(phits), ncol = length(ghits))
  rownames(gamma) <-phits
  colnames(gamma) <- ghits
  dim(gamma)
  sum(gamma)

  num.cells <- dim(alpha)[2]
  num.genes <- length(ghits)
  num.peaks <- length(phits)

  it <- 0
  start = c(sample(1:(num.peaks-1), size = 1), sample(1:(num.genes-1), size = 1))
  while (it < max.it) { # One Gibbs iteration
    it <- it + 1
    if (it %% 100000 == 0) {
      cat("\r", it)
    }

    # first chose a cell
    cell.id <- sample(1:num.cells, size = 1)

    # next for a chosen cell and defined gene chose a peak
    gene.name <- ghits[start[2]]
    if (gene.name == "gnoise") {
      opeaks <- phits[-num.peaks]
    } else {
      opeaks <- overlaps[overlaps$genes == gene.name, ]$peaks
    }
    pdist <- alpha[opeaks, cell.id]
    pdist <- pdist[pdist != 0]
    pdist <- cumsum(pdist / sum(pdist))
    if (length(pdist) == 0) {
      start[1] <- num.peaks
    } else {
      random.num <- runif(1, min = 0, max = 1)
      peak.name <- names(pdist)[Position(function(x) x > random.num, pdist)]
      start[1] <- Position(function(x) x == peak.name, phits)
    }

    # next for a chosen cell and defined peak chose a gene
    peak.name <- phits[start[1]]
    if (peak.name == "pnoise") {
      ogenes <- ghits[-num.genes]
    } else {
      ogenes <- overlaps[overlaps$peaks == peak.name, ]$genes
    }
    gdist <- beta[ogenes, cell.id]
    names(gdist) <- ogenes

    gdist <- gdist[gdist != 0]
    gdist <- cumsum(gdist / sum(gdist))
    if (length(gdist) == 0) {
      start[2] <- num.genes
    } else {
      random.num <- runif(1, min = 0, max = 1)
      gene.name <- names(gdist)[Position(function(x) x > random.num, gdist)]
      start[2] <- Position(function(x) x == gene.name, ghits)
    }

    gamma[start[1], start[2]] <-  gamma[start[1], start[2]] + 1
  }

  return (gamma)
}
