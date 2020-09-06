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
GetTauOverlaps <- function(peak.names, gene.names, gene.ranges, tau=100000) {
  pranges <- StringToGRanges(peak.names, sep = c("-", "-"))
  pranges <- Extend(x = pranges, upstream = tau, downstream = tau)

  overlaps <- data.frame(findOverlaps(pranges, gene.ranges))
  overlaps$queryHits <- peak.names[overlaps$queryHits]
  overlaps$subjectHits <- gene.ranges[overlaps$subjectHits]$gene_name

  overlaps <- unique(overlaps)
  overlaps <- overlaps[overlaps$subjectHits %in% gene.names, ]
  colnames(overlaps) <- c("peaks", "genes")
  return (overlaps)
}

#' @export
GetConnectedRegion <- function(overlaps, gname){
  ghits <- c(gname)
  num.ghits <- 1
  gchange <- TRUE

  num.phits <- 0
  qchange <- TRUE

  # recursively keep adding gene and peaks until
  # neither of them add new element.
  while (gchange || qchange) {
    phits <- unique(overlaps[overlaps$genes %in% ghits, ]$peaks)
    ghits <- unique(overlaps[overlaps$peaks %in% phits, ]$genes)
    new.num.phits <- length(phits)
    new.num.ghits <- length(ghits)

    if (new.num.ghits != num.ghits) {
      num.ghits <- new.num.ghits
      gchange = TRUE
    } else {
      gchange = FALSE
    }

    if (new.num.phits != num.phits) {
      num.phits <- new.num.phits
      qchange = TRUE
    } else {
      qchange = FALSE
    }
  }

  return(ghits)
}

#' @export
GetConnectedGenes <- function(overlaps, gene.name = NULL) {
  regions <- list()
  genes <- if(is.null(gene.name)) unique(overlaps$genes) else gene.name
  while(length(genes) > 0) {
    start.gene <- genes[1]
    gene.group <-GetConnectedRegion(overlaps, start.gene)
    regions <- append(regions, list(gene.group))

    genes <- setdiff(genes, gene.group)
  }

  regions
}

#' @export
DiffPlot <- function(gammas.1, gammas.2, gname) {
  idx <- which(unlist(lapply(gammas.1, function(gamma) { gname %in% colnames(gamma) })))
  diff <- gammas.1[[idx]] - gammas.2[[idx]]
  hm.df <- setNames(melt(diff), c('peaks', 'genes', 'diffs'))
  ggplot(hm.df, aes(genes, peaks)) +
    geom_tile(aes(fill=diffs), colour="white") +
    scale_fill_gradient(low="red", high="green")
}

#' @export
GenerateGamma <- function(signac.object, chr.name, keep.cells,
                          atac.assay = "macs.atac", rna.assay = "RNA",
                          tau = 100000, num.samples = 1000000,
                          gene.name = NULL, overlaps = NULL,
                          verbose = TRUE) {
  annotations <- signac.object[[atac.assay]]@annotation

  keep.peaks <- grep(paste0(chr.name, "-"), rownames(signac.object[[atac.assay]]))
  keep.genes <- unique(annotations[annotations@seqnames == chr.name]$gene_name)
  keep.genes <- intersect(keep.genes, rownames(signac.object[[rna.assay]]))

  atac.counts <- signac.object[[atac.assay]]@counts[keep.peaks, ]
  rna.counts <- signac.object[[rna.assay]]@counts[keep.genes, ]

  if (verbose) {
    print(paste0("Found ", length(keep.peaks), " peaks and ",
                 length(keep.genes), " genes."))
    print("Finding gene overlaps with peaks")
  }
  # extract overlapping peaks and genes within tau bound
  if (is.null(overlaps)) {
    overlaps <- GetTauOverlaps(peak.names = rownames(atac.counts),
                               gene.names = rownames(rna.counts),
                               gene.ranges = annotations,
                               tau)
  }

  if (verbose) {
    print(paste0("Found ", dim(overlaps)[[1]], " overlaps."))
    print("Finding connected gene groups")
  }
  # group genes which are overlapping into regions
  gene.groups <- GetConnectedGenes(overlaps, gene.name)

  if (verbose) {
    print(paste0("Found ", length(gene.groups), " gene groups."))
    print("Generating gammas")
  }
  # loop over gene groups and generate a gamma matrix
  # for each group
  gammas <- lapply(gene.groups, function(gene.hits){
    peak.hits <- unique(overlaps[overlaps$genes %in% gene.hits, ]$peaks)
    if (length(gene.hits) == 1 || length(peak.hits) == 1) {
      return (list())
    }

    alpha <- as.matrix(atac.counts[peak.hits, keep.cells])
    beta <- as.matrix(rna.counts[gene.hits, keep.cells])

    genes <- seq(length(gene.hits))
    names(genes) <- gene.hits

    peaks <- seq(length(peak.hits))
    names(peaks) <- peak.hits

    group.overlaps <- overlaps[overlaps$genes %in% gene.hits, ]
    group.overlaps$genes <- genes[group.overlaps$genes]
    group.overlaps$peaks <- peaks[group.overlaps$peaks]

    gamma <- getGamma(alpha, beta, group.overlaps, num.samples, FALSE)
    rownames(gamma) <- names(peaks)
    colnames(gamma) <- names(genes)

    gamma[, ] / num.samples
  })

  gammas
}

###########################
# pseudotime
###########################

#' @export
PrincipleTime <- function(cell.embeddings) {
  rna.princurve <- principal_curve(cell.embeddings)
  rna.pseudo.time <- rna.princurve$lambda
  rna.pseudo.time <- rna.pseudo.time / max(rna.pseudo.time)
  rna.pseudo.time
}

#' @export
GetSlingshotPseudotime <- function(counts) {
  # https://bioconductor.org/packages/release/bioc/vignettes/slingshot/inst/doc/vignette.html
  sim <- SingleCellExperiment(assays = List(counts = counts))
  geneFilter <- apply(assays(sim)$counts,1,function(x){
    sum(x >= 3) >= 10
  })
  sim <- sim[geneFilter, ]

  FQnorm <- function(counts){
    rk <- apply(counts,2,rank,ties.method='min')
    counts.sort <- apply(counts,2,sort)
    refdist <- apply(counts.sort,1,median)
    norm <- apply(rk,2,function(r){ refdist[r] })
    rownames(norm) <- rownames(counts)
    return(norm)
  }
  assays(sim)$norm <- FQnorm(assays(sim)$counts)

  pca <- prcomp(t(log1p(assays(sim)$norm)), scale. = FALSE)
  dm <- DiffusionMap(t(log1p(assays(sim)$norm)))

  rd1 <- pca$x[,1:2]
  rd2 <- cbind(DC1 = dm$DC1, DC2 = dm$DC2)
  reducedDims(sim) <- SimpleList(PCA = rd1, DiffMap = rd2)

  cl1 <- Mclust(rd1)$classification
  cl2 <- kmeans(rd1, centers = 5)$cluster
  colData(sim)$GMM <- cl1
  colData(sim)$kmeans <- cl2
  sim <- slingshot(sim, clusterLabels = 'GMM', reducedDim = 'PCA')
  #colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
  #plotcol <- colors[cut(sim$slingPseudotime_1, breaks=100)]
  #plot(reducedDims(sim)$PCA, col = plotcol, pch=16, asp = 1)
  #lines(SlingshotDataSet(sim), lwd=2, col='black')

  pseudotime <- sim$slingPseudotime_1
  pseudotime <- pseudotime / max(pseudotime)
  names(pseudotime) <- rownames(sim@colData)
  pseudotime
}
