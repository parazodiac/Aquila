library(Seurat)
library(Signac)

# The file is on O5
pbmc <- readRDS(file = "/home/srivastavaa/avi/tim/peak_test/pbmc.Rds")

# extract ATAC/RNA counts for all peaks/genes in chromosome 14
gname <- "CSGALNACT1"
chr.name <- "chr8"
counts <- SubsetCounts(chr.name, pbmc)
atac.counts <- counts[[1]]
rna.counts <- counts[[2]]
dim(atac.counts)
# 4530 11412
dim(rna.counts)
# 621 11412

# extract cells with naive and memory B cell annotation
naive.b.cells <- names(pbmc$predicted.id[pbmc$predicted.id  == "Naive B"])
mem.b.cells <- names(pbmc$predicted.id[pbmc$predicted.id  == "Memory B"])
paste(length(naive.b.cells), length(mem.b.cells))
# "434 516"

# extract overlapping peaks and genes within tau bound
# Note: This does not uses the counts just the annotations.
overlaps <- GetTauOverlaps(peak.names = rownames(atac.counts),
                           gene.names = rownames(rna.counts),
                           gene.ranges = pbmc@assays$macs.atac@annotation,
                           tau = 100000)

# chose a gene name to extract all the connected gene
# Note: The number of peak is 68 not 140
# I missed the unique in the line 35 resulting in
# redundant peak names but the slides were on the
# corrected counts, I forgot it's 68 not 140.
gene.hits <- GetConnectedRegion(overlaps, gname)
peak.hits <- unique(overlaps[overlaps$genes %in% gene.hits, ]$peaks)
paste(length(peak.hits), length(gene.hits))
#"68 7"

#keep.cells <- naive.b.cells
# generate gammas
gammas <- lapply(list(naive.b.cells, mem.b.cells), function(keep.cells){
  # extracting alpha,beta matrices for the relevant cell types.
  # and the peaks and genes within a query region
  # Note: I had to densify matrices because it's very slow to query
  # a sparse matrix
  alpha <- as.matrix(atac.counts[peak.hits, keep.cells])
  beta <- as.matrix(rna.counts[gene.hits, keep.cells])

  # reducing the max.it to 100k, 1M can be slow
  GibbsGamma(overlaps, alpha, beta, max.it = 100000)
})

gamma.naive <- gammas[[1]]
gamma.mem <- gammas[[2]]

# Marginalizing over P to show DE genes
dims = dim(gamma.mem)
df.genes <- data.frame(colSums(gamma.naive[-dims[1], ]), colSums(gamma.mem[-dims[1], ]))
rownames(df.genes) <- colnames(gamma.naive)
df.genes <- df.genes[-dims[2], ]
df.genes

# Marginalizing over G to show DE peaks
df.peaks <- data.frame(rowSums(gamma.naive[, -dims[2]]), rowSums(gamma.mem[, -dims[2]]))
rownames(df.peaks) <- rownames(gamma.naive)
df.peaks <- df.peaks[-dims[1], ]
colnames(df.peaks) <- c("naive", "mem")
df.peaks["diff"] <- df.peaks[, "mem"] - df.peaks[, "naive"]
tail(df.peaks[order(abs(df.peaks$diff)), ], 20)

#Magic test
pname <- "chr8-19816678-19818513"
df.naive <- gamma.naive[pname, -dims[2]]
df.mem <- gamma.mem[pname, -dims[2]]
df.naive <- df.naive[df.naive != 0]
df.mem <- df.mem[df.mem != 0]
df.naive
df.mem

# Signac Plot
CoveragePlot(
  object = pbmc,
  region = pname,
  links=FALSE,
  extend.upstream = 10000,
  extend.downstream = 10000,
  group.by = "predicted.id",
  features = gname
)

##############################
##############################
##### Helper Functions #######
##### Load them Before running
##### the above code.
###############################
##############################

SubsetCounts <- function(chr.name, object) {
  annotation <- pbmc@assays$macs.atac@annotation
  keep.gnames <- unique(annotation[annotation@seqnames == chr.name]$gene_name)
  keep.gnames <- intersect(keep.gnames, rownames(pbmc@assays$RNA@counts))
  keep.peaks <- grep(paste0(chr.name, "-"), rownames(pbmc@assays$macs.atac@counts))

  atac.counts <- pbmc@assays$macs.atac@counts[keep.peaks, ]
  rna.counts <- pbmc@assays$RNA@counts[keep.gnames, ]
  list(atac.counts, rna.counts)
}

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

GibbsGamma <- function(overlaps, alpha, beta, max.it = 1000000) {
  phits <- c(rownames(alpha), "pnoise")
  ghits <- c(rownames(beta), "gnoise")

  gamma <- matrix(0, nrow = length(phits), ncol = length(ghits))
  rownames(gamma) <-phits
  colnames(gamma) <- ghits
  cat("Generating gamma matrix of dimension", dim(gamma))

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
      # if there is nowhere to move assign noisy peak
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
      # if there is nowhere to move assign noisy gene
      start[2] <- num.genes
    } else {
      random.num <- runif(1, min = 0, max = 1)
      gene.name <- names(gdist)[Position(function(x) x > random.num, gdist)]
      start[2] <- Position(function(x) x == gene.name, ghits)
    }

    # Increment the count of the new destination state.
    gamma[start[1], start[2]] <-  gamma[start[1], start[2]] + 1
  }

  return (gamma)
}
