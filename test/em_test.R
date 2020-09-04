library(Seurat)
library(Signac)

library(devtools)
devtools::load_all()

# The file is on O5
pbmc <- readRDS(file = "/home/srivastavaa/avi/tim/peak_test/pbmc.Rds")

# extract ATAC/RNA counts for all peaks/genes in chromosome 14
gname <- "TCL1A"
chr.name <- "chr14"
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
overlaps <- GetTauOverlaps(peak.names = rownames(atac.counts),
                           gene.names = rownames(rna.counts),
                           gene.ranges = pbmc@assays$macs.atac@annotation,
                           tau = 100000)


regions <- GetConnectedRegions(overlaps)
gammas <- lapply(list(naive.b.cells, mem.b.cells), function(keep.cells){
  lapply(regions, function(gene.hits){
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

    roi.overlaps <- overlaps[overlaps$genes %in% gene.hits, ]
    roi.overlaps$genes <- genes[roi.overlaps$genes]
    roi.overlaps$peaks <- peaks[roi.overlaps$peaks]

    # reducing the max.it to 100k, 1M can be slow
    gamma <- getGamma(alpha, beta, roi.overlaps, 1000000, FALSE)
    rownames(gamma) <- names(peaks)
    colnames(gamma) <- names(genes)

    gamma
  })
})



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


######################################
######################################
#%%%%%%%% Motif function
######################################
######################################
library(TFBSTools)
library(JASPAR2018)
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg38)
features = StringToGRanges("chr14-95713613-95714401", sep = c("-", "-"))
pfm <- getMatrixSet(
  x = JASPAR2018,
  opts = list(species = 9606, all_versions = FALSE)
)
matched.motifs <- matchMotifs(pfm, features, genome = BSgenome.Hsapiens.UCSC.hg38)
