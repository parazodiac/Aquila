library(Seurat)
library(Signac)
library(princurve)

library(devtools)
devtools::load_all()

library(reshape2)
library(ggplot2)

library(SingleCellExperiment)
library(slingshot)
library(destiny, quietly = TRUE)
library(mclust, quietly = TRUE)
library(RColorBrewer)


# The file is on O5
pbmc <- readRDS(file = "/home/srivastavaa/avi/tim/peak_test/pbmc.Rds")

# extract cells with naive and memory B cell annotation
naive.b.cells <- names(pbmc$predicted.id[pbmc$predicted.id  == "Naive B"])
mem.b.cells <- names(pbmc$predicted.id[pbmc$predicted.id  == "Memory B"])
paste(length(naive.b.cells), length(mem.b.cells))
# "434 516"

naive.gammas <- GenerateGamma(signac.object = pbmc, chr.name = "chr14",
                              keep.cells = naive.b.cells)
mem.gammas <- GenerateGamma(signac.object = pbmc, chr.name = "chr14",
                            keep.cells = mem.b.cells)
DiffPlot(naive.gammas, mem.gammas, "TCL1A")

##############################################
signac.object <- pbmc
atac.assay <- "macs.atac"
chr.name <- "chr14"
rna.assay <- "RNA"
tau <- 100000

annotations <- signac.object[[atac.assay]]@annotation

keep.peaks <- grep(paste0(chr.name, "-"), rownames(signac.object[[atac.assay]]))
keep.genes <- unique(annotations[annotations@seqnames == chr.name]$gene_name)
keep.genes <- intersect(keep.genes, rownames(signac.object[[rna.assay]]))

atac.counts <- signac.object[[atac.assay]]@counts[keep.peaks, ]
rna.counts <- signac.object[[rna.assay]]@counts[keep.genes, ]

overlaps <- GetTauOverlaps(peak.names = rownames(atac.counts),
                           gene.names = rownames(rna.counts),
                           gene.ranges = annotations,
                           tau)

############################
## pseudotime
############################
time.method <- "slingshot"
if (time.method == "slingshot") {
  rna.counts <- pbmc@assays$RNA@counts[, c(naive.b.cells, mem.b.cells)]
  pseudotime <- GetSlingshotPseudotime(rna.counts)
  pseudotime <-pseudotime[order(pseudotime)]
} else {
  pseudotime <- PrincipleTime(pbmc@reductions$umap.rna@cell.embeddings[c(naive.b.cells, mem.b.cells), ])
  pseudotime <-pseudotime[order(pseudotime)]
}
length(pseudotime)

#############################

gammas <- lapply(seq(length(pseudotime)), function(cell.idx) {
  cat("\r", cell.idx)
  cell.name <- names(pseudotime)[cell.idx]
  GenerateGamma(signac.object = pbmc, chr.name = "chr14", keep.cells = cell.name,
                gene.name = "TCL1A", overlaps = overlaps, verbose=F)
})

vals <- lapply(as.vector(gammas), function(gamma) {
  #data.frame(gamma)["chr14-95156441-95158788", "DICER1"]
  data.frame(gamma)["chr14-95713613-95714401", "TCL1A"]
})
plot(unlist(vals))

diff <- naive.gammas[[40]] - mem.gammas[[40]]
which.max(diff[, "TCL1A"])












