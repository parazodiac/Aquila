library(Seurat)
library(EnsDb.Hsapiens.v86)
library(GenomeInfoDb)
library(ggplot2)
library(patchwork)
library(princurve)

#install.packages("~/avi/tim/signac-private/", repos = NULL, type = "source")
library(Signac)

#install.packages("~/avi/tim/Circus/", repos = NULL, type = "source")
devtools::load_all()
#library(Circus)

###########################
pbmc <- GetSignacObject(mrna.h5.file = "~/avi/tim/peak_test/filtered_feature_bc_matrix.h5",
                        peaks.fragments.file = "~/avi/tim/peak_test/atac_fragments.tsv.gz")

# filter cells
VlnPlot(object = pbmc, features = c('nCount_ATAC', 'nCount_RNA'), pt.size = 0.1, ncol = 2)
pbmc <- subset(x = pbmc, subset = nCount_ATAC < 100000 & nCount_RNA < 25000 &
                 nCount_ATAC > 1000 & nCount_RNA > 1000 )

pbmc <- ProcessATAC(pbmc, "ATAC")

# add new macs assay
pbmc <- AddAssayFromPeaks(pbmc,
                          "~/avi/tim/peak_test/peaks/bulk_peaks.narrowPeak",
                          "macs.atac")
pbmc <- ProcessATAC(pbmc, "macs.atac")

# process RNA
pbmc <- ProcessRNA(pbmc)
pbmc

###########################
# extract annotations from 10x v3 RNA
pbmc <- TransferRnaPBMC(pbmc)
DimPlot(pbmc, group.by = "predicted.id", label = T)

######################
# The sexy coverage plot

DefaultAssay(pbmc) <- "macs.atac"
CoveragePlot(
  object = pbmc,
  region = "BACH2",
  links=FALSE,
  extend.upstream = 50000,
  extend.downstream = 50000
)

########################
# August 24, 2020
#saveRDS(object = pbmc, file = "~/avi/tim/peak_test/pbmc.Rds")
pbmc <- readRDS(file = "~/avi/tim/peak_test/pbmc.Rds")
pbmc

#saveRDS(object = atac.counts, file = "~/avi/tim/peak_test/atac.Rds")
#saveRDS(object = rna.counts, file = "~/avi/tim/peak_test/rna.Rds")
#atac.counts <- readRDS(file = "~/avi/tim/peak_test/atac.Rds")
#rna.counts <- readRDS(file = "~/avi/tim/peak_test/rna.Rds")

## Validations
#test <- CreateSeuratObject(rna.counts[, combined.cells])
#pseudotime <- PrincipleTime(cell.embeddings = pbmc@reductions$umap.rna@cell.embeddings[combined.cells, ])
#test$pseudotime <- pseudotime
#test[["umap"]] <- CreateDimReducObject(
#  pbmc@reductions$umap.rna@cell.embeddings[combined.cells, ],
#  key = "UMAP"
#)
#
#test$pseudotime <- round(pseudotime, digits = 1)
#test$predicted.id <- pbmc$predicted.id[combined.cells]
#
#DimPlot(test, group.by = "pseudotime")
#DimPlot(test, group.by = "predicted.id", label = T)

# filtering things for chromosome 1
gname <- "TCL1A"
counts <- SubsetCounts("chr14", pbmc)

gname <- "OSBPL10"
counts <- SubsetCounts("chr3", pbmc)

gname <- "EEF1A1"
counts <- SubsetCounts("chr6", pbmc)

atac.counts <- counts[[1]]
rna.counts <- counts[[2]]
dim(atac.counts)
dim(rna.counts)
############################

mem.b.cells <- names(pbmc$predicted.id[pbmc$predicted.id  == "Memory B"])
naive.b.cells <- names(pbmc$predicted.id[pbmc$predicted.id  == "Naive B"])

plot(rna.counts[gname, mem.b.cells])
plot(rna.counts[gname, naive.b.cells])

markers <- FindMarkers(pbmc, assay = "SCT", ident.1 = mem.b.cells, ident.2 = naive.b.cells, only.pos = T)
head(markers, 20)

########################

combined.cells <- c(mem.b.cells, naive.b.cells)
length(combined.cells)

devtools::load_all()
modes <- list("dot", "dot_depth", "ks", "spearman")
gammas.test <- lapply(modes, function(use.mode) {
  print(use.mode)
  gammas <- lapply(list(combined.cells, mem.b.cells, naive.b.cells), function(keep.cells) {
    gamma <- SubsetGamma(pbmc, keep.cells, atac.counts, rna.counts, use.mode)
    print(dim(gamma))
    gamma[gamma$subjectHits == gname, ]
  })
  list(use.mode, gammas)
})

idx <- 2
mode <- gammas.test[[idx]][[1]]
gammas <-gammas.test[[idx]][[2]]
print(mode)

gammas[[1]]$type <- "combined"
gammas[[2]]$type <- "memory"
gammas[[3]]$type <- "naive"

ggplot(rbind(gammas[[1]], gammas[[2]], gammas[[3]]),
       aes(corr, fill=type)) + geom_density(alpha = 0.2) +
  ggtitle(paste(mode, gname, ks.test(gammas[[2]]$corr, gammas[[3]]$corr)[[1]]))

########################
annotation[annotation$gene_name == "EEF1A1", ]
gamma[gamma$corr == max(gamma$corr), ]
library(ggplot2)
gname <- "TCL1A"
pname <- ""
PlotATACnRNA(atac.counts, rna.counts, pname, gname, naive.b.cells)
PlotATACnRNA(atac.counts, rna.counts, pname, gname, mem.b.cells)

















