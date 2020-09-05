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
#%%%%%%%% Pseudotime
########################

# August 24, 2020
#saveRDS(object = pbmc, file = "~/avi/tim/peak_test/pbmc.Rds")
pbmc <- readRDS(file = "~/avi/tim/peak_test/pbmc.Rds")
pbmc

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


gname <- "TCL1A"#"chr14"
gname <- "OSBPL10"#"chr3"
gname <- "EEF1A1"#"chr6"

########################
#%%%%%%%% Find Marker
########################

mem.b.cells <- names(pbmc$predicted.id[pbmc$predicted.id  == "Memory B"])
naive.b.cells <- names(pbmc$predicted.id[pbmc$predicted.id  == "Naive B"])

plot(rna.counts[gname, mem.b.cells])
plot(rna.counts[gname, naive.b.cells])

markers <- FindMarkers(pbmc, assay = "SCT", ident.1 = mem.b.cells, ident.2 = naive.b.cells, only.pos = T)
head(markers, 20)

########################
#%%%%%%%% ATAC/ RNA plot
########################

library(ggplot2)
gname <- "TCL1A"
pname <- ""
PlotATACnRNA(atac.counts, rna.counts, pname, gname, naive.b.cells)
PlotATACnRNA(atac.counts, rna.counts, pname, gname, mem.b.cells)

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

















