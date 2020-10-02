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



#########################################
## make gif
#########################################
vals <- lapply(seq(length(gammas)), function(g.idx) {
  gamma <- data.frame(gammas[g.idx])
  df <- melt(as.matrix(gamma))
  colnames(df) <- c("peak", "gene", "mass")
  com.idx <- round(sum(df$mass * seq(dim(df)[1])))

  gamma[,] <- 0
  gamma[as.character(df[com.idx, ]$peak), as.character(df[com.idx, ]$gene)] <- 1
  p <- ggplot(melt(as.matrix(gamma)), aes(Var2, Var1)) +
    geom_tile(aes(fill=value), colour="white")

  ggsave(plot = p,
         filename = paste0("files/", g.idx, "_file.png"),
         device = "png")
})


imgs <- list.files("~/avi/tim/Circus/files", full.names = TRUE)
img_list <- lapply(imgs, image_read)

img_joined <- image_join(img_list)
img_animated <- image_animate(img_joined, fps = 2)
#img_animated

image_write(image = img_animated,
            path = "files/tx-sales.gif")



#########################################
## Gamma activation
#########################################
all_vals <- lapply(seq(dim(roi.comb)[1]), function(row.idx){
  print(row.idx)
  pname <- roi.comb[row.idx, ]$Var1
  gname <- roi.comb[row.idx, ]$Var2

  gammas <- lapply(seq(length(table(disc.pseudotime))), function(cell.idx) {
    cat("\r", cell.idx)
    cell.name <- names(disc.pseudotime[disc.pseudotime == ((cell.idx-1) / 10)])
    GenerateGamma(signac.object = pbmc, chr.name = chr.name, keep.cells = cell.name,
                  overlaps = overlaps, gene.name = gname, verbose=F)
  })

  vals <- unlist(lapply(gammas, function(gamma) {
    data.frame(gamma)[pname, gname]
  }))
  vals <- (vals - min(vals))
  vals <- vals / max(vals)
  vals
})

#plot(vals, main=gname, ylab="gamma", xlab="pseudotime")

df <- data.frame(t(data.frame(all_vals)))
colnames(df) <- paste0("t",seq(dim(df)[2]))
rownames(df) <- paste(roi.comb$Var1, roi.comb$Var2)

order <- c()
for (name in colnames(df)) {
  order <- c(order, rownames(df[df[, name] == 1.0, ]))
}
length(order)

hmap <- as.matrix(df[order, ])
heatmap(hmap, Rowv = NA, Colv = NA, xlab = "pseudotime", main = "gamma_activation")
df









