library(Seurat)
library(Signac)
library(princurve)

library(devtools)
devtools::load_all()

library(dplyr)
library(reshape2)
library(ggplot2)

library(SingleCellExperiment)
library(slingshot)
library(destiny, quietly = TRUE)
library(mclust, quietly = TRUE)
library(RColorBrewer)
library(magick)

# The file is on O5
pbmc <- readRDS(file = "/home/srivastavaa/avi/tim/peak_test/pbmc_tim_rep1/pbmc.Rds")
pbmc@assays$macs.atac@fragments[[1]]@path <- "/home/srivastavaa/avi/tim/peak_test/pbmc_tim_rep1/atac_fragments.tsv.gz"

# extract cells with naive and memory B cell annotation
naive.b.cells <- names(pbmc$predicted.id[pbmc$predicted.id  == "Naive B"])
mem.b.cells <- names(pbmc$predicted.id[pbmc$predicted.id  == "Memory B"])
paste(length(naive.b.cells), length(mem.b.cells))
# "434 516"

############################
## pseudotime
############################
time.method <- "slingshot"
time.method <- ""
if (time.method == "slingshot") {
  rna.counts <- pbmc@assays$RNA@counts[, c(naive.b.cells, mem.b.cells)]
  pseudotime <- GetSlingshotPseudotime(rna.counts)
} else {
  rna.embeddings <- pbmc@reductions$umap.rna@cell.embeddings[c(naive.b.cells, mem.b.cells), ]
  pseudotime <- PrincipleTime(rna.embeddings)
}
length(pseudotime)
disc.pseudotime <- round(pseudotime, digits = 1)
pbmc$pseudotime <- disc.pseudotime
DimPlot(subset(pbmc, cells=c(naive.b.cells, mem.b.cells)), group.by = "pseudotime")

##############################################
##############################################
##############################################

chr.name <- "chr14"
gname <- "TCL1A"
is.high <- TRUE
naive.gammas <- GenerateGamma(signac.object = pbmc, chr.name = chr.name,
                              keep.cells = sample(naive.b.cells, size = 217) )
mem.gammas <- GenerateGamma(signac.object = pbmc, chr.name = chr.name,
                            keep.cells = sample(mem.b.cells, size = 258))
DiffPlot(naive.gammas, mem.gammas, gname)
RelDiffPlot(naive.gammas, mem.gammas, gname)
FCPlot(naive.gammas, mem.gammas, gname)
FCPlot(mem.gammas, naive.gammas, gname)

##############################################
gname <- "TCL1A"
pname <- "chr14-95713613-95714401"

df <- data.frame(cbind(pbmc@assays$macs.atac@data[pname, names(disc.pseudotime)],
                       pbmc@assays$RNA@data[gname, names(disc.pseudotime)])
                 )
df["time"] <- disc.pseudotime[rownames(df)]
colnames(df) <- c(gsub("-", "_", pname), gname, "time")
dim(df)
head(df)

p1 <- ggplot(data=df, aes_string(x="time", y=gsub("-", "_", pname))) +
  geom_point()  + geom_smooth(method="gam")
p2 <- ggplot(data=df, aes_string(x="time", y=gname)) + geom_point() + geom_smooth(method="gam")
p1 / p2


gname <- "CLMN"
pnames <- overlaps[overlaps$genes == gname, ]$peaks
pnames <- rownames(naive.gammas[[idx]])
RegionPlot(pbmc, mem.gammas, gname, pnames)

##############################################
idx <- which(unlist(lapply(naive.gammas, function(gamma) { gname %in% colnames(gamma) })))
diff <- naive.gammas[[idx]] - mem.gammas[[idx]]
if (is.high) {
  pname <- names(which.max(diff[, gname]))
} else {
  pname <- names(which.min(diff[, gname]))
}

atac.assay <- "macs.atac"
rna.assay <- "RNA"
tau <- 100000

annotations <- pbmc[[atac.assay]]@annotation

keep.peaks <- grep(paste0(chr.name, "-"), rownames(pbmc[[atac.assay]]))
keep.genes <- unique(annotations[annotations@seqnames == chr.name]$gene_name)
keep.genes <- intersect(keep.genes, rownames(pbmc[[rna.assay]]))

atac.counts <- pbmc[[atac.assay]]@counts[keep.peaks, ]
rna.counts <- pbmc[[rna.assay]]@counts[keep.genes, ]

overlaps <- GetTauOverlaps(peak.names = rownames(atac.counts),
                           gene.names = rownames(rna.counts),
                           gene.ranges = annotations,
                           tau)

####################################

roi.comb <- tibble(Var1 = character(), Var2 = character(), value = numeric())
for(reg.id in  seq(length(naive.gammas))){
  if (length(naive.gammas[[reg.id]]) != 0) {
    df <- melt( naive.gammas[[reg.id]] - mem.gammas[[reg.id]] )
    roi.comb <- roi.comb %>% add_row(df[df$value > (max(df[df$value > 0, ]$value) * 0.50), ])
    roi.comb <- roi.comb %>% add_row(df[df$value < (min(df[df$value < 0, ]$value) * 0.50), ])
  }
}
head(roi.comb)
dim(roi.comb)
roi.comb[roi.comb$Var2 == "TCL1A", ]

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



















