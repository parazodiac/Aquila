library(Seurat)
library(Signac)
library(princurve)
devtools::load_all()

pbmc <- readRDS(file = "~/avi/tim/peak_test/pbmc.Rds")
counts <- SubsetCounts("chr14", pbmc)

atac.counts <- counts[[1]]
rna.counts <- counts[[2]]
dim(atac.counts)
dim(rna.counts)

mem.b.cells <- names(pbmc$predicted.id[pbmc$predicted.id  == "Memory B"])
naive.b.cells <- names(pbmc$predicted.id[pbmc$predicted.id  == "Naive B"])
paste(length(mem.b.cells), length(naive.b.cells))

######################################
# Question 3

gamma.naive <- GetCellTypeGamma(atac.counts, rna.counts, keep.cells = naive.b.cells)
gamma.mem <- GetCellTypeGamma(atac.counts, rna.counts, keep.cells = mem.b.cells)

dim(gamma.naive)
dims <- dim(gamma.naive)
sum(gamma.naive[-dims[1], -dims[2]])
sum(gamma.mem[-dims[1], -dims[2]])

df.peaks <- data.frame(rowSums(gamma.naive[, -dims[2]]), rowSums(gamma.mem[, -dims[2]]))
rownames(df.peaks) <- rownames(gamma.naive)
df.peaks <- df.peaks[-dims[1], ]
tail(df.peaks)

df.genes <- data.frame(colSums(gamma.naive[-dims[1], ]), colSums(gamma.mem[-dims[1], ]))
rownames(df.genes) <- colnames(gamma.naive)
df.genes <- df.genes[-dims[2], ]
df.genes

pname <- "chr14-95326969-95328278"
df.naive <- gamma.naive[pname, -dims[2]]
df.mem <- gamma.mem[pname, -dims[2]]
df.naive <- df.naive[df.naive != 0]
df.mem <- df.mem[df.mem != 0]
df.naive
df.mem

df <- data.frame(rowSums(gamma.naive[, -dims[2]]), rowSums(gamma.mem[, -dims[2]]))
colnames(df) <- c("naive", "mem")
df["diff"] <- df[, "mem"] - df[, "naive"]
df["diff"] <- df["diff"] / (df[, "mem"] + df[, "naive"])
df[order(df$diff), ]

######################################
# Question 2
gamma <- matrix(0.0, nrow = length(rownames(alpha)), ncol = length(rownames(beta)))
rownames(gamma) <- rownames(alpha)
colnames(gamma) <- rownames(beta)
dim(gamma)
sum(gamma)

done <- 0
for (peak in rownames(alpha)) {
  done <- done + 1
  if (done %% 1000 == 0) {
    cat("\r", done, "/", length(rownames(alpha)))
  }

  if (sum(alpha[peak, ]) == 0) { next }

  gene.eqclass <- overlaps[overlaps$queryHits == peak, ]$subjectHits
  if (sum(beta[gene.eqclass, ]) == 0) { next }

  cell.id <- 0
  for ( source.count in alpha[peak, ] ) {
    cell.id <- cell.id + 1
    if (source.count == 0) { next }

    dest.count <- beta[gene.eqclass, cell.id]
    dest.sum <- sum(dest.count)
    if (dest.sum == 0) { next }

    dest.prob <- dest.count / dest.sum
    dest.count <- dest.prob * source.count

    for (gname in names(dest.count)) {
      gamma[peak, gname] <- gamma[peak, gname] + dest.count[gname]
    }
  }
}

sum(alpha)
sum(gamma)

gamma <- gamma / sum(gamma)
sum(gamma)

######################################
#%%%%%%%% Motif function
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
######################################
