library(Seurat)
library(Signac)

library(devtools)
devtools::load_all()

library(reshape2)
library(ggplot2)

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

diff <- naive.gammas[[4]] - mem.gammas[[4]]
hm.df <- setNames(melt(diff), c('peaks', 'genes', 'diffs'))

ggplot(hm.df, aes(genes, peaks)) +
  geom_tile(aes(fill=diffs), colour="white") +
  scale_fill_gradient(low="red", high="green")

