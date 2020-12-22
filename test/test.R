library(Matrix)
library(aquila)
library(Seurat)

obj <- readRDS("/home/haoy/ASAP/moran.i.test.rds")
adt.set <-c("CD4-1","CD14", "CD34", "Cadherin")
adt.value <- obj[["ADT"]]@data[adt.set, ]
snn.graph <- obj[["wsnn"]]

snn.graph <- Matrix(snn.graph, sparse = TRUE)
adt.value <- Matrix(adt.value, sparse = TRUE)

df <- MoransI(adt.value, snn.graph, "/home/srivastavaa/")
df
