library(Seurat)
library(devtools)
load_all("~/avi/tim/signac-private/")

library(EnsDb.Hsapiens.v86)
library(GenomeInfoDb)
library(ggplot2)
library(patchwork)


annotations <- GetGRangesFromEnsDb(EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- "hg38"
annotations

# create Seurat object from gene counts
counts <- Read10X_h5("~/avi/tim/peak_test/filtered_feature_bc_matrix.h5")
pbmc <- CreateSeuratObject(counts = counts$`Gene Expression`, 
                           project = "coassay",
                           assay = "RNA")

# create chromatin assay from the peak counts
fragments <- "~/avi/tim/peak_test/atac_fragments.tsv.gz"
pbmc[["ATAC"]] <- CreateChromatinAssay(counts = counts$Peaks,
                                       min.cells = 5,
                                       sep = c(":", "-"),
                                       genome = "hg38",
                                       fragments = fragments,
                                       annotation = annotations)

# filter cells
pbmc <- subset(x = pbmc, subset = nCount_ATAC < 100000 & nCount_RNA < 25000 &
                 nCount_ATAC > 1000 & nCount_RNA > 1000 )

# Using Signac filters
DefaultAssay(pbmc) <- "ATAC"
pbmc <- NucleosomeSignal(pbmc)
pbmc <- TSSEnrichment(pbmc)

# perform SVD and generate UMAP for chromatin data
pbmc$blacklist <- FractionCountsInRegion(pbmc, regions = blacklist_hg38_unified)
pbmc <- FindTopFeatures(pbmc, min.cutoff = 10)
pbmc <- RunTFIDF(pbmc)
pbmc <- RunSVD(pbmc)
pbmc <- RunUMAP(pbmc, reduction = "lsi", dims = 2:40, reduction.name = "umap.atac")

## perform PCA and generate UMAP for RNA data
#DefaultAssay(pbmc) <- "RNA"
#pbmc <- SCTransform(pbmc, ncells = 5000)
#pbmc <- RunPCA(pbmc)
#pbmc <- RunUMAP(pbmc, dims = 1:40, reduction.name = "umap.rna")
#
## transfer labels
#pbmc_rna <- readRDS("~/pbmc_10k_v3_annotated.rds")
#pbmc_rna <- SCTransform(pbmc_rna, ncells = 5000)
#anchors <- FindTransferAnchors(reference = pbmc_rna, query = pbmc,
#                               normalization.method = "SCT")
#labels <- TransferData(anchorset = anchors, refdata = pbmc_rna$annotated,
#                       weight.reduction = pbmc[["pca"]])
#pbmc <- AddMetaData(pbmc, metadata = labels)
#
## write cell type classifications to file so we can split the bam
#predictions <- pbmc[[]]
#predictions <- predictions[, "predicted.id", drop = FALSE]
#
## remove whitespace
#predictions$predicted.id <- gsub(pattern = " ", replacement = "_", x = predictions$predicted.id)
#write.table(x = predictions,
#            file = "processed/coassay_celltypes_pbmc.tsv",
#            sep = "\t",
#            col.names = FALSE,
#            row.names = TRUE,
#            quote = FALSE)
#saveRDS(object = pbmc, file = "~/pbmc_multiomic.rds")

######################
# The sexy coverage plot

#library(Signac)
df <- read.table(
  file = "~/peak_test/bulk_peaks.narrowPeak",
  col.names = c("chr", "start", "end", "name", "score", "strand", "one", "two", "three", "four")
)
gr <- makeGRangesFromDataFrame(df = df)
p1 <- CoveragePlot(
  object = pbmc,
  region = "CD4",
  links=FALSE,
  extend.upstream = 50000,
  extend.downstream = 50000,
  peaks = gr
)
p1
######################

#pbmc <- readRDS("object/pbmc_multiomic.rds")
#DefaultAssay(pbmc) <- "ATAC"

########### pseudo bulk clustering
# quantify counts
new.counts <- FeatureMatrix(
  fragments = Fragments(pbmc),
  features = gr,
  cells = colnames(pbmc)
)

pbmc[["macs_pseudobulk"]] <- CreateChromatinAssay(
  counts = new.counts,
  genome = "hg38",
  fragments = Fragments(pbmc),
  annotation = Annotation(pbmc)
)

# Using Signac filters
DefaultAssay(pbmc) <- "macs_pseudobulk"
pbmc <- NucleosomeSignal(pbmc)
pbmc <- TSSEnrichment(pbmc)

# perform SVD and generate UMAP for chromatin data
pbmc$blacklist <- FractionCountsInRegion(pbmc, regions = blacklist_hg38_unified)
pbmc <- FindTopFeatures(pbmc, min.cutoff = 10)
pbmc <- RunTFIDF(pbmc)
pbmc <- RunSVD(pbmc)
pbmc <- RunUMAP(pbmc, reduction = "lsi", dims = 2:40, reduction.name = "umap.macs.bulk.atac")

########### cell type clustering
# read all peak files
peakfiles <- list.files(path = "/home/stuartt/data/10x_coassay/pbmc/", 
                        pattern = "*.narrowPeak", full.names = TRUE)
peaks.sep <- sapply(X = peakfiles, FUN = function(x){
  df <- read.table(
    file = x,
    col.names = c("chr", "start", "end", "name", "score", "dot", "v7", "v8", "v9", "v10")
  )
  gr <- makeGRangesFromDataFrame(df = df)
  return(gr)
})

unified <- Reduce(f = c, x = peaks.sep)
unified <- reduce(x = unified)
unified <- keepStandardChromosomes(unified, pruning.mode = "coarse")
sep.counts <- FeatureMatrix(
  fragments = Fragments(pbmc),
  features = unified,
  cells = colnames(pbmc)
)

pbmc[["macs_sep"]] <- CreateChromatinAssay(
  counts = sep.counts,
  genome = "hg38",
  fragments = Fragments(pbmc),
  annotation = Annotation(pbmc)
)

# Using Signac filters
DefaultAssay(pbmc) <- "macs_sep"
pbmc <- NucleosomeSignal(pbmc)
pbmc <- TSSEnrichment(pbmc)

# perform SVD and generate UMAP for chromatin data
#pbmc$blacklist <- FractionCountsInRegion(pbmc, regions = blacklist_hg38_unified)
pbmc <- FindTopFeatures(pbmc, min.cutoff = 10)
pbmc <- RunTFIDF(pbmc)
pbmc <- RunSVD(pbmc)
pbmc <- RunUMAP(pbmc, reduction = "lsi", dims = 2:40, reduction.name = "umap.macs.sep.atac")

pbmc
