# ---
# Description: Prepare data from Vieira Braga et al. 2019
# Author: Josh Peters
# ---

# [1] Basic setup ----

rm(list = ls())
library(future)
library(pbapply)
library(Seurat)
source("src/Functions.R")
source("src/Utilities.R")
source("src/Utilities.R")
stamp <- PrepEnv()
set <- "Braga_2019"

# [2] Load files ----

meta <- read_delim(gzfile("data/raw_data/GSE130148_barcodes_cell_types.txt.gz"),
  delim="\t", skip = 0, col_names = TRUE)
meta <- meta %>% janitor::clean_names()
meta <- meta[, c(1, 4, 6, 7, 8, 9)]
rownames(meta) <- meta$cell_barcode
colnames(meta)[grep("barcode", colnames(meta))] <- "cell_barcode"

raw_counts <- read_csv("data/raw_data/GSE130148_raw_counts.csv.gz")
raw_counts <- ConvertToSparseCounts(raw_counts)

# [3] Prepare object ----

object <- CreateSeuratObject(raw_counts, min.cells = 0, min.features = 100, meta.data = meta, project = set)
object <- AddExprMeta(object)

names(object[[]])
head(object[[]])

object$batch <- object$id
object$study <- set
object <- AddBatchFreq(object, "batch")
table(object$batch_frequency < 100)
object <- Scrub(object, batch.var = "batch")

names(object[[]])
head(object[[]])

genes <- ConvertGeneNames(genes = rownames(object), verbose = TRUE)
assertthat::assert_that(all(rownames(object) %in% rownames(genes)))
genes <- genes[rownames(object), ]
assertthat::assert_that(all.equal(rownames(object@assays$RNA@meta.features), rownames(genes)))
object@assays$RNA@meta.features <- cbind(object@assays$RNA@meta.features, genes)

assertthat::assert_that(!any(duplicated(genes$name)))
assertthat::assert_that(all.equal(rownames(object), genes$orig_name))
rownames(object@assays$RNA@counts) <- genes$name
rownames(object@assays$RNA@data) <- genes$name
object@assays$RNA@meta.features <- genes
object <- UpdateSeuratObject(object)

# Harmonizing donor, condition, sample_id, batch, celltype
colnames(object[[]])
head(object[[]])
colnames(object@meta.data)[grep("^id$", colnames(object@meta.data))] <- "sample_id"
colnames(object@meta.data)[grep("^orig_ident$", colnames(object@meta.data))] <- "donor"
object$condition <- "health"
object$disease_classification <- "health"

# [4] Write files ----
saveRDS(object, glue("data/objects/{set}_base_object.rds"))
Write10X(object, rootpath = glue("data/10x/{set}/"))
