# ---
# Description: Prepare data from Madissoon et al. 2019
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
set <- "Madissoon_2019"
Clean()

# [2] Load files ----

object <- readRDS("data/raw_data/lung_ts.rds")

meta <- object@meta.data
meta <- meta %>% janitor::clean_names()
head(meta)
meta$cell_barcode <- rownames(meta)
meta <- meta[, c(1, 2, 3, 4, 5, 6, 10, 11)]

# [3] Prepare object ----

object <- CreateSeuratObject(object@assays$RNA@counts, project = set, min.cells = 0, min.features = 100, meta.data = meta)
object$batch <- object$donor_time
object <- AddExprMeta(object)
names(object[[]])
head(object[[]])
object$study <- set
object$cell_barcode <- Cells(object)

object <- AddBatchFreq(object, batch = "batch")
table(object$batch_frequency < 100)
object <- Scrub(object, batch.var = "batch")

# Harmonizing cell_barcode, study, donor, condition, sample_id, batch, disease_classification, celltype
colnames(object[[]])
head(object[[]])
table(object$sample)
colnames(object@meta.data)[grep("^sample$", colnames(object@meta.data))] <- "sample_id"
colnames(object@meta.data)[grep("^patient$", colnames(object@meta.data))] <- "donor"
colnames(object@meta.data)[grep("^celltypes$", colnames(object@meta.data))] <- "celltype"
object$disease_classification <- "health"
object$condition <- "health"

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

# [4] Write files ----
saveRDS(object, glue("data/objects/{set}_base_object.rds"))
Write10X(object, rootpath = glue("data/10x/{set}/"))
