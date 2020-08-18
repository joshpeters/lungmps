# ---
# Description: Prepare data from Chua et al. 2020
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
set <- "Chua_2020"
Clean()

# [2] Load files ----

object <- readRDS("data/raw_data/covid_nbt_loc.rds")

# [3] Prepare object ----

names(object[[]])
head(object[[]])
object <- AddExprMeta(object)
object$percent.mt <- NULL
object$batch <- object$sample
object$study <- set
object$cell_barcode <- Cells(object)
object <- AddBatchFreq(object, "batch")
table(object$batch_frequency < 100)

# Harmonizing cell_barcode, study, donor, condition, sample_id, batch, disease_classification, celltype
colnames(object[[]])
head(object[[]])
table(object$infection)
table(object$patient)
table(object$severity)

colnames(object@meta.data)[grep("^patient$", colnames(object@meta.data))] <- "donor"
colnames(object@meta.data)[grep("^severity$", colnames(object@meta.data))] <- "condition"
colnames(object@meta.data)[grep("^sample$", colnames(object@meta.data))] <- "sample_id"

table(object$condition)
object$disease_classification <- "disease"
#object <- Scrub(object, batch.var = "batch")
table(object$batch)

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

object <- DietSeurat(object)
object@assays$integrated <- NULL

# [4] Write files ----
saveRDS(object, glue::glue("data/objects/{set}_base_object.rds"))
Write10X(object, rootpath = glue("data/10x/{set}/"))
