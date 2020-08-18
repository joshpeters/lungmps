# ---
# Description: Prepare data from Mayr et al. 2020
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
set <- "Mayr_2020"
Clean()

# [2] Load files ----

object <- readRDS("data/objects/Mayr_2020_base_object.rds")
counts <- object@assays$RNA@counts
meta <- object@meta.data
View(meta)
colnames(meta)

# [3] Prepare object ----

object <- CreateSeuratObject(counts = counts, min.cells = 0, min.features = 100, meta.data = meta, project = set)
object <- AddExprMeta(object)

names(object[[]])
head(object[[]])

# Harmonizing cell_barcode, study, donor, condition, sample_id, batch, disease_classification, celltype
object$qc_nmads3 <- NULL
object$qc_nmads5 <- NULL
object$percent_mito <- NULL
object$percent_ribo <- NULL
object$percent_hb <- NULL

colnames(object[[]])
head(object[[]])
object$study <- set
object$sample_id <- object$donor

table(object$condition)
colnames(object@meta.data)[grep("^health_status$", colnames(object@meta.data))] <- "condition"
colnames(object@meta.data)[grep("^patient_id$", colnames(object@meta.data))] <- "donor"
colnames(object@meta.data)[grep("^cell_type$", colnames(object@meta.data))] <- "celltype"
object$disease_classification <- "disease"
object$disease_classification[object$condition == "control donor"] <- "health"

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