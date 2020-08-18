# ---
# Description: Prepare data from Lambrechts et al. 2018
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
set <- "Lambrechts_2018"
Clean()

# [2] Load files ----

lfile <- loomR::connect(filename = "data/raw_data/Thienpont_Tumors_52k_v4_R_fixed.loom", mode = "r+")
matrix <- Matrix::Matrix(lfile$matrix[, ], sparse = TRUE)

lfile[["col_attrs"]]
lfile[["row_attrs"]]
genes <- lfile$row.attrs$Gene[]
attrs <- names(lfile$col.attrs)
meta <- lfile$get.attribute.df(MARGIN = 2)
attrs <-  c("CellFromTumor", "CellID", "ClusterID", "ClusterName", "PatientNumber")
meta <- lfile$get.attribute.df(attributes = attrs,
  MARGIN = 2)

meta <- meta %>% janitor::clean_names()
colnames(meta)[grep("cell_id", colnames(meta))] <- "cell_barcode"
rownames(meta) <- meta$cell_barcode

dim(matrix)
head(rownames(matrix))
head(matrix)
matrix <- Matrix::t(matrix)
rownames(matrix) <- genes
colnames(matrix) <- meta$cell_barcode
dim(matrix)
all.equal(rownames(meta), colnames(matrix))

# [3] Prepare object ----

object <- CreateSeuratObject(matrix, project = set, names.field = 2, names.delim = "_",
  min.cells = 0, min.features = 100, meta.data = meta)

object$cell_from_tumor <- as.logical(as.numeric(object$cell_from_tumor))
object$batch <- object$patient_number

object <- AddExprMeta(object)

names(object[[]])
head(object[[]])
object$study <- set

object <- AddBatchFreq(object, "batch")
table(object$batch_frequency < 100)
object <- Scrub(object, batch.var = "batch")

# Harmonizing cell_barcode, study, donor, condition, sample_id, batch, disease_classification, celltype
colnames(object[[]])
head(object[[]])
colnames(object@meta.data)[grep("cell_from_tumor", colnames(object@meta.data))] <- "condition"
colnames(object@meta.data)[grep("patient_number", colnames(object@meta.data))] <- "sample_id"
colnames(object@meta.data)[grep("cluster_name", colnames(object@meta.data))] <- "celltype"

colnames(object@meta.data)[grep("^disease$", colnames(object@meta.data))] <- "condition"
object$disease_classification <- "disease"
object$disease_classification[object$condition == FALSE] <- "health"

object$donor <- object$sample_id

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

# lambrechts <- object_ls[[object_names[3]]]
# View(lambrechts@meta.data)
# lambrechts$condition[lambrechts$condition == TRUE] <- "tumor"
# lambrechts$condition[lambrechts$condition == FALSE] <- "lung"
# table(lambrechts@meta.data[, c("condition", "disease_classification")])

# [4] Write files ----
saveRDS(object, glue("data/objects/{set}_base_object.rds"))
Write10X(object, rootpath = glue("data/10x/{set}/"))
