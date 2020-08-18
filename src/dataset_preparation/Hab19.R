# ---
# Description: Prepare data from Habermann et al. 2019
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
set <- "Habermann_2019"

# [2] Load files ----

barcodes <- read_tsv(gzfile("data/raw_data/GSE135893_barcodes.tsv.gz"), col_names = FALSE)
barcodes <- barcodes$X1

genes <- read_tsv(gzfile("data/raw_data/GSE135893_genes.tsv.gz"), col_names = FALSE)
genes <- genes$X1

meta <- read_csv(gzfile("data/raw_data/GSE135893_IPF_metadata.csv.gz"), col_names = TRUE)
length(unique(meta$Sample_Name))
colnames(meta)[1] <- "cell_barcode"
meta <- meta %>% janitor::clean_names()
View(meta)
rownames(meta) <- meta$cell_barcode

matrix <- Matrix::readMM(gzfile("data/raw_data/GSE135893_matrix.mtx.gz"))
colnames(matrix) <- barcodes
rownames(matrix) <- genes

# [3] Prepare object ----

object <- CreateSeuratObject(matrix, project = set, names.field = 1, names.delim = "_",
  min.cells = 0, min.features = 100, meta.data = meta)
object <- AddExprMeta(object)

names(object[[]])
head(object[[]])
table(is.na(object$cell_barcode))
object <- object[, !is.na(object$cell_barcode)]

object$batch <- object$sample_name
object$study <- set
object <- AddBatchFreq(object, "batch")
table(object$batch_frequency < 100)

# Harmonizing cell_barcode, study, donor, condition, sample_id, batch, disease_classification, celltype
colnames(object[[]])
head(object[[]])
colnames(object@meta.data)[grep("^sample_name$", colnames(object@meta.data))] <- "sample_id"
table(object$orig_ident)
colnames(object@meta.data)[grep("^orig_ident$", colnames(object@meta.data))] <- "donor"
table(object$diagnosis)
colnames(object@meta.data)[grep("^diagnosis$", colnames(object@meta.data))] <- "condition"
table(object$status)

cols_to_remove <- grep("n_count", colnames(object@meta.data))
cols_to_remove <- c(cols_to_remove, grep("n_feature", colnames(object@meta.data)))
cols_to_remove
object@meta.data[, cols_to_remove] <- NULL

table(object$condition)
object$disease_classification <- "disease"
object$disease_classification[object$condition == "Control"] <- "health"

object <- Scrub(object, batch.var = "batch")

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
