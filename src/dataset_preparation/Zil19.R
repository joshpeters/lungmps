# ---
# Description: Prepare data from Zilionis et al. 2019
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
set <- "Zilionis_2019"
Clean()

# [2] Load files ----

path <- "data/raw_data/GSE127465_RAW"
file_paths <- list.files(path, full.names = TRUE)
file_names <- list.files(path, full.names = FALSE)
file_paths

human_files <- file_paths[grep("GSM.*_human_", file_paths)]; human_files
sample_ids <- str_match(human_files, "^.*_(p[0-9][t|b][0-9])_.*$")
patient_ids <- str_match(human_files, "^.*_(p[0-9][t|b])[0-9]_.*$")

meta <- read_tsv(gzfile("data/raw_data/GSE127465_RAW/GSE127465_human_cell_metadata_54773x25.tsv.gz"))
meta <- janitor::clean_names(meta)

# load in counts
counts <- list()
for (i in seq_along(1:length(human_files))) {
  message(glue("Loading file #{i}: {human_files[i]}\nPatient data: {sample_ids[i, 2]}\n"))
  c <- data.table::fread(file = human_files[i], sep = "\t")
  sparse_c <- Matrix::Matrix(as.matrix(c[, -1]), sparse = TRUE)
  rownames(sparse_c) <- paste(c[[1]], sample_ids[i, 2], sep = "_")
  counts[[i]] <- sparse_c
}

merged_counts <- do.call(rbind, counts)
dim(merged_counts)
matrix <- Matrix::t(merged_counts)
dim(matrix)
matrix <- matrix[, colSums(matrix) > 100]
head(colnames(matrix))

meta$barcode_library <- paste(meta$barcode, meta$library, sep = "_")
head(meta$barcode_library)
rownames(meta) <- meta$barcode_library
length(unique(meta$barcode_library))

table(colnames(matrix) %in% meta$barcode_library)
table(meta$barcode_library %in% colnames(matrix))
reduced_matrix <- matrix[, meta$barcode_library]

# create object -----------------------------------------------------------

object <- CreateSeuratObject(reduced_matrix, project = set, names.field = 1, names.delim = "_",
  min.cells = 0, min.features = 100, meta.data = meta)
object <- AddExprMeta(object)

object$batch <- object$library
length(unique(object$batch))
table(object$batch)
unique(object$tissue)
object <- object[, WhichCells(object, expression = tissue == "tumor")]
object$cell_barcode <- Cells(object)

names(object[[]])
head(object[[]])
object$study <- set

object <- AddBatchFreq(object, "batch")
table(object$batch_frequency < 100)
object <- Scrub(object, batch.var = "batch")

# Harmonizing cell_barcode, study, donor, condition, sample_id, batch, disease_classification, celltype
colnames(object[[]])
head(object[[]])
colnames(object@meta.data)[grep("patient", colnames(object@meta.data))] <- "donor"
colnames(object@meta.data)[grep("major_cell_type", colnames(object@meta.data))] <- "celltype"
colnames(object@meta.data)[grep("tissue", colnames(object@meta.data))] <- "condition"
colnames(object@meta.data)[grep("library", colnames(object@meta.data))] <- "sample_id"

table(object$condition)
object$disease_classification <- "disease"

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