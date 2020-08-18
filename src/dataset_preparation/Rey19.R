# ---
# Description: Prepare data from Reyfman et al. 2019
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
set <- "Reyfman_2019"

# [2] Load files ----

path <- "data/raw_data/GSE122960_RAW"
file_paths <- list.files(path, full.names = TRUE)
file_names <- list.files(path, full.names = FALSE)
file_paths

filtered_files <- file_paths[grep("filtered", file_paths)]
sample_ids <- str_match(filtered_files, "GSM[0-9]{7}_(.*)_([0-9]{2})_.*")
sample_ids <- as.data.frame(sample_ids)
colnames(sample_ids) <- c("file", "donor_type", "donor_number")
sample_ids <- sample_ids[-5, ]
filtered_files <- filtered_files[-5]

counts <- list()
meta <- list()

for (i in seq_along(filtered_files)) {
  ui_todo("Loading file #{i}: {filtered_files[i]}\nPatient data: {sample_ids[i, 2]}, {sample_ids[i, 3]}\n")
  c <- Read10X_h5(filtered_files[i])
  colnames(c) <- paste(sample_ids[i, 2], sample_ids[i, 3], colnames(c), sep = "_")
  counts[[i]] <- c
  meta[[i]] <- data.frame(barcode = colnames(c), sample_id = rep(paste(sample_ids[i, 2], sample_ids[i, 3], sep = "_"), length(colnames(c))))
}

lapply(counts, dim)
reduced_counts <- lapply(counts, function(x) {
  x <- x[, Matrix::colSums(x) > 100]
})

matrix <- do.call(cbind, reduced_counts)
meta <- do.call(rbind, meta)
rownames(meta) <- meta$barcode
all.equal(rownames(meta), colnames(matrix))
head(colnames(matrix))

# [3] Prepare object ----

object <- CreateSeuratObject(matrix, project = set, names.field = 2, names.delim = "_",
  min.cells = 0, min.features = 100, meta.data = meta)
object <- AddExprMeta(object)

names(object[[]])
head(object[[]])

object$batch <- object$sample_id
object$study <- set
colnames(object@meta.data)[grep("barcode", colnames(object@meta.data))] <- "cell_barcode"
object <- AddBatchFreq(object, "batch")
table(object$batch_frequency < 100)
object <- Scrub(object, batch.var = "batch")

# Harmonizing cell_barcode, study, donor, condition, sample_id, batch, disease_classification, celltype
colnames(object[[]])
head(object[[]])
table(object$sample_id)
object$condition <- str_match(object$sample_id, "(.*)\\_([[:digit:]]{1,2})")[, 2, drop = TRUE]
object$donor_number <- str_match(object$sample_id, "(.*)\\_([[:digit:]]{1,2})")[, 3, drop = TRUE]
table(object$condition)
object$disease_classification <- "disease"
object$disease_classification[object$condition == "Donor"] <- "health"
object$donor <- object$batch

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
