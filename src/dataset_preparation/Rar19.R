# ---
# Description: Prepare data from Raredon et al. 2019
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
set <- "Raredon_2019"
Clean()

# [2] Load files ----

path <- "data/raw_data/GSE133747_RAW"
file_paths <- list.files(path, full.names = TRUE)
file_names <- list.files(path, full.names = FALSE)

sample_ids <- as.data.frame(str_match(file_paths,
  "^.*GSM[[:digit:]]{7}_(.{4,5})\\.?barcodes.*gz"))
sample_ids <- sample_ids[!is.na(sample_ids$V1), ]
sample_ids <- as.character(sample_ids$V2)
sample_ids[1] <- "Hum1"
sample_ids[2] <- "Hum2"

matrices <- file_paths[grep("^.*mtx.gz", file_paths)]
genes <- file_paths[grep("^.*genes.tsv.gz", file_paths)]
barcodes <- file_paths[grep("^.*barcodes.tsv.gz", file_paths)]

files <- data.frame(matrix = matrices, genes = genes, barcodes = barcodes)
files$sample_name <- sample_ids

# load in counts
counts <- list()
meta <- list()
for (i in seq_along(1:dim(files)[1])) {

  message(glue("Loading file #{i}"))
  matrix_file <- as.character(files[i, "matrix", drop = TRUE])
  message(glue("Loading matrix: {matrix_file}"))

  matrix <- as(Matrix::readMM(file = gzfile(as.character(files[i, "matrix", drop = TRUE]))), "dgCMatrix")
  genes <- read_tsv(gzfile(as.character(files[i, "genes", drop = TRUE])),
    col_names = TRUE)
  genes <- genes[, 2, drop = TRUE]
  barcodes <- read_tsv(gzfile(as.character(files[i, "barcodes", drop = TRUE])),
    col_names = TRUE)
  barcodes <- barcodes[, 2, drop = TRUE]

  rownames(matrix) <- genes
  colnames(matrix) <- barcodes

  message(glue("Matrix dimensions: {dim(matrix)[1]}x{dim(matrix)[2]}"))

  counts[[i]] <- matrix
  meta[[i]] <- data.frame(barcode = barcodes, sample_id = rep(files$sample_name[i], length(colnames(matrix))))

}

reduced_counts <- lapply(counts, function(x) {
  x <- x[, Matrix::colSums(x) > 100]
})

matrix <- do.call(cbind, reduced_counts)
meta <- do.call(rbind, meta)
colnames(meta) <- c("cell_barcode", "sample_id")
rownames(meta) <- meta$cell_barcode
meta <- meta[colnames(matrix), ]
rownames(meta) <- meta$cell_barcode
all.equal(rownames(meta), colnames(matrix))

# [3] Prepare object ----

object <- CreateSeuratObject(matrix, project = set, names.field = 1, names.delim = "__",
  min.cells = 0, min.features = 100, meta.data = meta)
object <- AddExprMeta(object)

names(object[[]])
head(object[[]])
object$study <- set
object$donor <- object$sample_id
object$batch <- object$sample_id
object$condition <- "health"
object$disease_classification <- "health"

object <- AddBatchFreq(object, "batch")
table(object$batch_frequency < 100)
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