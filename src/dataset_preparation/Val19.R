# ---
# Description: Prepare data from Valenzi et al. 2019
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
set <- "Valenzi_2019"
Clean()

# [2] Load files ----

path <- "data/raw_data/GSE128169_RAW"
file_paths <- list.files(path, full.names = TRUE)
file_names <- list.files(path, full.names = FALSE)

sample_ids <- as.data.frame(str_match(file_paths,
  "^.*GSM[[:digit:]]{7}_(SC[[:digit:]]{2,3})(SSC)(LOW|UP)_.*gz"))
sample_ids <- sample_ids[, -1]
sample_ids[] <- apply(sample_ids, 2, as.character)
sample_ids <- sample_ids[!duplicated(sample_ids), ]
sample_ids <- sample_ids[!is.na(sample_ids$V2), ]
colnames(sample_ids) <- c("sample_id", "condition", "location")
sample_ids$patient <- rep(c("5", "6", "10", "13"), each = 2)

matrices <- file_paths[grep("^.*mtx.gz", file_paths)]
genes <- file_paths[grep("^.*genes.tsv.gz", file_paths)]
barcodes <- file_paths[grep("^.*barcodes.tsv.gz", file_paths)]

files <- data.frame(matrix = matrices, genes = genes, barcodes = barcodes)
files <- cbind(files, sample_ids)

# load in counts
counts <- list()
meta <- list()

for (i in seq_along(1:dim(files)[1])) {

  message(glue("Loading file #{i}"))
  matrix_file <- as.character(files[i, "matrix", drop = TRUE])
  message(glue("Loading matrix: {matrix_file}"))

  matrix <- as(Matrix::readMM(file = gzfile(as.character(files[i, "matrix", drop = TRUE]))), "dgCMatrix")
  genes <- read_tsv(gzfile(as.character(files[i, "genes", drop = TRUE])),
    col_names = c("gene_id", "gene_symbol"), col_types = cols(gene_id = col_character(), gene_symbol = col_character()))
  barcodes <- read_tsv(gzfile(as.character(files[i, "barcodes", drop = TRUE])),
    col_names = c("barcode"), col_types = cols(barcode = col_character()))

  barcodes$barcode <- paste(barcodes$barcode, files$sample_id[i], sep = "_")
  rownames(matrix) <- genes$gene_symbol
  colnames(matrix) <- barcodes$barcode

  message(glue("Matrix dimensions: {dim(matrix)[1]}x{dim(matrix)[2]}"))

  counts[[i]] <- matrix
  length(colnames(matrix)) == length(barcodes$barcode)
  meta[[i]] <- data.frame(barcode = barcodes$barcode,
    sample_id = rep(files$sample_id[i], length(colnames(matrix))),
    condition = rep(files$condition[i], length(colnames(matrix))),
    location = rep(files$location[i], length(colnames(matrix))),
    patient = rep(files$patient[i], length(colnames(matrix))))
}

reduced_counts <- lapply(counts, function(x) {
  x <- x[, Matrix::colSums(x) > 100]
})

matrix <- do.call(cbind, reduced_counts)
head(colnames(matrix))
dim(matrix)

meta <- do.call(rbind, meta)
head(meta)
dim(meta)

colnames(meta) <- c("cell_barcode", "sample_id", "condition", "location", "patient")
rownames(meta) <- meta$cell_barcode
meta <- meta[colnames(matrix), ]
rownames(meta) <- meta$cell_barcode

dim(matrix)
dim(meta)

all.equal(rownames(meta), colnames(matrix))
table(is.na(meta))

# [3] Prepare object ----

object <- CreateSeuratObject(matrix, project = set, names.field = 2, names.delim = "__",
  min.cells = 0, min.features = 100, meta.data = meta)
object <- AddExprMeta(object)

names(object[[]])
head(object[[]])
object$batch <- paste0(object$sample_id, "_", object$location)
length(unique(object$batch))

object$study <- set
object$donor <- object$sample_id
object$sample_id <- object$batch
object$orig.ident <- object$donor

object <- AddBatchFreq(object, "batch")
table(object$batch_frequency < 100)
object <- Scrub(object, batch.var = "batch")

# Harmonizing cell_barcode, study, donor, condition, sample_id, batch, disease_classification, celltype
colnames(object[[]])
head(object[[]])
object$disease_classification <- "disease"
table(object$condition)

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
