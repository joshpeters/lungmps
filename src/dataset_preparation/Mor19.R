# ---
# Description: Prepare data from Morse et al. 2019
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
set <- "Morse_2019"

# [2] Load files ----

# grab file names
path <- "data/raw_data/GSE128033_RAW"
file_paths <- list.files(path, full.names = TRUE)
file_names <- list.files(path, full.names = FALSE)

sample_ids <- as.data.frame(str_match(file_paths,
  "^.*GSM[[:digit:]]{7}_(SC[[:digit:]]{2,3})(IPF|NOR|DNOR)(|UP|LOW|bal)_.*gz"))
samples <- data.frame(sample = sample_ids$V2, disease = sample_ids$V3, location = sample_ids$V4)
samples[] <- lapply(samples, as.character)
samples[samples == ""] <- "MIXED"
samples <- tidyr::unite(samples, unique_id, sample:location, na.rm = TRUE, remove = FALSE)
samples <- samples[seq(1, nrow(samples), 3), ]

matrices <- file_paths[grep("^.*mtx.gz", file_paths)]
genes <- file_paths[grep("^.*genes.tsv.gz", file_paths)]
barcodes <- file_paths[grep("^.*barcodes.tsv.gz", file_paths)]

files <- data.frame(matrix = matrices, genes = genes, barcodes = barcodes)
files <- cbind(files, samples)

# load in counts
counts <- list()
meta <- list()
files <- files[files$location != "bal", ]
for (i in seq_along(1:dim(files)[1])) {

  ui_todo("Loading file #{i}")
  matrix_file <- as.character(files[i, "matrix", drop = TRUE])

  matrix <- as(Matrix::readMM(file = gzfile(as.character(files[i, "matrix", drop = TRUE]))), "dgCMatrix")
  genes <- read_tsv(gzfile(as.character(files[i, "genes", drop = TRUE])),
    col_names = c("gene_id", "gene_symbol"), col_types = cols(gene_id = col_character(), gene_symbol = col_character()))
  barcodes <- read_tsv(gzfile(as.character(files[i, "barcodes", drop = TRUE])),
    col_names = c("barcode"), col_types = cols(barcode = col_character()))

  barcodes$barcode <- paste(barcodes$barcode, files$unique_id[i], sep = "_")
  rownames(matrix) <- genes$gene_symbol
  colnames(matrix) <- barcodes$barcode

  ui_info("Matrix dimensions: {dim(matrix)[1]}x{dim(matrix)[2]}")

  counts[[i]] <- matrix
  meta[[i]] <- data.frame(barcode = barcodes$barcode, sample_id = rep(files$unique_id[i], length(colnames(matrix))),
    donor = rep(files$sample[i], length(colnames(matrix))),
    disease = rep(files$disease[i], length(colnames(matrix))),
    location = rep(files$location[i], length(colnames(matrix))))
}

reduced_counts <- pblapply(counts, function(x) {
  x <- x[, Matrix::colSums(x) > 100]
})

matrix <- do.call(cbind, reduced_counts)
head(colnames(matrix))

meta <- do.call(rbind, meta)
rownames(meta) <- meta$barcode
meta <- meta[colnames(matrix), ]

dim(matrix)
dim(meta)
all.equal(rownames(meta), colnames(matrix))

# [3] Prepare object ----

object <- CreateSeuratObject(matrix, project = set, names.field = 2, names.delim = "_", min.cells = 0, min.features = 100, meta.data = meta)
object <- AddExprMeta(object)

names(object[[]])
head(object[[]])

object$batch <- object$sample_id
object$study <- set
colnames(object@meta.data)[grep("barcode", colnames(object@meta.data))] <- "cell_barcode"
object <- AddBatchFreq(object, "batch")
table(object$batch_frequency < 100)
object <- Scrub(object, batch.var = "batch")

names(object[[]])
head(object[[]])

object <- object[, Cells(object)[object$donor %in% c("SC56", "SC59", "SC155", "SC156",
  "SC87", "SC88", "SC93", "SC94", "SC153", "SC154")]]

# Harmonizing cell_barcode, study, donor, condition, sample_id, batch, disease_classification, celltype
colnames(object[[]])
head(object[[]])
colnames(object@meta.data)[grep("^disease$", colnames(object@meta.data))] <- "condition"
object$disease_classification <- "disease"
table(object$condition)
object$disease_classification[object$condition == "NOR"] <- "health"

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
