# ---
# Description: Prepare data from Laughney et al. 2020
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
set <- "Laughney_2020"
Clean()

# [2] Load files ----

file_paths <- list.files("data/raw_data/GSE123902_RAW", full.names = TRUE)
file_names <- list.files("data/raw_data/GSE123902_RAW", full.names = FALSE)
counts <- pbapply::pbsapply(file_paths, function(x) {

  mtx <- read_csv(gzfile(x), skip = 1,
    col_types = cols(X1 = col_character(), .default = col_integer()),
    trim_ws = TRUE, col_names = FALSE)
  head(mtx)
  barcodes <- mtx[, 1, drop = TRUE]
  length(barcodes)
  mtx <- Matrix::t(Matrix::Matrix(as.matrix(mtx[, -1]), sparse = TRUE))
  dim(mtx)

  genes <- read_csv(gzfile(x), n_max = 1,
    col_types = cols(.default = col_character()),
    trim_ws = TRUE, col_names = FALSE)
  genes <- as.character(genes[1, , drop = TRUE])[-1]
  length(genes)

  length(genes) == nrow(mtx)
  length(barcodes) == ncol(mtx)

  rownames(mtx) <- genes
  colnames(mtx) <- barcodes
  return(mtx)
})

sample_ids <- str_split(file_names, "_")
len <- sapply(sample_ids, length)
n <- max(len)
len <- n - len
sample_ids <- as.data.frame(mapply( function(x,y) c( x , rep( NA , y ) ) , sample_ids , len ))
sample_ids <- as.data.frame(t(sample_ids))
sample_ids <- sample_ids[, -c(2, 6)]
sample_ids[] <- lapply(sample_ids, as.character)
sample_ids[sample_ids$V4 == "PRIMARY", "V4"] <- "TUMOR"
sample_ids <- sample_ids[, -c(4)]
colnames(sample_ids) <- c("gsm", "sample_id", "condition")
sample_ids$sample_name <- paste0(sample_ids$sample_id, "_", sample_ids$condition)

# [3] Prepare object ----

objects <- pblapply(counts, function(x) {object <- CreateSeuratObject(x, project = set, min.cells = 0, min.features = 100); return(object)})
object <- merge(objects[[1]], objects[2:length(objects)], add.cell.ids = sample_ids$sample_name)
object <- AddExprMeta(object)

barcodes <- data.frame(barcodes = rownames(object@meta.data))
barcodes$sample_name <- str_match(barcodes$barcodes, "^(.{5,6}_(.*))_.*$")[, 2, drop = TRUE]
meta <- merge(x = barcodes, y = sample_ids, by = "sample_name", all.x = TRUE, all.y = FALSE)
nrow(meta) == nrow(barcodes)
rownames(meta) <- meta$barcodes
object <- AddMetaData(object, meta)

head(object[[]])
length(unique(object$orig.ident))

object$batch <- object$sample_name
length(unique(object$batch))

colnames(object@meta.data)[grep("barcode", colnames(object@meta.data))] <- "cell_barcode"
object$study <- set

object <- AddBatchFreq(object, "batch")
table(object$batch_frequency < 100)
object <- Scrub(object, batch.var = "batch")

names(object[[]])
head(object[[]])

table(object$condition)
object <- object[, Cells(object)[object$condition %in% c("NORMAL", "TUMOR")]]

# Harmonizing cell_barcode, study, donor, condition, sample_id, batch, disease_classification, celltype
colnames(object[[]])
head(object[[]])
colnames(object@meta.data)[grep("^sample_id$", colnames(object@meta.data))] <- "donor"
colnames(object@meta.data)[grep("^sample_name$", colnames(object@meta.data))] <- "sample_id"
object$disease_classification <- "disease"
table(object$condition)
object$disease_classification[object$condition == "NORMAL"] <- "health"

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