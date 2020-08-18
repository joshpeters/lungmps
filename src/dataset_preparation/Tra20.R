# ---
# Description: Prepare data from Travaglini et al. 2020
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
set <- "Travaglini_2020"
Clean()

# [2] Load files ----

load("data/raw_data/droplet_normal_lung_blood_seurat_ntiss10x.P1.anno.20191002.RC4.Robj")
load("data/raw_data/droplet_normal_lung_blood_seurat_ntiss10x.P3.anno.20191002.RC4.Robj")
load("data/raw_data/droplet_normal_lung_seurat_ntiss10x.P2.anno.20191002.RC4.Robj")

p1 <- UpdateSeuratObject(ntiss10x.P1.anno)
p2 <- UpdateSeuratObject(ntiss10x.P2.anno)
p3 <- UpdateSeuratObject(ntiss10x.P3.anno)

p1 <- DietSeurat(p1)
p2 <- DietSeurat(p2)
p3 <- DietSeurat(p3)

p <- merge(p1, list(p2, p3))
rm(list = c("ntiss10x.P1.anno", "ntiss10x.P2.anno", "ntiss10x.P3.anno", "p1", "p2", "p3"))
object <- p
rm(p)

# [3] Prepare object ----

object <- AddExprMeta(object)

names(object[[]])
head(object[[]])

table(object$tissue)
object <- object[, WhichCells(object, expression = tissue == "lung")]

object@meta.data[, "percent.ribo"] <- NULL
object@meta.data <- janitor::clean_names(object@meta.data)
colnames(object@meta.data)[c(1, 11, 12)] <- c("orig.ident", "nCount_RNA", "nFeature_RNA")

object$batch <- object$channel
object$study <- set

object$cell_barcode <- rownames(object[[]])
object <- AddBatchFreq(object, "batch")
object <- Scrub(object, batch.var = "batch")

# Harmonizing cell_barcode, study, donor, condition, sample_id, batch, disease_classification, celltype
colnames(object@meta.data)[grep("patient", colnames(object@meta.data))] <- "donor"
colnames(object@meta.data)[grep("sample", colnames(object@meta.data))] <- "sample_id"
colnames(object@meta.data)[grep("free_annotation", colnames(object@meta.data))] <- "celltype"
object$donor <- object$batch
object$condition <- "normal"
object$alt_id <- object$orig.ident
object$orig.ident <- object$batch

head(object[[]])
object$disease_classification <- "health"
table(object$batch)
table(object$sample_id)
object$study <- set

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
