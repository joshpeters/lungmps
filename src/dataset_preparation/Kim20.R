# ---
# Description: Prepare data from Kim et al. 2020
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
set <- "Kim_2020"
Clean()

# [2] Load files ----

meta <- read_tsv(gzfile("data/raw_data/GSE131907_Lung_Cancer_cell_annotation.txt.gz"))
counts <- read_tsv(gzfile("data/raw_data/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"))
genes <- counts[, 1, drop = TRUE]
head(genes)

meta_to_use <- meta %>% filter(Sample_Origin %in% c("nLung", "tLung", "tL/B"))
meta_to_use <- meta_to_use %>% janitor::clean_names()
colnames(meta_to_use) <- c("cell_barcode", "barcode", "sample_id", "sample_origin", "cell_type_prior", "cell_type", "cell_subtype")
barcodes_to_use <- meta_to_use %>% pull(cell_barcode)
rownames(meta_to_use) <- meta_to_use$cell_barcode

reduced_counts <- counts[, barcodes_to_use]
Check5(reduced_counts)
sparse_counts <- Matrix::Matrix(as.matrix(reduced_counts), sparse = TRUE)

head(rownames(sparse_counts))
head(colnames(sparse_counts))
rownames(sparse_counts) <- genes

dim(meta_to_use)
dim(sparse_counts)
table(meta_to_use$cell_barcode %in% colnames(sparse_counts))

# counts <- Read10X(data.dir = "/Users/jpeters/data/sets/KimLee_2020/processed/raw_10x/")
# meta <- read_tsv(gzfile("/Users/jpeters/data/sets/KimLee_2020/processed/raw_10x/meta.tsv.gz"))
# dim(meta)
# dim(counts)
# head(meta)
# rownames(meta) <- meta$cell_barcode

# [3] Prepare object ----

object <- CreateSeuratObject(sparse_counts, project = set, meta.data = meta_to_use)
#object <- CreateSeuratObject(counts, project = set, meta.data = meta, min.cells = 0, min.features = 100)
# colnames(object@meta.data)[grep("^percent_mito$", colnames(object@meta.data))] <- "percent_mit"
# colnames(object@meta.data)[grep("^percent_ribo$", colnames(object@meta.data))] <- "percent_rib"
# colnames(object@meta.data)[grep("^percent_hb$", colnames(object@meta.data))] <- "percent_hgb"
colnames(object@meta.data)[grep("^cell_subtype$", colnames(object@meta.data))] <- "celltype"

object <- AddExprMeta(object)
object$batch <- object$sample_id
object$study <- set
object <- AddBatchFreq(object, "batch")
object <- Scrub(object, batch.var = "batch")

# Harmonizing cell_barcode, study, donor, condition, sample_id, batch, disease_classification, celltype
colnames(object[[]])
head(object[[]])

table(object$sample_origin)
colnames(object@meta.data)[grep("^sample_origin$", colnames(object@meta.data))] <- "condition"
object$disease_classification <- "disease"
object$disease_classification[object$condition == "nLung"] <- "health"

sample_ids <- unique(object$sample_id)
sample_ids <- sample_ids[grep(pattern = "LUNG", sample_ids)]
object <- object[, Cells(object)[object$sample_id %in% sample_ids]]
object
object$donor <- gsub("LUNG_", "", object$sample_id)

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

# # [4] Prepare objects ----
#
# object_ls <- list(object)
# names(object_ls) <- set
# set <- SelfName(set)
#
# hgnc <- LoadHGNC()
# genes <- lapply(set, function(x, object_ls, hgnc) {
#
#   ui_todo("Processing {x}")
#   genes <- ConvertGeneNames(genes = rownames(object_ls[[x]]), verbose = TRUE, hgnc = hgnc)
#   assertthat::assert_that(all(rownames(object_ls[[x]]) %in% rownames(genes)))
#   genes <- genes[rownames(object_ls[[x]]), ]
#   return(genes)
#
# }, object_ls = object_ls, hgnc = hgnc)
#
# object_ls <- lapply(set, function(x) {
#
#   assertthat::assert_that(!any(duplicated(genes[[x]]$name)))
#   assertthat::assert_that(all.equal(rownames(object_ls[[x]]), genes[[x]]$orig_name))
#   rownames(object_ls[[x]]@assays$RNA@counts) <- genes[[x]]$name
#   rownames(object_ls[[x]]@assays$RNA@data) <- genes[[x]]$name
#   object_ls[[x]]@assays$RNA@meta.features <- genes[[x]]
#   object_ls[[x]] <- UpdateSeuratObject(object_ls[[x]])
#   return(object_ls[[x]])
#
# })
#
# i = 1
# View(object_ls[[i]]@meta.data)
# object_ls[[i]]$donor <- object_ls[[i]]$sample_id
# object_ls[[i]]$orig.ident <- object_ls[[i]]$batch
# object_ls[[i]]$condition <- object_ls[[i]]$sample_origin
# object_ls[[i]]$disease_classification <- "disease"
# object_ls[[i]]$disease_classification[object_ls[[i]]$sample_origin == "nLung"] <- "health"
# table(object_ls[[i]]$disease_classification)
# table(object_ls[[i]]$sample_origin)
#
# split_object_ls <- lapply(set, function(x) {
#   split <- Seurat::SplitObject(object_ls[[x]], split.by = "disease_classification")
#   return(split)
# })
#
# split_object_ls <- unlist(split_object_ls)
# names(split_object_ls) <- gsub("\\.", "_", names(split_object_ls))
# object_names <- names(split_object_ls)
# names(object_names) <- object_names
#
# walk(object_names, ~ {
#   ui_info("\nSaving {.x} ({round(as.numeric(object.size(split_object_ls[[.x]])/1E6), 2)} mb)")
#   saveRDS(Seurat::DietSeurat(split_object_ls[[.x]]),
#     file = glue::glue("/Users/jpeters/Projects/lungmps/data/objects/{.x}_base.rds"))
#   Clean()
#   ui_done(("\nCompleted {.x}"))
# })
#
# # [5] Process object ----
#
# object_ls <- pblapply(X = object_names, function(x) {
#   object_ls[[x]] <- FilterCells(
#     object_ls[[x]],
#     nmads = 5,
#     variables = c("nCount_RNA", "nFeature_RNA", "percent_mito", "percent_ribo"),
#     batch = "batch",
#     percent_mito_cutoff = 20,
#     percent_ribo_cutoff = 50,
#     percent_hb_cutoff = 5,
#     percent_hsp_cutoff = 5,
#     nCount_RNA_cutoff = 200,
#     nFeature_RNA_cutoff = 100,
#     remove = FALSE,
#     qc_column = "passing_qc")
#   Clean()
#   return(object_ls[[x]])
# })
#
# object_ls <- lapply(X = object_names, function(x) {
#   assertthat::assert_that(all.equal(colnames(object_ls[[x]]), rownames(object_ls[[x]][[]])))
#   keep_cells <- WhichCells(object_ls[[x]], expression = passing_qc == TRUE)
#   ui_info("Number of cells pre-filtering: {ncol(object_ls[[x]])}")
#   object_ls[[x]]<- object_ls[[x]][, keep_cells]
#   ui_info("Number of cells post-filtering: {ncol(object_ls[[x]])}")
#   Clean()
#   return(object_ls[[x]])
# })
#
# object_ls <- pblapply(X = object_names, function(x, object_ls) {
#   object_ls[[x]] <- BaseProcess(object = object_ls[[x]], nfeatures = 4000)
#   object_ls[[x]] <- HarmonyCorrection(object = object_ls[[x]], batch = "batch")
#   Clean()
#   return(object_ls[[x]])
# }, object_ls = object_ls)
#
# object_ls <- pblapply(X = object_names, function(x) {
#   dims <- intrinsicDimension::maxLikGlobalDimEst(object_ls[[x]]@reductions$harmony@cell.embeddings, k = 10)
#   dims <- ceiling(dims$dim.est)
#   ui_done("\nUsing {ui_field(dims)} dims for post-Harmony functions for dataset {x}")
#   object_ls[[x]]$postharmony_dims_used <- dims
#   Clean()
#   return(object_ls[[x]])
# })
# ui_info("Average number of post-Harmony dimensions: {median(unlist(lapply(object_ls, function(x) return(x$postharmony_dims_used))))}")
#
# object_ls <- pblapply(X = object_names, function(x) {
#   object_ls[[x]] <- CalculatePCHeuristics(object_ls[[x]], reduction = "harmony", store.name = "harmony_heuristics")
#   Clean()
#   return(object_ls[[x]])
# })
#
# object_ls <- pblapply(X = object_names, function(x) {
#   dims <- 20
#   method <- "harmony"
#   object_ls[[x]] <- RunUMAP(object = object_ls[[x]], reduction = method,
#     dims = 1:dims, seed.use = 1, reduction.name = glue("{method}_umap"),
#     reduction.key = glue("{substr(method, 1, 1)}UMAP_"))
#   Clean()
#   return(object_ls[[x]])
# })
#
# pblapply(X = object_names, function(x) {
#   ui_info("\nSaving {x} ({round(as.numeric(object.size(object_ls[[x]])/1E6), 2)} mb)")
#   saveRDS(object_ls[[x]],
#     file = glue::glue("/Users/jpeters/Projects/lungmps/data/objects/{x}_annotated.rds"))
#   Clean()
#   ui_done(("\nCompleted {x}"))
# })
#
# DimPlot(object_ls[[1]], reduction = "harmony_umap", group.by = "sample_id")
#
# verbose = TRUE
# path = "data/objects"
# pattern = "KimLee_2020_(disease|health)_annotated.rds"
# object_files <- list.files(path = path, full.names = TRUE)
# object_files_to_load <- object_files[grep(pattern = pattern, x = object_files)]
# object_names <- stringr::str_match(object_files_to_load, pattern = paste0("objects\\/(.*)","_annotated.rds"))[, 2]
# names(object_names) <- object_names
# names(object_files_to_load) <- object_names
# object_ls <- future.apply::future_lapply(object_names, function(x) {
#   if (verbose) usethis::ui_info("\nLoading {x} ({round(file.info(object_files_to_load[x])$size/1e6, 2)} mb)")
#   obj <- readRDS(file = object_files_to_load[x])
#   if (verbose) usethis::ui_done(("\nCompleted {x}"))
#   Clean()
#   return(obj)
# })
# names(object_ls) <- object_names
#
#
# object_ls <- pblapply(X = object_names, function(x) {
#   dims <- 20
#   object_ls[[x]] <- FindNeighbors(object_ls[[x]], reduction = "harmony", dims = 1:dims,
#     k.param = 20, prune.SNN = 1/15, compute.SNN = TRUE, graph.name = glue("harmony_snn"), verbose = TRUE)
#   object_ls[[x]] <- BaseCluster(object = object_ls[[x]], dims = 20, verbose = TRUE,
#     res.start = 1E-5, res.end = 1, num.res = 30, num.clusters.lower = 5, num.clusters.upper = 30, mod.similarity = 0.995)
#   snn <- igraph::graph_from_adjacency_matrix(object_ls[[x]]@graphs$harmony_snn, mode = "directed", weighted = TRUE)
#   walktrap_cluster <- igraph::cluster_walktrap(graph = snn, steps = 4)
#   object_ls[[x]]$base_walktrap <- walktrap_cluster$membership
#   Clean()
#   return(object_ls[[x]])
# })
#
# minimum.cells <- 10
# cluster <- "base_leiden"
# object_ls <- pblapply(X = object_names, function(x) {
#   ui_info("Processing {x}...")
#   print(table(object_ls[[x]]$base_leiden))
#   object_ls[[x]] <- SetIdent(object_ls[[x]], value = cluster)
#   object_ls[[x]] <- BuildClusterTree(object_ls[[x]], reorder = TRUE, reorder.numeric = TRUE)
#   object_ls[[x]]$base_clustering <- object_ls[[x]]@active.ident
#   return(object_ls[[x]])
# })
#
# sheets <- readxl::excel_sheets(path = "/Users/jpeters/Projects/lungmps/data/raw/travaglini_2020_markers.xlsx")
# sheets <- sheets[-grep("SS2|SS", sheets)]
# degs <- pblapply(sheets, function(x) {
#   degs <- readxl::read_xlsx(path = "/Users/jpeters/Projects/lungmps/data/raw/travaglini_2020_markers.xlsx",
#     skip = 1, sheet = x)
#   degs <- degs %>% filter(avg_logFC > log(2) & p_val_adj < 1E-10) %>% top_n(20, avg_logFC) %>% pull(Gene)
#   return(degs)
# })
# deg_names <- pblapply(sheets, function(x) {
#   deg_names <- readxl::read_xlsx(path = "/Users/jpeters/Projects/lungmps/data/raw/travaglini_2020_markers.xlsx",
#     n_max = 1, sheet = x, col_names = FALSE)
#   deg_names <- as.character(deg_names[1, 1])
#   return(deg_names)
# })
# deg_names <- snakecase::to_snake_case(unlist(deg_names))
# deg_names <- paste0(deg_names, "_score")
# names(degs) <- deg_names
# type_df <- read.csv("data/processed/tra20_type.csv")
#
# object_ls <- pblapply(X = object_names, function(x) {
#   object_ls[[x]] <- ScoreLineages(object = object_ls[[x]], grouping_var = "base_clustering", genes = degs)
#   Clean()
#   return(object_ls[[x]])
# })
#
# deg_names_df <- data.frame(index = seq(1:length(deg_names)), name = deg_names)
# deg_names_df <- merge(deg_names_df, type_df, by.x = "name", by.y = "score_name", sort = FALSE)
# score_names_touse <- as.character(deg_names_df$name[deg_names_df$usage == 1])
# deg_names_touse <- deg_names_df[deg_names_df$usage == 1, ]
# deg_names_touse$index <- seq(1:nrow(deg_names_touse))
#
# score_labels <- gsub("_", " ", score_names_touse)
# score_labels <- snakecase::to_title_case(score_labels)
# score_labels <- gsub("Score", " ", score_labels)
# score_labels <- str_trim(score_labels)
#
# ncol(object_ls[[1]])
# length(unique(object_ls[[1]]$cell_barcode))
#
# object_ls <- pblapply(X = object_names, function(x) {
#
#   dotplot <- DotPlot_ol(object_ls[[x]], features = score_names_touse,
#     cols = rev(colorspace::sequential_hcl(2, "Blues")),
#     dot.min = 0.05, dot.scale = 4, scale.by = "radius") +
#     scale_x_discrete(labels = rev(score_labels)) +
#     labs(x = "Expression Score", y = "Cluster") +
#     theme(
#       plot.subtitle = element_text(face = "bold"),
#       legend.position = "none",
#       axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 6),
#       axis.text.y = element_text(size = 10),
#       axis.title.x = element_text(face = "bold", margin = margin(t = 2, unit = "pt")),
#       axis.title.y = element_text(face = "bold", margin = margin(r = 12, unit = "pt")),
#       legend.title = element_text(vjust = 1, hjust = 0),
#       panel.grid.major = element_line(color = "Gray90", size = 0.5),
#       plot.margin = unit(rep(18, 4), "pt"),
#       axis.line = element_blank(),
#       panel.border = element_rect(size = 0.5, color = "black")
#     )
#   cowplot::save_plot(plot = dotplot, filename = glue("plots/scores_{x}.png"), base_height = 6, base_asp = 1.5)
#
#   max_scores <- apply(object_ls[[x]][[]][, score_names_touse, drop = TRUE], 1, which.max)
#   max_scores <- ConvertNamedVecToDF(max_scores)
#   max_scores <- merge(max_scores, deg_names_touse, by.x = "value", by.y = "index", sort = FALSE)
#   max_scores <- max_scores %>% select(name.x, name.y, jmp_class, jmp_subclass, jmp_type)
#   colnames(max_scores) <- c("cell_barcode", "score_name", "jmp_class", "jmp_subclass", "jmp_type")
#   rownames(max_scores) <- max_scores$cell_barcode
#   max_scores <- max_scores[, -1]
#   object_ls[[x]] <- AddMetaData(object_ls[[x]], max_scores)
#
#   cluster_barcodes <- object_ls[[x]][[]] %>% select(cell_barcode, base_clustering)
#   cluster_assignments <- object_ls[[x]][[]] %>% select(base_clustering, jmp_class)
#   cluster_assignments <- cluster_assignments %>%
#     group_by(base_clustering, jmp_class) %>% summarize(n = n())
#
#   cols <- ggthemes::colorblind_pal()(8)[c(2, 3, 4, 7)]
#   names(cols) <- levels(cluster_assignments$jmp_class)
#   count_plot <- ggplot(cluster_assignments, aes(x = base_clustering, y = n, fill = jmp_class)) +
#     geom_col(position = "stack", width = 1, color = "black") +
#     scale_y_continuous(expand = c(0, 0), limits = c(0, max(cluster_assignments$n) + 200)) +
#     scale_fill_manual(values = cols, name = "Classes", labels = c("Endothelial", "Epithelial", "Immune", "Stromal")) +
#     labs(x = "Cluster", y = "Number of cells") +
#     theme_classic(base_size = 14) +
#     ggeasy::easy_all_text_color("black") +
#     theme(
#       axis.title.y = element_text(angle = 90, vjust = 0.5, margin = margin(0, 12, 0, 0), face = "bold"),
#       axis.title.x = element_text(margin = margin(12, 0, 0, 0), face = "bold"),
#       axis.text.y = element_text(margin = margin(0, 4, 0, 0)),
#       axis.line = element_blank(),
#       plot.margin = margin(4, 4, 4, 4, unit = "mm"),
#       panel.border = element_rect(size = 1, color = "black", fill = "transparent"),
#       legend.title = element_text(size = 12, face = "bold"),
#       legend.text = element_text(size = 10),
#       legend.position = "bottom",
#       legend.justification = "left",
#       legend.key.size = unit(1, "line")
#     ) + coord_flip()
#   cowplot::save_plot(plot = count_plot, filename = glue("plots/cluster_class_counts_{x}.png"), base_height = 6, base_asp = 1)
#   cluster_assignments <- cluster_assignments %>%
#     top_n(n = 1, wt = n)
#   cluster_assignments <- merge(cluster_barcodes, cluster_assignments, .by = "base_clustering", sort = FALSE)
#   rownames(cluster_assignments ) <- cluster_assignments$cell_barcode
#   cluster_assignments  <- cluster_assignments [, -1]
#   object_ls[[x]] <- AddMetaData(object_ls[[x]], cluster_assignments[, "jmp_class", drop = FALSE])
#
#   subcluster_barcodes <- object_ls[[x]][[]] %>% select(cell_barcode, base_clustering)
#   subcluster_assignments <- object_ls[[x]][[]] %>% select(base_clustering, jmp_subclass)
#   subcluster_assignments <- subcluster_assignments %>%
#     group_by(base_clustering, jmp_subclass) %>% summarize(n = n()) %>%
#     top_n(n = 1, wt = n)
#   subcluster_assignments <- merge(cluster_barcodes, subcluster_assignments, .by = "base_clustering", sort = FALSE)
#   length(unique(subcluster_assignments$cell_barcode)) == length(subcluster_assignments$cell_barcode)
#   subcluster_assignments <- subcluster_assignments[!duplicated(subcluster_assignments$cell_barcode), ]
#   rownames(subcluster_assignments) <- subcluster_assignments$cell_barcode
#   subcluster_assignments  <- subcluster_assignments [, -1]
#   object_ls[[x]] <- AddMetaData(object_ls[[x]], subcluster_assignments[, "jmp_subclass", drop = FALSE])
#   table(object_ls[[x]]$jmp_subclass)
#
#   Clean()
#   return(object_ls[[x]])
# })
#
# pblapply(X = object_names, function(x) {
#   if ("celltype" %in% colnames(object_ls[[x]])) {
#     plots <- DimPlot(object_ls[[x]], reduction = "harmony_umap",
#       label = TRUE, group.by = c("base_clustering", "celltype", "jmp_class", "jmp_subclass"), combine = FALSE, label.size = 4, repel = TRUE)
#     plots <- lapply(plots, function(x) {
#       x + colorspace::scale_color_discrete_qualitative("Dark 2") +
#         labs(x = "UMAP 1", y = "UMAP 2") +
#         theme_void(base_size = 14) +
#         theme(
#           axis.title.x = element_text(hjust = 0, face = "plain", color = "gray60"),
#           axis.title.y = element_text(hjust = 0, angle = 90, face = "plain", color = "gray60"),
#           panel.background = element_rect(fill = "transparent", color = "transparent"), # bg of the panel
#           plot.background = element_rect(fill = "transparent", color = "transparent"), # bg of the plot
#           panel.grid.major = element_blank(), # get rid of major grid
#           panel.grid.minor = element_blank(), # get rid of minor grid
#           legend.background = element_rect(fill = "transparent"), # get rid of legend bg
#           legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
#           legend.position = "none"
#         )
#     })
#     plotgrid <- cowplot::plot_grid(plotlist = plots, align = c("v"), ncol = 4)
#     cowplot::save_plot(filename = glue("plots/{x}_cellclass.png"), plot = plotgrid,
#       base_height = 5, base_asp = 4)
#   } else {
#     plots <- DimPlot(object_ls[[x]], reduction = "harmony_umap",
#       label = TRUE, group.by = c("base_clustering", "jmp_class", "jmp_subclass"), combine = FALSE, label.size = 4, repel = TRUE)
#     plots <- lapply(plots, function(x) {
#       x + colorspace::scale_color_discrete_qualitative("Dark 2") +
#         labs(x = "UMAP 1", y = "UMAP 2") +
#         theme_void(base_size = 14) +
#         theme(
#           axis.title.x = element_text(hjust = 0, face = "plain", color = "gray60"),
#           axis.title.y = element_text(hjust = 0, angle = 90, face = "plain", color = "gray60"),
#           panel.background = element_rect(fill = "transparent", color = "transparent"), # bg of the panel
#           plot.background = element_rect(fill = "transparent", color = "transparent"), # bg of the plot
#           panel.grid.major = element_blank(), # get rid of major grid
#           panel.grid.minor = element_blank(), # get rid of minor grid
#           legend.background = element_rect(fill = "transparent"), # get rid of legend bg
#           legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
#           legend.position = "none"
#         )
#     })
#     plotgrid <- cowplot::plot_grid(plotlist = plots, align = c("v"), ncol = 3)
#     cowplot::save_plot(filename = glue("plots/{x}_cellclass.png"), plot = plotgrid,
#       base_height = 5, base_asp = 3)
#   }
#   Clean()
# })
#
# object_ls <- pblapply(X = object_names, function(x) {
#   object_ls[[x]] <- LabelCyclingCells(object = object_ls[[x]], grouping_var = "base_clustering")
#   Clean()
#   return(object_ls[[x]])
# })
# lapply(object_ls, function(x) table(x$cc_diff_label))
#
# lm22_expr <- readxl::read_xlsx(path = "/Users/jpeters/Projects/lungmps/data/raw/lm22_expr.xlsx")
# is_expr <- readxl::read_xlsx(path = "/Users/jpeters/Projects/lungmps/data/raw/immunoStates_expr.xlsx")
#
# is_expr <- janitor::clean_names(is_expr)
# colnames(is_expr)
# ref_genes <- is_expr$gene
# is_expr <- is_expr[, -c(1)]
# rownames(is_expr) <- ref_genes
#
# lm22_expr <- janitor::clean_names(lm22_expr)
# colnames(lm22_expr)
# ref_genes <- lm22_expr$gene_symbol
# lm22_expr <- lm22_expr[, -c(1:3)]
# rownames(lm22_expr) <- ref_genes
#
# object_ls <- pblapply(X = object_names, function(x) {
#   lls <- LLClassifier(object_ls[[x]], is_expr)
#   object_ls[[x]]@misc$is_llpp <- lls$pp
#   calls <- ConvertNamedVecToDF(lls$call)
#   calls <- calls %>% select(value)
#   object_ls[[x]] <- AddMetaData(object_ls[[x]], calls, col.name = "IS")
#
#   lls <- LLClassifier(object_ls[[x]], lm22_expr)
#   object_ls[[x]]@misc$is_llpp <- lls$pp
#   calls <- ConvertNamedVecToDF(lls$call)
#   calls <- calls %>% select(value)
#   object_ls[[x]] <- AddMetaData(object_ls[[x]], calls, col.name = "LM22")
#
#   Clean()
#   return(object_ls[[x]])
# })
#
# object_ls <- pblapply(X = object_names, function(x) {
#   object_ls[[x]] <- SetIdent(object = object_ls[[x]], value = "base_clustering")
#   global_markers <- presto::wilcoxauc(object_ls[[x]], assay = "data", seurat_assay = "RNA")
#   object_ls[[x]]@misc$global_markers <- global_markers
#   Clean()
#   return(object_ls[[x]])
# })
#
# pblapply(X = object_names, function(x) {
#   ui_info("\nSaving {x} ({round(as.numeric(object.size(object_ls[[x]])/1E6), 2)} mb)")
#   saveRDS(object_ls[[x]],
#     file = glue::glue("/Users/jpeters/Projects/lungmps/data/objects/{x}_annotated.rds"))
#   Clean()
#   ui_done(("\nCompleted {x}"))
# })
#
# object_ls <- pblapply(X = object_names, function(x) {
#
#   ui_info("Classifying {x}")
#   nonmp <- as.character(type_df$score_name[type_df$mnp == 0])
#   nonmp <- nonmp[nonmp %in% colnames(object_ls[[x]]@meta.data)]
#   mp <- as.character(type_df$score_name[type_df$mnp == 1])
#   assertthat::assert_that(all(nonmp %in% colnames(object_ls[[x]]@meta.data)))
#   mnp_diff_score <- rowSums(object_ls[[x]]@meta.data[, mp]) - rowSums(object_ls[[x]]@meta.data[, nonmp])
#   assertthat::assert_that(all.equal(names(mnp_diff_score), rownames(object_ls[[x]]@meta.data)))
#   object_ls[[x]]$mp_diff_score <- mnp_diff_score
#
#   mp_scores <- object_ls[[x]][[]] %>% select(base_clustering, mp_diff_score)
#   mm <- mixtools::normalmixEM(mp_scores$mp_diff_score, k = 2)
#   plot(mm, 2)
#   x_data <- with(mm, seq(min(x), max(x), len = 1000))
#   pars <- with(mm, data.frame(comp = colnames(posterior), mu, sigma, lambda))
#   em_df <- data.frame(x = rep(x_data, each = nrow(pars)), pars)
#   em_df$y <- with(em_df, lambda * dnorm(x, mean = mu, sd = sigma))
#   min <- which.min(mm$mu)
#   threshold <- mm$mu[min] + 2.5*mm$sigma[min]
#   if (threshold >= 5) threshold <- 0
#   ui_info("Threshold: {threshold}")
#
#   mp_clusters <- as.character(mp_scores %>%
#       group_by(base_clustering) %>%
#       summarize(median_mp_score = median(mp_diff_score)) %>%
#       filter(median_mp_score > threshold & median_mp_score > 0) %>%
#       pull(base_clustering))
#   nonmp_clusters <- setdiff(unique(object_ls[[x]]$base_clustering), mp_clusters)
#
#   cc_props <- as.matrix(table(object_ls[[x]]$base_clustering, object_ls[[x]]$cc_diff_label))
#   cc_props <- data.frame(cluster = rownames(cc_props), cycling = cc_props[, 1, drop = TRUE], noncycling = cc_props[, 2, drop = TRUE])
#   cc_props <- cc_props %>% group_by(cluster) %>% mutate(props = cycling/(cycling+noncycling))
#   cycling_cluster <- as.character(cc_props %>% ungroup %>% top_n(1, props) %>% pull(cluster))
#
#   nonmp_clusters <- setdiff(nonmp_clusters, cycling_cluster)
#   mp_clusters <- setdiff(mp_clusters, cycling_cluster)
#
#   counts <- object_ls[[x]]@meta.data %>%
#     select(base_clustering, jmp_type) %>%
#     filter(base_clustering %in% mp_clusters) %>%
#     group_by(base_clustering, jmp_type) %>%
#     summarize(n = n()) %>%
#     mutate(prop = n/sum(n)) %>%
#     filter(prop > 0.1 & prop < 0.50)
#   contam_clusters <- unique(as.character(counts$base_clustering[!counts$jmp_type %in% c("monocyte", "dc", "macrophage")]))
#   counts <- object_ls[[x]]@meta.data %>%
#     select(base_clustering, jmp_type) %>%
#     filter(base_clustering %in% mp_clusters) %>%
#     group_by(base_clustering, jmp_type) %>%
#     summarize(n = n()) %>%
#     mutate(prop = n/sum(n)) %>%
#     filter(prop >= 0.50)
#   highcontam_clusters <- unique(as.character(counts$base_clustering[!counts$jmp_type %in% c("monocyte", "dc", "macrophage")]))
#   mp_clusters <- setdiff(mp_clusters, contam_clusters)
#   mp_clusters <- setdiff(mp_clusters, highcontam_clusters)
#
#   object_ls[[x]]$base_labels <- NA
#   object_ls[[x]]$base_labels[object_ls[[x]]$base_clustering %in% mp_clusters] <- "mnp"
#   object_ls[[x]]$base_labels[object_ls[[x]]$base_clustering %in% contam_clusters] <- "prelim_mnp"
#   object_ls[[x]]$base_labels[object_ls[[x]]$base_clustering %in% c(nonmp_clusters, highcontam_clusters)] <- "non_mnp"
#   object_ls[[x]]$base_labels[object_ls[[x]]$base_clustering %in% cycling_cluster] <- "cycling"
#
#   cc_scores <- object_ls[[x]][[]] %>% select(cell_barcode, base_clustering, mp_diff_score) %>% filter(base_clustering == cycling_cluster)
#   mp_cycling <- cc_scores %>%
#     group_by(base_clustering) %>%
#     filter(mp_diff_score > threshold) %>% pull(cell_barcode)
#
#   object_ls[[x]]$mp <- FALSE
#   object_ls[[x]]$mp[object_ls[[x]]$base_clustering %in% c(mp_clusters, contam_clusters)] <- TRUE
#   object_ls[[x]]$mp[mp_cycling] <- TRUE
#   ui_info("Found {sum(object_ls[[x]]$mp)} MNP cells")
#
#   a <- DimPlotPlus(object_ls[[x]], reduction = "harmony_umap", group.by = "base_clustering",
#     label = TRUE, label.size = 6, repel = TRUE) +
#     colorspace::scale_color_discrete_qualitative("Dark 3", name = "Cluster") +
#     labs(x = "UMAP 1", y = "UMAP 2") +
#     theme_void(base_size = 14) +
#     theme(
#       axis.title.x = element_text(hjust = 0, face = "plain", color = "gray60"),
#       axis.title.y = element_text(hjust = 0, angle = 90, face = "plain", color = "gray60"),
#       panel.background = element_rect(fill = "transparent", color = "transparent"), # bg of the panel
#       plot.background = element_rect(fill = "transparent", color = "transparent"), # bg of the plot
#       panel.grid.major = element_blank(), # get rid of major grid
#       panel.grid.minor = element_blank(), # get rid of minor grid
#       legend.background = element_rect(fill = "transparent"), # get rid of legend bg
#       legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
#       legend.position = "none",
#       plot.margin = margin(6, 6, 6, 6, unit = "mm")
#     )
#
#   b <- ggplot(data.frame(x = mm$x), aes(x, y = ..density..)) +
#     geom_histogram(fill = NA, color = "black") +
#     geom_polygon(data = em_df, aes(x, y, fill = comp), color = "black", alpha = 0.5) +
#     geom_vline(xintercept = threshold, color = "gray40", linetype = "dashed") +
#     colorspace::scale_fill_discrete_sequential("Light Grays", name = "Component Means") +
#       #labels = format(x = em_df$mu, digits = 2)) +
#     scale_x_continuous(expand = c(0, 0.05)) +
#     scale_y_continuous(expand = c(0, 0)) +
#     labs(x = "MP Score", y = "Density") +
#     theme_classic(base_size = 14) +
#     ggeasy::easy_all_text_color(color = "black") +
#     theme(
#       axis.title.y = element_text(angle = 90, vjust = 0.5, margin = margin(0, 12, 0, 0), face = "bold"),
#       axis.title.x = element_text(margin = margin(12, 0, 0, 0), face = "bold"),
#       axis.text.y = element_text(margin = margin(0, 4, 0, 0)),
#       axis.line = element_blank(),
#       panel.border = element_rect(size = 1, color = "black", fill = "transparent"),
#       legend.title = element_text(size = 12, face = "bold"),
#       legend.text = element_text(size = 10),
#       legend.justification = "left",
#       legend.position = "bottom",
#       legend.key.size = unit(1, "line"),
#       plot.margin = margin(6, 6, 6, 6, unit = "mm"))
#
#   c <- ggplot(mp_scores, aes(x = mp_diff_score, y = base_clustering, fill = base_clustering)) +
#     ggridges::geom_density_ridges(alpha = 0.8) +
#     geom_vline(xintercept = threshold, color = "gray40", linetype = "dashed") +
#     colorspace::scale_fill_discrete_qualitative("Dark 2", name = "Cluster") +
#     scale_x_continuous(expand = c(0, 0.05)) +
#     scale_y_discrete(expand = c(0, 1)) +
#     labs(x = "MP Score", y = "Cluster") +
#     theme_classic(base_size = 14) +
#     ggeasy::easy_all_text_color(color = "black") +
#     theme(
#       axis.title.y = element_text(angle = 90, vjust = 0.5, margin = margin(0, 12, 0, 0), face = "bold"),
#       axis.title.x = element_text(margin = margin(12, 0, 0, 0), face = "bold"),
#       axis.text.y = element_text(margin = margin(0, 4, 0, 0)),
#       axis.line = element_blank(),
#       panel.border = element_rect(size = 1, color = "black", fill = "transparent"),
#       legend.position = "none",
#       plot.margin = margin(6, 6, 6, 6, unit = "mm"))
#
#   d <- DimPlotPlus(object_ls[[x]], reduction = "harmony_umap", group.by = "mp",
#     label = TRUE, label.size = 6, repel = TRUE) +
#     colorspace::scale_color_discrete_qualitative("Dark 3", name = "Filtered as MP") +
#     labs(x = "UMAP 1", y = "UMAP 2",
#       caption = glue::glue("MP assignment ({sum(object_ls[[x]]$mp)} MNP cells assigned, {round(sum(object_ls[[x]]$mp)/ncol(object_ls[[x]])*100, 2)}% of total)")) +
#     theme_void(base_size = 14) +
#     theme(
#       axis.title.x = element_text(hjust = 0, face = "plain", color = "gray60"),
#       axis.title.y = element_text(hjust = 0, angle = 90, face = "plain", color = "gray60"),
#       plot.caption = element_text(size = 10, hjust = 1),
#       legend.position = "none",
#       plot.margin = margin(6, 6, 6, 6, unit = "mm"),
#     )
#
#   full <- cowplot::plot_grid(a, b, c, d, ncol = 2, labels = "auto")
#   cowplot::save_plot(filename = glue::glue("/Users/jpeters/Projects/lungmps/plots/{x}_mnp_classification.png"),
#     plot = full, base_height = 10, base_asp = 1.2, device = "png")
#
#   return(object_ls[[x]])
# })
#
# pblapply(X = object_names, function(x) {
#   ui_info("\nSaving {x} ({round(as.numeric(object.size(object_ls[[x]])/1E6), 2)} mb)")
#   saveRDS(object_ls[[x]],
#     file = glue::glue("/Users/jpeters/Projects/lungmps/data/objects/{x}_annotated.rds"))
#   Clean()
#   ui_done(("\nCompleted {x}"))
# })
#
# # [6] Filter MPs ----
#
# object_ls <- LoadObjectList("data/objects", glue("{set}.*_annotated.rds"), "_annotated.rds", TRUE)
# object_ls <- pblapply(X = object_names, function(x) {
#   ui_done("Subsetting {sum(object_ls[[x]]$mp)} MNPs from {x}")
#   mnp_object <- object_ls[[x]][, WhichCells(object_ls[[x]], expression = mp == TRUE)]
#   mnp_object<- DietSeurat(object = mnp_object)
#   Clean()
#   return(mnp_object)
# })
#
# ui_info("Total MNPs: {sum(unlist(lapply(object_ls, ncol)))}")
#
# object_ls <- pblapply(X = object_names, function(x) {
#   object_ls[[x]] <- BaseProcess(object = object_ls[[x]], nfeatures = 4000)
#   object_ls[[x]] <- HarmonyCorrection(object = object_ls[[x]], batch = "batch")
#   Clean()
#   return(object_ls[[x]])
# })
#
# object_ls <- pblapply(X = object_names, function(x) {
#   object_ls[[x]] <- CalculatePCHeuristics(object_ls[[x]], reduction = "harmony", store.name = "harmony_heuristics")
#   dims <- ceiling(unique(object_ls[[x]]@misc$harmony_heuristics$tp))
#   ui_done("\nUsing {ui_field(dims)} dims for post-Harmony functions for dataset {x}")
#   object_ls[[x]]$postharmony_dims_used <- dims
#   Clean()
#   return(object_ls[[x]])
# })
# ui_info("Average number of post-Harmony dimensions: {median(unlist(lapply(object_ls, function(x) return(unique(x$postharmony_dims_used)))))}")
#
# object_ls <- pblapply(X = object_names, function(x) {
#   dims <- 20
#   #dims <- unique(object_ls[[x]]$postharmony_dims_used)
#   method <- "harmony"
#   object_ls[[x]] <- RunUMAP(object = object_ls[[x]], reduction = method,
#     dims = 1:dims, seed.use = 1, reduction.name = glue("{method}_umap"),
#     reduction.key = glue("{substr(method, 1, 1)}UMAP_"))
#   Clean()
#   return(object_ls[[x]])
# })
#
# object_ls <- pblapply(X = object_names, function(x) {
#   dims <- 20
#   object_ls[[x]] <- FindNeighbors(object_ls[[x]], reduction = "harmony", dims = 1:dims,
#     k.param = 20, prune.SNN = 1/15, compute.SNN = TRUE, graph.name = glue("harmony_snn"), verbose = TRUE)
#   object_ls[[x]] <- BaseCluster(object = object_ls[[x]], dims = 20, verbose = TRUE,
#     res.start = 1E-5, res.end = 1, num.res = 30, num.clusters.lower = 5, num.clusters.upper = 30, mod.similarity = 0.995)
#   snn <- igraph::graph_from_adjacency_matrix(object_ls[[x]]@graphs$harmony_snn, mode = "directed", weighted = TRUE)
#   walktrap_cluster <- igraph::cluster_walktrap(graph = snn, steps = 4)
#   object_ls[[x]]$base_walktrap <- walktrap_cluster$membership
#   Clean()
#   return(object_ls[[x]])
# })
#
# minimum.cells <- 10
# cluster <- "base_leiden"
# object_ls <- pblapply(X = object_names, function(x) {
#   ui_info("Processing {x}...")
#   print(table(object_ls[[x]]$base_leiden))
#   object_ls[[x]] <- SetIdent(object_ls[[x]], value = cluster)
#   object_ls[[x]] <- BuildClusterTree(object_ls[[x]], reorder = TRUE, reorder.numeric = TRUE)
#   object_ls[[x]]$base_clustering <- object_ls[[x]]@active.ident
#   return(object_ls[[x]])
# })
#
# # Find markers for walktrap clusters
# object_ls <- pblapply(X = object_names, function(x) {
#
#   object_ls[[x]] <- SetIdent(object = object_ls[[x]], value = "base_walktrap")
#   walktrap_markers <- presto::wilcoxauc(object_ls[[x]], assay = "data", seurat_assay = "RNA")
#   object_ls[[x]]@misc$walktrap_markers <- walktrap_markers
#
#   object_ls[[x]] <- SetIdent(object = object_ls[[x]], value = "base_leiden")
#   leiden_markers <- presto::wilcoxauc(object_ls[[x]], assay = "data", seurat_assay = "RNA")
#   object_ls[[x]]@misc$leiden_markers <- leiden_markers
#
#   Clean()
#   return(object_ls[[x]])
# })
#
# all_genes <- unlist(lapply(object_ls, rownames))
# ig <- all_genes[grep("^IG[H/L/V/K].*$", all_genes)]
# tcrb <- all_genes[grep("^TRBV.*$", all_genes)]
# tcrg <- all_genes[grep("^TRGV.*$", all_genes)]
# tcrc <- all_genes[grep("^TRGC.*$", all_genes)]
# tcrj <- all_genes[grep("^TRGJ.*$", all_genes)]
# tcra <- all_genes[grep("^TRA[C|V].*$", all_genes)]
# tcrs <- c(tcrb, tcrg, tcrc, tcrj, tcra)
#
# blacklist_features <- readRDS("/Users/jpeters/Projects/lungmps/data/processed/blacklist_features.rds")
# nonmp_markers <- readRDS("/Users/jpeters/Projects/lungmps/data/processed/nonMP_markers.rds")
# nonmp_markers <- c(nonmp_markers$feature, "NAMPT", "CSF3R", "CD79A", "CD79B", "MS4A1", "IGKC", ig, tcrs)
# nonmp_markers <- nonmp_markers[!(nonmp_markers %in% c("A2M", "GZMB", "RGS1"))]
#
# total_mp_scores <- pblapply(object_names, function(x) {
#   return(object_ls[[x]][[]] %>% select(base_leiden, mp_diff_score, cell_barcode))
# })
# total_mp_scores <- Reduce(rbind, total_mp_scores)
# total_threshold <- median(total_mp_scores$mp_diff_score) - 3*mad(total_mp_scores$mp_diff_score)
# total_threshold <- 0
# ui_done("Total threshold for MP score: {total_threshold}")
# mean(total_mp_scores$mp_diff_score) - 2*sd(total_mp_scores$mp_diff_score)
# hist(total_mp_scores$mp_diff_score)
#
# blacklist_hits <- pblapply(X = object_names, function(x) {
#
#   ui_info("Printing figure for {x}")
#   object_ls[[x]]$base_leiden <- as.factor(object_ls[[x]]$base_leiden)
#   object_ls[[x]]$base_leiden <- factor(object_ls[[x]]$base_leiden, levels = levels(object_ls[[x]]$base_leiden)[order(as.numeric(levels(object_ls[[x]]$base_leiden)))] )
#
#   mp_scores <- object_ls[[x]][[]] %>% select(base_leiden, mp_diff_score)
#
#   plot_markers <- object_ls[[x]]@misc$leiden_markers %>%
#     dplyr::filter(padj < 1E-10 & logFC > log(1.25) & auc > 0.5) %>%
#     group_by(group) %>%
#     arrange(as.numeric(group)) %>%
#     top_n(10, logFC)
#   blacklist_hits <- plot_markers[(plot_markers$feature %in% c(nonmp_markers)), ]
#   blacklist_hits <- names(table(blacklist_hits$group))[table(blacklist_hits$group) >= 2]
#
#   nonmp_clusters <- as.character(mp_scores %>%
#       group_by(base_leiden) %>%
#       summarize(median_mp_score = median(mp_diff_score)) %>%
#       filter(median_mp_score < 0) %>%
#       pull(base_leiden))
#   nonmp_cells <- colnames(object_ls[[x]])[as.character(object_ls[[x]]$base_leiden) %in% unique(c(blacklist_hits, nonmp_clusters))]
#   nonmp_cells <- unique(c(nonmp_cells, rownames(mp_scores)[mp_scores$mp_diff_score < total_threshold]))
#   ui_info("Identified {length(nonmp_cells)} contaminating cell profiles")
#
#   genes <- unique(plot_markers %>%
#       group_by(group) %>%
#       top_n(5, logFC) %>% pull(feature))
#
#   cluster_umap <- DimPlotPlus(object_ls[[x]], reduction = "harmony_umap", group.by = "base_leiden",
#     label = TRUE, label.size = 6, repel = TRUE) +
#     colorspace::scale_color_discrete_qualitative("Dark 3", name = "Cluster") +
#     labs(x = "UMAP 1", y = "UMAP 2") +
#     theme_void(base_size = 14) +
#     theme(
#       axis.title.x = element_text(hjust = 0, face = "plain", color = "gray60"),
#       axis.title.y = element_text(hjust = 0, angle = 90, face = "plain", color = "gray60"),
#       panel.background = element_rect(fill = "transparent", color = "transparent"), # bg of the panel
#       plot.background = element_rect(fill = "transparent", color = "transparent"), # bg of the plot
#       panel.grid.major = element_blank(), # get rid of major grid
#       panel.grid.minor = element_blank(), # get rid of minor grid
#       legend.background = element_rect(fill = "transparent"), # get rid of legend bg
#       legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
#       legend.position = "none",
#       plot.margin = margin(6, 6, 6, 6, unit = "mm")
#     )
#
#   contam_umap <- DimPlotPlus(object_ls[[x]], reduction = "harmony_umap", group.by = "base_leiden", pt.size = 1,
#     cells.highlight = nonmp_cells) +
#     scale_color_manual(values = c("gray90", "black"),
#       name = "Potential\ncontamination", labels = c("MP", paste(blacklist_hits, collapse = ", "))) +
#     labs(x = "UMAP 1", y = "UMAP 2") +
#     theme_void(base_size = 14) +
#     theme(
#       axis.title.x = element_text(hjust = 0, face = "plain", color = "gray60"),
#       axis.title.y = element_text(hjust = 0, angle = 90, face = "plain", color = "gray60"),
#       panel.background = element_rect(fill = "transparent", color = "transparent"), # bg of the panel
#       plot.background = element_rect(fill = "transparent", color = "transparent"), # bg of the plot
#       panel.grid.major = element_blank(), # get rid of major grid
#       panel.grid.minor = element_blank(), # get rid of minor grid
#       legend.background = element_rect(fill = "transparent"), # get rid of legend bg
#       legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
#       legend.position = "none",
#       plot.margin = margin(6, 6, 6, 6, unit = "mm")
#     )
#
#   scores_plot <- DotPlot_ol(object_ls[[x]], features = score_names_touse,
#     cols = rev(colorspace::sequential_hcl(2, "Blues")),
#     dot.min = 0.05, dot.scale = 4, scale.by = "radius", group.by = "base_leiden") +
#     scale_x_discrete(labels = rev(score_labels)) +
#     labs(x = "Expression Score", y = "Cluster") +
#     theme(
#       plot.subtitle = element_text(face = "bold"),
#       legend.position = "none",
#       axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 6),
#       axis.text.y = element_text(size = 10),
#       axis.title.x = element_text(face = "bold", margin = margin(t = 2, unit = "pt")),
#       axis.title.y = element_text(face = "bold", margin = margin(r = 12, unit = "pt")),
#       legend.title = element_text(vjust = 1, hjust = 0),
#       panel.grid.major = element_line(color = "Gray90", size = 0.5),
#       plot.margin = unit(rep(18, 4), "pt"),
#       axis.line = element_blank(),
#       panel.border = element_rect(size = 0.5, color = "black")
#     )
#
#   markers_plot <- DotPlot_ol(object_ls[[x]], features = genes,
#     cols = rev(colorspace::sequential_hcl(2, "Blues")),
#     dot.min = 0.05, dot.scale = 6, scale.by = "radius", group.by = "base_leiden") +
#     labs(x = "Gene markers (top 5 ~ logFC)", y = "Cluster") +
#     theme(
#       plot.subtitle = element_text(face = "bold"),
#       legend.position = "none",
#       axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 6),
#       axis.text.y = element_text(size = 10),
#       axis.title.x = element_text(face = "bold", margin = margin(t = 2, unit = "pt")),
#       axis.title.y = element_text(face = "bold", margin = margin(r = 12, unit = "pt")),
#       legend.title = element_text(vjust = 1, hjust = 0),
#       panel.grid.major = element_line(color = "Gray90", size = 0.5),
#       plot.margin = unit(rep(18, 4), "pt"),
#       axis.line = element_blank(),
#       panel.border = element_rect(size = 0.5, color = "black")
#     )
#
#   top <- cowplot::plot_grid(cluster_umap, contam_umap, align = "hv")
#   full <- cowplot::plot_grid(top, scores_plot, markers_plot, ncol = 1)
#   cowplot::save_plot(filename = glue("/Users/jpeters/Projects/lungmps/plots/{x}_mp_filtering.png"),
#     plot = full, base_height = 12, base_asp = round(8.5/11, 2))
#
#   return(nonmp_cells)
# })
#
# remove <- list(c(9), c(12))
# names(remove) <- object_names
#
# nonremove <- list(NA, c(16))
# names(nonremove) <- object_names
#
# blacklist_hits <- pblapply(X = object_names, function(x) {
#
#   ui_info("Printing figure for {x}")
#   mp_scores <- object_ls[[x]][[]] %>% select(base_leiden, mp_diff_score)
#   nonmp_clusters <- as.character(remove[[x]])
#   nonmp_cells <- colnames(object_ls[[x]])[as.character(object_ls[[x]]$base_leiden) %in% nonmp_clusters]
#   nonmp_cells <- unique(c(nonmp_cells, rownames(mp_scores)[mp_scores$mp_diff_score < total_threshold]))
#   if (!is.na(nonremove[[x]])) nonmp_cells <- setdiff(nonmp_cells, colnames(object_ls[[x]])[as.character(object_ls[[x]]$base_leiden) %in% nonremove[[x]]])
#   ui_info("Identified {length(nonmp_cells)} contaminating cell profiles")
#
#   cluster_umap <- DimPlotPlus(object_ls[[x]], reduction = "harmony_umap", group.by = "base_leiden",
#     label = TRUE, label.size = 6, repel = TRUE) +
#     colorspace::scale_color_discrete_qualitative("Dark 3", name = "Cluster") +
#     labs(x = "UMAP 1", y = "UMAP 2") +
#     theme_void(base_size = 14) +
#     theme(
#       axis.title.x = element_text(hjust = 0, face = "plain", color = "gray60"),
#       axis.title.y = element_text(hjust = 0, angle = 90, face = "plain", color = "gray60"),
#       panel.background = element_rect(fill = "transparent", color = "transparent"), # bg of the panel
#       plot.background = element_rect(fill = "transparent", color = "transparent"), # bg of the plot
#       panel.grid.major = element_blank(), # get rid of major grid
#       panel.grid.minor = element_blank(), # get rid of minor grid
#       legend.background = element_rect(fill = "transparent"), # get rid of legend bg
#       legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
#       legend.position = "none",
#       plot.margin = margin(6, 6, 6, 6, unit = "mm")
#     )
#
#   contam_umap <- DimPlotPlus(object_ls[[x]], reduction = "harmony_umap", group.by = "base_leiden", pt.size = 1,
#     cells.highlight = nonmp_cells) +
#     scale_color_manual(values = c("gray90", "black"),
#       name = "Potential\ncontamination", labels = c("MP", paste(blacklist_hits, collapse = ", "))) +
#     labs(x = "UMAP 1", y = "UMAP 2") +
#     theme_void(base_size = 14) +
#     theme(
#       axis.title.x = element_text(hjust = 0, face = "plain", color = "gray60"),
#       axis.title.y = element_text(hjust = 0, angle = 90, face = "plain", color = "gray60"),
#       panel.background = element_rect(fill = "transparent", color = "transparent"), # bg of the panel
#       plot.background = element_rect(fill = "transparent", color = "transparent"), # bg of the plot
#       panel.grid.major = element_blank(), # get rid of major grid
#       panel.grid.minor = element_blank(), # get rid of minor grid
#       legend.background = element_rect(fill = "transparent"), # get rid of legend bg
#       legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
#       legend.position = "none",
#       plot.margin = margin(6, 6, 6, 6, unit = "mm")
#     )
#
#   title <- cowplot::ggdraw() +
#     cowplot::draw_label(
#       glue("{length(nonmp_cells)} contaminating cell profiles \nCluster(s) {paste(nonmp_clusters, collapse = ', ')}"),
#       fontface = 'bold',
#       x = 0,
#       hjust = 0
#     ) +
#     theme(
#       plot.margin = margin(0, 0, 0, 10)
#     )
#
#   removal <- cowplot::plot_grid(cluster_umap, contam_umap, align = "hv")
#   removal <- cowplot::plot_grid(title, removal, ncol = 1, rel_heights = c(0.1, 1))
#   cowplot::save_plot(filename = glue("/Users/jpeters/Projects/lungmps/plots/{x}_mp_filtered.png"),
#     plot = removal, base_height = 6, base_asp = 2)
#   return(nonmp_cells)
# })
#
# # [7] Confirm filtering ----
#
# object_ls <- pblapply(X = object_names, function(x) {
#   ui_done("Subsetting {length(setdiff(Cells(object_ls[[x]]), blacklist_hits[[x]]))} MPs from {x}")
#   mp_object <- object_ls[[x]][, setdiff(Cells(object_ls[[x]]), blacklist_hits[[x]])]
#   mp_object<- DietSeurat(object = mp_object)
#   Clean()
#   return(mp_object)
# })
#
# object_ls <- pblapply(X = object_names, function(x) {
#   object_ls[[x]] <- BaseProcess(object = object_ls[[x]], nfeatures = 4000)
#   object_ls[[x]] <- HarmonyCorrection(object = object_ls[[x]], batch = "batch")
#   Clean()
#   return(object_ls[[x]])
# })
#
# object_ls <- pblapply(X = object_names, function(x) {
#   object_ls[[x]] <- CalculatePCHeuristics(object_ls[[x]], reduction = "harmony", store.name = "harmony_heuristics")
#   dims <- ceiling(unique(object_ls[[x]]@misc$harmony_heuristics$tp))
#   ui_done("\nUsing {ui_field(dims)} dims for post-Harmony functions for dataset {x}")
#   object_ls[[x]]$postharmony_dims_used <- dims
#   Clean()
#   return(object_ls[[x]])
# })
# ui_info("Average number of post-Harmony dimensions: {median(unlist(lapply(object_ls, function(x) return(unique(x$postharmony_dims_used)))))}")
#
# object_ls <- pblapply(X = object_names, function(x) {
#   dims <- 20
#   #dims <- unique(object_ls[[x]]$postharmony_dims_used)
#   method <- "harmony"
#   object_ls[[x]] <- RunUMAP(object = object_ls[[x]], reduction = method,
#     dims = 1:dims, seed.use = 1, reduction.name = glue("{method}_umap"),
#     reduction.key = glue("{substr(method, 1, 1)}UMAP_"))
#   Clean()
#   return(object_ls[[x]])
# })
#
# object_ls <- pblapply(X = object_names, function(x) {
#   dims <- 20
#   object_ls[[x]] <- FindNeighbors(object_ls[[x]], reduction = "harmony", dims = 1:dims,
#     k.param = 20, prune.SNN = 1/15, compute.SNN = TRUE, graph.name = glue("harmony_snn"), verbose = TRUE)
#   object_ls[[x]] <- BaseCluster(object = object_ls[[x]], dims = 20, verbose = TRUE,
#     res.start = 1E-5, res.end = 1, num.res = 30, num.clusters.lower = 5, num.clusters.upper = 30, mod.similarity = 0.995)
#   snn <- igraph::graph_from_adjacency_matrix(object_ls[[x]]@graphs$harmony_snn, mode = "directed", weighted = TRUE)
#   walktrap_cluster <- igraph::cluster_walktrap(graph = snn, steps = 4)
#   object_ls[[x]]$base_walktrap <- walktrap_cluster$membership
#   Clean()
#   return(object_ls[[x]])
# })
#
# minimum.cells <- 10
# cluster <- "base_leiden"
# object_ls <- pblapply(X = object_names, function(x) {
#   ui_info("Processing {x}...")
#   print(table(object_ls[[x]]$base_leiden))
#   object_ls[[x]] <- SetIdent(object_ls[[x]], value = cluster)
#   object_ls[[x]] <- BuildClusterTree(object_ls[[x]], reorder = TRUE, reorder.numeric = TRUE)
#   object_ls[[x]]$base_clustering <- object_ls[[x]]@active.ident
#   return(object_ls[[x]])
# })
#
# object_ls <- pblapply(X = object_names, function(x) {
#
#   object_ls[[x]] <- SetIdent(object = object_ls[[x]], value = "base_walktrap")
#   walktrap_markers <- presto::wilcoxauc(object_ls[[x]], assay = "data", seurat_assay = "RNA")
#   object_ls[[x]]@misc$walktrap_markers <- walktrap_markers
#
#   object_ls[[x]] <- SetIdent(object = object_ls[[x]], value = "base_leiden")
#   leiden_markers <- presto::wilcoxauc(object_ls[[x]], assay = "data", seurat_assay = "RNA")
#   object_ls[[x]]@misc$leiden_markers <- leiden_markers
#
#   Clean()
#   return(object_ls[[x]])
# })
#
# total_mp_scores <- pblapply(object_names, function(x) {
#   return(object_ls[[x]][[]] %>% select(base_leiden, mp_diff_score, cell_barcode))
# })
# total_mp_scores <- Reduce(rbind, total_mp_scores)
# total_threshold <- median(total_mp_scores$mp_diff_score) - 3*mad(total_mp_scores$mp_diff_score)
# total_threshold <- 0
# mean(total_mp_scores$mp_diff_score) - 2*sd(total_mp_scores$mp_diff_score)
# hist(total_mp_scores$mp_diff_score)
#
# blacklist_hits <- pblapply(X = object_names, function(x) {
#
#   ui_info("Printing figure for {x}")
#   object_ls[[x]]$base_leiden <- as.factor(object_ls[[x]]$base_leiden)
#   object_ls[[x]]$base_leiden <- factor(object_ls[[x]]$base_leiden, levels = levels(object_ls[[x]]$base_leiden)[order(as.numeric(levels(object_ls[[x]]$base_leiden)))])
#   mp_scores <- object_ls[[x]][[]] %>% select(base_leiden, mp_diff_score)
#
#   plot_markers <- object_ls[[x]]@misc$leiden_markers %>%
#     dplyr::filter(padj < 1E-10 & logFC > log(1.25) & auc > 0.5) %>%
#     group_by(group) %>%
#     arrange(as.numeric(group)) %>%
#     top_n(10, logFC)
#   blacklist_hits <- plot_markers[(plot_markers$feature %in% c(nonmp_markers)), ]
#   blacklist_hits <- names(table(blacklist_hits$group))[table(blacklist_hits$group) >= 2]
#
#   nonmp_clusters <- as.character(mp_scores %>%
#       group_by(base_leiden) %>%
#       summarize(median_mp_score = median(mp_diff_score)) %>%
#       filter(median_mp_score < 0) %>%
#       pull(base_leiden))
#   nonmp_cells <- colnames(object_ls[[x]])[as.character(object_ls[[x]]$base_leiden) %in% unique(c(blacklist_hits, nonmp_clusters))]
#   nonmp_cells <- unique(c(nonmp_cells, rownames(mp_scores)[mp_scores$mp_diff_score < total_threshold]))
#   ui_info("Identified {length(nonmp_cells)} contaminating cell profiles")
#
#   genes <- unique(plot_markers %>%
#       group_by(group) %>%
#       top_n(5, logFC) %>% pull(feature))
#
#   cluster_umap <- DimPlotPlus(object_ls[[x]], reduction = "harmony_umap", group.by = "base_leiden",
#     label = TRUE, label.size = 6, repel = TRUE) +
#     colorspace::scale_color_discrete_qualitative("Dark 3", name = "Cluster") +
#     labs(x = "UMAP 1", y = "UMAP 2") +
#     theme_void(base_size = 14) +
#     theme(
#       axis.title.x = element_text(hjust = 0, face = "plain", color = "gray60"),
#       axis.title.y = element_text(hjust = 0, angle = 90, face = "plain", color = "gray60"),
#       panel.background = element_rect(fill = "transparent", color = "transparent"), # bg of the panel
#       plot.background = element_rect(fill = "transparent", color = "transparent"), # bg of the plot
#       panel.grid.major = element_blank(), # get rid of major grid
#       panel.grid.minor = element_blank(), # get rid of minor grid
#       legend.background = element_rect(fill = "transparent"), # get rid of legend bg
#       legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
#       legend.position = "none",
#       plot.margin = margin(6, 6, 6, 6, unit = "mm")
#     )
#
#   contam_umap <- DimPlotPlus(object_ls[[x]], reduction = "harmony_umap", group.by = "base_leiden", pt.size = 1,
#     cells.highlight = nonmp_cells) +
#     scale_color_manual(values = c("gray90", "black"),
#       name = "Potential\ncontamination", labels = c("MP", paste(blacklist_hits, collapse = ", "))) +
#     labs(x = "UMAP 1", y = "UMAP 2") +
#     theme_void(base_size = 14) +
#     theme(
#       axis.title.x = element_text(hjust = 0, face = "plain", color = "gray60"),
#       axis.title.y = element_text(hjust = 0, angle = 90, face = "plain", color = "gray60"),
#       panel.background = element_rect(fill = "transparent", color = "transparent"), # bg of the panel
#       plot.background = element_rect(fill = "transparent", color = "transparent"), # bg of the plot
#       panel.grid.major = element_blank(), # get rid of major grid
#       panel.grid.minor = element_blank(), # get rid of minor grid
#       legend.background = element_rect(fill = "transparent"), # get rid of legend bg
#       legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
#       legend.position = "none",
#       plot.margin = margin(6, 6, 6, 6, unit = "mm")
#     )
#
#   scores_plot <- DotPlot_ol(object_ls[[x]], features = score_names_touse,
#     cols = rev(colorspace::sequential_hcl(2, "Blues")),
#     dot.min = 0.05, dot.scale = 4, scale.by = "radius", group.by = "base_leiden") +
#     scale_x_discrete(labels = rev(score_labels)) +
#     labs(x = "Expression Score", y = "Cluster") +
#     theme(
#       plot.subtitle = element_text(face = "bold"),
#       legend.position = "none",
#       axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 6),
#       axis.text.y = element_text(size = 10),
#       axis.title.x = element_text(face = "bold", margin = margin(t = 2, unit = "pt")),
#       axis.title.y = element_text(face = "bold", margin = margin(r = 12, unit = "pt")),
#       legend.title = element_text(vjust = 1, hjust = 0),
#       panel.grid.major = element_line(color = "Gray90", size = 0.5),
#       plot.margin = unit(rep(18, 4), "pt"),
#       axis.line = element_blank(),
#       panel.border = element_rect(size = 0.5, color = "black")
#     )
#
#   markers_plot <- DotPlot_ol(object_ls[[x]], features = genes,
#     cols = rev(colorspace::sequential_hcl(2, "Blues")),
#     dot.min = 0.05, dot.scale = 6, scale.by = "radius", group.by = "base_leiden") +
#     labs(x = "Gene markers (top 5 ~ logFC)", y = "Cluster") +
#     theme(
#       plot.subtitle = element_text(face = "bold"),
#       legend.position = "none",
#       axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 6),
#       axis.text.y = element_text(size = 10),
#       axis.title.x = element_text(face = "bold", margin = margin(t = 2, unit = "pt")),
#       axis.title.y = element_text(face = "bold", margin = margin(r = 12, unit = "pt")),
#       legend.title = element_text(vjust = 1, hjust = 0),
#       panel.grid.major = element_line(color = "Gray90", size = 0.5),
#       plot.margin = unit(rep(18, 4), "pt"),
#       axis.line = element_blank(),
#       panel.border = element_rect(size = 0.5, color = "black")
#     )
#
#   top <- cowplot::plot_grid(cluster_umap, contam_umap, align = "hv")
#   full <- cowplot::plot_grid(top, scores_plot, markers_plot, ncol = 1)
#   cowplot::save_plot(filename = glue("/Users/jpeters/Projects/lungmps/plots/{x}_mp_filtering_second.png"),
#     plot = full, base_height = 12, base_asp = round(8.5/11, 2))
#   return(nonmp_cells)
# })
#
# # [8] Finalize filtering and save ----
#
# object_ls <- pblapply(X = object_names, function(x) {
#
#   results <- IterativeCluster(object_ls[[x]], )
#   plotgrid <- cowplot::plot_grid(plotlist = results$plots, align = "hv", labels = "auto")
#   cowplot::save_plot(filename = glue("/Users/jpeters/Projects/lungmps/plots/{x}_mnp_cluster_merges_{stamp[1]}.png"),
#     plot = plotgrid, base_height = 12, base_asp = round(16/9, 2))
#   object_ls[[x]] <- results$object
#   object_ls[[x]][["merged_leiden"]] <- Seurat::Idents(object_ls[[x]])
#   markers <- IdentifyQualityMarkers(object_ls[[x]], tlogFC = log(1.1))
#   object_ls[[x]]@misc$merged_clusters <- markers
#
#   res_df <- object_ls[[x]]@misc$leiden_results$second_results
#   num_current <- length(unique(Idents(object_ls[[x]])))
#   res <- res_df$resolution_parameter[which.min(abs(res_df$cluster_count - num_current))]
#   g <- igraph::graph_from_adjacency_matrix(adjmatrix = object_ls[[x]]@graphs$harmony_snn, mode = "directed", weighted = TRUE)
#   final_result <- leidenbase::leiden_find_partition(g, partition_type = "CPMVertexPartition", seed = 1,
#     resolution_parameter = res, num_iter = 30, verbose = TRUE)
#   out_result <- list(membership = final_result[['membership']],
#     modularity = final_result[['modularity']])
#   names(out_result$membership) <- colnames(object_ls[[x]]@graphs$harmony_snn)
#   if (any(table(out_result$membership) < 5)) {
#     ids <- GroupSingletons(ids = out_result$membership, SNN = object_ls[[x]]@graphs$harmony_snn,
#       min.size = 5, group.singletons = TRUE, verbose = TRUE)
#     object_ls[[x]] <- AddMetaData(object_ls[[x]], ids, col.name = "equiv_leiden")
#   } else {
#     object_ls[[x]] <- AddMetaData(object_ls[[x]], out_result$membership, col.name = "equiv_leiden")
#   }
#   object_ls[[x]] <- SetIdent(object_ls[[x]], value = "equiv_leiden")
#   object_ls[[x]]@misc$leiden_mergedequiv_markers <- IdentifyQualityMarkers(object_ls[[x]], tlogFC = log(1.1))
#
#   a <- PlotAnnotatedClusters(object_ls[[x]], group.by = "merged_leiden")
#   b <- PlotAnnotatedClusters(object_ls[[x]], group.by = "equiv_leiden")
#   title <- cowplot::ggdraw() +
#     cowplot::draw_label(
#       glue::glue("Clustering results for {x} (ARI = {round(mclust::adjustedRandIndex(as.character(object_ls[[x]]$merged_leiden), as.character(object_ls[[x]]$equiv_leiden)), 3)})"),
#       fontface = 'bold', x = 0, hjust = 0) +
#     theme(plot.margin = margin(0, 0, 0, 12, unit = "pt"))
#   plotgrid <- cowplot::plot_grid(title, cowplot::plot_grid(a, b, labels = "auto", align = "hv"), ncol = 1, rel_heights = c(0.05, 1))
#   cowplot::save_plot(filename = glue("/Users/jpeters/Projects/lungmps/plots/{x}_mnp_cluster_results_{stamp[1]}.png"),
#     plot = plotgrid, base_height = 8, base_asp = 2)
#   Clean()
#   return(object_ls[[x]])
# })
#
# pblapply(X = object_names, function(x) {
#   ui_info("\nSaving {x} ({round(as.numeric(object.size(object_ls[[x]])/1E6), 2)} mb)")
#   saveRDS(object_ls[[x]],
#     file = glue::glue("/Users/jpeters/Projects/lungmps/data/objects/{x}_mnp_final.rds"))
#   Clean()
#   ui_done(("\nCompleted {x}"))
# })
#
#
#
#
#
