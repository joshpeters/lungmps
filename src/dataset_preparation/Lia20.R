# ---
# Description: Prepare data from Liao et al. 2020
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
set <- "Liao_2020"
Clean()

# [2] Load files ----

path <- "data/raw_data/GSE145926_RAW/"
meta <- read_tsv("data/raw_data/GSE145926_RAW/GSE145926_series_matrix.txt")
file_paths <- list.files(path, full.names = TRUE)
file_names <- list.files(path, full.names = FALSE)
file_paths <- file_paths[grep(".h5$", file_paths)]
file_names <- file_names[grep(".h5$", file_names)]

h5s <- map(file_paths, Read10X_h5)
file_names <- str_match(file_names, "(GSM.*)_(C[[:digit:]]{2,3})_.*\\.h5$")
meta <- merge(file_names, meta, by.x = "V3", by.y = "sample_id", sort = FALSE)
meta$V2 <- NULL
colnames(meta)[c(1,2)] <- c("sample_id", "filename")
names(h5s) <- meta$sample_id

# [3] Prepare object ----

objects <- imap(h5s, ~ {
  object_meta <- meta[meta$sample_id == .y, ] %>% slice(rep(1:n(), times = ncol(.x)))
  rownames(object_meta) <- colnames(.x)
  object <- CreateSeuratObject(.x, project = set, meta.data = object_meta, min.cells = 0, min.features = 100)
  return(object)
})

object <- merge(objects[[1]], objects[c(2:length(objects))])
dim(object)

Check5(object[[]])
colnames(object[[]])
object@meta.data <- janitor::clean_names(object@meta.data)
colnames(object[[]])
colnames(object@meta.data)[c(1, 2, 3)] <- c("orig.ident", "nCount_RNA", "nFeature_RNA")
colnames(object[[]])
object <- AddExprMeta(object)

names(object[[]])
head(object[[]])
object$batch <- object$sample_id
object$study <- set

object$cell_barcode <- rownames(object[[]])
object <- AddBatchFreq(object, "batch")
object <- Scrub(object, batch.var = "batch")

# Harmonizing cell_barcode, study, donor, condition, sample_id, batch, disease_classification, celltype
colnames(object[[]])
head(object[[]])
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
# table(object_ls[[i]]$batch)
# table(object_ls[[i]]$condition)
# table(object_ls[[i]]$sample_id)
# object_ls[[i]]$donor <- object_ls[[i]]$sample_id
# object_ls[[i]]$orig.ident <- object_ls[[i]]$batch
# object_ls[[i]]$disease_classification <- "disease"
# object_ls[[i]]$disease_classification[object_ls[[i]]$condition == "Healthy"] <- "health"
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
# object_ls <- split_object_ls
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
# object_names
# DimPlot(object_ls[[1]], reduction = "harmony_umap", group.by = "condition")
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
#   b <- ggplot(data.frame(x = mm$x), aes(x, y=..density..)) +
#     geom_histogram(fill = NA, color = "black") +
#     geom_polygon(data = em_df, aes(x, y, fill = comp), color = "black", alpha = 0.5) +
#     geom_vline(xintercept = threshold, color = "gray40", linetype = "dashed") +
#     colorspace::scale_fill_discrete_sequential("Light Grays", name = "Component Means",
#       labels = format(em_df$mu, digits = 2)) +
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
# object_names <- names(object_ls)
# object_names <- SelfName(object_names)
#
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
# ui_done("Total threshold for MP score: {total_threshold}")
# total_threshold <- 0
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
# remove <- list(c(17, 21, 22, 23), c(6, 17))
# names(remove) <- object_names
#
# nonremove <- list(c(18), c(13))
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
#   cowplot::save_plot(filename = glue("/Users/jpeters/Projects/lungmps/plots/{x}_mp_filtering_fourth.png"),
#     plot = full, base_height = 12, base_asp = round(8.5/11, 2))
#   return(nonmp_cells)
# })
#
# DimPlot(object_ls[[2]], group.by = "base_leiden", label = TRUE)
# VlnPlot(object_ls[[2]], c("FCGR3A", "NAMPT", "CD14", "CSF3R"), group.by = "base_leiden", pt.size = 0.1)
# DotPlot(object_ls[[2]], features = c("FCGR3A", "NAMPT", "CD14", "CSF3R"), group.by = "base_leiden")
# # remove <- list(c(16, 17, 13), c(2, 13, 12))
# # names(remove) <- object_names
# remove <- list(c(7), c(NA))
# names(remove) <- object_names
#
# nonremove <- list(c(NA), c(NA))
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
#   cowplot::save_plot(filename = glue("/Users/jpeters/Projects/lungmps/plots/{x}_mp_second_filtered.png"),
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
#   cowplot::save_plot(filename = glue("/Users/jpeters/Projects/lungmps/plots/{x}_mp_filtering_5th.png"),
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
# # [9] Apply scores to clusters ----
#
# object_ls <- LoadObjectList("data/objects", "Bost_2020_(disease|health)_mnp_final.rds", "_mnp_final.rds", TRUE)
# object_names <- names(object_ls)
# object_names <- SelfName(object_names)
#
# i = 2
# #FeaturePlot(object_ls[[i] ], features = c("percent_ribo"))
# cluster_umap <- DimPlotPlus(object_ls[[i]], reduction = "harmony_umap", group.by = "merged_leiden",
#   label = TRUE, label.size = 6, repel = TRUE) +
#   colorspace::scale_color_discrete_qualitative("Dark 3", name = "Cluster") +
#   labs(x = "UMAP 1", y = "UMAP 2") +
#   GeneralTheme(base_size = 14) +
#   theme(
#     axis.text = element_blank(),
#     axis.ticks = element_blank(),
#     axis.title.x = element_text(hjust = 0, face = "bold", color = "black", margin = margin(10,0,0,0)),
#     axis.title.y = element_text(hjust = 0, angle = 90, face = "bold", color = "black", margin = margin(0,10,0,0)),
#     panel.background = element_rect(fill = "transparent", color = "black", size = 1), # bg of the panel
#     plot.background = element_rect(fill = "transparent", color = "transparent"), # bg of the plot
#     panel.grid.major = element_blank(), # get rid of major grid
#     panel.grid.minor = element_blank(), # get rid of minor grid
#     legend.background = element_rect(fill = "transparent"), # get rid of legend bg
#     legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
#     legend.position = "none",
#     plot.margin = margin(6, 6, 6, 6, unit = "mm")
#   )
# cluster_umap
# #SavePlot("plots/bost2020_umap_mergedclusters.png", plot = cluster_umap)
#
# library(ggdendro)
# library(dendextend)
# PlotClusterTree(object_ls[[i]])
# tree <- as.dendrogram(object_ls[[i]]@tools$BuildClusterTree)
# ddata <- dendro_data(tree, type = "rectangle")
# p <- ggplot(segment(ddata)) +
#   geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
#   geom_text(data = ddata$labels,
#     aes(x = x, y = y, label = label), size = 4, vjust = 0.5, hjust = 1) +
#   ylim(c(-2, round(max(ddata$segments$yend))+2)) +
#   coord_flip() +
#   theme_dendro() +
#   theme(plot.margin = unit(c(0.01, 0.01, 0.01, 0.01), "cm"))
# #umap_tree <- cowplot::plot_grid(cluster_umap, p, rel_widths = c(1, 0.5))
# #umap_tree
# #SavePlot("plots/bost2020_disease_umaptree.png", base_asp = 1.5, plot = umap_tree)
#
# genes <- unique(object_ls[[i]]@misc$merged_clusters %>%
#     group_by(group) %>%
#     top_n(2, logFC) %>% select(feature, group))
# markers_plot <- DotPlot_ol(object_ls[[i]], features = genes$feature,
#   cols = rev(colorspace::sequential_hcl(2, "Grays")),
#   dot.min = 0.05, dot.scale = 4, scale.by = "radius", group.by = "merged_leiden") +
#   labs(x = "", y = "Cluster") +
#   theme(
#     plot.subtitle = element_text(face = "bold"),
#     legend.position = "none",
#     axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 8),
#     axis.text.y = element_text(size = 10),
#     axis.title.x = element_text(face = "bold", margin = margin(t = 2, unit = "pt")),
#     axis.title.y = element_text(face = "bold", margin = margin(r = 12, unit = "pt")),
#     legend.title = element_text(vjust = 1, hjust = 0),
#     panel.grid.major = element_line(color = "Gray90", size = 0.5),
#     plot.margin = unit(rep(18, 4), "pt"),
#     axis.line = element_blank(),
#     panel.border = element_rect(size = 0.5, color = "black")
#   )
# #markers_plot
# #SavePlot("plots/bost2020_disease_markers.png", base_asp = 2, base_height= 4, plot = markers_plot)
#
# total <- cowplot::plot_grid(cluster_umap, markers_plot, NULL, p, rel_widths = c(1, 2, 0, 0.75), nrow = 1, align = "h")
# total_disease <- total
# SavePlot("plots/bost2020d_summary.pdf", base_asp = 6.5/2, base_height = 4, plot = total_disease)
# total_health <- total
# SavePlot("plots/bost2020h_summary.pdf", base_asp = 6.5/2, base_height = 4, plot = total_health)
#
#
# programs <- readRDS("/Users/jpeters/Projects/lungmps/data/processed/consensusGEPs_q95_k23.rds")
# programs_ls <- split(programs$gene, programs$partition)
# programs_ls <- programs_ls[sapply(programs_ls, length) >= 20]
# names(programs_ls) <- paste0("P", "_", names(programs_ls))
# programs <- names(programs_ls)
# programs <- SelfName(programs)
#
# markers <- readRDS("data/processed/stmmarkers_consensus_20200511.rds")
# markers_ls <- split(markers$feature, markers$partition)
# markers_ls <- markers_ls[sapply(markers_ls, length) >= 20]
# names(markers_ls) <- paste0("M", "_", names(markers_ls))
# markers <- names(markers_ls)
# markers <- SelfName(markers)
#
# # generate all scores across all programs
# group <- "merged_leiden"
# object_ls <- lapply(object_names, function(x) {
#   object_ls[[x]]$set <- x
#   return(object_ls[[x]])
# })
#
# merged_object <- merge(object_ls[[1]], object_ls[[2]], add.cell.ids = c("d", "h"), merge.data = TRUE)
# obj = merged_object
#
# genesets <- programs_ls
# names <- programs
# group <- "merged_leiden"
#
# ui_todo("Preparing gene sets...")
# genesets <- map(names, ~ {
#   genesets[[.x]] <- genesets[[.x]][genesets[[.x]] %in% rownames(obj)]
#   return(genesets[[.x]])
# }, obj = obj)
# percents <- map_dfr(.x = names, ~ {
#   count_subset <- GetAssayData(obj, slot = "counts")[genesets[[.x]], ]
#   percents <- Matrix::colSums(count_subset > 0)
#   percents <- percents/nrow(count_subset)
#   percents <- data.frame(percents = percents, geneset = .x, cell_barcode = colnames(obj))
#   return(percents)
# }, obj = obj)
# ui_done("Gene sets prepared")
# assertthat::assert_that(all(lapply(genesets, length) > 0))
#
# ui_todo("Scoring object...")
# obj <- AddScore(object = obj, features = genesets, nbin = 30, ctrl = 30, name = "Geneset", seed = 1)
# assertthat::assert_that(length(grep("Geneset[[:digit:]]", colnames(obj@meta.data))) == length(genesets))
# colnames(obj@meta.data)[grep("Geneset[[:digit:]]", colnames(obj@meta.data))] <- names
# ui_done("Object scored")
#
# df <- obj@meta.data[, c("nCount_RNA", "disease_classification", "set", group, "batch", names)]
# df <- df %>% rownames_to_column("cell_barcode")
# df <- df %>% gather(key = "program", value = "score", names)
# df <- cbind(df, percents[, c("percents", "geneset")])
# df$source <- set_source
# df
# colnames(df)[5] <- "group"
# colnames(df)[10] <- "percents_geneset"
# colnames(df)[7] <- "geneset"
#
# marker_scores <- df
# program_scores <- df
# all_scores <- rbind(marker_scores, program_scores)
# all_scores$raw_cell_barcode <- gsub("^(d|h)_", "", all_scores$cell_barcode)
# conditions <- obj@meta.data %>% select(cell_barcode, condition)
# all_scores <- merge(all_scores, conditions, by.x = "raw_cell_barcode", by.y = "cell_barcode", all.x = TRUE)
# head(all_scores)
#
# sum_mscores <- all_scores %>%
#   group_by(condition, geneset) %>%
#   filter(geneset %in% markers) %>%
#   summarize(score = mean(score), percent = median(percents))
# # plot markers list first
# mscore_levels <- sum_mscores %>% filter(condition == "Severe") %>% arrange(desc(score)) %>% pull(geneset)
# sum_mscores$geneset <- factor(sum_mscores$geneset, levels = mscore_levels)
# reordered_marker_names <- filtered_markers_ls[mscore_levels]
# #unique(scores$program)[order(as.numeric(gsub("M_", "", unique(scores$program))))]
# #sum_scores$geneset <- factor(sum_scores$geneset, levels = levels_ordered)
# mscores_dotplot <- ggplot(sum_mscores, aes(x = geneset, y = condition, fill = score, size = percent)) +
#   geom_point(shape = 21, color = "black") +
#   scale_x_discrete(labels = paste(gsub("M_", "", levels(sum_mscores$geneset)), reordered_marker_names, sep = " ")) +
#   #scale_y_discrete(labels = gsub("M_", "", levels_ordered)) +
#   scale_size_continuous(range = c(2, 8), name = "Average \n% detected", breaks = c(0.25, 0.50, 0.75)) +
#   colorspace::scale_fill_continuous_diverging(palette = "Blue-Red 3", name = "Average\nScore",
#     limits = c(-0.5, 0.5), breaks = c(-0.25, 0, 0.25)) +
#   labs(x = "Marker geneset", y = "Condition") +
#   GeneralTheme(18) + theme(
#     axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
#     legend.box = "bottom",
#     legend.direction = "vertical",
#     panel.grid.major = element_line(size = 0.5, color = "gray90"),
#     panel.grid.minor = element_line(size = 0.5, color = "gray90")
#   )
# mscores_dotplot
#
# sum_pscores <- all_scores %>%
#   group_by(condition, geneset) %>%
#   filter(geneset %in% programs) %>%
#   summarize(score = mean(score), percent = median(percents))
# pscore_levels <- sum_pscores %>% filter(condition == "Severe") %>% arrange(desc(score)) %>% pull(geneset)
# sum_pscores$geneset <- factor(sum_pscores$geneset, levels = pscore_levels)
# reordered_program_names <- program_gene_names[pscore_levels]
# #unique(scores$program)[order(as.numeric(gsub("M_", "", unique(scores$program))))]
# #sum_scores$geneset <- factor(sum_scores$geneset, levels = levels_ordered)
# pscores_dotplot <- ggplot(sum_pscores, aes(x = geneset, y = condition, fill = score, size = percent)) +
#   geom_point(shape = 21, color = "black") +
#   scale_x_discrete(labels = paste(gsub("P_", "", levels(sum_pscores$geneset)), reordered_program_names, sep = " ")) +
#   scale_size_continuous(range = c(2, 8), name = "Average \n% detected", breaks = c(0.25, 0.50, 0.75)) +
#   colorspace::scale_fill_continuous_diverging(palette = "Blue-Red 3", name = "Average\nScore",
#     limits = c(-0.5, 0.5), breaks = c(-0.25, 0, 0.25)) +
#   labs(x = "Program geneset", y = "Condition") +
#   GeneralTheme(18) + theme(
#     axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
#     legend.position = "none",
#     legend.box = "top",
#     legend.direction = "vertical",
#     panel.grid.major = element_line(size = 0.5, color = "gray90"),
#     panel.grid.minor = element_line(size = 0.5, color = "gray90")
#   )
# pscores_dotplot
#
# ab <- cowplot::plot_grid(mscores_dotplot, pscores_dotplot, nrow = 2, rel_widths = c(1, 1), hjust = -1, align = "hv", axis = "r")
# ab
# SavePlot(filename = "plots/bost_scores_{stamp$date}.pdf", plot = ab, base_asp = 1.75, base_height = 1.9*3)
#
# # all_scores_plot <- all_scores
# # all_scores_plot$score[all_scores_plot$score > 1] <- 1
# # scoresdist_plot <- ggplot(all_scores_plot, aes(x = score, y = geneset, fill = condition)) +
# #   ggridges::geom_density_ridges(color = "black", alpha = 0.6) +
# #   #scale_fill_manual(values = c("#0072B2", "#D55E00"), labels = c("Health", "Disease"), name = "Condition") +
# #   #scale_y_discrete(labels = gsub("M_", "", levels_ordered)) +
# #   labs(x = "Score", y = "Marker gene set", title = "Consensus marker gene sets' scores", subtitle = "COVID-19 BAL samples") +
# #   ggridges::theme_ridges(center_axis_labels = TRUE) +
# #   theme(
# #     axis.title = element_text(face = "bold"),
# #     legend.title = element_text(face = "bold", size = 12),
# #     legend.text = element_text(size = 12),
# #     legend.position = "right",
# #     legend.justification = "top")
# # scoresdist_plot
# # SavePlot("plots/bost2020_marker_scores_dist.png", base_asp = 1, plot = scoresdist_plot)
#
# all_scores_plot <- all_scores %>% filter(geneset %in% c("P_3", "P_9")) %>% select(-percents, -percents_geneset)
# all_scores_plot <- all_scores_plot %>% spread(key = geneset, value = score)
# shiftplot <- ggplot(all_scores_plot, aes(x = P_3, y = P_9, color = condition)) +
#   geom_vline(xintercept = 0, color = "Gray80", linetype = "dotted", size = 1) +
#   geom_hline(yintercept = 0, color = "Gray80", linetype = "dotted", size = 1) +
#   geom_density_2d(size = 1) +
#   scale_color_manual(values = c("#009E73", "#0072B2", "#D55E00"), labels = c("Healthy", "Mild", "Severe"), name = "Condition") +
#   labs(x = "Program 3 Score", y = "Program 9 Score") +
#   GeneralTheme(base_size = 18)
# SavePlot("plots/bost2020_program39_shift.pdf", base_asp = 1.2, plot = shiftplot)
#
# # umap_coords <- object_ls[[1]]@reductions$harmony_umap@cell.embeddings
# # umapscores <- merge(all_scores %>% filter(disease_classification == "disease"), umap_coords, by.x = "raw_cell_barcode", by.y = 0, all.x = TRUE)
# # featureplot <- ggplot(umapscores %>% filter(geneset == "M_4") %>% arrange(score), aes(x = hUMAP_1, y = hUMAP_2, fill = score)) +
# #   geom_point(shape = 21, size = 1, stroke = 0.1) +
# #   #colorspace::scale_fill_continuous_diverging("Blue-Red 3", limits = c(-1, 3), name = "Score") +
# #   labs(x = "UMAP1", y = "UMAP2") +
# #   PointTheme(14) + theme(
# #     legend.title = element_text(face = "bold"),
# #     panel.grid = element_blank(),
# #     panel.grid.major = element_blank(),
# #     panel.grid.minor = element_blank(),
# #     axis.ticks = element_blank(),
# #     axis.text = element_blank())
# # featureplot
# # SavePlot("plots/bost2020_disease_m8_umap.png", base_asp = 1.2, base_height = 4, plot = featureplot)
#
#
# percent_above <- map_dfr(unique(all_scores$geneset), .f = ~ {
#   test_scores <- all_scores %>% filter(geneset == .x)
#   mm <- mixtools::normalmixEM(test_scores$score, k = 2, maxit = 10000, maxrestarts = 30)
#   min_mm <- which.min(mm$mu)
#   threshold <- mm$mu[min_mm] + mm$sigma[min_mm]*2
#   percent_above <- test_scores %>% group_by(batch) %>% summarize(n = n(), nabove = sum(score >= threshold), pabove = nabove/n*100, condition = unique(condition), nCount_RNA = mean(nCount_RNA))
#   percent_above$geneset <- .x
#   return(percent_above)
# })
#
# use_genesets <- c("M_1", "M_4", "M_5", "M_6", "M_8", "M_10", "M_11", "M_17", "P_1", "P_5", "P_7", "P_9", "P_10", "P_15", "P_16", "P_22")
# fpercent_above <- percent_above %>% filter(geneset %in% use_genesets)
# mlabels <- filtered_markers_ls[names(filtered_markers_ls) %in% use_genesets]
# plabels <- program_gene_names[names(program_gene_names) %in% use_genesets]
# labels <- c(mlabels, plabels)
# labels <- labels[use_genesets]
# labels <- paste0(gsub("_", "", use_genesets), " ", labels)
# labels
#
# fpercent_above$geneset <- factor(fpercent_above$geneset,
#   levels =  c("M_1", "M_4", "M_5", "M_6", "M_8", "M_10", "M_11", "M_17", "P_1", "P_5", "P_7", "P_9", "P_10", "P_15", "P_16", "P_22"),
#   labels = labels)
# fpercent_above$geneset
#
# pabove_plot <- ggplot(fpercent_above,
#   aes(x = condition, y = pabove, fill = condition)) +
#   geom_boxplot(alpha = 0.5) +
#   geom_point(shape = 21) +
#   guides(fill = FALSE) +
#   labs(x = "Condition", y = "% positive") +
#   scale_y_continuous(labels = function(x) paste0(x, "%")) +
#   scale_fill_manual(values = c("#009E73", "#0072B2", "#D55E00"), labels = c("Healthy", "Mild", "Severe"), name = "Condition") +
#   GeneralTheme(18) +
#   facet_wrap(~ geneset, nrow = 2) + theme(
#     strip.background = element_rect(size = 1),
#     strip.text = element_text(size = 14),
#     axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
#   )
# pabove_plot
# SavePlot(plot = pabove_plot, filename = "plots/bost_percentabove.pdf", base_height = 1.9*3, base_asp = 4/1.9)
#
# # Working -----------------------------------------------------------------
# # geneset <- geneset[geneset %in% rownames(object)]
# # count_subset <- GetAssayData(object, slot = "counts")[geneset, ]
# # percents <- Matrix::colSums(count_subset > 0)
# # percents <- percents/nrow(count_subset)
# # object <- AddModuleScore(object, features = list(geneset), name = "Geneset", nbin = 30, ctrl = 100)
# # if (any(grepl(name, colnames(object@meta.data)))) {
# #   ui_info("{name} already present in object, replacing")
# #   object@meta.data[grep(name, colnames(object@meta.data))] <- NULL
# # }
# # colnames(object@meta.data)[grep("Geneset1", colnames(object@meta.data))] <- name
# #
# # df <- object@meta.data %>% select(!!sym(group), batch, !!sym(name), nCount_RNA, set, condition)
# # df <- df %>% rownames_to_column("cell_barcode")
# # df <- df %>% gather(key = "program", value = "score", -!!sym(group), -cell_barcode, -batch, -nCount_RNA, -set, -condition)
# # df$percent <- percents
# #
# # ggplot(df, aes(x = score, fill = condition)) +
# #   geom_density() +
# #   GeneralTheme(14)
# #
# # ggplot(df, aes(x = score, y = nCount_RNA, color = set)) +
# #   geom_density2d() +
# #   GeneralTheme(14)
# #
# # ggplot(df, aes(x = condition, y = score, fill = batch)) +
# #   #geom_jitter(width = 0.25, alpha = 0.2) +
# #   geom_violin(alpha = 0.5, width = 0.5) +
# #   GeneralTheme(14)
# #
# # hist(df$score)
# # multimode::modetest(df$score, method = "SI")
# #
# # library(mclust)
# # x.gmm = mclust::Mclust(df$score, G = 2)
# # summary(x.gmm)
# # x.gmm.1 = Mclust(df$score, G = 1)
# # logLik(x.gmm.1)
# # ll <- logLik(x.gmm) - logLik(x.gmm.1)
# # 1 - pchisq(as.numeric(ll), df = 3)
# # min <- which.min(x.gmm$parameters$mean)
# # threshold <- x.gmm$parameters$mean[min] + 3 * sqrt(x.gmm$parameters$variance$sigmasq[1])
# #
# # gated_pops <- df %>% group_by(condition) %>% summarize(n = n(), n_above = sum(score > threshold), frac_above = n_above/n)
# # gated_pops$plot <- round(gated_pops$frac_above*100, 2)
# #
# # #raw_scores <- map_dfr(.x = markers, ~ ScoreObject(object = merged_object, geneset = markers_ls[[.x]], name = .x, group = group))
# #
# # mmodel <- lme4::lmer(formula = score ~ condition + (1 | batch), data = df %>% filter(condition %in% c("Healthy", "Severe")))
# # lmtest::lrtest(mmodel, 1)[2, 5]
# #
# # t.test(df$score[df$condition == "Healthy"], df$score[df$condition == "Severe"])
# # t.test(df$score[df$condition == "Healthy"], df$score[df$condition == "Mild"])
# # t.test(df$score[df$condition == "Mild"], df$score[df$condition == "Severe"])
# #
# # wilcox.test(df$score[df$condition == "Healthy"], df$score[df$condition == "Severe"])
# # wilcox.test(df$score[df$condition == "Healthy"], df$score[df$condition == "Mild"])
# # wilcox.test(df$score[df$condition == "Mild"], df$score[df$condition == "Severe"])
#
# # scaled average expression score (not centered)
# # norm_expr <- GetAssayData(merged_object, slot = "data")
# # norm_expr <- norm_expr[geneset, ]
# # scaled_expr <- apply(norm_expr, 1, scale, center = FALSE)
# # scaled_expr <- t(scaled_expr)
# # sum_expr <- rowSums(scaled_expr)
# # hist(sum_expr)
# #
# # object <- merged_object
# # pool <- NULL
# # pool <- pool %||% rownames(x = object)
# # assay.data <-  GetAssayData(object, slot = "data")
# # nbin <- 30
# # ctrl <- 100
# # features <- list(geneset)
# # cluster.length = length(features)
# #
# # # bin data
# # data.avg <- Matrix::rowMeans(x = assay.data[pool, ])
# # data.avg <- data.avg[order(data.avg)]
# # data.cut <- cut_number(x = data.avg + rnorm(n = length(data.avg))/1e30, n = nbin, labels = FALSE, right = FALSE)
# # names(x = data.cut) <- names(x = data.avg)
# # ctrl.use <- vector(mode = "list", length = cluster.length)
# #
# # # set control genes
# # for (i in 1:cluster.length) {
# #   features.use <- features[[i]]
# #   for (j in 1:length(x = features.use)) {
# #     ctrl.use[[i]] <- c(
# #       ctrl.use[[i]],
# #       names(x = sample(
# #         x = data.cut[which(x = data.cut == data.cut[features.use[j]])],
# #         size = ctrl,
# #         replace = FALSE
# #       ))
# #     )
# #   }
# # }
# #
# # ctrl.use <- lapply(X = ctrl.use, FUN = unique)
# # ctrl.scores <- matrix(
# #   data = numeric(length = 1L),
# #   nrow = length(x = ctrl.use),
# #   ncol = ncol(x = object)
# # )
# #
# # for (i in 1:length(ctrl.use)) {
# #   features.use <- ctrl.use[[i]]
# #   ctrl_data <- assay.data[features.use, ]
# #   ctrl_data <- t(apply(ctrl_data, 1, scale, center = FALSE))
# #   dim(ctrl_data)
# #   ctrl.scores[i, ] <- Matrix::colMeans(x = ctrl_data)
# # }
# #
# # features.scores <- matrix(
# #   data = numeric(length = 1L),
# #   nrow = cluster.length,
# #   ncol = ncol(x = object)
# # )
# # for (i in 1:cluster.length) {
# #   features.use <- features[[i]]
# #   feature_data <- assay.data[features.use, ]
# #   feature_data <- t(apply(feature_data, 1, scale, center = FALSE))
# #   dim(feature_data)
# #   features.scores[i, ] <- Matrix::colMeans(x = feature_data)
# # }
# #
# # features.scores.use <- features.scores - ctrl.scores
# #
# # rownames(x = features.scores.use) <- paste0(name, 1:cluster.length)
# # features.scores.use <- as.data.frame(x = t(x = features.scores.use))
# # rownames(x = features.scores.use) <- colnames(x = object)
# # hist(df$score)
#
# # test_scores <- AddScore(object, list(geneset), nbin = 30, ctrl = 100)
# # hist(test_scores$Cluster1)
# # hist(features.scores.use)
#
# markers_ls[[1]]
#
# df$disease_classification
# m1_scores <- ggplot(marker_scores %>% filter(geneset == "M_1"), aes(x = disease_classification, y = score, fill = batch)) +
#   #geom_hline(yintercept = threshold, linetype = "dotted", color = "gray60", size = 1) +
#   geom_jitter(alpha = 0.1, shape = 21, color = "transparent", position = position_jitterdodge(jitter.width = 0.1, dodge.width = 1)) +
#   geom_violin(alpha = 0.8, width = 1, position = position_dodge(width = 1), draw_quantiles = c(0.5)) +
#   colorspace::scale_fill_discrete_qualitative("Dark 3") +
#   #annotate(geom = "text", x = 1, y = 2, label = glue("{gated_pops$plot[1]}%"), size = 6) +
#   #annotate(geom = "text", x = 2, y = 2, label = glue("{gated_pops$plot[2]}%"), size = 6) +
#   #annotate(geom = "text", x = 3, y = 2, label = glue("{gated_pops$plot[3]}%"), size = 6) +
#   ylim(c(-0.5, 2)) +
#   guides(fill = FALSE) +
#   #labs(x = "Condition (split by batch)", y = "Score", title = "M1 Score", subtitle = glue("Mixed effects (Healthy vs. Severe), {scales::scientific(lmtest::lrtest(mmodel, 1)[2, 5], digits = 2)}")) +
#   GeneralTheme(14)
# m1_scores
#
# SavePlot(plot = m1_scores, filename = "plots/bost2020_m1scores.png", base_height = 6, base_asp = 1.2)
#
#
# scores <- lapply(raw_scores, function(x) {
#   x$source <- gsub("(_disease)|(_health)", "", x$set)
#   x$category <- str_match(pattern = "(_disease|_health)", string = x$set)[, 2]
#   x$category <- gsub("_", "", x$category)
#   x$category <- factor(x$category, levels = c("health", "disease"))
#   x$response <- as.numeric(x$category)-1
#   x <- x %>% mutate(compare = ifelse(length(unique(x$category)) == 2, TRUE, FALSE))
#   x$set_cluster <- paste0(x$set, "_", x$merged_leiden)
#   x$umis <- (x$nCount_RNA - mean(x$nCount_RNA))/sd(x$nCount_RNA)
#   return(x)
# })
#
# scores <- bind_rows(scores)
# sum_scores <- scores %>% group_by(category, program) %>% summarize(score = mean(score), percent = median(percent))
# sum_scores$program <- factor(sum_scores$program, levels = levels_ordered)
# sumscores_plot <- ggplot(sum_scores, aes(x = category, y = program, fill = score, size = percent)) +
#   geom_point(shape = 21, color = "black") +
#   colorspace::scale_fill_continuous_sequential("Blues", rev = 1, name = "Score") +
#   #scale_x_discrete(labels = c("Health", "Disease")) +
#   scale_y_discrete(labels = gsub("M_", "", levels_ordered)) +
#   labs(x = "Classification", y = "Marker gene set") +
#   scale_size_continuous(range = c(2, 6), name = "Average \n% detected") +
#   PointTheme(14) +
#   theme(
#     legend.title = element_text(face = "bold")
#   )
# sumscores_plot
# SavePlot("plots/bost2020_marker_scores_dot.png", base_asp = 0.8, plot = sumscores_plot)
#
# LoadColors()
# levels_ordered <- unique(scores$program)[order(as.numeric(gsub("M_", "", unique(scores$program))))]
# scores$program <- factor(scores$program, levels = levels_ordered)
# scoresdist_plot <- ggplot(scores, aes(x = score, y = program, fill = category)) +
#   ggridges::geom_density_ridges(color = "black", alpha = 0.6) +
#   scale_fill_manual(values = c("#0072B2", "#D55E00"), labels = c("Health", "Disease"), name = "Condition") +
#   scale_y_discrete(labels = gsub("M_", "", levels_ordered)) +
#   labs(x = "Score", y = "Marker gene set", title = "Consensus marker gene sets' scores", subtitle = "COVID-19 BAL samples") +
#   ggridges::theme_ridges(center_axis_labels = TRUE) +
#   theme(
#     axis.title = element_text(face = "bold"),
#     legend.title = element_text(face = "bold", size = 12),
#     legend.text = element_text(size = 12),
#     legend.position = "right",
#     legend.justification = "top")
# scoresdist_plot
# SavePlot("plots/bost2020_marker_scores_dist.png", base_asp = 1, plot = scoresdist_plot)
#
#   ggplot(scores, aes(x = nCount_RNA, fill = category)) +
#   geom_density() +
#   scale_x_log10() +
#   DensityTheme(18)
#
# umap_coords <- object_ls[[2]]@reductions$harmony_umap@cell.embeddings
# scores <- merge(scores, umap_coords, by.x = "cell_barcode", by.y = 0, all.x = TRUE)
#
# featureplot <- ggplot(scores %>% filter(category == "disease", program == "M_4") %>% arrange(score), aes(x = hUMAP_1.x, y = hUMAP_2.x, fill = score)) +
#   geom_point(shape = 21, size = 1, stroke = 0.1) +
#   colorspace::scale_fill_continuous_diverging("Blue-Red 3", limits = c(-1, 2.5), name = "Score") +
#   labs(x = "UMAP1", y = "UMAP2") +
#   PointTheme(14) + theme(
#     legend.title = element_text(face = "bold"),
#     panel.grid = element_blank(),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     axis.ticks = element_blank(),
#     axis.text = element_blank())
# SavePlot("plots/bost2020_disease_m8_umap.png", base_asp = 1.2, base_height = 4, plot = featureplot)
# featureplot <- ggplot(scores %>% filter(category == "health", program == "M_8") %>% arrange(score),
#   aes(x = hUMAP_1.y, y = hUMAP_2.y, fill = score)) +
#   geom_point(shape = 21, size = 1, stroke = 0.1) +
#   colorspace::scale_fill_continuous_diverging("Blue-Red 3", limits = c(-1, 2.5), name = "Score") +
#   labs(x = "UMAP1", y = "UMAP2") +
#   PointTheme(14) + theme(
#     legend.title = element_text(face = "bold"),
#     panel.grid = element_blank(),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     axis.ticks = element_blank(),
#     axis.text = element_blank())
# SavePlot("plots/bost2020_health_m8_umap.png", base_asp = 1.2, base_height = 4, plot = featureplot)
#
# geneset = "M_17"
# num_top = 10
# tpscores <- bind_rows(scores) %>%
#   filter(program == geneset)
# passing_clusters <- tpscores %>%
#   group_by(set, merged_leiden) %>%
#   summarize(mean = median(score), set_cluster = unique(set_cluster)) %>%
#   top_n(num_top, mean) %>%
#   select(set_cluster)
# tpscores <- tpscores %>% mutate(plot = ifelse(set_cluster %in% passing_clusters$set_cluster, TRUE, FALSE))
#
# dt <- diptest::dip.test(tpscores$score)
# dt$p.value
# if (dt$p.value >= 0.5) {
#   median_score <- median(tpscores$score)
#   mad_score <- mad(tpscores$score)
#   threshold <- median_score + 1.5*mad_score
# } else {
#   mm <- mixtools::normalmixEM(tpscores$score, k = 2)
#   plot(mm, 2)
#   min_mm <- which.min(mm$mu)
#   threshold <- mm$mu[min_mm] + mm$sigma[min_mm]*2
# }
#
# sum_scores <- tpscores %>%
#   group_by(set, program) %>%
#   summarize(n = n(), mean_score = mean(score),
#     n_above = sum(score > threshold)/n)
#
# a <- ggplot(tpscores, aes(x = score, y = category, fill = category)) +
#   ggridges::geom_density_ridges(alpha = 0.3, color = "black") +
#   geom_vline(xintercept = threshold, size = 1, linetype = "dotted") +
#   scale_fill_manual(values = c("Gray60", "#D55E00")) +
#   guides(fill = FALSE) +
#   scale_y_discrete(expand = expansion(mult = c(0.25, 1), add = c(0, 0.5)),
#     labels = snakecase::to_title_case(rev(unique(as.character(tpscores$category))))) +
#   annotate(geom = "text", label = glue("{round(sum_scores$n_above[1], 2)}% above threshold"),
#     y = 2.5, x = 1.25, fontface = "bold") +
#   annotate(geom = "text", label = glue("{round(sum_scores$n_above[2], 2)}% above threshold"),
#     y = 1.5, x = 1.25, fontface = "bold") +
#   labs(x = "Score", y = "", subtitle = "% cells above global threshold (1.5 MADs)") +
#   DensityTheme(14)
# a
#
# b <- ggplot(tpscores %>% filter(plot == TRUE),
#   aes(x = nCount_RNA, y = score , color = category)) +
#   geom_point(alpha = 0.5, size = 0.2) +
#   geom_density2d(size = 0.5) +
#   scale_x_continuous(trans = scales::log10_trans(),
#     breaks = scales::trans_breaks("log10", function(x) 10^x),
#     labels = scales::trans_format("log10", scales::math_format(10^.x))) +
#   scale_color_manual(values = c("Gray60", "#D55E00")) +
#   labs(x = "log10(Library size)", y = "Score", subtitle = "Gene set score vs. library size (top 2 clusters)") +
#   ViolinTheme() + NoLegend()
# b <- ggExtra::ggMarginal(b, tpscores %>% filter(plot == TRUE),
#   x = score, y = nCount_RNA, type = "density", groupFill = TRUE)
# b
#
# numeric_flvls <- gsub(glue("{dataset}_(health|disease)_"), "", unique(tpscores$set_cluster))
# lvls <- unique(tpscores$set_cluster)[order(as.numeric(numeric_flvls))]
# lvls <- c(lvls[grep("disease", lvls)], lvls[grep("health", lvls)])
# tpscores$set_cluster <- factor(tpscores$set_cluster, levels = rev(lvls))
# lvl_labels <- gsub(glue("{dataset}_(health|disease)_"), "", rev(lvls))
# #lvl_labels <- snakecase::to_title_case(gsub("_", "", lvl_labels))
#
# c <- ggplot(tpscores,
#   aes(x = score, y = set_cluster, fill = category)) +
#   #geom_jitter(shape = 21, alpha = 0.25, width = 0.1, size = 0.2) +
#   geom_boxplot(alpha = 0.75, width = 0.75, position = position_dodge(width = 0.5), outlier.color = "Gray60", outlier.alpha = 0.25) +
#   scale_fill_manual(values = c("Gray60", "#d55e00")) +
#   scale_y_discrete(labels = lvl_labels) +
#   scale_x_continuous(limits = c(min(tpscores$score),
#     ceiling(max(tpscores$score)))) +
#   # annotate(geom = "text", label = glue("P = {format.pval(unique(kpis$plogit))}, \u03b2 = {round(unique(kpis$blogit), 2)}"), x = "health",
#   #   y = ceiling(max(scores_to_plot$score)), hjust = 0) +
#   labs(y = "Cluster", x = "Score", subtitle = c("Cluster-specific expression of gene set")) +
#   ViolinTheme(base_size = 14) + NoLegend()
# c
#
# title <- cowplot::ggdraw() +
#   cowplot::draw_label(glue("Marker {gsub('M_', '', geneset)}, {gsub('_', ' ', dataset)}"), fontface = 'bold', x = 0, hjust = 0) +
#   theme(plot.margin = margin(0, 0, 0, 24))
#
# grid <- cowplot::plot_grid(title, cowplot::plot_grid(a, b, c, nrow = 1), nrow = 2, rel_heights = c(0.1, 1))
# SavePlot("plots/{dataset}_{geneset}.png", grid, base_asp = 3, base_height = 6)
#
#
# lisi_res <- lisi::compute_lisi(object_ls[[1]]@reductions$harmony@cell.embeddings, object_ls[[1]]@meta.data,
#   c("merged_leiden"))
# lisi_res
#
# object_ls[[1]]@reductions$harmony_umap@cell.embeddings %>%
#   cbind(lisi_res) %>%
#   dplyr::sample_frac(1L, FALSE) %>%
#   tidyr::gather(key, lisi_value, merged_leiden) %>%
#   ggplot(aes(hUMAP_1, hUMAP_2, color = lisi_value)) + geom_point(shape = 21) +
#   facet_wrap(~key)
#
# all.equal(rownames(object_ls[[1]]@meta.data), rownames(lisi_res))
# object_ls[[1]]$merged_leiden_LISI <- lisi_res$merged_leiden
# VlnPlot(object_ls[[1]], features = "merged_leiden_LISI", group.by = "merged_leiden", pt.size = 0.01) + scale_y_continuous(limits = c(0, 10))
# object_ls[[1]][[]] %>% group_by(merged_leiden) %>% summarize(median = median(merged_leiden_LISI))
#
