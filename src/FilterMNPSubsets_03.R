# ---
# Description: Filter MNP subsets
# Author: Josh Peters
# ---

# [1] Basic setup ----

# load libraries
rm(list = ls())
library(future)
library(pbapply)
library(Seurat)
source("src/Functions.R")
source("src/Utilities.R")
source("src/Plotting.R")
stamp <- PrepEnv()

# load metadata and prepare for pipeline
suffix <- "mnp.rds"
dir <- "data/objects"

load_files <- list.files(path = dir, pattern = suffix, full.names = TRUE)
ids <- str_match(string = load_files, pattern = glue("(data/objects/)(.*_(health|disease))_{suffix}"))[, 3, drop = TRUE]
ui_todo("Loading in {length(load_files)} files")
names(load_files) <- ids

# [2] Load supporting data ----
degs <- readRDS(file = "data/interim/tra20_degs.rds")
deg_names_touse <- readRDS(file = "data/interim/tra20_degs_names.rds")
type_df <- readRDS(file = "data/interim/tra20_type_df.rds")
score_labels <- readRDS(file = "data/interim/tra20_score_labels.rds")
score_names_touse <- readRDS(file = "data/interim/tra20_score_names.rds")
is_expr <- readRDS(file = "data/interim/is_expr.rds")
lm22_expr <- readRDS(file = "data/interim/lm22_expr.rds")

# [3] Extract markers ----
object_ls <- LoadObjectList("data/objects/", file_pattern = "_mnp.rds", verbose = TRUE)
names(object_ls) <- ids
ids <- SelfName(ids)
Clean()

# run once ---
all_genes <- lapply(object_ls, rownames)
saveRDS(all_genes, file = "data/annotated_all_genes.rds")
all_genes <- unique(unname(unlist(all_genes)))
ig <- all_genes[grep("^IG[H/L/V/K].*$", all_genes)]
tcrb <- all_genes[grep("^TRBV.*$", all_genes)]
tcrg <- all_genes[grep("^TRGV.*$", all_genes)]
tcrc <- all_genes[grep("^TRGC.*$", all_genes)]
tcrj <- all_genes[grep("^TRGJ.*$", all_genes)]
tcra <- all_genes[grep("^TRA[C|V].*$", all_genes)]
tcrs <- c(tcrb, tcrg, tcrc, tcrj, tcra)

auc_threshold <- 0.95
nonMP_markers <- pbapply::pblapply(X = ids, FUN = function(x) {
  nonMP_clusters <- as.character(object_ls[[x]][[]] %>%
      group_by(base_clustering) %>%
      summarize(num_mnp = sum(mp)) %>%
      filter(num_mnp == 0) %>%
      pull(base_clustering))
  nonMP_markers <- object_ls[[x]]@misc$global_markers %>%
    group_by(group) %>%
    mutate(set = x) %>%
    filter(logFC > log(1.5) & auc > auc_threshold & padj < 1E-5) %>%
    filter(group %in% nonMP_clusters) %>%
    arrange(desc(logFC))
  return(nonMP_markers)
})

nonMP_markers_df <- bind_rows(nonMP_markers)
nonMP_markers_df <- nonMP_markers_df %>%
  group_by(feature) %>%
  mutate(n = n(), median_auc = median(auc)) %>%
  filter(n >= 3) %>%
  distinct(feature, .keep_all = TRUE)
saveRDS(nonMP_markers_df, glue("data/nonMPmarkers_auc{auc_threshold}.rds"))

gene_markers <- pblapply(list.files(path = "data/blacklist_genes/", full.names = TRUE), FUN = readRDS)
blacklist_genes <- bind_rows(gene_markers)
blacklist_genes <- blacklist_genes %>%
  group_by(feature) %>%
  mutate(n = n()) %>%
  filter(n >= 2) %>%
  distinct(feature, .keep_all = TRUE)
saveRDS(blacklist_genes, glue("data/blacklist_genes_auc0.90.rds"))

blacklist_features <- readRDS("data/blacklist_genes_auc0.90.rds")
nonmp_markers <- readRDS("data/nonMPmarkers_auc0.95.rds")
nonmp_markers <- c(nonmp_markers$feature, "NAMPT", "CSF3R", "CD79A", "CD79B", "MS4A1", "IGKC", ig, tcrs)
nonmp_markers <- nonmp_markers[!(nonmp_markers %in% c("A2M", "GZMB", "RGS1"))]
saveRDS(nonmp_markers, "data/nonmp_markers_total.rds")

# [4] Identify poor quality cells ----
total_mp_scores <- pblapply(ids, function(x) {
  return(object_ls[[x]][[]] %>% select(base_leiden, mp_diff_score, cell_barcode))
})
total_mp_scores <- Reduce(rbind, total_mp_scores)
total_threshold <- median(total_mp_scores$mp_diff_score) - 2.5*mad(total_mp_scores$mp_diff_score); total_threshold
hist(total_mp_scores$mp_diff_score)

nonmp_markers <- readRDS("data/nonmp_markers_total.rds")
blacklist_hits <- pblapply(X = ids, function(x) {

  ui_info("Printing figure for {x}")
  object_ls[[x]]$base_leiden <- as.factor(object_ls[[x]]$base_leiden)
  object_ls[[x]]$base_leiden <- factor(object_ls[[x]]$base_leiden, levels = levels(object_ls[[x]]$base_leiden)[order(as.numeric(levels(object_ls[[x]]$base_leiden)))] )

  mp_scores <- object_ls[[x]][[]] %>% select(base_leiden, mp_diff_score)

  plot_markers <- object_ls[[x]]@misc$leiden_markers %>%
    dplyr::filter(padj < 1E-10 & logFC > log(1.25) & auc > 0.5) %>%
    group_by(group) %>%
    arrange(as.numeric(group)) %>%
    top_n(10, logFC)
  blacklist_hits <- plot_markers[(plot_markers$feature %in% c(nonmp_markers)), ]
  blacklist_hits <- names(table(blacklist_hits$group))[table(blacklist_hits$group) >= 2]

  nonmp_clusters <- as.character(mp_scores %>%
      group_by(base_leiden) %>%
      summarize(median_mp_score = median(mp_diff_score)) %>%
      filter(median_mp_score < 0) %>%
      pull(base_leiden))
  nonmp_cells <- colnames(object_ls[[x]])[as.character(object_ls[[x]]$base_leiden) %in% unique(c(blacklist_hits, nonmp_clusters))]
  nonmp_cells <- unique(c(nonmp_cells, rownames(mp_scores)[mp_scores$mp_diff_score < total_threshold]))
  ui_info("Identified {length(nonmp_cells)} contaminating cell profiles")

  genes <- unique(plot_markers %>%
      group_by(group) %>%
      top_n(5, logFC) %>% pull(feature))

  cluster_umap <- DimPlotPlus(object_ls[[x]], reduction = "harmony_umap", group.by = "base_leiden",
    label = TRUE, label.size = 6, repel = TRUE) +
    colorspace::scale_color_discrete_qualitative("Dark 3", name = "Cluster") +
    labs(x = "UMAP 1", y = "UMAP 2") +
    theme_void(base_size = 14) +
    theme(
      axis.title.x = element_text(hjust = 0, face = "plain", color = "gray60"),
      axis.title.y = element_text(hjust = 0, angle = 90, face = "plain", color = "gray60"),
      panel.background = element_rect(fill = "transparent", color = "transparent"), # bg of the panel
      plot.background = element_rect(fill = "transparent", color = "transparent"), # bg of the plot
      panel.grid.major = element_blank(), # get rid of major grid
      panel.grid.minor = element_blank(), # get rid of minor grid
      legend.background = element_rect(fill = "transparent"), # get rid of legend bg
      legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
      legend.position = "none",
      plot.margin = margin(6, 6, 6, 6, unit = "mm")
    )

  contam_umap <- DimPlotPlus(object_ls[[x]], reduction = "harmony_umap", group.by = "base_leiden", pt.size = 1,
    cells.highlight = nonmp_cells) +
    scale_color_manual(values = c("gray90", "black"),
      name = "Potential\ncontamination", labels = c("MP", paste(blacklist_hits, collapse = ", "))) +
    labs(x = "UMAP 1", y = "UMAP 2") +
    theme_void(base_size = 14) +
    theme(
      axis.title.x = element_text(hjust = 0, face = "plain", color = "gray60"),
      axis.title.y = element_text(hjust = 0, angle = 90, face = "plain", color = "gray60"),
      panel.background = element_rect(fill = "transparent", color = "transparent"), # bg of the panel
      plot.background = element_rect(fill = "transparent", color = "transparent"), # bg of the plot
      panel.grid.major = element_blank(), # get rid of major grid
      panel.grid.minor = element_blank(), # get rid of minor grid
      legend.background = element_rect(fill = "transparent"), # get rid of legend bg
      legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
      legend.position = "none",
      plot.margin = margin(6, 6, 6, 6, unit = "mm")
    )

  scores_plot <- DotPlot_ol(object_ls[[x]], features = score_names_touse,
    cols = rev(colorspace::sequential_hcl(2, "Blues")),
    dot.min = 0.05, dot.scale = 4, scale.by = "radius", group.by = "base_leiden") +
    scale_x_discrete(labels = rev(score_labels)) +
    labs(x = "Expression Score", y = "Cluster") +
    theme(
      plot.subtitle = element_text(face = "bold"),
      legend.position = "none",
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 6),
      axis.text.y = element_text(size = 10),
      axis.title.x = element_text(face = "bold", margin = margin(t = 2, unit = "pt")),
      axis.title.y = element_text(face = "bold", margin = margin(r = 12, unit = "pt")),
      legend.title = element_text(vjust = 1, hjust = 0),
      panel.grid.major = element_line(color = "Gray90", size = 0.5),
      plot.margin = unit(rep(18, 4), "pt"),
      axis.line = element_blank(),
      panel.border = element_rect(size = 0.5, color = "black")
    )

  markers_plot <- DotPlot_ol(object_ls[[x]], features = genes,
    cols = rev(colorspace::sequential_hcl(2, "Blues")),
    dot.min = 0.05, dot.scale = 6, scale.by = "radius", group.by = "base_leiden") +
    labs(x = "Gene markers (top 5 ~ logFC)", y = "Cluster") +
    theme(
      plot.subtitle = element_text(face = "bold"),
      legend.position = "none",
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 6),
      axis.text.y = element_text(size = 10),
      axis.title.x = element_text(face = "bold", margin = margin(t = 2, unit = "pt")),
      axis.title.y = element_text(face = "bold", margin = margin(r = 12, unit = "pt")),
      legend.title = element_text(vjust = 1, hjust = 0),
      panel.grid.major = element_line(color = "Gray90", size = 0.5),
      plot.margin = unit(rep(18, 4), "pt"),
      axis.line = element_blank(),
      panel.border = element_rect(size = 0.5, color = "black")
    )

  top <- cowplot::plot_grid(cluster_umap, contam_umap, align = "hv")
  full <- cowplot::plot_grid(top, scores_plot, markers_plot, ncol = 1)
  cowplot::save_plot(filename = glue("plots/{x}_mnp_filtering.png"),
    plot = full, base_height = 12, base_asp = round(8.5/11, 2))

  return(nonmp_cells)
})

ids[12]
VlnPlot(object_ls[[12]], c("FCGR3A", "NAMPT", "CD14", "CSF3R"), group.by = "base_leiden")

# manually identify problematic clusters using metrics
remove <- list(
  c(6, 13, 14), # braga_health
  c(),
  c(),
  c(11),
  c(10), # lambrechts_health
  c(11),
  c(8, 16, 19), # laughney_health
  c(3, 12),
  c(11, 12, 13, 14), # mayr_disease
  c(8, 9),
  c(3, 16, 20, 21, 22), # morse disease
  c(13),
  c(11, 13, 17, 18, 19), # raredon_health
  c(16, 17, 19),
  c(11), # reyfman_health
  c(19, 20, 23, 24),
  c(19),
  c(3, 4, 14, 17) # zilionis_disease
)
names(remove) <- ids

# print results of manual ID

blacklist_hits <- pblapply(X = ids, function(x) {

  ui_info("Printing figure for {x}")

  mp_scores <- object_ls[[x]][[]] %>% select(base_leiden, mp_diff_score)
  nonmp_clusters <- as.character(remove[[x]])
  nonmp_cells <- colnames(object_ls[[x]])[as.character(object_ls[[x]]$base_leiden) %in% nonmp_clusters]
  nonmp_cells <- unique(c(nonmp_cells, rownames(mp_scores)[mp_scores$mp_diff_score < total_threshold]))
  ui_info("Identified {length(nonmp_cells)} contaminating cell profiles")

  cluster_umap <- DimPlotPlus(object_ls[[x]], reduction = "harmony_umap", group.by = "base_leiden",
    label = TRUE, label.size = 6, repel = TRUE) +
    colorspace::scale_color_discrete_qualitative("Dark 3", name = "Cluster") +
    labs(x = "UMAP 1", y = "UMAP 2") +
    theme_void(base_size = 14) +
    theme(
      axis.title.x = element_text(hjust = 0, face = "plain", color = "gray60"),
      axis.title.y = element_text(hjust = 0, angle = 90, face = "plain", color = "gray60"),
      panel.background = element_rect(fill = "transparent", color = "transparent"), # bg of the panel
      plot.background = element_rect(fill = "transparent", color = "transparent"), # bg of the plot
      panel.grid.major = element_blank(), # get rid of major grid
      panel.grid.minor = element_blank(), # get rid of minor grid
      legend.background = element_rect(fill = "transparent"), # get rid of legend bg
      legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
      legend.position = "none",
      plot.margin = margin(6, 6, 6, 6, unit = "mm")
    )

  contam_umap <- DimPlotPlus(object_ls[[x]], reduction = "harmony_umap", group.by = "base_leiden", pt.size = 1,
    cells.highlight = nonmp_cells) +
    scale_color_manual(values = c("gray90", "black"),
      name = "Potential\ncontamination", labels = c("MP", paste(blacklist_hits, collapse = ", "))) +
    labs(x = "UMAP 1", y = "UMAP 2") +
    theme_void(base_size = 14) +
    theme(
      axis.title.x = element_text(hjust = 0, face = "plain", color = "gray60"),
      axis.title.y = element_text(hjust = 0, angle = 90, face = "plain", color = "gray60"),
      panel.background = element_rect(fill = "transparent", color = "transparent"), # bg of the panel
      plot.background = element_rect(fill = "transparent", color = "transparent"), # bg of the plot
      panel.grid.major = element_blank(), # get rid of major grid
      panel.grid.minor = element_blank(), # get rid of minor grid
      legend.background = element_rect(fill = "transparent"), # get rid of legend bg
      legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
      legend.position = "none",
      plot.margin = margin(6, 6, 6, 6, unit = "mm")
    )

  title <- cowplot::ggdraw() +
    cowplot::draw_label(
      glue("{length(nonmp_cells)} contaminating cell profiles \nCluster(s) {paste(nonmp_clusters, collapse = ', ')}"),
      fontface = 'bold',
      x = 0,
      hjust = 0
    ) +
    theme(
      plot.margin = margin(0, 0, 0, 10)
    )

  removal <- cowplot::plot_grid(cluster_umap, contam_umap, align = "hv")
  removal <- cowplot::plot_grid(title, removal, ncol = 1, rel_heights = c(0.1, 1))
  cowplot::save_plot(filename = glue("plots/{x}_mp_filtered.png"),
    plot = removal, base_height = 6, base_asp = 2)
  return(nonmp_cells)
})

saveRDS(blacklist_hits, glue("data/interim/removed_cells_round1.rds"))
first_round <- object_ls

# [5] Reprocess filtered cells ----

object_ls <- pblapply(X = ids, function(x) {
  ui_done("Subsetting {length(setdiff(Cells(object_ls[[x]]), blacklist_hits[[x]]))} MPs from {x}")
  mp_object <- object_ls[[x]][, setdiff(Cells(object_ls[[x]]), blacklist_hits[[x]])]
  mp_object <- DietSeurat(object = mp_object)
  Clean()
  return(mp_object)
})

object_ls <- pblapply(X = ids, function(x) {
  object_ls[[x]] <- BaseProcess(object = object_ls[[x]], nfeatures = 4000)
  object_ls[[x]] <- HarmonyCorrection(object = object_ls[[x]], batch = "batch")
  Clean()
  return(object_ls[[x]])
})

object_ls <- pblapply(X = ids, function(x) {
  object_ls[[x]] <- CalculatePCHeuristics(object_ls[[x]], reduction = "harmony", store.name = "harmony_heuristics")
  dims <- ceiling(unique(object_ls[[x]]@misc$harmony_heuristics$tp))
  ui_done("\nUsing {ui_field(dims)} dims for post-Harmony functions for dataset {x}")
  object_ls[[x]]$postharmony_dims_used <- dims
  Clean()
  return(object_ls[[x]])
})
ui_info("Average number of post-Harmony dimensions: {median(unlist(lapply(object_ls, function(x) return(unique(x$postharmony_dims_used)))))}")

object_ls <- pblapply(X = ids, function(x) {
  dims <- 20
  method <- "harmony"
  object_ls[[x]] <- RunUMAP(object = object_ls[[x]], reduction = method,
    dims = 1:dims, seed.use = 1, reduction.name = glue("{method}_umap"),
    reduction.key = glue("{substr(method, 1, 1)}UMAP_"))
  Clean()
  return(object_ls[[x]])
})

object_ls <- pblapply(X = ids, function(x) {
  dims <- 20
  object_ls[[x]] <- FindNeighbors(object_ls[[x]], reduction = "harmony", dims = 1:dims, k.param = 20, prune.SNN = 1/15,
    compute.SNN = TRUE, graph.name = c("harmony_nn", "harmony_snn"), verbose = TRUE)
  object_ls[[x]] <- BaseCluster(object = object_ls[[x]], dims = 20, verbose = TRUE,
    res.start = 1E-5, res.end = 1, num.res = 30, num.clusters.lower = 5, num.clusters.upper = 30, mod.similarity = 0.995)
  snn <- igraph::graph_from_adjacency_matrix(object_ls[[x]]@graphs$harmony_snn, mode = "directed", weighted = TRUE)
  walktrap_cluster <- igraph::cluster_walktrap(graph = snn, steps = 4)
  object_ls[[x]]$base_walktrap <- walktrap_cluster$membership
  Clean()
  return(object_ls[[x]])
})

minimum.cells <- 20
cluster <- "base_leiden"
object_ls <- pblapply(X = ids, function(x) {
  ui_info("Processing {x}...")
  print(table(object_ls[[x]]$base_leiden))
  object_ls[[x]] <- SetIdent(object_ls[[x]], value = cluster)
  object_ls[[x]] <- BuildClusterTree(object_ls[[x]], reorder = TRUE, reorder.numeric = TRUE)
  object_ls[[x]]$base_clustering <- object_ls[[x]]@active.ident
  return(object_ls[[x]])
})

object_ls <- pblapply(X = ids, function(x) {

  object_ls[[x]] <- SetIdent(object = object_ls[[x]], value = "base_walktrap")
  walktrap_markers <- presto::wilcoxauc(object_ls[[x]], assay = "data", seurat_assay = "RNA")
  object_ls[[x]]@misc$walktrap_markers <- walktrap_markers

  object_ls[[x]] <- SetIdent(object = object_ls[[x]], value = "base_leiden")
  leiden_markers <- presto::wilcoxauc(object_ls[[x]], assay = "data", seurat_assay = "RNA")
  object_ls[[x]]@misc$leiden_markers <- leiden_markers

  Clean()
  return(object_ls[[x]])
})

total_mp_scores <- pblapply(ids, function(x) {
  return(object_ls[[x]][[]] %>% select(base_leiden, mp_diff_score, cell_barcode))
})
total_mp_scores <- Reduce(rbind, total_mp_scores)
total_threshold <- median(total_mp_scores$mp_diff_score) - 2.5*mad(total_mp_scores$mp_diff_score); total_threshold
hist(total_mp_scores$mp_diff_score)

blacklist_hits <- pblapply(X = ids, function(x) {

  ui_info("Printing figure for {x}")
  object_ls[[x]]$base_leiden <- as.factor(object_ls[[x]]$base_leiden)
  object_ls[[x]]$base_leiden <- factor(object_ls[[x]]$base_leiden, levels = levels(object_ls[[x]]$base_leiden)[order(as.numeric(levels(object_ls[[x]]$base_leiden)))])
  mp_scores <- object_ls[[x]][[]] %>% select(base_leiden, mp_diff_score)

  plot_markers <- object_ls[[x]]@misc$leiden_markers %>%
    dplyr::filter(padj < 1E-10 & logFC > log(1.25) & auc > 0.5) %>%
    group_by(group) %>%
    arrange(as.numeric(group)) %>%
    top_n(10, logFC)
  blacklist_hits <- plot_markers[(plot_markers$feature %in% c(nonmp_markers)), ]
  blacklist_hits <- names(table(blacklist_hits$group))[table(blacklist_hits$group) >= 2]

  nonmp_clusters <- as.character(mp_scores %>%
      group_by(base_leiden) %>%
      summarize(median_mp_score = median(mp_diff_score)) %>%
      filter(median_mp_score < 0) %>%
      pull(base_leiden))
  nonmp_cells <- colnames(object_ls[[x]])[as.character(object_ls[[x]]$base_leiden) %in% unique(c(blacklist_hits, nonmp_clusters))]
  nonmp_cells <- unique(c(nonmp_cells, rownames(mp_scores)[mp_scores$mp_diff_score < total_threshold]))
  ui_info("Identified {length(nonmp_cells)} contaminating cell profiles")

  genes <- unique(plot_markers %>%
      group_by(group) %>%
      top_n(5, logFC) %>% pull(feature))

  cluster_umap <- DimPlotPlus(object_ls[[x]], reduction = "harmony_umap", group.by = "base_leiden",
    label = TRUE, label.size = 6, repel = TRUE) +
    colorspace::scale_color_discrete_qualitative("Dark 3", name = "Cluster") +
    labs(x = "UMAP 1", y = "UMAP 2") +
    theme_void(base_size = 14) +
    theme(
      axis.title.x = element_text(hjust = 0, face = "plain", color = "gray60"),
      axis.title.y = element_text(hjust = 0, angle = 90, face = "plain", color = "gray60"),
      panel.background = element_rect(fill = "transparent", color = "transparent"), # bg of the panel
      plot.background = element_rect(fill = "transparent", color = "transparent"), # bg of the plot
      panel.grid.major = element_blank(), # get rid of major grid
      panel.grid.minor = element_blank(), # get rid of minor grid
      legend.background = element_rect(fill = "transparent"), # get rid of legend bg
      legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
      legend.position = "none",
      plot.margin = margin(6, 6, 6, 6, unit = "mm")
    )

  contam_umap <- DimPlotPlus(object_ls[[x]], reduction = "harmony_umap", group.by = "base_leiden", pt.size = 1,
    cells.highlight = nonmp_cells) +
    scale_color_manual(values = c("gray90", "black"),
      name = "Potential\ncontamination", labels = c("MP", paste(blacklist_hits, collapse = ", "))) +
    labs(x = "UMAP 1", y = "UMAP 2") +
    theme_void(base_size = 14) +
    theme(
      axis.title.x = element_text(hjust = 0, face = "plain", color = "gray60"),
      axis.title.y = element_text(hjust = 0, angle = 90, face = "plain", color = "gray60"),
      panel.background = element_rect(fill = "transparent", color = "transparent"), # bg of the panel
      plot.background = element_rect(fill = "transparent", color = "transparent"), # bg of the plot
      panel.grid.major = element_blank(), # get rid of major grid
      panel.grid.minor = element_blank(), # get rid of minor grid
      legend.background = element_rect(fill = "transparent"), # get rid of legend bg
      legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
      legend.position = "none",
      plot.margin = margin(6, 6, 6, 6, unit = "mm")
    )

  scores_plot <- DotPlot_ol(object_ls[[x]], features = score_names_touse,
    cols = rev(colorspace::sequential_hcl(2, "Blues")),
    dot.min = 0.05, dot.scale = 4, scale.by = "radius", group.by = "base_leiden") +
    scale_x_discrete(labels = rev(score_labels)) +
    labs(x = "Expression Score", y = "Cluster") +
    theme(
      plot.subtitle = element_text(face = "bold"),
      legend.position = "none",
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 6),
      axis.text.y = element_text(size = 10),
      axis.title.x = element_text(face = "bold", margin = margin(t = 2, unit = "pt")),
      axis.title.y = element_text(face = "bold", margin = margin(r = 12, unit = "pt")),
      legend.title = element_text(vjust = 1, hjust = 0),
      panel.grid.major = element_line(color = "Gray90", size = 0.5),
      plot.margin = unit(rep(18, 4), "pt"),
      axis.line = element_blank(),
      panel.border = element_rect(size = 0.5, color = "black")
    )

  markers_plot <- DotPlot_ol(object_ls[[x]], features = genes,
    cols = rev(colorspace::sequential_hcl(2, "Blues")),
    dot.min = 0.05, dot.scale = 6, scale.by = "radius", group.by = "base_leiden") +
    labs(x = "Gene markers (top 5 ~ logFC)", y = "Cluster") +
    theme(
      plot.subtitle = element_text(face = "bold"),
      legend.position = "none",
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 6),
      axis.text.y = element_text(size = 10),
      axis.title.x = element_text(face = "bold", margin = margin(t = 2, unit = "pt")),
      axis.title.y = element_text(face = "bold", margin = margin(r = 12, unit = "pt")),
      legend.title = element_text(vjust = 1, hjust = 0),
      panel.grid.major = element_line(color = "Gray90", size = 0.5),
      plot.margin = unit(rep(18, 4), "pt"),
      axis.line = element_blank(),
      panel.border = element_rect(size = 0.5, color = "black")
    )

  top <- cowplot::plot_grid(cluster_umap, contam_umap, align = "hv")
  full <- cowplot::plot_grid(top, scores_plot, markers_plot, ncol = 1)
  cowplot::save_plot(filename = glue("plots/{x}_mp_confirmation.png"),
    plot = full, base_height = 12, base_asp = round(8.5/11, 2))
  return(nonmp_cells)
})

# confirm any remaining questionable cells - NOT REQUIRED
# working_set <- ids[13]
# plot_markers <- object_ls[[working_set]]@misc$leiden_markers %>%
#   dplyr::filter(padj < 1E-10 & logFC > log(1.25) & auc > 0.5) %>%
#   group_by(group) %>%
#   arrange(as.numeric(group)) %>%
#   top_n(10, logFC)
# VlnPlot(object_ls[[working_set]], c("FCGR3A", "NAMPT", "CD14"), group.by = "base_leiden")
# View(plot_markers)
# nonmp_cells <- colnames(object_ls[[working_set]])[as.character(object_ls[[working_set]]$base_leiden) %in% c("15")]
# object_ls[[working_set]] <-  object_ls[[working_set]][, setdiff(Cells(object_ls[[working_set]]), nonmp_cells)]

pblapply(X = ids, function(x) {
  ui_info("\nSaving {x} ({round(as.numeric(object.size(object_ls[[x]])/1E6), 2)} mb)")
  saveRDS(object_ls[[x]],
    file = glue::glue("data/objects/{x}_mnp_final.rds"))
  Clean()
  ui_done(("\nCompleted {x}"))
})
