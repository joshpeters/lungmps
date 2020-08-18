# ---
# Description: Extract MNP subsets from annotated objects
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
suffix <- "annotated.rds"
dir <- "data/objects"

load_files <- list.files(path = dir, pattern = suffix, full.names = TRUE)
ids <- str_match(string = load_files, pattern = "(data/objects/)(.*_(health|disease))_annotated.rds")[, 3, drop = TRUE]
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

# [3] Annotate objects ----
# for each object, load, process, save
walk(.x = ids, .f = ~ {
  usethis::ui_info("Loading {.x} ({round(file.info(load_files[.x])$size/1e6, 2)} mb)")
  obj <- readRDS(file = load_files[.x])
  usethis::ui_done(("\nCompleted loading {.x}"))
  Clean()

  ui_info("Classifying {.x}")
  nonmp <- as.character(type_df$score_name[type_df$mnp == 0])
  nonmp <- nonmp[nonmp %in% colnames(obj@meta.data)]
  mp <- as.character(type_df$score_name[type_df$mnp == 1])
  assertthat::assert_that(all(nonmp %in% colnames(obj@meta.data)))
  mnp_diff_score <- rowSums(obj@meta.data[, mp]) - rowSums(obj@meta.data[, nonmp])
  assertthat::assert_that(all.equal(names(mnp_diff_score), rownames(obj@meta.data)))
  obj$mp_diff_score <- mnp_diff_score

  mp_scores <- obj[[]] %>% select(base_clustering, mp_diff_score)
  mm <- mixtools::normalmixEM(mp_scores$mp_diff_score, k = 2)
  x_data <- with(mm, seq(min(x), max(x), len = 1000))
  pars <- with(mm, data.frame(comp = colnames(posterior), mu, sigma, lambda))
  em_df <- data.frame(x = rep(x_data, each = nrow(pars)), pars)
  em_df$y <- em_df$lambda * dnorm(em_df$x, mean = em_df$mu, sd = em_df$sigma)
  min <- which.min(mm$mu)

  threshold <- mm$mu[min] + 2.5*mm$sigma[min]
  ui_info("Threshold: {threshold}")
  if (threshold >= 3) { threshold <- 3}

  mp_clusters <- as.character(mp_scores %>%
      group_by(base_clustering) %>%
      summarize(median_mp_score = median(mp_diff_score)) %>%
      filter(median_mp_score > threshold & median_mp_score > 0) %>%
      pull(base_clustering))
  nonmp_clusters <- setdiff(unique(obj$base_clustering), mp_clusters)

  cc_props <- as.matrix(table(obj$base_clustering, obj$cc_diff_label))
  cc_props <- data.frame(cluster = rownames(cc_props), cycling = cc_props[, 1, drop = TRUE], noncycling = cc_props[, 2, drop = TRUE])
  cc_props <- cc_props %>% group_by(cluster) %>% mutate(props = cycling/(cycling+noncycling))
  cycling_cluster <- as.character(cc_props %>% ungroup %>% filter(props >= 0.25) %>% pull(cluster))

  nonmp_clusters <- setdiff(nonmp_clusters, cycling_cluster)
  mp_clusters <- setdiff(mp_clusters, cycling_cluster)

  counts <- obj@meta.data %>%
    select(base_clustering, jmp_type) %>%
    filter(base_clustering %in% mp_clusters) %>%
    group_by(base_clustering, jmp_type) %>%
    summarize(n = n()) %>%
    mutate(prop = n/sum(n)) %>%
    filter(prop > 0.1 & prop < 0.50)
  contam_clusters <- unique(as.character(counts$base_clustering[!counts$jmp_type %in% c("monocyte", "dc", "macrophage")]))
  counts <- obj@meta.data %>%
    select(base_clustering, jmp_type) %>%
    filter(base_clustering %in% mp_clusters) %>%
    group_by(base_clustering, jmp_type) %>%
    summarize(n = n()) %>%
    mutate(prop = n/sum(n)) %>%
    filter(prop >= 0.50)
  highcontam_clusters <- unique(as.character(counts$base_clustering[!counts$jmp_type %in% c("monocyte", "dc", "macrophage")]))
  mp_clusters <- setdiff(mp_clusters, contam_clusters)
  mp_clusters <- setdiff(mp_clusters, highcontam_clusters)

  obj$base_labels <- NA
  obj$base_labels[obj$base_clustering %in% mp_clusters] <- "mnp"
  obj$base_labels[obj$base_clustering %in% contam_clusters] <- "prelim_mnp"
  obj$base_labels[obj$base_clustering %in% c(nonmp_clusters, highcontam_clusters)] <- "non_mnp"
  obj$base_labels[obj$base_clustering %in% cycling_cluster] <- "cycling"

  cc_scores <- obj[[]] %>% select(cell_barcode, base_clustering, mp_diff_score) %>% filter(base_clustering == cycling_cluster)
  mp_cycling <- cc_scores %>%
    group_by(base_clustering) %>%
    filter(mp_diff_score > threshold) %>% pull(cell_barcode)

  length(mp_cycling)
  obj$mp <- FALSE
  obj$mp[obj$base_clustering %in% c(mp_clusters, contam_clusters)] <- TRUE
  obj$mp[mp_cycling] <- TRUE
  ui_info("Found {sum(obj$mp)} MNP cells")

  a <- DimPlotPlus(obj, reduction = "harmony_umap", group.by = "base_clustering",
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

  b <- ggplot(data.frame(x = mm$x), aes(x, y = ..density..)) +
    geom_histogram(fill = NA, color = "black") +
    geom_polygon(data = em_df, mapping = aes(x = x, y = y, fill = comp), color = "black", alpha = 0.5) +
    geom_vline(xintercept = threshold, color = "gray40", linetype = "dashed") +
    colorspace::scale_fill_discrete_sequential(palette = "Light Grays", name = "Component Means") +
    guides(fill = FALSE) +
    scale_x_continuous(expand = c(0, 0.05)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = "MP Score", y = "Density") +
    theme_classic(base_size = 14) +
    ggeasy::easy_all_text_color(color = "black") +
    theme(
      axis.title.y = element_text(angle = 90, vjust = 0.5, margin = margin(0, 12, 0, 0), face = "bold"),
      axis.title.x = element_text(margin = margin(12, 0, 0, 0), face = "bold"),
      axis.text.y = element_text(margin = margin(0, 4, 0, 0)),
      axis.line = element_blank(),
      panel.border = element_rect(size = 1, color = "black", fill = "transparent"),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 10),
      legend.justification = "left",
      legend.position = "bottom",
      legend.key.size = unit(1, "line"),
      plot.margin = margin(6, 6, 6, 6, unit = "mm"))

  c <- ggplot(mp_scores, aes(x = mp_diff_score, y = base_clustering, fill = base_clustering)) +
    ggridges::geom_density_ridges(alpha = 0.8) +
    geom_vline(xintercept = threshold, color = "gray40", linetype = "dashed") +
    colorspace::scale_fill_discrete_qualitative("Dark 2", name = "Cluster") +
    scale_x_continuous(expand = c(0, 0.05)) +
    scale_y_discrete(expand = c(0, 1)) +
    labs(x = "MP Score", y = "Cluster") +
    theme_classic(base_size = 14) +
    ggeasy::easy_all_text_color(color = "black") +
    theme(
      axis.title.y = element_text(angle = 90, vjust = 0.5, margin = margin(0, 12, 0, 0), face = "bold"),
      axis.title.x = element_text(margin = margin(12, 0, 0, 0), face = "bold"),
      axis.text.y = element_text(margin = margin(0, 4, 0, 0)),
      axis.line = element_blank(),
      panel.border = element_rect(size = 1, color = "black", fill = "transparent"),
      legend.position = "none",
      plot.margin = margin(6, 6, 6, 6, unit = "mm"))

  d <- DimPlotPlus(obj, reduction = "harmony_umap", group.by = "mp",
    label = TRUE, label.size = 6, repel = TRUE) +
    colorspace::scale_color_discrete_qualitative("Dark 3", name = "Filtered as MP") +
    labs(x = "UMAP 1", y = "UMAP 2",
      caption = glue::glue("MP assignment ({sum(obj$mp)} MNP cells assigned, {round(sum(obj$mp)/ncol(obj)*100, 2)}% of total)")) +
    theme_void(base_size = 14) +
    theme(
      axis.title.x = element_text(hjust = 0, face = "plain", color = "gray60"),
      axis.title.y = element_text(hjust = 0, angle = 90, face = "plain", color = "gray60"),
      plot.caption = element_text(size = 10, hjust = 1),
      legend.position = "none",
      plot.margin = margin(6, 6, 6, 6, unit = "mm"),
    )

  full <- cowplot::plot_grid(a, b, c, d, ncol = 2, labels = "auto")
  cowplot::save_plot(filename = glue::glue("plots/{.x}_mnp_classification_r1.png"),
    plot = full, base_height = 10, base_asp = 1.2, device = "png")

  # generate data for blacklist genes
  auc_threshold <- 0.9

  mnp_clusters <- unique(obj$base_clustering[obj$base_labels == "mnp"])
  nonmnp_clusters <- unique(obj$base_clustering[obj$base_labels == "non_mnp"])

  bl_genes <- obj@misc$global_markers %>%
    group_by(group) %>%
    mutate(set = .x) %>%
    filter(logFC > log(1.5) & auc > auc_threshold & padj < 1E-5) %>%
    filter(group %in% nonmnp_clusters) %>%
    arrange(desc(logFC))
  saveRDS(bl_genes, file = glue("data/blacklist_genes/{.x}_blacklist_genes.rds"))

  # save object
  ui_todo("Saving object...")
  ui_info("\nSaving {.x} ({round(as.numeric(object.size(obj)/1E6), 2)} mb)")
  saveRDS(obj,
    file = glue::glue("data/objects/{.x}_annotated.rds"))
  Clean()
  ui_done(("\nCompleted {.x}"))

  # generate MNP object
  ui_done("Subsetting {sum(obj$mp)} MNPs from {.x}")
  obj <- obj[, WhichCells(obj, expression = mp == TRUE)]
  obj <- DietSeurat(object = obj)
  Clean()

  # base process
  ui_todo("Base processing...")
  obj <- BaseProcess(object = obj, nfeatures = 3000)
  obj$batch_LISI <- lisi::compute_lisi(obj@reductions$pca@cell.embeddings, obj@meta.data, "batch")$batch
  obj <- HarmonyCorrection(object = obj, batch = "batch")
  obj <- CalculatePCHeuristics(obj, reduction = "harmony", store.name = "harmony_heuristics", force_tp = TRUE)

  # umap and clustering
  ui_todo("UMAP and clustering...")
  dims <- 20
  method <- "harmony"
  obj <- RunUMAP(object = obj, reduction = method,
    dims = 1:dims, seed.use = 1, reduction.name = glue("{method}_umap"),
    reduction.key = glue("{substr(method, 1, 1)}UMAP_"))
  Clean()
  obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:dims, k.param = 20, prune.SNN = 1/15,
    compute.SNN = TRUE, graph.name = c("harmony_nn", "harmony_snn"), verbose = TRUE)
  obj <- BaseCluster(object = obj, dims = 20, verbose = TRUE, graph.name = "harmony_snn",
    res.start = 1E-5, res.end = 1, num.res = 30, num.clusters.lower = 5, num.clusters.upper = 30, mod.similarity = 0.995)
  snn <- igraph::graph_from_adjacency_matrix(obj@graphs$harmony_snn, mode = "directed", weighted = TRUE)
  walktrap_cluster <- igraph::cluster_walktrap(graph = snn, steps = 4)
  obj$base_walktrap <- walktrap_cluster$membership
  Clean()
  minimum.cells <- 10
  cluster <- "base_leiden"
  print(table(obj$base_leiden))
  obj <- SetIdent(obj, value = cluster)
  obj <- BuildClusterTree(obj, reorder = TRUE, reorder.numeric = TRUE)
  obj$base_clustering <- obj@active.ident

  obj <- SetIdent(object = obj, value = "base_walktrap")
  walktrap_markers <- presto::wilcoxauc(obj, assay = "data", seurat_assay = "RNA")
  obj@misc$walktrap_markers <- walktrap_markers

  obj <- SetIdent(object = obj, value = "base_leiden")
  leiden_markers <- presto::wilcoxauc(obj, assay = "data", seurat_assay = "RNA")
  obj@misc$leiden_markers <- leiden_markers

  ui_todo("Saving object...")
  ui_info("\nSaving {.x} ({round(as.numeric(object.size(obj)/1E6), 2)} mb)")
  saveRDS(obj,
    file = glue::glue("data/objects/{.x}_mnp.rds"))
  Clean()
  ui_done(("\nCompleted {.x}"))

})

