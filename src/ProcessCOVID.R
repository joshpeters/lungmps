# ---
# Description: Process test COVID object
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

# [2] Load data ----

full_obj <- readRDS("data/objects/Liao_2020_base_object.rds")
.x = "Liao_2020"

degs <- readRDS(file = "data/interim/tra20_degs.rds")
deg_names_touse <- readRDS(file = "data/interim/tra20_degs_names.rds")
type_df <- readRDS(file = "data/interim/tra20_type_df.rds")
score_labels <- readRDS(file = "data/interim/tra20_score_labels.rds")
score_names_touse <- readRDS(file = "data/interim/tra20_score_names.rds")
is_expr <- readRDS(file = "data/interim/is_expr.rds")
lm22_expr <- readRDS(file = "data/interim/lm22_expr.rds")

# [3] Initial processing ----

full_obj <- SetIdent(full_obj, value = "batch")
"disease_classification" %in% colnames(full_obj[[]])

if (length(unique(full_obj$disease_classification)) == 1) {
  split_obj <- list(full_obj)
  names(split_obj) <- unique(full_obj$disease_classification)
} else {
  full_obj$disease_classification <- as.factor(full_obj$disease_classification)
  split_obj <- Seurat::SplitObject(full_obj, split.by = "disease_classification")
}

ui_info("Processing {names(split_obj)}")
for (i in names(split_obj)) {
  ui_info("Processing {i} object")
  print(split_obj[[i]])
  obj <- split_obj[[i]]

  if (any(table(obj$batch) < 500)) {
    ui_info("Removing low frequency batches")
    passing_batches <- names(table(obj$batch))[table(obj$batch) >= 500]
    obj <- obj[, obj$batch %in% passing_batches]
  } else {
    ui_info("No low frequency batches")
  }

  # filter object
  obj <- FilterCells(
    obj,
    nmads = 5,
    variables = c("nCount_RNA", "nFeature_RNA", "percent_mit", "percent_rib"),
    batch = "batch",
    percent_mito_cutoff = 20,
    percent_ribo_cutoff = 50,
    percent_hb_cutoff = 5,
    percent_hsp_cutoff = 5,
    nCount_RNA_cutoff = 200,
    nFeature_RNA_cutoff = 100,
    remove = FALSE,
    qc_column = "passing_qc")

  assertthat::assert_that(all.equal(colnames(obj), rownames(obj[[]])))
  keep_cells <- WhichCells(obj, expression = passing_qc == TRUE)
  ui_info("Number of cells pre-filtering: {ncol(obj)}")
  obj<- obj[, keep_cells]
  ui_info("Number of cells post-filtering: {ncol(obj)}")
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
    res.start = 1E-5, res.end = 1, num.res = 30, num.clusters.lower = 5, num.clusters.upper = 30, mod.similarity = 1)
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

  # downstream scoring and annotation
  ui_todo("Downstream scoring...")
  obj <- ScoreLineages(object = obj, grouping_var = "base_clustering", genes = degs)
  dotplot <- DotPlot_ol(obj, features = score_names_touse,
    cols = rev(colorspace::sequential_hcl(2, "Blues")),
    dot.min = 0.05, dot.scale = 4, scale.by = "radius") +
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
  cowplot::save_plot(plot = dotplot, filename = glue("plots/global_tra20_scores_{.x}_{i}.png"), base_height = 6, base_asp = 1.5)

  max_scores <- apply(obj[[]][, score_names_touse, drop = TRUE], 1, which.max)
  max_scores <- ConvertNamedVecToDF(max_scores)
  max_scores <- merge(max_scores, deg_names_touse, by.x = "value", by.y = "index", sort = FALSE)
  max_scores <- max_scores %>% select(name.x, name.y, jmp_class, jmp_subclass, jmp_type)
  colnames(max_scores) <- c("cell_barcode", "score_name", "jmp_class", "jmp_subclass", "jmp_type")
  rownames(max_scores) <- max_scores$cell_barcode
  max_scores <- max_scores[, -1]
  obj <- AddMetaData(obj, max_scores)

  cluster_barcodes <- obj[[]] %>% select(cell_barcode, base_clustering)
  cluster_assignments <- obj[[]] %>% select(base_clustering, jmp_class)
  cluster_assignments <- cluster_assignments %>%
    group_by(base_clustering, jmp_class) %>% summarize(n = n())

  cols <- ggthemes::colorblind_pal()(8)[c(2, 3, 4, 7)]
  names(cols) <- levels(cluster_assignments$jmp_class)
  count_plot <- ggplot(cluster_assignments, aes(x = base_clustering, y = n, fill = jmp_class)) +
    geom_col(position = "stack", width = 1, color = "black") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, max(cluster_assignments$n) + 200)) +
    scale_fill_manual(values = cols, name = "Classes", labels = c("Endothelial", "Epithelial", "Immune", "Stromal")) +
    labs(x = "Cluster", y = "Number of cells") +
    theme_classic(base_size = 14) +
    ggeasy::easy_all_text_color("black") +
    theme(
      axis.title.y = element_text(angle = 90, vjust = 0.5, margin = margin(0, 12, 0, 0), face = "bold"),
      axis.title.x = element_text(margin = margin(12, 0, 0, 0), face = "bold"),
      axis.text.y = element_text(margin = margin(0, 4, 0, 0)),
      axis.line = element_blank(),
      plot.margin = margin(4, 4, 4, 4, unit = "mm"),
      panel.border = element_rect(size = 1, color = "black", fill = "transparent"),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 10),
      legend.position = "bottom",
      legend.justification = "left",
      legend.key.size = unit(1, "line")
    ) + coord_flip()
  cowplot::save_plot(plot = count_plot, filename = glue("plots/global_cluster_class_counts_{.x}_{i}.png"), base_height = 6, base_asp = 1)

  cluster_assignments <- cluster_assignments %>%
    top_n(n = 1, wt = n)
  cluster_assignments <- merge(cluster_barcodes, cluster_assignments, .by = "base_clustering", sort = FALSE)
  rownames(cluster_assignments ) <- cluster_assignments$cell_barcode
  cluster_assignments  <- cluster_assignments [, -1]
  obj <- AddMetaData(obj, cluster_assignments[, "jmp_class", drop = FALSE])

  subcluster_barcodes <- obj[[]] %>% select(cell_barcode, base_clustering)
  subcluster_assignments <- obj[[]] %>% select(base_clustering, jmp_subclass)
  subcluster_assignments <- subcluster_assignments %>%
    group_by(base_clustering, jmp_subclass) %>% summarize(n = n()) %>%
    top_n(n = 1, wt = n)
  subcluster_assignments <- merge(cluster_barcodes, subcluster_assignments, .by = "base_clustering", sort = FALSE)
  rownames(subcluster_assignments) <- subcluster_assignments$cell_barcode
  subcluster_assignments  <- subcluster_assignments [, -1]
  obj <- AddMetaData(obj, subcluster_assignments[, "jmp_subclass", drop = FALSE])

  if ("celltype" %in% colnames(obj[[]])) {
    plots <- DimPlot(obj, reduction = "harmony_umap",
      label = TRUE, group.by = c("base_clustering", "celltype", "jmp_class", "jmp_subclass"), combine = FALSE, label.size = 4, repel = TRUE)
    plots <- lapply(plots, function(x) {
      x + colorspace::scale_color_discrete_qualitative("Dark 2") +
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
          legend.position = "none"
        )
    })
    plotgrid <- cowplot::plot_grid(plotlist = plots, align = c("v"), ncol = 4)
    cowplot::save_plot(filename = glue("plots/umap_cellclass_{.x}_{i}.png"), plot = plotgrid,
      base_height = 5, base_asp = 4)
  } else {
    plots <- DimPlot(obj, reduction = "harmony_umap",
      label = TRUE, group.by = c("base_clustering", "jmp_class", "jmp_subclass"), combine = FALSE, label.size = 4, repel = TRUE)
    plots <- lapply(plots, function(x) {
      x + colorspace::scale_color_discrete_qualitative("Dark 2") +
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
          legend.position = "none"
        )
    })
    plotgrid <- cowplot::plot_grid(plotlist = plots, align = c("v"), ncol = 3)
    cowplot::save_plot(filename = glue("plots/umap_cellclass_{.x}_{i}.png"), plot = plotgrid,
      base_height = 5, base_asp = 3)
  }

  Clean()
  obj <- LabelCyclingCells(object = obj, grouping_var = "base_clustering")

  ui_todo("Classification...")
  lls <- LLClassifier(obj, is_expr)
  obj@misc$is_llpp <- lls$pp
  calls <- ConvertNamedVecToDF(lls$call)
  calls <- calls %>% select(value)
  obj <- AddMetaData(obj, calls, col.name = "IS")
  lls <- LLClassifier(obj, lm22_expr)
  obj@misc$is_llpp <- lls$pp
  calls <- ConvertNamedVecToDF(lls$call)
  calls <- calls %>% select(value)
  obj <- AddMetaData(obj, calls, col.name = "LM22")
  Clean()

  obj <- SetIdent(object = obj, value = "base_clustering")
  global_markers <- presto::wilcoxauc(obj, assay = "data", seurat_assay = "RNA")
  obj@misc$global_markers <- global_markers

  ui_todo("Saving object...")
  ui_info("\nSaving {.x} ({round(as.numeric(object.size(obj)/1E6), 2)} mb)")
  saveRDS(obj,
    file = glue::glue("data/objects/{.x}_{i}_annotated.rds"))
  Clean()
  ui_done(("\nCompleted {.x}"))

}

# [4] Extract subsets ----
suffix <- "Liao_2020_.*annotated.rds"
dir <- "data/objects"
load_files <- list.files(path = dir, pattern = suffix, full.names = TRUE)
ids <- str_match(string = load_files, pattern = "(data/objects/)(Liao_2020_(health|disease))_annotated.rds")[, 3, drop = TRUE]
ui_todo("Loading in {length(load_files)} files")
names(load_files) <- ids
ids <- SelfName(ids)

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

# [5] Filter MNP subsets ----

suffix <- "Liao_2020_.*mnp.rds"
dir <- "data/objects"
load_files <- list.files(path = dir, pattern = suffix, full.names = TRUE)
ids <- str_match(string = load_files, pattern = "(data/objects/)(.*(health|disease))_mnp.rds")[, 3, drop = TRUE]
ui_todo("Loading in {length(load_files)} files")
names(load_files) <- ids
ids <- SelfName(ids)
ids

object_ls <- future.apply::future_lapply(ids, function(x) {
  usethis::ui_info("\nLoading {x} ({round(file.info(load_files[x])$size/1e6, 2)} mb)")
  obj <- readRDS(file = load_files[x])
  usethis::ui_done(("\nCompleted {x}"))
  Clean()
  return(obj)
})
names(object_ls) <- ids

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

remove <- list(c(21, 25, 26, 27, 28, 29), c(4, 5, 16, 18, 21))
names(remove) <- ids

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

saveRDS(blacklist_hits, glue("data/interim/liao2020_removed_cells_round1.rds"))
first_round <- object_ls

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
  object_ls[[x]] <- BaseCluster(object = object_ls[[x]], dims = 20, verbose = TRUE, graph.name = "harmony_snn",
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

# [6] Define clusters ----
ids
pblapply(X = ids, function(x) {

  object <- object_ls[[x]]

  # perform clustering
  object <- FindNeighbors(object = object, reduction = "harmony", graph.name = c("harmony_nn", "harmony_snn"),
    dims = 1:20, verbose= TRUE, compute.SNN = TRUE, force.recalc = TRUE)
  results <- IterativeCluster(object, graph_name = "harmony_snn")
  plotgrid <- cowplot::plot_grid(plotlist = results$plots, align = "hv", labels = "auto")
  cowplot::save_plot(filename = glue("plots/{x}_mnp_cluster_merges_RNAsnn.png"),
    plot = plotgrid, base_height = 12, base_asp = round(16/9, 2))
  object <- results$object
  object[["merged_leiden"]] <- Seurat::Idents(object)
  markers <- IdentifyQualityMarkers(object, tlogFC = log(1.1))
  object@misc$merged_clusters <- markers

  res_df <- object@misc$leiden_results$second_results
  num_current <- length(unique(Idents(object)))
  res <- res_df$resolution_parameter[which.min(abs(res_df$cluster_count - num_current))]
  g <- igraph::graph_from_adjacency_matrix(adjmatrix = object@graphs$harmony_snn, mode = "directed", weighted = TRUE)
  final_result <- leidenbase::leiden_find_partition(g, partition_type = "CPMVertexPartition", seed = 1,
    resolution_parameter = res, num_iter = 30, verbose = TRUE)
  out_result <- list(membership = final_result[["membership"]],
    modularity = final_result[["modularity"]])
  names(out_result$membership) <- colnames(object@graphs$harmony_snn)
  if (any(table(out_result$membership) < 5)) {
    ids <- GroupSingletons(ids = out_result$membership, SNN = object@graphs$harmony_snn,
      min.size = 5, group.singletons = TRUE, verbose = TRUE)
    object <- AddMetaData(object, ids, col.name = "equiv_leiden")
  } else {
    object <- AddMetaData(object, out_result$membership, col.name = "equiv_leiden")
  }
  object <- SetIdent(object, value = "equiv_leiden")
  object@misc$leiden_mergedequiv_markers <- IdentifyQualityMarkers(object, tlogFC = log(1.1))

  a <- PlotAnnotatedClusters(object, group.by = "merged_leiden")
  b <- PlotAnnotatedClusters(object, group.by = "equiv_leiden")
  title <- cowplot::ggdraw() +
    cowplot::draw_label(
      glue::glue("Clustering results for {x} (ARI = {round(mclust::adjustedRandIndex(as.character(object$merged_leiden), as.character(object$equiv_leiden)), 3)})"),
      fontface = "bold", x = 0, hjust = 0) +
    theme(plot.margin = margin(0, 0, 0, 12, unit = "pt"))
  plotgrid <- cowplot::plot_grid(title, cowplot::plot_grid(a, b, labels = "auto", align = "hv"), ncol = 1, rel_heights = c(0.05, 1))
  cowplot::save_plot(filename = glue("plots/{x}_mnp_cluster_results.png"),
    plot = plotgrid, base_height = 8, base_asp = 2)
  Clean()

  # resave object
  ui_todo("Saving object...")
  ui_info("\nSaving {x} ({round(as.numeric(object.size(object)/1E6), 2)} mb)")
  saveRDS(object, file = glue::glue("data/objects/{x}_mnp_final.rds"))
  Clean()
  ui_done(("\nCompleted {x}"))

})

# [7] Save summary plots ----

suffix <- "Liao_2020_.*mnp_final.rds"
dir <- "data/objects"
load_files <- list.files(path = dir, pattern = suffix, full.names = TRUE)
load_files <- load_files[c(1,2)]
ids <- str_match(string = load_files, pattern = "(data/objects/)(.*(health|disease))_mnp_final.rds")[, 3, drop = TRUE]
ui_todo("Loading in {length(load_files)} files")
names(load_files) <- ids
ids <- SelfName(ids)

object_ls <- future.apply::future_lapply(ids, function(x) {
  usethis::ui_info("\nLoading {x} ({round(file.info(load_files[x])$size/1e6, 2)} mb)")
  obj <- readRDS(file = load_files[x])
  usethis::ui_done(("\nCompleted {x}"))
  Clean()
  return(obj)
})
names(object_ls) <- ids

i = 2
#FeaturePlot(object_ls[[i] ], features = c("percent_ribo"))
cluster_umap <- DimPlotPlus(object_ls[[i]], reduction = "harmony_umap", group.by = "merged_leiden",
  label = TRUE, label.size = 6, repel = TRUE) +
  colorspace::scale_color_discrete_qualitative("Dark 3", name = "Cluster") +
  labs(x = "UMAP 1", y = "UMAP 2") +
  GeneralTheme(base_size = 14) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_text(hjust = 0, face = "bold", color = "black", margin = margin(10,0,0,0)),
    axis.title.y = element_text(hjust = 0, angle = 90, face = "bold", color = "black", margin = margin(0,10,0,0)),
    panel.background = element_rect(fill = "transparent", color = "black", size = 1), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = "transparent"), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
    legend.position = "none",
    plot.margin = margin(6, 6, 6, 6, unit = "mm")
  )
cluster_umap
#SavePlot("plots/bost2020_umap_mergedclusters.png", plot = cluster_umap)

library(ggdendro)
library(dendextend)
PlotClusterTree(object_ls[[i]])
tree <- as.dendrogram(object_ls[[i]]@tools$BuildClusterTree)
ddata <- dendro_data(tree, type = "rectangle")
p <- ggplot(segment(ddata)) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_text(data = ddata$labels,
    aes(x = x, y = y, label = label), size = 4, vjust = 0.5, hjust = 1) +
  ylim(c(-2, round(max(ddata$segments$yend))+2)) +
  coord_flip() +
  theme_dendro() +
  theme(plot.margin = unit(c(0.01, 0.01, 0.01, 0.01), "cm"))
#umap_tree <- cowplot::plot_grid(cluster_umap, p, rel_widths = c(1, 0.5))
#umap_tree
#SavePlot("plots/bost2020_disease_umaptree.png", base_asp = 1.5, plot = umap_tree)

genes <- unique(object_ls[[i]]@misc$merged_clusters %>%
    group_by(group) %>%
    top_n(2, logFC) %>% select(feature, group))
markers_plot <- DotPlot_ol(object_ls[[i]], features = genes$feature,
  cols = rev(colorspace::sequential_hcl(2, "Grays")),
  dot.min = 0.05, dot.scale = 4, scale.by = "radius", group.by = "merged_leiden") +
  labs(x = "", y = "Cluster") +
  theme(
    plot.subtitle = element_text(face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 8),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(face = "bold", margin = margin(t = 2, unit = "pt")),
    axis.title.y = element_text(face = "bold", margin = margin(r = 12, unit = "pt")),
    legend.title = element_text(vjust = 1, hjust = 0),
    panel.grid.major = element_line(color = "Gray90", size = 0.5),
    plot.margin = unit(rep(18, 4), "pt"),
    axis.line = element_blank(),
    panel.border = element_rect(size = 0.5, color = "black")
  )
#markers_plot
#SavePlot("plots/bost2020_disease_markers.png", base_asp = 2, base_height= 4, plot = markers_plot)

total <- cowplot::plot_grid(cluster_umap, markers_plot, NULL, p, rel_widths = c(1, 2, 0, 0.75), nrow = 1, align = "h")
total_disease <- total
SavePlot("plots/bost2020d_summary.pdf", base_asp = 6.5/2, base_height = 4, plot = total_disease)
total_health <- total
SavePlot("plots/bost2020h_summary.pdf", base_asp = 6.5/2, base_height = 4, plot = total_health)

# [8] Score the dataset ----

group <- "merged_leiden"
object_ls <- lapply(ids, function(x) {
  object_ls[[x]]$set <- x
  return(object_ls[[x]])
})

obj <- merge(object_ls[[1]], object_ls[[2]], add.cell.ids = c("d", "h"), merge.data = TRUE)

# SCORE PROGRAMS
genesets <- programs_ls
names <- programs
group <- "merged_leiden"

ui_todo("Preparing gene sets...")
genesets <- map(names, ~ {
  genesets[[.x]] <- genesets[[.x]][genesets[[.x]] %in% rownames(obj)]
  return(genesets[[.x]])
}, obj = obj)

 percents <- map_dfr(.x = names, ~ {
  count_subset <- GetAssayData(obj, slot = "counts")[genesets[[.x]], ]
  percents <- Matrix::colSums(count_subset > 0)
  percents <- percents/nrow(count_subset)
  percents <- data.frame(percents = percents, geneset = .x, cell_barcode = colnames(obj))
  return(percents)
}, obj = obj)
ui_done("Gene sets prepared")
assertthat::assert_that(all(lapply(genesets, length) > 0))

ui_todo("Scoring object...")
genesets # check genesets before running
obj <- AddScore(object = obj, features = genesets, nbin = 30, ctrl = 30, name = "Geneset", seed = 1)
assertthat::assert_that(length(grep("Geneset[[:digit:]]", colnames(obj@meta.data))) == length(genesets))
colnames(obj@meta.data)[grep("Geneset[[:digit:]]", colnames(obj@meta.data))] <- names
ui_done("Object scored")

df <- obj@meta.data[, c("nCount_RNA", "disease_classification", "set", group, "batch", names)]
df <- df %>% rownames_to_column("cell_barcode")
df <- df %>% gather(key = "program", value = "score", names)
df <- cbind(df, percents[, c("percents", "geneset")])
colnames(df)[5] <- "group"
colnames(df)[10] <- "percents_geneset"
colnames(df)[7] <- "geneset"
program_scores <- df

# SCORE MARKERS
genesets <- markers_ls
names <- markers
group <- "merged_leiden"

ui_todo("Preparing gene sets...")
genesets <- map(names, ~ {
  genesets[[.x]] <- genesets[[.x]][genesets[[.x]] %in% rownames(obj)]
  return(genesets[[.x]])
}, obj = obj)

percents <- map_dfr(.x = names, ~ {
  count_subset <- GetAssayData(obj, slot = "counts")[genesets[[.x]], ]
  percents <- Matrix::colSums(count_subset > 0)
  percents <- percents/nrow(count_subset)
  percents <- data.frame(percents = percents, geneset = .x, cell_barcode = colnames(obj))
  return(percents)
}, obj = obj)
ui_done("Gene sets prepared")
assertthat::assert_that(all(lapply(genesets, length) > 0))

ui_todo("Scoring object...")
genesets # check genesets before running
obj <- AddScore(object = obj, features = genesets, nbin = 30, ctrl = 30, name = "Geneset", seed = 1)
assertthat::assert_that(length(grep("Geneset[[:digit:]]", colnames(obj@meta.data))) == length(genesets))
colnames(obj@meta.data)[grep("Geneset[[:digit:]]", colnames(obj@meta.data))] <- names
ui_done("Object scored")

df <- obj@meta.data[, c("nCount_RNA", "disease_classification", "set", group, "batch", names)]
df <- df %>% rownames_to_column("cell_barcode")
df <- df %>% gather(key = "program", value = "score", names)
df <- cbind(df, percents[, c("percents", "geneset")])
colnames(df)[5] <- "group"
colnames(df)[10] <- "percents_geneset"
colnames(df)[7] <- "geneset"
marker_scores <- df

# COMBINE SCORES
all_scores <- rbind(marker_scores, program_scores)
all_scores$raw_cell_barcode <- gsub("^(d|h)_", "", all_scores$cell_barcode)
conditions <- obj@meta.data %>% select(cell_barcode, condition)
all_scores <- merge(all_scores, conditions, by.x = "raw_cell_barcode", by.y = "cell_barcode", all.x = TRUE)
head(all_scores)

sum_mscores <- all_scores %>%
  group_by(condition, geneset) %>%
  filter(geneset %in% markers) %>%
  summarize(score = mean(score), percent = median(percents))

# plot markers list first
mscore_levels <- sum_mscores %>% filter(condition == "severe") %>% arrange(desc(score)) %>% pull(geneset)
sum_mscores$geneset <- factor(sum_mscores$geneset, levels = mscore_levels)
reordered_marker_names <- markers_ls[mscore_levels]
#unique(scores$program)[order(as.numeric(gsub("M_", "", unique(scores$program))))]
#sum_scores$geneset <- factor(sum_scores$geneset, levels = levels_ordered)
mscores_dotplot <- ggplot(sum_mscores, aes(x = geneset, y = condition, fill = score, size = percent)) +
  geom_point(shape = 21, color = "black") +
  #scale_x_discrete(labels = paste(gsub("M_", "", levels(sum_mscores$geneset)), reordered_marker_names, sep = " ")) +
  scale_y_discrete(labels = rev(c("Severe", "Mild", "Control"))) +
  scale_size_continuous(range = c(2, 8), name = "Average \n% detected", breaks = c(0.25, 0.50, 0.75)) +
  colorspace::scale_fill_continuous_diverging(palette = "Blue-Red 3", name = "Average\nScore",
    limits = c(-0.5, 0.5), breaks = c(-0.25, 0, 0.25)) +
  labs(x = "Marker geneset", y = "Condition") +
  GeneralTheme(18) + theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.box = "bottom",
    legend.direction = "vertical",
    panel.grid.major = element_line(size = 0.5, color = "gray90"),
    panel.grid.minor = element_line(size = 0.5, color = "gray90")
  )
mscores_dotplot
SavePlot(filename = "plots/liao_marker_scores.pdf",
  plot = mscores_dotplot, base_asp = 1.75, base_height = 1.9*3)


sum_pscores <- all_scores %>%
  group_by(condition, geneset) %>%
  filter(geneset %in% programs) %>%
  summarize(score = mean(score), percent = median(percents))
pscore_levels <- sum_pscores %>% filter(condition == "severe") %>% arrange(desc(score)) %>% pull(geneset)
sum_pscores$geneset <- factor(sum_pscores$geneset, levels = pscore_levels)
#reordered_program_names <- program_gene_names[pscore_levels]
#unique(scores$program)[order(as.numeric(gsub("M_", "", unique(scores$program))))]
#sum_scores$geneset <- factor(sum_scores$geneset, levels = levels_ordered)
pscores_dotplot <- ggplot(sum_pscores, aes(x = geneset, y = condition, fill = score, size = percent)) +
  geom_point(shape = 21, color = "black") +
  #scale_x_discrete(labels = paste(gsub("P_", "", levels(sum_pscores$geneset)), reordered_program_names, sep = " ")) +
  scale_size_continuous(range = c(2, 8), name = "Average \n% detected", breaks = c(0.25, 0.50, 0.75)) +
  colorspace::scale_fill_continuous_diverging(palette = "Blue-Red 3", name = "Average\nScore",
    limits = c(-0.5, 0.5), breaks = c(-0.25, 0, 0.25)) +
  labs(x = "Program geneset", y = "Condition") +
  GeneralTheme(18) + theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "none",
    legend.box = "bottom",
    legend.direction = "vertical",
    panel.grid.major = element_line(size = 0.5, color = "gray90"),
    panel.grid.minor = element_line(size = 0.5, color = "gray90")
  )
pscores_dotplot

library(patchwork)
ab <- mscores_dotplot / pscores_dotplot + plot_layout(guides = "collect")
SavePlot(filename = "plots/liao_scores.pdf", plot = ab, base_asp = 1.75, base_height = 1.9*3)

# [9] Plot shift ----

all_scores_plot <- all_scores %>% filter(geneset %in% c("P1", "P8")) %>% select(-percents, -percents_geneset)
all_scores_plot <- all_scores_plot %>% spread(key = geneset, value = score)
shiftplot <- ggplot(all_scores_plot, aes(x = P1, y = P8, color = condition)) +
  geom_vline(xintercept = 0, color = "Gray80", linetype = "dotted", size = 1) +
  geom_hline(yintercept = 0, color = "Gray80", linetype = "dotted", size = 1) +
  scale_y_continuous(breaks = c(seq(-0.25, 0.5, 0.25))) +
  geom_density_2d(size = 1) +
  scale_color_manual(values = c("#009E73", "#0072B2", "#D55E00"), labels = c("Healthy", "Mild", "Severe"), name = "Condition") +
  labs(x = "P1 Score", y = "P8 Score") +
  GeneralTheme(base_size = 18)
shiftplot
SavePlot("plots/bost2020_program18_shift.png", base_asp = 1.2, plot = shiftplot)

# [10] Plot percent above ----

percent_above <- map_dfr(unique(all_scores$geneset), .f = ~ {
  test_scores <- all_scores %>% filter(geneset == .x)
  mm <- mixtools::normalmixEM(test_scores$score, k = 2, maxit = 10000, maxrestarts = 30)
  min_mm <- which.min(mm$mu)
  threshold <- mm$mu[min_mm] + mm$sigma[min_mm]*2
  percent_above <- test_scores %>% group_by(batch) %>% summarize(n = n(), nabove = sum(score >= threshold), pabove = nabove/n*100, condition = unique(condition), nCount_RNA = mean(nCount_RNA))
  percent_above$geneset <- .x
  return(percent_above)
})

use_genesets <- c("M3", "M4", "M7", "M12", "M5", "M13", "M18", "P11", "P13", "P19", "P3", "P1", "P8", "P9")
fpercent_above <- percent_above %>% filter(geneset %in% use_genesets)
mlabels <- unlist(marker_gene_name)[names(unlist(marker_gene_name)) %in% use_genesets]
plabels <- program_gene_names[names(program_gene_names) %in% use_genesets]
labels <- c(mlabels, plabels)
labels <- labels[use_genesets]
labels <- paste0(gsub("_", "", use_genesets), " ", labels)
labels

fpercent_above$geneset <- factor(fpercent_above$geneset,
  levels =  stringr::str_sort(use_genesets, numeric = TRUE))
fpercent_above$geneset

pabove_plot <- ggplot(fpercent_above,
  aes(x = condition, y = pabove, fill = condition)) +
  geom_boxplot(alpha = 0.5) +
  geom_point(shape = 21) +
  guides(fill = FALSE) +
  labs(x = "Condition", y = "% above threshold") +
  scale_x_discrete(labels = c("Control", "Mild", "Severe")) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  scale_fill_manual(values = c("#009E73", "#0072B2", "#D55E00"), labels = c("Healthy", "Mild", "Severe"), name = "Condition") +
  GeneralTheme(18) +
  facet_wrap(~ geneset, nrow = 2) + theme(
    strip.background = element_rect(size = 1),
    strip.text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )
pabove_plot
SavePlot(plot = pabove_plot, filename = "plots/liao_percentabove.pdf", base_height = 1.9*3, base_asp = 4/1.9)
SavePlot(plot = pabove_plot, filename = "plots/liao_percentabove.png", base_height = 1.9*3, base_asp = 4/1.9)
