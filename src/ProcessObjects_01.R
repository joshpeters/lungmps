# ---
# Description: Process objects
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
sets <- readxl::read_excel("data/study_metadata.xlsx")
ids <- sets$ID[sets$Train == TRUE]
ids <- SelfName(ids)
ids <- ids[order(ids)]
suffix <- "base_object.rds"
dir <- "data/objects"

load_files <- list.files(path = dir, pattern = suffix, full.names = TRUE)
load_files <- load_files[grep(paste0(ids, collapse = "|"), x = load_files)]
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

# [3] Process objects ----
# for each object, load, process, save
ids
walk(.x = ids, .f = ~ {
  usethis::ui_info("Loading {.x} ({round(file.info(load_files[.x])$size/1e6, 2)} mb)")
  full_obj <- readRDS(file = load_files[.x])
  usethis::ui_done(("\nCompleted loading {.x}"))
  Clean()

  full_obj <- SetIdent(full_obj, value = "batch")
  "disease_classification" %in% colnames(full_obj[[]])
  if (length(unique(full_obj$disease_classification)) == 1) {
    split_obj <- list(full_obj)
    names(split_obj) <- unique(full_obj$disease_classification)
  } else {
    full_obj$disease_classification <- as.factor(full_obj$disease_classification)
    split_obj <- Seurat::SplitObject(full_obj, split.by = "disease_classification")
  }


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
    obj <- BaseCluster(object = obj, graph.name = "harmony_snn", dims = 20, verbose = TRUE,
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

})
