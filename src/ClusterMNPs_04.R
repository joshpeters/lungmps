# ---
# Description: Determine downstream clustering
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
suffix <- "mnp_final.rds"
dir <- "data/objects"

load_files <- list.files(path = dir, pattern = suffix, full.names = TRUE)
ids <- str_match(string = load_files, pattern = glue("(data/objects/)(.*_(health|disease))_{suffix}"))[, 3, drop = TRUE]
ui_todo("Loading in {length(load_files)} files")
names(load_files) <- ids
ids <- SelfName(ids)

# [2] Load and determine most-optimal clustering ----
print(ids)
walk(.x = ids, .f = ~ {

  # load object
  usethis::ui_info("Loading {.x} ({round(file.info(load_files[.x])$size/1e6, 2)} mb)")
  object <- readRDS(file = load_files[.x])
  usethis::ui_done(("\nCompleted loading {.x}"))
  Clean()

  # perform clustering
  object <- FindNeighbors(object = object, reduction = "harmony", graph.name = c("harmony_nn", "harmony_snn"),
    dims = 1:20, verbose= TRUE, compute.SNN = TRUE, force.recalc = TRUE)
  results <- IterativeCluster(object, graph_name = "harmony_snn")
  plotgrid <- cowplot::plot_grid(plotlist = results$plots, align = "hv", labels = "auto")
  cowplot::save_plot(filename = glue("plots/{.x}_mnp_cluster_merges_RNAsnn.png"),
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
      glue::glue("Clustering results for {.x} (ARI = {round(mclust::adjustedRandIndex(as.character(object$merged_leiden), as.character(object$equiv_leiden)), 3)})"),
      fontface = "bold", x = 0, hjust = 0) +
    theme(plot.margin = margin(0, 0, 0, 12, unit = "pt"))
  plotgrid <- cowplot::plot_grid(title, cowplot::plot_grid(a, b, labels = "auto", align = "hv"), ncol = 1, rel_heights = c(0.05, 1))
  cowplot::save_plot(filename = glue("plots/{.x}_mnp_cluster_results.png"),
    plot = plotgrid, base_height = 8, base_asp = 2)
  Clean()

  # resave object
  ui_todo("Saving object...")
  ui_info("\nSaving {.x} ({round(as.numeric(object.size(object)/1E6), 2)} mb)")
  saveRDS(object, file = glue::glue("data/objects/{.x}_mnp_final.rds"))
  Clean()
  ui_done(("\nCompleted {.x}"))

})
