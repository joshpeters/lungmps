# ---
# Description: Generate consensus cell states
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

library(ComplexHeatmap)

# load data
suffix <- "mnp_final.rds"
dir <- "data/objects"

load_files <- list.files(path = dir, pattern = suffix, full.names = TRUE)
ids <- str_match(string = load_files, pattern = glue("(data/objects/)(.*_(health|disease))_{suffix}"))[, 3, drop = TRUE]
ui_todo("Loading in {length(load_files)} files")
names(load_files) <- ids
ids <- SelfName(ids)

object_ls <- pblapply(ids, function(x) { return(readRDS(load_files[x])) })

# [2] Define markers ----

# define stm markers
stm_markers <- map_df(ids, .f = ~ {
  markers <- object_ls[[.x]]@misc$merged_clusters
  markers$set <- .x
  return(markers)
})

# define broad mwu markers
pb <- progress_estimated(length(ids))
mwu_markers <- map_df(.x = ids, .f = ~ {
  object_ls[[.x]] <- SetIdent(object_ls[[.x]], value = "merged_leiden")
  markers <- presto::wilcoxauc(object_ls[[.x]], assay = "data", seurat_assay = "RNA")
  markers$set <- .x
  Clean()
  pb$tick()$print()
  return(markers)
})

# double check stm markers
pb <- progress_estimated(length(ids))
stm_markers <- map_df(.x = ids, .f = ~ {
  object_ls[[.x]] <- SetIdent(object_ls[[.x]], value = "merged_leiden")
  markers <- IdentifyQualityMarkers(object_ls[[.x]], tlogFC = log(1.1))
  markers$set <- .x
  Clean()
  pb$tick()$print()
  return(markers)
})

saveRDS(stm_markers, file = glue("data/interim/ccs_stmmarkers.rds"))
saveRDS(mwu_markers, file = glue("data/interim/ccs_mwumarkers.rds"))

write_csv(stm_markers, path = "data/publication/TableS4.csv")

# [3] Pull and filter markers ----

markers = stm_markers
object_names = ids
k = 21
logFC_threshold = log(1.25)
auc_threshold = 0.6
padj_threshold = 1E-3
multi_marker_threshold = 0.6
num_genes = 50
min_degs = 5
min_jacc_remove = 0.8

# pull markers
pb <- progress_estimated(length(object_names))
def_markers <- map(.x = object_names, .f = ~ {
  fil_markers <- markers[markers$set == .x, ] %>%
    dplyr::filter(logFC >= logFC_threshold & auc >= auc_threshold & padj <= padj_threshold)
  num_group <- length(unique(fil_markers$group))
  bad_markers <- fil_markers %>% dplyr::group_by(feature) %>% dplyr::summarize(frac = n()/num_group) %>%
    filter(frac >= multi_marker_threshold) %>% dplyr::pull(feature)
  fil_markers <- fil_markers %>% dplyr::filter(!feature %in% bad_markers)
  fil_markers <- fil_markers %>% dplyr::group_by(group) %>% dplyr::top_n(num_genes, auc)
  pb$tick()$print()
  return(fil_markers)
}, markers = markers)

# split each dataset into a list of genes
markers_ls <- map(.x = def_markers, .f = ~ {
  marker_ls <- split(.x$feature, .x$group)
  return(marker_ls)
})
num_clusters <- unlist(lapply(markers_ls, length))
summary(num_clusters)
markers_ls <- unlist(markers_ls, recursive = FALSE)
length(markers_ls)

# remove clusters with less than 5 markers
markers_ls <- markers_ls[sapply(markers_ls, length) >=  min_degs]
length(markers_ls)

# tabulate and calculate jaccard distance
markers_totab <- markers_ls
markers_tab <- qdapTools::mtabulate(markers_totab)
jacc_dist <- (ade4::dist.binary(markers_tab, method = 1, diag = FALSE, upper = FALSE))^2
jacc_mtx <- as.matrix(jacc_dist)
diag(jacc_mtx) <- 1
org_jacc_mtx <- jacc_mtx

# filter poor performing matrices, lists with no matches above 0.1
thresholds <- (rowSums(jacc_mtx < min_jacc_remove))
remove <- names(thresholds[thresholds <= 0])
remove <- !(names(markers_totab) %in% remove)
markers_totab <- markers_totab[remove]
markers_tab <- qdapTools::mtabulate(markers_totab)
jacc_dist <- (ade4::dist.binary(markers_tab, method = 1, diag = FALSE, upper = FALSE))^2
jacc_mtx <- as.matrix(jacc_dist)
diag(jacc_mtx) <- 1

# [4] Assess clustering ----

# ensure ward method is best
m <- c( "average", "single", "complete", "ward")
names(m) <- c( "average", "single", "complete", "ward")
map_dbl(m, ~ { cluster::agnes(jacc_dist, method = .x)$ac })

# perform clustering
hc <- hclust(jacc_dist, method = "ward.D2")
pheatmap::pheatmap(
  jacc_dist,
  cluster_rows = hc,
  cluster_cols = hc,
  show_rownames = FALSE,
  show_colnames = FALSE,
  color = viridis::inferno(10)
)

# generate heuristic plots
wss <- factoextra::fviz_nbclust(jacc_mtx, FUN = factoextra::hcut, method = "wss", k.max = 35)
silo <- factoextra::fviz_nbclust(jacc_mtx, FUN = factoextra::hcut, method = "silhouette", k.max = 35)
stats <- data.frame(k = wss$data$clusters, wss = wss$data$y, silo = silo$data$y)
stats <- stats %>% mutate(rank_wss = dense_rank(wss), rank_silo = dense_rank(dplyr::desc(silo)))
stats <- stats %>% mutate(total_rank = rank_wss + rank_silo)
a <- ggplot(stats, aes(x = k, y = wss)) +
  geom_point(shape = 21, fill = "gray90", color = "black", size = 3, alpha = 1) +
  scale_x_discrete(breaks = c(1, seq(5, 30, 5))) +
  colorspace::scale_fill_discrete_sequential("Grays", name = "Statistic") +
  labs(x = "k", y = "Within sum-of-squares") +
  GeneralTheme(14)
b <- ggplot(stats, aes(x = k, y = silo)) +
  geom_point(shape = 21, fill = "gray90", color = "black", size = 3, alpha = 1) +
  scale_x_discrete(breaks = c(1, seq(5, 30, 5))) +
  colorspace::scale_fill_discrete_sequential("Grays", name = "Statistic") +
  labs(x = "k", y = "Silhouette Coefficient") +
  GeneralTheme(14)
heuristics_plot <- cowplot::plot_grid(a, b)
SavePlot(filename = "plots/ccs_clustering_heuristics.png", plot = heuristics_plot,
  base_height = 4, base_asp = 2)

# determine bootstrapping
cboots <- fpc::clusterboot(data = jacc_dist, B = 1000, distances = TRUE, clustermethod = fpc::hclustCBI, method = "ward.D2", k = 21)
cboot_results <- cboots$result$partition
table(cboot_results)
View(cboots$bootresult)
plot(cboots$bootmean)
cboots_results <- data.frame(partition = 1:21, bootmean = cboots$bootmean)
c <- ggplot(cboots_results, aes(x = partition, y = bootmean)) +
  scale_x_continuous(breaks = seq(1, 21, 2)) +
  scale_y_continuous(limits = c(0, 1)) +
  geom_point(shape = 21, fill = "gray90", color = "black", size = 3, alpha = 1) +
  labs(x = "Partition", y = "Mean bootstrap similarities") +
  GeneralTheme(14)
SavePlot(filename = "plots/ccs_clustering_bootstraps.png", plot = c, base_height = 4, base_asp = 1)

# define groups and generate cluster annotation dataframe
clustering <- cutree(hc, k = 21)
mclust::adjustedRandIndex(cboot_results, clustering)

anno <- data.frame(clustering)
anno$id <- rownames(anno)
anno$id <- factor(anno$id, levels = hc$labels[hc$order])

# add dataset information to anno
dataset_names <- stringr::str_match(pattern = "^(.*)[[:punct:]](\\d{1,2})$", string = anno$id)
anno$dataset <- dataset_names[, 2, drop = TRUE]
anno$dataset <- as.character(anno$dataset)
anno$orig_cluster <- dataset_names[, 3, drop = TRUE]
colnames(anno) <- c("Agglo_Cluster", "ID", "Dataset", "Orig_Cluster")

# add dynamic tree information
dtree <- dynamicTreeCut::cutreeDynamic(hc,
  minClusterSize = 5,
  deepSplit = 1,
  distM = jacc_mtx
)
table(dtree)
dtree <- ConvertNamedVecToDF(dtree)
colnames(dtree) <- c("name", "Dynamic_Tree")
dtree$clusters <- hc$labels
anno <- merge(anno, dtree %>% select(-name), by.x = 0, by.y = "clusters", all.x = TRUE)
mclust::adjustedRandIndex(anno$Agglo_Cluster, anno$Dynamic_Tree)
stm_results <- list(orig_matrix = org_jacc_mtx, markers = def_markers, matrix = jacc_mtx, hclust = hc, anno = anno)
saveRDS(stm_results, file = glue("data/interim/stm_results.rds"))

# [5] Heatmap and similarity plotting ----

PlotHeatmap(stm_results$matrix, stm_results$hclust, stm_results$hclust, anno = stm_results$anno,
  filename = glue("plots/markers_simmtx_stm_21_fullscale.pdf"))

jacc_mtx <- stm_results$matrix
rownames(stm_results$anno) <- stm_results$anno$Row.names
clusters <- split(stm_results$anno[, c("Dataset", "ID", "Orig_Cluster", "Agglo_Cluster")], stm_results$anno$Agglo_Cluster)
cluster_stats <- map_df(.x = names(clusters), .f = ~ {
  subset_jacc_mtx <- jacc_mtx[rownames(clusters[[.x]]), rownames(clusters[[.x]])]
  subset_jacc_mtx <- subset_jacc_mtx[lower.tri(subset_jacc_mtx)]
  avg_in <- mean(subset_jacc_mtx)
  med_in <- median(subset_jacc_mtx)
  sd_in <- sd(subset_jacc_mtx)
  mad_in <- mad(subset_jacc_mtx)

  subrow_jmtx <- jacc_mtx[rownames(clusters[[.x]]), !(rownames(jacc_mtx)  %in% rownames(clusters[[.x]]))]
  subcol_jmtx <- jacc_mtx[!(rownames(jacc_mtx)  %in% rownames(clusters[[.x]])), rownames(clusters[[.x]])]
  out_jmtx <- c(subrow_jmtx, subcol_jmtx)
  avg_out <- mean(out_jmtx)
  med_out <- median(out_jmtx)
  sd_out <- sd(out_jmtx)
  mad_out <- mad(out_jmtx)

  return(data.frame(avg_in = avg_in, med_in = med_in, sd_in = sd_in, mad_in = mad_in,
    avg_out = avg_out, med_out = med_out, sd_out = sd_out, mad_out = mad_out, cluster = .x))
})

LoadColors()
cluster_stats$cluster <- as.numeric(cluster_stats$cluster)
d <- ggplot(cluster_stats) +
  geom_point(mapping = aes(x = avg_in, y = cluster), shape = 21, fill = "#009E73", alpha = 0.6, color = "black", size = 3) +
  geom_point(mapping = aes(x = avg_out, y = cluster), shape = 21, fill = "#D55E00", alpha = 0.6, color = "black", size = 3) +
  geom_segment(mapping = aes(y = cluster, yend = cluster, xend = avg_out, x = avg_in)) +
  scale_y_continuous(breaks = seq(1, 21, 2)) +
  scale_x_continuous(limits = c(0.5, 1.0), breaks = seq(0.6, 1, 0.2)) +
  labs(x = "Partition", y = "Mean bootstrap similarities") +
  GeneralTheme(14) +
  theme(panel.grid.major.y = element_line(size = 1, color = "gray95", linetype = "dotted"))
d
SavePlot(filename = "plots/ccs_clustering_similarities.png", plot = d, base_height = 4, base_asp = 0.8)

# [6] Plot cluster diversity ----

df_sum <- stm_results$anno %>% group_by(Agglo_Cluster) %>% summarize(n = n())
dcasted <- reshape2::dcast(Dataset ~ Agglo_Cluster, data = stm_results$anno, fun.aggregate = length, value.var = "ID")
shannon <- vegan::diversity(dcasted[, -1], MARGIN = 2, index = "shannon")
simpson <- vegan::diversity(dcasted[, -1], MARGIN = 2, index = "simpson")
brillouin_v <- apply(dcasted[, -1], 2, Brillouin)

stats <- data.frame(cluster = names(shannon), Shannon = shannon, Simpson = simpson, Brillouin = brillouin_v)
stats <- gather(stats, key = "key", value = "value", -cluster)
stats$cluster <- factor(stats$cluster, levels = as.character(unique(sort(as.numeric(stats$cluster)))))
stats <- merge(stats, df_sum, by.x = "cluster", by.y = "Agglo_Cluster.x")
statplot <- ggplot(stats %>% filter(key == "Shannon"),
  aes(x = cluster, y = value)) +
  geom_point(shape = 21, fill = "#0072B2", color = "black", size = 3, alpha = 1) +
  labs(x = "Partition", y = "Shannon\nDiversity") +
  GeneralTheme(14)
statplot

dataset_labels <- gsub("_", " ", ids)
dataset_labels <- snakecase::to_title_case(dataset_labels)
dataset_labels

dcasted_gather <- gather(dcasted, key = "cluster", value = "sum", -Dataset)
dcasted_gather$cluster <- factor(dcasted_gather$cluster, levels = as.character(unique(sort(as.numeric(dcasted_gather$cluster)))))
div_tile <- ggplot(dcasted_gather %>% filter(cluster %in% unique(stm_markers$partition)), aes(x = cluster, y = Dataset, fill = sum)) +
  geom_tile(color = "black", size = 0.5) +
  scale_y_discrete(labels = dataset_labels, expand = c(0, 0)) +
  colorspace::scale_fill_continuous_sequential("Blues", name = "# datasets", breaks = c(0, 1), limits = c(0, 1),
    na.value = "#56B4E9", expand = c(0, 0)) +
  labs(x = "Partition", y = "Dataset") +
  TileTheme() + NoLegend()
div_tile
div_tile <- ggplot(dcasted_gather, aes(x = cluster, y = Dataset, fill = sum)) +
  geom_tile(color = "black", size = 0.5) +
  scale_y_discrete(labels = dataset_labels, expand = c(0, 0)) +
  scale_fill_gradient(low = "white", high = "#0072B2",  name = "# datasets", breaks = c(0, 1), limits = c(0, 1)) +
  # colorspace::scale_fill_continuous_sequential("Blues", name = "# datasets", breaks = c(0, 1), limits = c(0, 1),
  #   na.value = "#56B4E9", expand = c(0, 0)) +
  labs(x = "Partition", y = "Dataset") +
  TileTheme() + NoLegend()
div_tile
div_plot <- cowplot::plot_grid(statplot, div_tile, align = c("v"), axis = "lr", nrow = 2, rel_heights = c(0.5, 1))
SavePlot(filename = glue("plots/stmmarker_diversity_k21.pdf"), plot = div_plot, base_asp = 1, base_height = 6)
SavePlot(filename = glue("plots/stmmarker_diversity_k21.png"), plot = div_plot, base_asp = 1, base_height = 6)

# [7] Collate markers ----

clusters <- split(stm_results$anno[, c("Dataset", "ID", "Orig_Cluster", "Agglo_Cluster")], stm_results$anno$Agglo_Cluster)
markers <- stm_results$markers

stm_markers <- CollateMarkers(clusters, markers)
saveRDS(stm_markers, file = glue("data/interim/ccs_markers.rds"))

stm_markers_tosave <- stm_markers %>% select(feature, n_sets, median_auc, median_logFC, partition)
colnames(stm_markers_tosave) <- c("Feature", "Number of datasets", "Median auROC", "Median log FC", "Partition")
write_csv(stm_markers_tosave, path = "data/publication/TableS5.csv")

# check for similarities
stm_comp <- map(unique(stm_markers$partition), ~ {
  return(stm_markers %>% filter(partition == .x) %>% top_n(100, median_auc) %>% pull(feature))
})

gom.obj <- GeneOverlap::newGOM(
  gsetA = stm_comp,
  gsetB = stm_comp,
  genome.size = 3000)
jacc <- GeneOverlap::getMatrix(gom.obj, "Jaccard")
diag(jacc) <- 0
which(jacc > 0.1, arr.ind = TRUE)
image(jacc)
ddf <- as.data.frame(jacc)
ddf$c1 <- rownames(ddf)
ddf <- ddf %>% gather(key = "c2", value = "x", -c1)
ddf$c2 <- gsub("V", "", ddf$c2)
ddf$c1 <- factor(ddf$c1, levels = as.character(unique(ddf$c1)))
ddf$c2 <- factor(ddf$c2, levels = as.character(unique(ddf$c2)))
ddf$x <- 1-ddf$x
ddf_filtered <- ddf %>% group_by(c1) %>% filter(x < 1.0) %>% top_n(-1, x) %>% arrange(x)

# [8] Collate AUC values ----

top_genes <- lapply(object_names, function(x) {
  return(stm_results$markers[[x]] %>% pull(feature))
})
all_genes <- sort(unique(unlist(top_genes)))
all_genes <- data.frame(feature = all_genes, stringsAsFactors = FALSE)

# pull variable of interest from each marker file and merge
voi <- "auc"
fill_value <- 0

voi_mtxs <- pblapply(object_names, function(a) {
  ui_info("Generating matrices for {a}")
  df <- stm_results$markers[[a]] %>% mutate(group_set = paste0(group, "_", set))
  df <- df %>% ungroup() %>% select(feature, !!sym(voi), group_set)
  sdf <- split(df[, c("feature", voi)], df[, "group_set"])
  voi_values <- lapply(sdf, function(b) {
    merged_to_all <- merge(all_genes, b, by = "feature", all.x = TRUE)
    assertthat::assert_that(all.equal(merged_to_all$feature, all_genes$feature))
    merged_to_all[is.na(merged_to_all)] <- fill_value
    return_values <- merged_to_all %>% pull(!!sym(voi))
    #return_values <- -log10(return_values)
    return(return_values)
  })
  voi_values <- Reduce(cbind, voi_values)
  rownames(voi_values) <- all_genes$feature
  colnames(voi_values) <- names(sdf)
  return(voi_values)
})
voi_full <- Reduce(cbind, voi_mtxs)
Check5(voi_full)
dim(voi_full)
saveRDS(voi_full, file = "data/interim/AUCmatrix_allmarkers.rds")

stm_markers_filtered <- stm_markers %>% filter(!(partition %in% c(6, 14, 19, 21)))
stm_markers_filtered$partition <- as.factor(stm_markers_filtered$partition)
cluster_order <- as.character(unique(stm_results$anno$Agglo_Cluster[order(stm_results$anno$ID)]))
cluster_order <- cluster_order[cluster_order %in% unique(stm_markers_filtered$partition)]
stm_markers_filtered$partition <- factor(stm_markers_filtered$partition, levels = cluster_order)
ordered_markers_df <- stm_markers_filtered %>% arrange(partition, desc(median_auc))
ordered_markers_df$partition

sub_voi <- voi_full[unique(ordered_markers_df$feature), ]
sub_anno <- stm_results$anno
sub_anno$name <- paste0(sub_anno$Orig_Cluster, "_", sub_anno$Dataset)
sub_anno <- sub_anno[sub_anno$Agglo_Cluster %in% unique(stm_markers_filtered$partition), ]
sub_anno$Agglo_Cluster <- factor(sub_anno$Agglo_Cluster, levels = cluster_order)
sub_anno$name[order(sub_anno$Agglo_Cluster)]
sub_voi <- sub_voi[, sub_anno$name[order(sub_anno$Agglo_Cluster)]]
dim(sub_voi)

# set matrix colors
colfunc <- colorRampPalette(c("white", "#009E73"))
voi_colors <- circlize::colorRamp2(c(0, seq(0.5, 1, 0.05)), colfunc(12))

final_markers <- stm_markers_filtered %>%
  filter(!(feature %in% stm_markers_filtered$feature[grep("(MT-|RPS|RPL)", stm_markers_filtered$feature)])) %>%
  group_by(feature) %>%
  top_n(1, median_auc) %>%
  group_by(partition) %>%
  top_n(3, median_auc)
genes_to_label <- final_markers$feature
plot_indices <- match(genes_to_label, rownames(sub_voi))

genes_anno <- rowAnnotation(
  Genes = anno_mark(
    at = plot_indices,
    labels = genes_to_label,
    side = "left",
    labels_gp = gpar(fontsize = 8),
    link_width = unit(15, "mm"),
    padding = unit(2, "mm")
  )
)

voi <- Heatmap(sub_voi,
  border = TRUE,
  col = voi_colors,
  cluster_rows = FALSE, cluster_columns = FALSE,
  show_column_names = FALSE, show_column_dend = FALSE,
  column_title = "Gene AUC across clusters", column_title_side = "bottom",
  show_row_names = FALSE, show_row_dend = FALSE,
  left_annotation = genes_anno,
  heatmap_legend_param = list(
    title = "AUC",
    at = c(0, 0.25, 0.5, 0.75, 1.0),
    border = "black",
    legend_height = unit(4, "cm"),
    legend_width = unit(1, "cm")
  ),
  top_annotation = HeatmapAnnotation(
    Cluster = anno_block(
      gp = gpar(fill = "Gray80"),
      labels = as.character(cluster_order),
      labels_gp = gpar(fontsize = 8, col = "black", fontface = "bold")
    ),
    show_legend = c(FALSE)
  ),
  column_split = sub_anno$Agglo_Cluster[order(sub_anno$Agglo_Cluster)],
  height = 4,
  use_raster = TRUE, raster_quality = 1
)

dev.off()
#png(glue("plots/stmmarkers_auc_k21.png"), width = 10, height = 6, units = "in", res = 400)
pdf(glue("plots/stmmarkers_auc_k21.pdf"), width = 10, height = 6)
draw(voi,
  heatmap_legend_side = "right",
  use_raster = TRUE, raster_quality = 5
)
dev.off()

input_corr_mtx <- t(sub_voi)
input_corr_mtx <- merge(input_corr_mtx, sub_anno %>% select(name, Agglo_Cluster), by.x = 0, by.y = "name")
input_corr_mtx <- input_corr_mtx %>% group_by(Agglo_Cluster) %>% select(-Row.names) %>% summarize_all(median)
rownames(input_corr_mtx) <- input_corr_mtx[, 1, drop = TRUE]
corr_mtx <- cor(t(input_corr_mtx[, -1]))
rownames(corr_mtx) <- rownames(input_corr_mtx)
colnames(corr_mtx) <- rownames(input_corr_mtx)

map <- ComplexHeatmap::pheatmap(corr_mtx, border_color = NA, color = rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")),
  show_colnames = FALSE, fontsize = 14, treeheight_col = FALSE)
pdf(glue("plots/ccs_auc_corr.pdf"), width = 6, height = 5)
draw(map,
  heatmap_legend_side = "right",
  use_raster = TRUE, raster_quality = 5
)
dev.off()

dim(voi_full)
umap <- uwot::umap(t(voi_full), n_neighbors = 10, spread = 5, min_dist = 0.01)
umap <- as.data.frame(umap)
umap$sets <- colnames(voi_full)
anno$name <- paste0(anno$Orig_Cluster, "_", anno$Dataset)
umap_df <- merge(umap, anno, by.x = "sets", by.y = "name", all.x = TRUE)
umap_df$C1 <- factor(umap_df$Agglo_Cluster, levels = as.character(unique(anno$Agglo_Cluster[order(anno$ID)])))
centroids <- umap_df %>% group_by(C1) %>% summarize(mean_v1 = mean(V1), mean_v2 = mean(V2))
auc_umap <- ggplot(umap_df, aes(x = V1, y = V2, fill = C1)) +
  geom_point(size = 3, shape = 21, color = "black", alpha = 0.8, stroke = 0.5) +
  ggrepel::geom_label_repel(data = centroids, aes(x = mean_v1, y = mean_v2, label = C1), fill = "white", fontface = "bold", color = "black") +
  colorspace::scale_fill_discrete_qualitative("Dark 2", name = "Partitions") +
  labs(x = "UMAP 1", y = "UMAP 2") +
  UMAPTheme(14)
auc_umap
cowplot::save_plot(filename = "plots/ccs_umap_all.png", plot = auc_umap, base_asp = 1.2, base_height = 6)

dim(sub_voi)
View(sub_voi)
dim(input_corr_mtx)
sub_voi_withmedian <- cbind(sub_voi, t(input_corr_mtx[, -1]))
sub_voi_withmedian <- t(sub_voi_withmedian)
umap <- uwot::umap(X = sub_voi_withmedian, n_neighbors = 10, spread = 5, min_dist = 0.01)
umap <- as.data.frame(umap)
umap$sets <- c(colnames(sub_voi), input_corr_mtx$Agglo_Cluster)
umap$consensus <- FALSE
umap$consensus[c(139:155)] <- TRUE
auc_umap <- ggplot(umap, aes(x = V1, y = V2, fill = consensus)) +
  geom_point(size = 3, shape = 21, color = "black", alpha = 0.8, stroke = 0.5) +
  scale_fill_manual(values = c("gray80", "#009E73"), name = "Consensus") +
  labs(x = "UMAP 1", y = "UMAP 2") +
  UMAPTheme(14)
auc_umap
cowplot::save_plot(filename = "plots/ccs_umap_consensus.png", plot = auc_umap, base_asp = 1.2, base_height = 6)

# [9] Explore clusters that are poorly defined ----
markers_totab <- markers_ls
markers_tab <- qdapTools::mtabulate(markers_totab)
jacc_dist <- (ade4::dist.binary(markers_tab, method = 1, diag = FALSE, upper = FALSE))^2
jacc_mtx <- as.matrix(jacc_dist)
diag(jacc_mtx) <- 1
org_jacc_mtx <- jacc_mtx

dim(jacc_mtx)
missed_markers <- rownames(jacc_mtx[rowSums(jacc_mtx < 0.8) == 0, ])
passing_markers <- rownames(jacc_mtx[rowSums(jacc_mtx <= 0.8) > 0, ])
mean(sapply(markers_ls[missed_markers], length))
mean(sapply(markers_ls[passing_markers], length))

# [10] Enrichment ----

stm_comp <- map(unique(stm_markers$partition), ~ {
  return(stm_markers %>% filter(partition == .x & !(feature %in% stm_markers_filtered$feature[grep("(MT-|RPS|RPL)", stm_markers_filtered$feature)])) %>% top_n(100, median_auc) %>% pull(feature))
})

names(stm_comp) <- as.character(seq(1:length(stm_comp)))
stm_comp_filtered <- stm_comp[-c(6, 14, 19, 21)]

library(enrichR)
dbs <- listEnrichrDbs()
dbs <- c("GO_Biological_Process_2018")
dbs <- c("GO_Molecular_Function_2018")
dbs <- c("KEGG_2019_Human")
dbs <- c("BioCarta_2016")
dbs <- c("Reactome_2016")

dbs <- c("TRRUST_Transcription_Factors_2019")
dbs <- c("TRANSFAC_and_JASPAR_PWMs")
#dbs <- c("ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X", "ARCHS4_TFs_Coexp")
#dbs <- c("LINCS_L1000_Ligand_Perturbations_up")

enr <- map_dfr(names(stm_comp_filtered), .f = ~ {
  enr <- enrichr(stm_comp_filtered[[.x]], databases = dbs)
  enr <- lapply(names(enr), function(x) {
    enr[[x]]$db <- x
    return(enr[[x]])
  })
  enr <- bind_rows(enr)
  enr$geneset <- .x
  return(enr)
})
tophits <- enr %>% group_by(geneset) %>% filter(Adjusted.P.value < 0.1) %>% top_n(-5, Adjusted.P.value)

order <- as.character(unique(anno$Agglo_Cluster[order(anno$ID)]))
order <- order[!(order %in% c(6, 14, 19, 21))]
tophits$log10padj <- -log10(tophits$Adjusted.P.value)
tophits$log10padj[tophits$log10padj > 50] <- 50
tophits <- tophits %>% arrange(desc(log10padj))
tophits$geneset <- factor(tophits$geneset,
  levels = order)
#tophits$short_term <- gsub("( human tf ARCHS4 coexpression)|( ENCODE)|( CHEA)", "", tophits$Term)
#reordered_program_names <- program_gene_names[as.character(levels(tophits$geneset))]
# enr_toplot <- enr %>% filter(!is.na(factor) & neglog10padj > 1) %>%
#   group_by(factor) %>% top_n(1, neglog10padj) %>% arrange(desc(factor))
enrplot <- ggplot(tophits, aes(x = log10padj, y = geneset, label = Term)) +
  geom_label(hjust = 0, size = 4) +
  scale_x_continuous(limits = c(0, 60), breaks = c(0, 10, 25, 50)) +
  scale_y_discrete(labels = paste(gsub("M", "", levels(tophits$geneset)), order, sep = " ")) +
  geom_vline(xintercept = -log10(0.1), color = "darkred", linetype = "dotted") +
  labs(x = "-log10(adj. P value)", y = "Program geneset", title = "ARCHS4 TF Coexpression") +
  GeneralTheme(18)
enrplot
SavePlot(plot = enrplot, filename = "plots/program_enr_archs4tfcoexp_{stamp$date}.pdf", base_height = 6, base_asp = 1)

enrplot <- ggplot(tophits, aes(x = geneset, y = Term, fill = log10padj)) +
  geom_tile(color = "black", size = 0.5) +
  scale_fill_gradient(low = "white", high = "#0072B2", name = "-log10\nadj. P value") +
  labs(x = "Consensus cell states", y = "", title = "Reactome") +
  GeneralTheme(14) + theme(legend.position = "right") +
  theme(axis.text.y = element_blank(),
    panel.grid.major = element_line(size = 0.5, color = "gray90"),
    panel.grid.minor = element_line(size = 0.5, color = "gray90"))
SavePlot(plot = enrplot, filename = "plots/ccs_reactome.png", base_height = 6, base_asp = 0.8)

enrplot <- ggplot(tophits, aes(x = geneset, y = Term, fill = log10padj)) +
  geom_tile(color = "black", size = 0.5) +
  scale_fill_gradient(low = "white", high = "#0072B2", name = "-log10\nadj. P value") +
  labs(x = "Consensus cell states", y = "", title = "Reactome") +
  GeneralTheme(14) + theme(legend.position = "right")
SavePlot(plot = enrplot, filename = "plots/ccs_reactome.pdf", base_height = 6, base_asp = 2.5)

# [11] Calculate scmap similarities ----
library(SingleCellExperiment)
library(scmap)

bad_genes <- readRDS("data/nonmp_markers_total.rds")
names(markers)
tlogFC = log(1.25)
tauc = 0.6
tpadj = 1E-5
top_markers <- pblapply(ids, function(x, markers) {
  tm <- markers[[x]] %>% filter(logFC > tlogFC & auc > tauc & padj < tpadj)
  num_group <- length(unique(tm$group))
  bad_markers <- tm %>% group_by(feature) %>% summarize(frac = n()/num_group) %>%
    filter(frac >= 0.6) %>% pull(feature)
  tm <- tm %>% filter(!feature %in% bad_markers)
  tm <- tm %>% group_by(group) %>% top_n(100, auc)
  return(tm)
}, markers = markers)
top_markers <- bind_rows(top_markers)
markers_sum <- top_markers %>% group_by(feature) %>% summarize(n = n(), median_auc = median(auc))
markers_final <- markers_sum %>% filter(n >= 4 & median_auc > 0.6) %>% arrange(desc(median_auc))
features_to_use <- markers_final$feature[1:1000]

sce_ls <- pblapply(ids, function(x) {
  ui_todo("Preparing {x}")
  sce <- as.SingleCellExperiment(object_ls[[x]])
  rowData(sce)$feature_symbol <- rownames(sce)
  scmap_features <- features_to_use[features_to_use %in% rownames(sce)]
  rowData(sce)$scmap_features <- rowData(sce)$feature_symbol %in% scmap_features
  ui_todo("Indexing clusters for {x}")
  sce <- indexCluster(sce, cluster_col = "merged_leiden")
  Clean()
  return(sce)
})

# rip indices
cluster_indices <- pblapply(ids, function(x) {
  cluster_index <- metadata(sce_ls[[x]])$scmap_cluster_index
  return(cluster_index)
})
idx <- Reduce(rbind, cluster_indices)

# define labels
labels <- pblapply(ids, function(x) {
  num_clusters <- nrow(cluster_indices[[x]])
  labels <- rep(x, num_clusters)
  return(labels)
})
labels <- unlist(labels)

# cross-classify
#test_names <- sample(18, 1)
mappings <- pblapply(ids, function(x) {
  mapped <- ClusterClassifier(object = sce_ls[[x]], indices = cluster_indices, threshold = 0.3)
  Clean()
  return(mapped)
})

# save mappings
saveRDS(mappings, file = "data/interim/scmap_03threshold_ranked1000aucmarkers.rds")

sims <- pblapply(ids, function(x) {
  labels <- mappings[[x]]$scmap_cluster_labs
  rownames(labels) <- colnames(sce_ls[[x]])
  dfs <- lapply(colnames(labels), function(y) {
    df <- data.frame(cells = colnames(object_ls[[x]]),
      assignments = labels[, y],
      clusters = object_ls[[x]]$merged_leiden)
    byquery <- df %>%
      group_by(clusters, assignments) %>%
      summarize(nref = n()) %>%
      group_by(clusters) %>%
      mutate(prop = nref/sum(nref), ntotal = sum(nref))
    byquery$sqrtn <- sqrt(byquery$nref)
    byquery$assignments <- as.character(byquery$assignments)
    byquery$assignments[byquery$assignments == "unassigned"] <- "UA"
    byquery$prop[byquery$prop < 0.01 | byquery$nref < 5] <- 0
    byquery$ref <- y
    byquery$query <- x
    byquery$ref_id <- paste0(byquery$ref, "_", byquery$assignments)
    byquery$query_id <- paste0(byquery$query, "_", byquery$clusters)
    return(byquery)
  })
  return(dfs)
})
sims <- Reduce(rbind, sims)
sims <- Reduce(rbind, sims)
saveRDS(sims, file = "data/interim/scmap_03threshold_ranked1000aucmarkers_similaritiesdf.rds")

scmap_sims <- sims
# for each group of clusters, find average similarity
sims <- imap_dfr(.x = clusters, .f = ~ {
  .x$ID <- gsub("\\.", "_", .x$ID)
  sims_pull <- scmap_sims %>% filter(ref_id %in% .x$ID & query_id %in% .x$ID)
  sims_result <- data.frame(cluster = .y, mean_sim = mean(sims_pull$prop),
    median_sim = median(sims_pull$prop), mad_sim = mad(sims_pull$prop), sd_sim = sd(sims_pull$prop))
  return(sims_result)
})

# random
sets <- unique(scmap_sims$ref_id)
sets <- sets[-(grep("_UA$", sets))]
iterations <- 1:30
samples <- sapply(clusters, nrow)
samples <- samples[samples >= 4]
sampling <- lapply(samples, function(y) {
  sampled10_sims <- map_df(iterations, ~ {
    sampled <- sample(sets, size = y)
    sims_pull <- scmap_sims %>% filter(ref_id %in% sampled & query_id %in% sampled)
    sims_result <- data.frame(cluster = glue("r{.x}"), mean_sim = mean(sims_pull$prop),
      median_sim = median(sims_pull$prop), mad_sim = mad(sims_pull$prop), sd_sim = sd(sims_pull$prop))
    sims_result$type <- "Random"
    sims_result$n <- y
    return(sims_result)
  })
})

ggthemes::colorblind_pal()(8)
scales::show_col(ggthemes::colorblind_pal()(8))
sims$type <- "Test"
sampling <- bind_rows(sampling)
all_sims <- bind_rows(sims, sampling)
sim_plot <- ggplot(all_sims, aes(x = median_sim, y = ..scaled.., fill = type)) +
  geom_density(alpha = 0.5, adjust = 0.5, color = "black") +
  labs(x = "Median scmap similarity", y = "Scaled density") +
  scale_fill_manual(values = c("Gray60", "#D55E00"), name = "Comparison\nlists") +
  GeneralTheme(14) + theme(legend.position = "right")
sim_plot <- ggplot(all_sims, aes(x = type, y = median_sim)) +
  geom_point(color = "black", alpha = 0.2) +
  geom_boxplot(aes(fill = type), alpha = 0.8, width = 0.5) +
  scale_x_discrete(labels = c("Shuffled", "Partitions")) +
  scale_fill_manual(values = c("Gray60", "#D55E00"), name = "Comparison\nlists") +
  labs(x = "", y = "Median similarity") +
  GeneralTheme(14) + theme(legend.position = "none")
cowplot::save_plot(filename = glue("plots/ccs_scmapsim.png"), plot = sim_plot, base_asp = 1, base_height = 4)
