# ---
# Description: Cluster factors (LIGER results) using Jaccard distance
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

# [2] Load data ----
ws <- readRDS(file = "data/interim/median_W.rds")
object_names <- names(median_w)
object_names <- SelfName(object_names)

# [3] Cluster ----

# determine genes for each factor, top 95% percentile ~ loadings
fs <- pblapply(object_names, function(x) {
  fs <- apply(ws[[x]], 1, function(y) {
    y <-  y[y > 0]
    return(names(y)[y > quantile(y, 0.95, na.rm = TRUE)])
  })
})
any(sapply(fs, function(x) sapply(x, length)) < 10)

# unlist factors and add names
factors_ls <- unlist(unname(fs), recursive = FALSE, use.names = TRUE)
names(factors_ls) <- paste0(seq(1:20), "_", rep(object_names, each = 20))

# calculate distances
factors_ls_iter <- factors_ls
factors_tab <- qdapTools::mtabulate(factors_ls_iter)
d <- (ade4::dist.binary(factors_tab, method = 1, diag = FALSE, upper = FALSE))^2
d_mtx <- as.matrix(d)
diag(d_mtx) <- 1
orig_mtx <- d_mtx
length(factors_ls_iter)

# assess matrix
hc <- hclust(d, method = "ward.D2")
pheatmap::pheatmap(
  d_mtx,
  cluster_rows = hc,
  cluster_cols = hc,
  show_rownames = FALSE,
  show_colnames = FALSE,
  color = viridis::inferno(10)
)

# filter poor performing matrices
thresholds <- rowSums(d_mtx < 0.80)
keep <- names(thresholds[thresholds > 0])
factors_ls_iter <- factors_ls[keep]
factors_tab <- qdapTools::mtabulate(factors_ls_iter)
d <- (ade4::dist.binary(factors_tab, method = 1, diag = FALSE, upper = FALSE))^2
d_mtx <- as.matrix(d)
diag(d_mtx) <- 1

# reassess matrix
hc <- hclust(d, method = "ward.D2")
pheatmap::pheatmap(
  d_mtx,
  cluster_rows = hc,
  cluster_cols = hc,
  show_rownames = FALSE,
  show_colnames = FALSE,
  color = viridis::inferno(10)
)

# ensure ward method is best
m <- c( "average", "single", "complete", "ward")
names(m) <- c( "average", "single", "complete", "ward")
map_dbl(m, ~ { cluster::agnes(d_mtx, method = .x)$ac })

# generate heuristic plots
wss <- factoextra::fviz_nbclust(d_mtx, FUN = factoextra::hcut, method = "wss", k.max = 35)
silo <- factoextra::fviz_nbclust(d_mtx, FUN = factoextra::hcut, method = "silhouette", k.max = 35)
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
SavePlot(filename = "plots/cgp_clustering_heuristics.png", plot = heuristics_plot,
  base_height = 4, base_asp = 2)

# determine bootstrapping
cboots <- fpc::clusterboot(data = d_mtx, B = 1000, distances = TRUE, clustermethod = fpc::hclustCBI, method = "ward.D2", k = 23)
cboot_results <- cboots$result$partition
cboots_results <- data.frame(partition = 1:23, bootmean = cboots$bootmean)
c <- ggplot(cboots_results, aes(x = partition, y = bootmean)) +
  scale_x_continuous(breaks = seq(1, 21, 2)) +
  scale_y_continuous(limits = c(0, 1)) +
  geom_point(shape = 21, fill = "gray90", color = "black", size = 3, alpha = 1) +
  labs(x = "Partition", y = "Mean bootstrap similarities") +
  GeneralTheme(14)
c
SavePlot(filename = "plots/cgp_clustering_bootstraps.png", plot = c, base_height = 4, base_asp = 1)

# [4] Gather and save results ---

clustering <- cutree(hc, k = 23)
mclust::adjustedRandIndex(cboot_results, clustering)

anno <- data.frame(clustering)
anno$id <- rownames(anno)
anno$id <- factor(anno$id, levels = hc$labels[hc$order])

# add dataset information to anno
dataset_names <- stringr::str_match(pattern = "^(\\d{1,2})[[:punct:]](.*)$", string = anno$id)
anno$dataset <- dataset_names[, 3, drop = TRUE]
anno$dataset <- as.character(anno$dataset)
anno$orig_factor <- dataset_names[, 2, drop = TRUE]
colnames(anno) <- c("Agglo_Cluster", "ID", "Dataset", "Orig_Factor")

# add dynamic tree information
dtree <- dynamicTreeCut::cutreeDynamic(hc,
  minClusterSize = 5,
  deepSplit = 1,
  distM = d_mtx
)
dtree <- ConvertNamedVecToDF(dtree)
colnames(dtree) <- c("name", "Dynamic_Tree")
dtree$factors <- hc$labels
anno <- merge(anno, dtree %>% select(-name), by.x = "ID", by.y = "factors", all.x = TRUE)
mclust::adjustedRandIndex(anno$Agglo_Cluster, anno$Dynamic_Tree)
cgp_results <- list(orig_matrix = orig_mtx, factors = factors_ls, matrix = d_mtx, hclust = hc, anno = anno)
saveRDS(cgp_results, file = "data/interim/cgp_results.rds")

cgp_results <- readRDS("data/interim/cgp_results.rds")
mclust::adjustedRandIndex(cgp_results$anno$Agglo_Cluster, cgp_results$anno$Dynamic_Tree)
cboots <- fpc::clusterboot(data = cgp_results$matrix, B = 1000, distances = TRUE, clustermethod = fpc::hclustCBI, method = "ward.D2", k = 23)
cboot_results <- cboots$result$partition
cboot_df <- qdapTools::list2df(cboot_results)
check_cboot <- merge(cgp_results$anno, cboot_df, by.x = "ID", by.y = "X2")
mclust::adjustedRandIndex(check_cboot$X1, check_cboot$Agglo_Cluster)

# [5] Distance-based clustering ----

# create base matrix to merge all factor loadings
ws_genes <- lapply(ws, colnames) # genes within each W are not the same
genes <- Reduce(union, ws_genes)
intersect_genes <- Reduce(intersect, ws_genes)

collated_factors <- Matrix.utils::rBind.fill(ws)
dim(collated_factors)
rownames(collated_factors)
colnames(collated_factors)
nas <- which(is.na(collated_factors), arr.ind = TRUE)
rownames(collated_factors)[nas[, 1]]
colnames(collated_factors)[nas[, 2]]
collated_factors[is.na(collated_factors)] <- 0

# find compiled factors partitions
rownames(reduced_coll_factors) <- paste(rep(seq(1:20), times = 18), rep(object_names, each = 20), sep = "_")
coln <- colnames(reduced_coll_factors)
rown <- rownames(reduced_coll_factors)
compiled_W <- as.matrix(reduced_coll_factors)

# calculate UMAP
reduced_compiled_W <- compiled_W[, intersect_genes]
dim(reduced_compiled_W)
#reduced_compiled_W <- reduced_compiled_W[names(factors_ls_iter), ]
dist_W <- MatrixCosSim(reduced_compiled_W, reduced_compiled_W)
factor_umap <- uwot::umap(dist_W, n_neighbors = 20, spread = 5, min_dist = 0.5)
factor_umap <- as.data.frame(factor_umap)
factor_umap$dataset <- rep(object_names, each = 20)
factor_umap$factor <- rep(seq(1:20), times = 18)
factor_umap$id <- rownames(dist_W)
colnames(factor_umap) <- c("UMAP1", "UMAP2", "dataset", "factor_k", "id")

plot <- PlotFactorUMAP(factor_umap, "dataset", label = FALSE); plot
SavePlot(filename = glue("plots/factor_bydataset_umap.pdf"), plot = plot,
  base_height = 6, base_asp = 1)

# kmeans on umap coords
kmeans_umap <- kmeans(factor_umap[, c(1, 2)], centers = 23)
factor_umap$umap_kmeans <- as.character(kmeans_umap$cluster)

kmeans_W <- kmeans(compiled_W, centers = 23, nstart = 30, iter.max = 100)
factor_umap$w_kmeans <- as.character(kmeans_W$cluster)

sc <- kernlab::specc(compiled_W, centers = 23)
factor_umap$w_spec <- as.character(sc@.Data)

sc <- kernlab::specc(as.matrix(factor_umap[, c(1, 2)]), centers = 23, iterations = 1000)
factor_umap$umap_spec <- as.character(sc@.Data)

cluster_comps <- combn(colnames(factor_umap)[6:9], 2)
cluster_comps
map_dbl(1:ncol(cluster_comps), ~ {
  mclust::adjustedRandIndex(
    factor_umap[, cluster_comps[1, .x], drop = TRUE],
    factor_umap[, cluster_comps[2, .x], drop = TRUE])
}, cluster_comps = cluster_comps, factor_umap = factor_umap)

plot <- PlotFactorUMAP(factor_umap, fill = "w_kmeans"); plot
plot <- PlotFactorUMAP(factor_umap, fill = "umap_kmeans"); plot
plot <- PlotFactorUMAP(factor_umap, fill = "w_spec"); plot
plot <- PlotFactorUMAP(factor_umap, fill = "umap_spec"); plot

# [6] Finalize clusterings ----

PlotHeatmap(matrix = d_mtx, row_cluster = hc, column_cluster = hc, anno = anno,
  filename = glue("plots/cgp_23_fullscale.pdf"))
partition_order <- as.character(unique(anno$Agglo_Cluster[order(anno$ID)]))
partition_order <- partition_order[!(partition_order %in% c("21", "22", "23"))]
reduced_compiled_W <- compiled_W[, intersect_genes]
dim(reduced_compiled_W)
reduced_compiled_W <- reduced_compiled_W[names(factors_ls_iter), ]
dist_W <- MatrixCosSim(reduced_compiled_W, reduced_compiled_W)
factor_umap <- uwot::umap(dist_W, n_neighbors = 20, spread = 5, min_dist = 0.5)
factor_umap <- as.data.frame(factor_umap)
factor_umap$id <- rownames(dist_W)
dataset_names <- stringr::str_match(pattern = "^(\\d{1,2})[[:punct:]](.*)$", string = factor_umap$id)
factor_umap$dataset <- dataset_names[, 2, drop = TRUE]
factor_umap$dataset <- as.character(factor_umap$dataset)
factor_umap$orig_cluster <- dataset_names[, 3, drop = TRUE]
colnames(factor_umap) <- c("UMAP1", "UMAP2", "id", "factor_k", "dataset")

plot <- PlotFactorUMAP(factor_umap, "dataset", label = FALSE); plot
SavePlot(filename = glue("plots/factor_bydataset_umap_subset.pdf"), plot = plot,
  base_height = 6, base_asp = 1)

anno$id <- paste0(anno$Dataset, "_", anno$Orig_Factor)
factor_umap <- merge(factor_umap, anno, by.x = "id", by.y = "ID", all.x = TRUE)
factor_umap$Agglo_Cluster <- as.character(factor_umap$Agglo_Cluster)
#factor_umap$group <- fct_explicit_na(factor_umap$group, na_level = "(Filtered)")
plot <- PlotFactorUMAP(factor_umap, fill = "Agglo_Cluster"); plot
SavePlot(filename = "plots/factorumap_byjaccalgo_plot_subset.pdf", plot = plot,
  base_height = 6, base_asp = 1)

partitions <- split(as.character(anno$ID), as.character(anno$Agglo_Cluster))
sum_partitions <- map(partitions, ~ {
  p <- as.character(.x)
  num_p <- length(p)
  pfactors <- factors_tab[p, ]
  pfactors <- colSums(pfactors)
  pfactors <- pfactors[pfactors > 0]
  return(names(pfactors[pfactors > max(c(floor(num_p*0.33), 2))]))
})
cgp <- qdapTools::list2df(sum_partitions)
cgp <- cgp %>% filter(!(X2 %in% c("21", "22", "23")))

# [7] Plot cluster similarities ----

jacc_mtx <- d_mtx
rownames(anno) <- anno$ID
partitions_split <- split(anno[, c("Dataset", "ID", "Orig_Factor", "Agglo_Cluster")], anno$Agglo_Cluster)
partition_stats <- map_df(.x = names(partitions_split), .f = ~ {
  subset_jacc_mtx <- jacc_mtx[rownames(partitions_split[[.x]]), rownames(partitions_split[[.x]])]
  subset_jacc_mtx <- subset_jacc_mtx[lower.tri(subset_jacc_mtx)]
  avg_in <- mean(subset_jacc_mtx)
  med_in <- median(subset_jacc_mtx)
  sd_in <- sd(subset_jacc_mtx)
  mad_in <- mad(subset_jacc_mtx)

  subrow_jmtx <- jacc_mtx[rownames(partitions_split[[.x]]), !(rownames(jacc_mtx)  %in% rownames(partitions_split[[.x]]))]
  subcol_jmtx <- jacc_mtx[!(rownames(jacc_mtx)  %in% rownames(partitions_split[[.x]])), rownames(partitions_split[[.x]])]
  out_jmtx <- c(subrow_jmtx, subcol_jmtx)
  avg_out <- mean(out_jmtx)
  med_out <- median(out_jmtx)
  sd_out <- sd(out_jmtx)
  mad_out <- mad(out_jmtx)

  return(data.frame(avg_in = avg_in, med_in = med_in, sd_in = sd_in, mad_in = mad_in,
    avg_out = avg_out, med_out = med_out, sd_out = sd_out, mad_out = mad_out, cluster = .x))
})

partition_stats$cluster <- as.numeric(partition_stats$cluster)
d <- ggplot(partition_stats) +
  geom_point(mapping = aes(x = avg_in, y = cluster), shape = 21, fill = "#009E73", alpha = 0.6, color = "black", size = 3) +
  geom_point(mapping = aes(x = avg_out, y = cluster), shape = 21, fill = "#D55E00", alpha = 0.6, color = "black", size = 3) +
  geom_segment(mapping = aes(y = cluster, yend = cluster, xend = avg_out, x = avg_in)) +
  scale_y_continuous(breaks = seq(1, 21, 2)) +
  scale_x_continuous(limits = c(0, 1.0), breaks = seq(0, 1, 0.5)) +
  labs(x = "Partition", y = "Mean bootstrap similarities") +
  GeneralTheme(14) +
  theme(panel.grid.major.y = element_line(size = 1, color = "gray95", linetype = "dotted"))
d
SavePlot(filename = "plots/cgp_clustering_similarities.png", plot = d, base_height = 4, base_asp = 0.8)

# [8] Calculate diversity of consensus factors ----

group_sums <- factor_umap %>% filter(!is.na(Agglo_Cluster)) %>% group_by(Agglo_Cluster, dataset) %>% summarize(n = n())
dcasted <- reshape2::dcast(dataset ~ Agglo_Cluster, data = group_sums, length)
shannon <- vegan::diversity(dcasted[, -1], MARGIN = 2, index = "shannon")

stats <- data.frame(cluster = names(shannon), Shannon = shannon) #Simpson = simpson, Brillouin = brillouin_v)
stats <- gather(stats, key = "key", value = "value", -cluster)
stats$cluster <- factor(stats$cluster, levels = as.character(unique(sort(as.numeric(stats$cluster)))))
statplot <- ggplot(stats,
  aes(x = cluster, y = value)) +
  geom_point(shape = 21, fill = "#0072B2", color = "black", size = 3, alpha = 1) +
  labs(x = "Partition", y = "Shannon\nDiversity") +
  GeneralTheme(14)
statplot

dataset_labels <- gsub("_", " ", object_names)
dataset_labels <- snakecase::to_title_case(dataset_labels)

dcasted_gather <- gather(dcasted, key = "cluster", value = "sum", -dataset)
dcasted_gather$cluster <- factor(dcasted_gather$cluster, levels = as.character(unique(sort(as.numeric(dcasted_gather$cluster)))))
div_tile <- ggplot(dcasted_gather, aes(x = cluster, y = dataset, fill = sum)) +
  geom_tile(color = "black", size = 0.5) +
  scale_y_discrete(labels = dataset_labels, expand = c(0, 0)) +
  scale_fill_gradient(low = "white", high = "#0072B2",  name = "# datasets", breaks = c(0, 1), limits = c(0, 1)) +
  labs(x = "Partition", y = "Dataset") +
  TileTheme() + NoLegend()

div_plot <- cowplot::plot_grid(statplot, div_tile, align = c("v"), axis = "lr", nrow = 2, rel_heights = c(0.5, 1))
div_plot
SavePlot(filename = glue("plots/cgp_diversity_23.pdf"), plot = div_plot, base_asp = 1.1, base_height = 5)

# check for partition overlap
partitions_tab <- qdapTools::mtabulate(sum_partitions)
d <- (ade4::dist.binary(partitions_tab, method = 1, diag = FALSE, upper = FALSE))^2
cons_dmtx <- as.matrix(d)
range(cons_dmtx)
ddf <- as.data.frame(cons_dmtx)
pheatmap::pheatmap(cons_dmtx)

ddf$c1 <- rownames(ddf)
ddf <- ddf %>% gather(key = "c2", value = "x", -c1)
ddf$c1 <- factor(ddf$c1, levels = sort(as.numeric(unique(ddf$c1))))
ddf$c2 <- factor(ddf$c2, levels = sort(as.numeric(unique(ddf$c2))))
blues <- colorspace::sequential_hcl("Blues 2", n = length(c(0.1, seq(0.6, 1, 0.05))))
blues[length(blues)] <- "#FFFFFF"
colors <- circlize::colorRamp2(c(0.1, seq(0.6, 1, 0.05)), blues)
partition_comp_plot <- ggplot(ddf, aes(x = c1, c2, fill = x)) +
  geom_tile(color = "black", size = 0.5) +
  scale_x_discrete(breaks = c(1, seq(5, 30, 5), 23)) +
  scale_y_discrete(breaks = c(1, seq(5, 30, 5), 23)) +
  scale_color_manual(values = c("black", "#0072B2"), name = "Maximum\nmatch") +
  scale_fill_gradient2(high = "#FFFFFF", name = "Jaccard\ndistance", midpoint = 0.5,
    mid = colorspace::sequential_hcl(n = 10, palette = "Blues 2")[2],
    low = colorspace::sequential_hcl(n = 10, palette = "Blues 2")[1]) +
  guides(size = FALSE) +
  labs(x = "Partitions defined by unique DEGs", y = "Paritions defined by non-unique DEGs") +
  TileTheme(base_size = 18)
# geom_tile(color = "black", size = 0.5) +
# scale_x_discrete(breaks = c(1, seq(5, 30, 5), 33)) +
# scale_y_discrete(breaks = c(1, seq(5, 30, 5), 33)) +
# scale_fill_gradientn(colors = blues, name = "Jaccard\nDistance", breaks = c(0.25, 0.5, 0.75, 1), limits = c(0, 1)) +
# labs(x = "Factor partition", y = "Factor partition") +
# theme_classic(base_size = 18) +
# ggeasy::easy_all_text_color("black") +
# theme(
#   axis.title.y = element_text(angle = 90, vjust = 0.5, margin = margin(0, 8, 0, 0), face = "bold"),
#   axis.title.x = element_text(face = "bold"),
#   axis.text.y = element_text(margin = margin(0, 4, 0, 0)),
#   axis.ticks = element_blank(),
#   axis.line = element_blank(),
#   panel.border = element_blank(),
#   #panel.border = element_rect(size = 1, color = "black", fill = "transparent"),
#   legend.text = element_text(size = 10),
#   legend.justification = "top",
#   legend.key.size = unit(1, "line")
# )
partition_comp_plot
cowplot::save_plot(filename = glue("plots/factor_{k}_partition_comp.png"), plot = partition_comp_plot,
  base_height = 8, base_asp = 1.2)

# [7] Compile consensus loadings ----

# pull variable of interest from each marker file and merge
medW <- map_dfr(partitions, function(x) {
  subset_W <- reduced_compiled_W[x, ]
  sum_W <- matrixStats::colMedians(as.matrix(subset_W))
  return(sum_W)
})
rownames(medW) <- intersect_genes
medW <- as.matrix(medW)
genes_to_use <- intersect(rownames(medW), unique(cgp$X1))
medW <- medW[genes_to_use, ]
map <- ComplexHeatmap::pheatmap(cor(medW, method = "pearson"), border_color = NA, color = rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")),
  show_colnames = FALSE, fontsize = 14, treeheight_col = FALSE)
map
pdf(glue("plots/cgp_loading_corr.pdf"), width = 6, height = 5)
draw(map,
  heatmap_legend_side = "right",
  use_raster = TRUE, raster_quality = 5
)
dev.off()

# [8] Plot median loadings ----
# set matrix colors

colfunc <- colorRampPalette(c("white", "#009E73"))
voi_colors <- circlize::colorRamp2(seq(0, 1, 0.1), colfunc(11))

medW <- medW[, match(as.numeric(partition_order), colnames(medW))]
#medW <- medW[match(unique(cgp$X1), rownames(medW)), ]
colnames(medW) <- gsub("K_", "", colnames(medW))

top_genes <- genes_to_label <- apply(medW, 2, function(x) {
  names(sort(x, decreasing = TRUE))[1:20]
})

genes_to_label <- apply(medW, 2, function(x) {
  names(sort(x, decreasing = TRUE))[1:3]
})
genes_to_label[, 5] <- c("IRF7", "TCF4", "PLAC8")
genes_to_label[, 6] <- c("INHBA", "CD101", "FFAR4")

plot_indices <- match(genes_to_label, rownames(medW))
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

voi <- Heatmap(
  medW,
  border = TRUE,
  col = voi_colors,
  cluster_rows = TRUE, cluster_columns = FALSE,
  show_column_names = TRUE, show_column_dend = FALSE,
  column_names_side = "top",column_names_rot = 0, column_names_centered = TRUE,
  column_title = "", column_title_side = "bottom",
  show_row_names = FALSE, show_row_dend = FALSE,
  left_annotation = genes_anno,
  #top_annotation = cluster,
  # top_annotation = HeatmapAnnotation(
  #   Cluster = anno_block(
  #     gp = gpar(fill = colorspace::qualitative_hcl("Dark 3", n = k)),
  #     labels = levels(gene_umap$cluster),
  #     labels_gp = gpar(fontsize = 8, col = "black", fontface = "bold")
  #   )),
  #column_split = 15,
  #row_title = NULL, #row_gap = unit(0, "mm"), column_gap = unit(0, "mm"),
  heatmap_legend_param = list(
    title = "Median\nNormalized\nloading\n ",
    at = c(0, 0.25, 0.5, 0.75, 1.0),
    border = "black",
    legend_height = unit(4, "cm"),
    legend_width = unit(1, "cm")
    #direction = "horizontal",
    #title_position = "lefttop"
  ),
  height = 4,
  use_raster = TRUE, raster_quality = 1
)
voi

dev.off()
png(glue("plots/cgp_23_top3.png"), width = 10, height = 6, units = "in", res = 400)
pdf(glue("plots/cgp_23_top3.pdf"), width = 8, height = 8)
draw(voi,
  heatmap_legend_side = "right",
  use_raster = FALSE, raster_quality = 10
)
dev.off()

# [9] Save files ----

colnames(cgp) <- c("Gene", "CGP")
saveRDS(cgp, file = "data/interim/CGPs_k23.rds")
write_csv(x = cgp, path = "data/publication/TableS7.csv")

dim(medW)
saveRDS(medW, file = "data/interim/CGPs_k23_loadings.rds")

factors_df <- qdapTools::list2df(factors_ls)
colnames(factors_df) <- c("Gene", "Program_Dataset")
write_csv(x = factors_df, path = "data/publication/TableS6.csv")

# [10] Compare factors to markers ----

stm_results <- readRDS("data/interim/ccs_markers.rds")
stm_results <- stm_results %>% filter(!(partition %in% c(6, 14, 19, 21)))

gom <- GeneOverlap::newGOM(
  gsetA = split(cgp$Gene, cgp$CGP),
  gsetB = split(stm_results$feature, stm_results$partition),
  genome.size = 1E4
)
odds <- GeneOverlap::getMatrix(gom, "odds.ratio")
jacc <- GeneOverlap::getMatrix(gom, "Jaccard")
pvals <- GeneOverlap::getMatrix(gom, "pval")
df <- ConvertGOMatrix(pvals, jacc, odds,
  comp1 = "Programs", comp2 = "Markers")

df$set1 <- factor(df$set1, levels = as.character(unique(sort(as.numeric(df$set1)))))
df$set2 <- factor(df$set2, levels = as.character(unique(sort(as.numeric(df$set2)))))
partition_comp_plot <- ggplot(df, aes(x = set1, y = set2, size = log_adj_pvals, fill = jacc)) +
  geom_point(shape = 21, stroke = 0.8) +
  scale_size_continuous(range = c(2, 8), name = "-Log10(\nFDR-adj.\nP value)") +
  colorspace::scale_fill_continuous_sequential("Blues", rev = 1, name = "Jaccard\nIndex") +
  labs(x = "Gene expression programs", y = "Gene markers") +
  GeneralTheme(14)
partition_comp_plot
cowplot::save_plot(filename = "plots/stm_factor_comp.png", plot = partition_comp_plot,
  base_height = 6, base_asp = 1.2)

# [11] Enrichment ----
names(sum_partitions)[13]
enr_partitions <- sum_partitions[-c(14:16)]
enr_partitions <- map(enr_partitions, ~ {
  return(.x[!grepl("(MT-|RPS|RPL)", .x)])
})

library(enrichR)
dbs <- listEnrichrDbs()
dbs <- c("Reactome_2016")

enr <- map_dfr(names(enr_partitions), .f = ~ {
  enr <- enrichr(enr_partitions[[.x]], databases = dbs)
  enr <- lapply(names(enr), function(x) {
    enr[[x]]$db <- x
    return(enr[[x]])
  })
  enr <- bind_rows(enr)
  enr$geneset <- .x
  return(enr)
})
tophits <- enr %>% group_by(geneset) %>% filter(Adjusted.P.value < 0.1) %>% top_n(-5, Adjusted.P.value)

tophits$log10padj <- -log10(tophits$Adjusted.P.value)
tophits$log10padj[tophits$log10padj > 50] <- 50
tophits <- tophits %>% arrange(desc(log10padj))
tophits$geneset <- factor(tophits$geneset,
  levels = partition_order)

enrplot <- ggplot(tophits, aes(x = geneset, y = Term, fill = log10padj)) +
  geom_tile(color = "black", size = 0.5) +
  scale_fill_gradient(low = "white", high = "#0072B2", name = "-log10\nadj. P value") +
  labs(x = "Consensus gene programs", y = "", title = "Reactome") +
  GeneralTheme(14) + theme(legend.position = "right") +
  theme(axis.text.y = element_blank(),
    panel.grid.major = element_line(size = 0.5, color = "gray90"),
    panel.grid.minor = element_line(size = 0.5, color = "gray90"))
enrplot
SavePlot(plot = enrplot, filename = "plots/cgp_reactome.png", base_height = 6, base_asp = 1)

enrplot <- ggplot(tophits, aes(x = geneset, y = Term, fill = log10padj)) +
  geom_tile(color = "black", size = 0.5) +
  scale_fill_gradient(low = "white", high = "#0072B2", name = "-log10\nadj. P value") +
  labs(x = "Consensus cell states", y = "", title = "Reactome") +
  GeneralTheme(14) + theme(legend.position = "right")
SavePlot(plot = enrplot, filename = "plots/cgp_reactome.pdf", base_height = 6, base_asp = 2.5)




