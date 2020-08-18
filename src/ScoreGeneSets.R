# ---
# Description: Score and plot gene set trends across datasets
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

sets <- readxl::read_excel("data/study_metadata.xlsx")
ids <- sets$ID[sets$Train == TRUE]
ids <- SelfName(ids)
ids <- ids[order(ids)]

# [2] Load objects, markers, and programs ----

object_ls <- LoadObjectList(dir.path = "/Users/jpeters/Projects/lungmps/data/objects",
  file_pattern = "_mnp_final.rds$", name_pattern = "_mnp_final.rds$", exclude = "Liao_2020", verbose = 1)
object_names <- names(object_ls)
names(object_names) <- object_names
object_ls <- lapply(object_names, function(x) {
  object_ls[[x]]$set <- x
  return(object_ls[[x]])
})

cdf <- data.frame(set = object_names)
cdf$source <- gsub("(_disease)|(_health)", "", cdf$set)
cdf$category <- str_match(pattern = "(_disease|_health)", string = cdf$set)[, 2]
cdf$category <- gsub("_", "", cdf$category)
cdf$response <- as.numeric(factor(cdf$category, levels = c("health", "disease")))-1
cdf <- cdf %>% group_by(source) %>% mutate(compare = ifelse(length(unique(category)) == 2, TRUE, FALSE))
cdf$use <- TRUE
cdf$use[cdf$source %in% c("Bost_2020", "KimLee_2020")] <- FALSE
cdf$compare[cdf$source %in% c("Bost_2020", "KimLee_2020")] <- FALSE

markers <- readRDS("data/interim/ccs_markers.rds")
markers_ls <- split(markers$feature, markers$partition)

#! write function to remove elements matching a string vector
orig_length <- length(markers_ls)
remove <- as.character(c(6, 14, 19, 21))
remove <- paste0(remove, collapse = "|")
remove <- paste0("(", remove, ")")
stopifnot(any(grepl(remove, names(markers_ls))))
markers_ls <- markers_ls[!grepl(remove, names(markers_ls))]
stopifnot(orig_length == length(markers_ls))

#! write function to reorder factor numerically
names(markers_ls)[order(as.character(names(markers_ls)))]
paste(sort(as.integer(names(markers_ls))))
gtools::mixedsort(names(markers_ls))
markers_ls <- markers_ls[str_sort(names(markers_ls), numeric = TRUE)]
markers_ls <- markers_ls[sapply(markers_ls, length) >= 20]
names(markers_ls) <- paste0("M", names(markers_ls))
markers <- names(markers_ls)
markers <- SelfName(markers)

programs <- readRDS("data/interim/CGPs_k23.rds")
head(programs)
programs_ls <- split(programs$Gene, programs$CGP)
programs_ls <- programs_ls[sapply(programs_ls, length) >= 20]
programs_ls <- programs_ls[str_sort(names(programs_ls), numeric = TRUE)]
names(programs_ls) <- paste0("P", names(programs_ls))
programs <- names(programs_ls)
programs <- SelfName(programs)

all_markers <- unlist(markers_ls, use.names = FALSE)
all_markers <- names(table(all_markers))[(table(all_markers) == 1)]
marker_gene_name <- map(markers_ls, ~ { .x[(.x %in% all_markers)][1] })

all_programs <- unlist(programs_ls, use.names = FALSE)
all_programs <- names(table(all_programs))[(table(all_programs) == 1)]

loadings <- readRDS("data/interim/CGPs_k23_loadings.rds")
dim(loadings)
rownames(loadings)
loadings <- loadings[intersect(all_programs, rownames(loadings)), ]

program_gene_names <- c()
for (i in names(programs_ls)) {
  prog_num <- gsub("P", "", i)
  subset_genes <- programs_ls[[i]][programs_ls[[i]] %in% rownames(loadings)]
  subset_loadings <- loadings[subset_genes, prog_num]
  print(i)
  print(names(sort(subset_loadings, decreasing = TRUE))[1:10])
  program_gene_names <- c(program_gene_names, names(sort(subset_loadings, decreasing = TRUE))[1])
}
program_gene_names
names(program_gene_names) <- names(programs_ls)
program_gene_names[3] <- "INHBA"
program_gene_names[4] <- "TCF4"
program_gene_names[6] <- "STMN1"
program_gene_names[7] <- "LILRB2"
program_gene_names[12] <- "HLA-DQB2"
program_gene_names[15] <- "MARCO"
program_gene_names[20] <- "CPVL"

# [3] Score ----

sources <- unique(cdf$source[cdf$compare])
sources <- SelfName(sources)

scaled_marker_scores <- pbapply::pblapply(sources, function(x) {
  df <- ScoreObjectScaled(object_ls = object_ls, set_df = cdf, set_source = x,
    genesets = markers_ls, names = markers, group = "merged_leiden")
  return(df)
})
scaled_marker_scores <- bind_rows(scaled_marker_scores)
saveRDS(object = scaled_marker_scores, file = "data/interim/scaled_marker_scores.rds")
scaled_marker_scores <- readRDS("data/interim/scaled_marker_scores.rds")

scaled_program_scores <- pbapply::pblapply(sources, function(x) {
  df <- ScoreObjectScaled(object_ls = object_ls, set_df = cdf, set_source = x,
    genesets = programs_ls, names = programs, group = "merged_leiden")
  return(df)
})
scaled_program_scores <- bind_rows(scaled_program_scores)
saveRDS(object = scaled_program_scores, file = "data/interim/scaled_program_scores.rds")
scaled_program_scores <- readRDS("data/interim/scaled_program_scores.rds")

# [4] Build score models ----

sources
markers
iterations <- expand.grid(sources, markers)
nrow(iterations)
head(iterations)

StopFutureLapply()
StartFutureLapply()
scaled_scores <- scaled_marker_scores
scaled_scores$disease_var <- 0
scaled_scores$disease_var[scaled_scores$disease_classification == "disease"] <- 1
scaled_scores$set_cluster <- paste0(scaled_scores$set, "_", scaled_scores$group)
unique(scaled_marker_scores$set[is.nan(scaled_scores$score)])
head(scaled_scores)

results <- furrr::future_map(1:nrow(iterations), TestScores,
  iterations = iterations, cdf = cdf, scaled_scores = scaled_scores)
StopFutureLapply()

results <- lapply(results, function(x) {
  x$model$data <- NULL
  attr(x$model$terms, ".Environment") <- NULL
  x$model$modelStruct <- NULL
  return(x)
})

marker_results <- results
names(marker_results) <- paste(iterations$Var1, iterations$Var2, sep = "_")
marker_results <- saveRDS(object = marker_results, file =  "data/interim/models_markers_linear_top2.rds")
marker_results <- readRDS(file =  "data/interim/models_markers_linear_top2.rds")
results <- marker_results

iterations <- expand.grid(sources, markers)
names(results) <- paste(iterations$Var1, iterations$Var2, sep = "_")
Clean()

estimates <- map_dbl(marker_results, ~ { unique(.x$model_tidy$estimate[.x$model_tidy$term == "disease_var"]) })
model_results <- map_dfr(marker_results, "pvalue")
model_results$estimate <- estimates

model_results$bh_plinear <- p.adjust(model_results$plinear, method = "fdr")
model_results$nl_linearadjp <- -log10(model_results$bh_plinear)
model_results$nl_linearadjp[model_results$nl_linearadjp > -log10(2.2E-16)] <- 16
model_results$signif <- ifelse(model_results$bh_plinear <= 0.05, TRUE, FALSE)
marker_model_plot <- ggplot(model_results, aes(x = program, y = source, fill = estimate, size = nl_linearadjp)) +
  geom_point(shape = 21) +
  geom_point(data = model_results %>% filter(signif == TRUE),
    aes(x = program, y = source), shape = 4, stroke = 1, color = "black", size = 1, fill = "transparent") +
  scale_x_discrete(labels = paste(levels(model_results$program), marker_gene_name, sep = " ")) +
  scale_y_discrete(labels = gsub("_", " ", unique(model_results$source))) +
  scale_size_continuous(range = c(2, 8), name = "-log10(\n P value)") +
  colorspace::scale_fill_continuous_diverging(palette = "Blue-Red 3", name = "Model coefficient") +
  labs(x = "Program", y = "Dataset", title = "Consensus marker gene sets") +
  GeneralTheme(18) + theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.box = "top",
    legend.direction = "vertical",
    panel.grid.major = element_line(size = 0.5, color = "gray90"),
    panel.grid.minor = element_line(size = 0.5, color = "gray90")
  )

marker_model_plot
SavePlot(filename = "plots/linreg_markers_top2.pdf", plot = marker_model_plot, base_asp = 1.657895*2, base_height = 1.9*2.5)
SavePlot(filename = "plots/linreg_markers_top2.png", plot = marker_model_plot, base_asp = 1.657895*2, base_height = 1.9*2.5)

sources
programs
iterations <- expand.grid(sources, programs)
nrow(iterations)
head(iterations)

StopFutureLapply()
StartFutureLapply()
scaled_scores <- scaled_program_scores
scaled_scores$disease_var <- 0
scaled_scores$disease_var[scaled_scores$disease_classification == "disease"] <- 1
scaled_scores$set_cluster <- paste0(scaled_scores$set, "_", scaled_scores$group)
head(scaled_scores)

results <- furrr::future_map(1:nrow(iterations), TestScores,
  iterations = iterations, cdf = cdf, scaled_scores = scaled_scores)
StopFutureLapply()

results <- lapply(results, function(x) {
  x$model$data <- NULL
  attr(x$model$terms, ".Environment") <- NULL
  x$model$modelStruct <- NULL
  return(x)
})

program_results <- results
names(program_results) <- paste(iterations$Var1, iterations$Var2, sep = "_")
program_results  <- saveRDS(object = program_results, file =  "data/interim/models_programs_linear_top2.rds")
program_results <- readRDS(file =  "data/interim/models_programs_linear_top2.rds")
names(program_results)
results <- program_results

names(results) <- paste(iterations$Var1, iterations$Var2, sep = "_")
Clean()

estimates <- map_dbl(program_results, ~ { unique(.x$model_tidy$estimate[.x$model_tidy$term == "disease_var"]) })
model_results <- map_dfr(program_results, "pvalue")
model_results$estimate <- estimates

model_results$bh_plinear <- p.adjust(model_results$plinear, method = "fdr")
model_results$nl_linearadjp <- -log10(model_results$bh_plinear)
model_results$nl_linearadjp[model_results$nl_linearadjp > -log10(2.2E-16)] <- 16
model_results$signif <- ifelse(model_results$bh_plinear <= 0.05, TRUE, FALSE)
program_model_plot <- ggplot(model_results, aes(x = program, y = source, fill = estimate, size = nl_linearadjp)) +
  geom_point(shape = 21) +
  geom_point(data = model_results %>% filter(signif == TRUE),
    aes(x = program, y = source), shape = 4, stroke = 1, color = "black", size = 1, fill = "transparent") +
  scale_x_discrete(labels = paste(levels(model_results$program), program_gene_names[levels(model_results$program)], sep = " ")) +
  scale_y_discrete(labels = gsub("_", " ", unique(model_results$source))) +
  scale_size_continuous(range = c(3, 8), name = "-log10(\n P value)") +
  colorspace::scale_fill_continuous_diverging(palette = "Blue-Red 3", name = "Model coefficient") +
  labs(x = "Program", y = "Dataset", title = "Consensus program gene sets") +
  GeneralTheme(18) + theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.box = "top",
    legend.direction = "vertical",
    panel.grid.major = element_line(size = 0.5, color = "gray90"),
    panel.grid.minor = element_line(size = 0.5, color = "gray90")
  )
program_model_plot
SavePlot(filename = "plots/linreg_programs_top2.pdf", plot = program_model_plot, base_asp = 1.657895*2, base_height = 1.9*2.5)
SavePlot(filename = "plots/linreg_programs_top2.png", plot = program_model_plot, base_asp = 1.657895*2, base_height = 1.9*2.5)

# plot together using patchwork
model_plots <- (marker_model_plot + labs(x = "Consensus markers", title = "")) / (program_model_plot + labs(x = "Consensus programs", title = "")) +
  plot_layout(guides = "collect") &
  colorspace::scale_fill_continuous_diverging(palette = "Blue-Red 3", name = "Model coefficient", limits = c(-2, 2), breaks = c(-1, 0, 1)) &
  scale_size_continuous(range = c(3, 8), name = "-log10(\n P value)", limits = c(0, 20), breaks = c(2, 6, 10, 14)) &
  theme(legend.justification = "left", legend.box.just = "top")
ggsave(filename = "plots/linreg_model_plots.pdf", plot = model_plots, height = 1, width = 1.67, scale = 8, dpi = 400, units = "in")
ggsave(filename = "plots/linreg_model_plots.png", plot = model_plots, height = 1, width = 1.67, scale = 8, dpi = 400, units = "in")

# [5] Plot marker trends ----
num_top = 2
geneset_use = "M7"
head(scaled_marker_scores)
plot_scores <- scaled_marker_scores %>% filter(geneset == geneset_use)
plot_scores$set_cluster <- paste0(plot_scores$set, "_", plot_scores$group)
passing_clusters <- plot_scores %>%
  group_by(set, group) %>%
  summarize(mean = median(score), set_cluster = unique(set_cluster)) %>%
  top_n(num_top, mean) %>%
  select(set_cluster)
plot_scores <- plot_scores %>% filter(set_cluster %in% passing_clusters$set_cluster)
plot_scores <- plot_scores %>% filter(source %in% c("Habermann_2019", "Morse_2019", "Reyfman_2019"))
score_vlnplot <- ggplot(plot_scores,
  aes(x = source, y = score, fill = disease_classification)) +
  geom_jitter(shape = 21, alpha = 0.25, size = 0.2, position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.1)) +
  #geom_boxplot(alpha = 0.75, width = 0.75, position = position_dodge(width = 0.5), outlier.color = "Gray60", outlier.alpha = 0.25) +
  geom_violin(alpha = 0.6, width = 0.75, position = position_dodge(width = 0.5), draw_quantiles = c(0.25, 0.5, 0.75)) +
  scale_fill_manual(values = c("#d55e00", "#0072b2"), name = "Condition", labels = c("Disease", "Health")) +
  scale_x_discrete(labels = gsub("_", " ", unique(plot_scores$source))) +
  scale_y_continuous(limits = c(min(plot_scores$score),
    ceiling(max(plot_scores$score)))) +
  # annotate(geom = "text", label = glue("P = {format.pval(unique(kpis$plogit))}, \u03b2 = {round(unique(kpis$blogit), 2)}"), x = "health",
  #   y = ceiling(max(scores_to_plot$score)), hjust = 0) +
  labs(x = "Dataset", y = "M10 Score") +
  GeneralTheme(18) #+ theme(axis.title.y = element_text(angle = 0, hjust = 0.5, vjust = 1))
score_vlnplot
#SavePlot(filename = "plots/m5_violin_plot.xpdf", plot = m5_violin_plot, base_height = 1.9*2, base_asp = 3.15/1.9)

# check statistics if needed
unique(plot_scores$source)
source_test <- unique(plot_scores$source)[3]
wilcox.test(plot_scores$score[plot_scores$source == source_test & plot_scores$disease_classification == "disease"],
  plot_scores$score[plot_scores$source == source_test & plot_scores$disease_classification == "health"])

names(marker_results) # check names
stats <- imap(marker_results, ~ {
  df <- .x$stats
  df$source_geneset <- .y
  return(df)
})
stats <- bind_rows(stats)
subset_stats <- stats %>% filter(source_geneset %in% c("Habermann_2019_M10", "Morse_2019_M10", "Reyfman_2019_M10"))
percentile_plot <- ggplot(subset_stats,
  aes(x = source_geneset, y = diff)) +
  geom_jitter(alpha = 0.1, width = 0.15) +
  #tidybayes::stat_eye(side = "right") +
  #tidybayes::stat_dots(scale = 1, width = 1, alpha = 0.8, size = 20) +
  #geom_point() +
  tidybayes::stat_pointinterval(position = position_nudge(x = 0.2)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "Gray90", size = 1) +
  scale_color_brewer(name = "Median QI") +
  scale_y_continuous(limits = c(-0.1, 1.5)) +
  # annotate(geom = "text", label = glue("P = {format.pval(unique(kpis$plogit))}, \u03b2 = {round(unique(kpis$blogit), 2)}"), x = "health",
  #   y = ceiling(max(scores_to_plot$score)), hjust = 0) +
  labs(x = "", y = "90th quantile\n difference") +
  GeneralTheme(18) +
  guides(color = guide_legend(override.aes = list(size = 10))) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks = element_blank()
    #axis.title.y = element_text(angle = 0, hjust = 0.5, vjust = 1)
  )
percentile_plot

pw_b <- (percentile_plot + theme(plot.margin = margin(0,0,0,0))) / score_vlnplot + plot_layout(heights = c(1, 2))
SavePlot(filename = "plots/m10_plot.pdf", pw_b, base_height = 1.9*4, base_asp = 3.15/1.9)
SavePlot(filename = "plots/m10_plot.png", pw_b, base_height = 1.9*4, base_asp = 3.15/1.9)

# [6] Plot program trends ----
num_top = 2
geneset_use = "P9"
head(scaled_program_scores)
plot_scores <- scaled_program_scores %>% filter(geneset == geneset_use)
plot_scores$set_cluster <- paste0(plot_scores$set, "_", plot_scores$group)
passing_clusters <- plot_scores %>%
  group_by(set, group) %>%
  summarize(mean = median(score), set_cluster = unique(set_cluster)) %>%
  top_n(num_top, mean) %>%
  select(set_cluster)
plot_scores <- plot_scores %>% filter(set_cluster %in% passing_clusters$set_cluster)
#plot_scores <- plot_scores %>% filter(source %in% c("Habermann_2019", "Morse_2019", "Reyfman_2019"))
score_vlnplot <- ggplot(plot_scores,
  aes(x = source, y = score, fill = disease_classification)) +
  geom_jitter(shape = 21, alpha = 0.25, size = 0.2, position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.1)) +
  #geom_boxplot(alpha = 0.75, width = 0.75, position = position_dodge(width = 0.5), outlier.color = "Gray60", outlier.alpha = 0.25) +
  geom_violin(alpha = 0.6, width = 0.75, position = position_dodge(width = 0.5), draw_quantiles = c(0.25, 0.5, 0.75)) +
  scale_fill_manual(values = c("#d55e00", "#0072b2"), name = "Condition", labels = c("Disease", "Health")) +
  scale_x_discrete(labels = gsub("_", " ", unique(plot_scores$source))) +
  scale_y_continuous(limits = c(min(plot_scores$score),
    ceiling(max(plot_scores$score)))) +
  # annotate(geom = "text", label = glue("P = {format.pval(unique(kpis$plogit))}, \u03b2 = {round(unique(kpis$blogit), 2)}"), x = "health",
  #   y = ceiling(max(scores_to_plot$score)), hjust = 0) +
  labs(x = "Dataset", y = "P9 Score") +
  GeneralTheme(18) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)
  )
score_vlnplot
#SavePlot(filename = "plots/m5_violin_plot.xpdf", plot = m5_violin_plot, base_height = 1.9*2, base_asp = 3.15/1.9)

# check statistics if needed
unique(plot_scores$source)
source_test <- unique(plot_scores$source)[1]
wilcox.test(plot_scores$score[plot_scores$source == source_test & plot_scores$disease_classification == "disease"],
  plot_scores$score[plot_scores$source == source_test & plot_scores$disease_classification == "health"])
mean(plot_scores$score[plot_scores$source == source_test & plot_scores$disease_classification == "disease"])
mean(plot_scores$score[plot_scores$source == source_test & plot_scores$disease_classification == "health"])

names(program_results) # check names
stats <- imap(program_results, ~ {
  df <- .x$stats
  df$source_geneset <- .y
  return(df)
})
stats <- bind_rows(stats)
subset_stats <- stats %>% filter(source_geneset %in% paste0(sources, "_P9"))
percentile_plot <- ggplot(subset_stats,
  aes(x = source_geneset, y = diff)) +
  geom_jitter(alpha = 0.1, width = 0.15) +
  #tidybayes::stat_eye(side = "right") +
  #tidybayes::stat_dots(scale = 1, width = 1, alpha = 0.8, size = 20) +
  #geom_point() +
  tidybayes::stat_pointinterval(position = position_nudge(x = 0.2)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "Gray90", size = 1) +
  scale_color_brewer(name = "Median QI") +
  scale_y_continuous(limits = c(-2, 2)) +
  # annotate(geom = "text", label = glue("P = {format.pval(unique(kpis$plogit))}, \u03b2 = {round(unique(kpis$blogit), 2)}"), x = "health",
  #   y = ceiling(max(scores_to_plot$score)), hjust = 0) +
  labs(x = "", y = "90th quantile\n difference") +
  GeneralTheme(18) +
  guides(color = guide_legend(override.aes = list(size = 10))) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks = element_blank()
    #axis.title.y = element_text(angle = 0, hjust = 0.5, vjust = 1)
  )
percentile_plot

pw_d <- (percentile_plot + theme(plot.margin = margin(0,0,0,0))) / score_vlnplot + plot_layout(heights = c(1, 2))
pw_d
SavePlot(filename = "plots/p9_plot.pdf", pw_d, base_height = 1.9*4, base_asp = 3.15/1.9)
SavePlot(filename = "plots/p9_plot.png", pw_d, base_height = 1.9*4, base_asp = 3.15/1.9)

# [7] Identify small population ----

plotlist <- list()
for (i in unique(cdf$set[cdf$compare == TRUE])) {
  ui_info("Plotting {i}")
  object <- object_ls[[i]]
  # extract DC likelihood and merge with cross-presenting DC score
  # plot 2d density of that
  dc_ll <- colSums(object@misc$is_llpp[17:18, ])
  dc_ll[is.nan(dc_ll)] <- 0
  dc_ll <- qdapTools::vect2df(dc_ll, col1 = "cell_barcode", col2 = "dc_ll", order = FALSE)
  ms <- scaled_program_scores %>% filter(set == i & geneset == "P17")
  ms$cell_barcode <- gsub("(d|h)_", "", ms$cell_barcode)
  ms <- merge(ms, dc_ll, all.x = TRUE, by = "cell_barcode")
  msf <- ms %>% filter(dc_ll > 0.0001)
  msf$density <- GetDensity(msf$dc_ll, msf$score, n = 200)
  pd <- ggplot(msf, aes(x = dc_ll, y = score, color = density)) +
    geom_point() +
    colorspace::scale_color_continuous_sequential("Blues 3", rev = 0) +
    GeneralTheme(18) +
    NoLegend() +
    labs(x  = "DC Likelihood", y = "Score", title = snakecase::to_title_case(gsub("_", " ", i))) +
    theme(
      panel.grid.major = element_line(color = "Gray95", size = 0.5)
    )
  # add score
  msa <- ms %>% select(score)
  colnames(msa) <- c("meta_score_to_plot")
  rownames(msa) <- ms$cell_barcode
  object <- AddMetaData(object, metadata = msa)
  # make hexbin
  object <- schex::make_hexbin(object, nbins = 48, dimension_reduction = "harmony_umap")
  # plot hexbin
  hxmeta <- schex::plot_hexbin_meta(object, col = "meta_score_to_plot", action = "median")
  hxmeta <- hxmeta + colorspace::scale_fill_continuous_sequential("Blues 3", name = "Score") +
    GeneralTheme(18) +
    labs(title = NULL, x = "UMAP1", y = "UMAP2") +
    theme(axis.ticks = element_blank(), axis.text = element_blank())
  grid <- cowplot::plot_grid(pd, hxmeta, nrow = 1, rel_widths = c(0.8, 1), align = "h")
  plotlist[[i]] <- grid
}
all_grids <- cowplot::plot_grid(plotlist = plotlist, nrow = 6, ncol = 2)
SavePlot(plot = all_grids, filename = "plots/cpdcs_fig4.pdf", base_asp = 0.8, base_height = 24)
SavePlot(plot = all_grids, filename = "plots/cpdcs_fig4.png", base_asp = 0.8, base_height = 24)

# [8] Program correlation plot ----

cors <- pbapply::pblapply(cdf$set[cdf$compare], function(x) {
  ui_info("Plotting {x}")
  object <- object_ls[[x]]
  # extract DC likelihood and merge with cross-presenting DC score
  # plot 2d density of that
  mon_ll <- colSums(object@misc$is_llpp[13, , drop = FALSE])
  mac_ll <- colSums(object@misc$is_llpp[14:16, , drop = FALSE])
  mon_ll[is.nan(mon_ll)] <- 0
  mac_ll[is.nan(mac_ll)] <- 0
  mon_ll <- qdapTools::vect2df(mon_ll, col1 = "cell_barcode", col2 = "mon_ll", order = FALSE)
  mac_ll <- qdapTools::vect2df(mac_ll, col1 = "cell_barcode", col2 = "mac_ll", order = FALSE)
  dc_ll <- colSums(object@misc$is_llpp[17:18, ])
  dc_ll[is.nan(dc_ll)] <- 0
  dc_ll <- qdapTools::vect2df(dc_ll, col1 = "cell_barcode", col2 = "dc_ll", order = FALSE)
  lls <- cbind(mon_ll, mac_ll %>% select(mac_ll), dc_ll %>% select(dc_ll))

  ms <- scaled_program_scores %>% filter(set == x) %>% select(-percents, -percents_geneset)
  ms$cell_barcode <- gsub("(d|h)_", "", ms$cell_barcode)
  ms <- merge(ms, lls, all.x = TRUE, by = "cell_barcode")
  cors <- ms %>%
    group_by(geneset) %>%
    summarize(
      moncor = cor(score, mon_ll),
      maccor = cor(score, mac_ll),
      dccor = cor(score, dc_ll))
  norm_cors <- as.data.frame(apply(cors %>% select(-geneset), 2, function(x) {
    return((1-0.01)*((x - min(x, na.rm = TRUE)) /(max(x, na.rm = TRUE) - min(x, na.rm = TRUE))) + 0.01)
  }))
  colnames(norm_cors) <- paste0("norm_", colnames(cors %>% select(-geneset)))
  cors <- cbind(cors, norm_cors)
  cors$set <- x
  return(cors)
})
cors <- bind_rows(cors)

cors$geneset<- factor(cors$geneset,
  levels = unique(cors$geneset[order(as.numeric(gsub("P", "", cors$geneset)))]))
sum_cors <- cors %>% group_by(geneset) %>% summarize_all(mean)
a <- ggplot(cors, aes(x = moncor, y = maccor, color = geneset)) +
  geom_point(alpha = 0.5) +
  geom_point(data = sum_cors, shape = 21, size = 4, aes(x = moncor, y = maccor, fill = geneset), color = "black") +
  GeneralTheme(18) +
  scale_x_continuous(limits = c(-0.75, 0.75), breaks = c(-0.5, 0, 0.5), labels = c("-0.5", "0", "0.5")) +
  scale_y_continuous(limits = c(-0.75, 0.75), breaks = c(-0.5, 0, 0.5), labels = c("-0.5", "0", "0.5")) +
  colorspace::scale_fill_discrete_qualitative("Dark 3", labels = gsub("_", "", levels(cors$geneset)), name = "Program") +
  colorspace::scale_color_discrete_qualitative("Dark 3", labels = gsub("_", "", levels(cors$geneset)), name = "Program") +
  guides(color = FALSE, fill = FALSE) +
  labs(x = "Monocyte Correlation", y = "Macrophage Correlation") +
  theme(
    axis.text.x = element_text(color = "black", size = 14),
  )
b <- ggplot(cors, aes(x = moncor, y = dccor, color = geneset)) +
  geom_point(alpha = 0.5) +
  geom_point(data = sum_cors, shape = 21, size = 4, aes(x = moncor, y = dccor, fill = geneset), color = "black") +
  GeneralTheme(18) +
  scale_x_continuous(limits = c(-0.75, 0.75), breaks = c(-0.5, 0, 0.5), labels = c("-0.5", "0", "0.5")) +
  scale_y_continuous(limits = c(-0.75, 0.75), breaks = c(-0.5, 0, 0.5), labels = c("-0.5", "0", "0.5")) +
  colorspace::scale_fill_discrete_qualitative("Dark 3", labels = gsub("_", "", levels(cors$geneset)), name = "Program") +
  colorspace::scale_color_discrete_qualitative("Dark 3", labels = gsub("_", "", levels(cors$geneset)), name = "Program") +
  guides(color = FALSE, fill = FALSE) +
  labs(x = "Monocyte Correlation", y = "DC Correlation") +
  theme(
    axis.text.x = element_text(color = "black", size = 14)
  )
legend <-  cowplot::get_legend(a + guides(fill = guide_legend(nrow = 2)) + theme(legend.position = "bottom"))
ab <- cowplot::plot_grid(a, b, nrow = 1, rel_widths = c(1, 1), hjust = -1, align = "vh")
corr_plot <- cowplot::plot_grid(ab, legend, rel_heights = c(1, 0.1), nrow = 2, ncol = 1)
corr_plot
SavePlot(filename = "plots/programcorr.pdf", plot = corr_plot, base_asp = 2, base_height = 1.9*3)

# [9] Broad correlation plot ----

x = cdf$set[1]
lls <- pbapply::pblapply(cdf$set, function(x) {
  ui_info("Extracting from {x}...")
  object <- object_ls[[x]]

  # extract DC likelihood and merge with cross-presenting DC score
  # plot 2d density of that
  mon_ll <- colSums(object@misc$is_llpp[13, , drop = FALSE])
  mac_ll <- colSums(object@misc$is_llpp[14:16, , drop = FALSE])
  mon_ll[is.nan(mon_ll)] <- 0
  mac_ll[is.nan(mac_ll)] <- 0
  mon_ll <- qdapTools::vect2df(mon_ll, col1 = "cell_barcode", col2 = "mon_ll", order = FALSE)
  mac_ll <- qdapTools::vect2df(mac_ll, col1 = "cell_barcode", col2 = "mac_ll", order = FALSE)
  dc_ll <- colSums(object@misc$is_llpp[17:18, ])
  dc_ll[is.nan(dc_ll)] <- 0
  dc_ll <- qdapTools::vect2df(dc_ll, col1 = "cell_barcode", col2 = "dc_ll", order = FALSE)
  lls <- cbind(mon_ll, mac_ll %>% select(mac_ll), dc_ll %>% select(dc_ll))

  df <- object@meta.data[, c("nCount_RNA", "disease_classification", "set", "merged_leiden", "batch")]
  df <- df %>% rownames_to_column("cell_barcode")
  df <- merge(df, lls, by = "cell_barcode", all.x = TRUE)

  sum_lls <- df %>%
    group_by(merged_leiden) %>%
    summarize(
      moncor = mean(mon_ll),
      maccor = mean(mac_ll),
      dccor = mean(dc_ll))

  sum_lls$set <- x
  return(sum_lls)
})
lls <- bind_rows(lls)

library(patchwork)
a <- ggplot(lls, aes(x = moncor, y = maccor, fill = set)) +
  geom_point(alpha = 0.5, size = 3, shape = 21, color = "black") +
  GeneralTheme(14) +
  colorspace::scale_fill_discrete_qualitative("Dark 3") +
  guides(color = FALSE, fill = FALSE) +
  labs(x = "", y = "Macrophage Likelihood")
b <- ggplot(lls, aes(x = moncor, y = dccor, fill = set)) +
  geom_point(alpha = 0.5, size = 3, shape = 21, color = "black") +
  GeneralTheme(14) +
  colorspace::scale_fill_discrete_qualitative("Dark 3") +
  guides(color = FALSE, fill = FALSE) +
  labs(x = "Monocyte Likelihood", y = "DC Likelihood")
c <- ggplot(lls, aes(x = maccor, y = dccor, fill = set)) +
  geom_point(alpha = 0.5, size = 3, shape = 21, color = "black") +
  GeneralTheme(14) +
  colorspace::scale_fill_discrete_qualitative("Dark 3") +
  guides(color = FALSE, fill = FALSE) +
  labs(x = "Macrophage Likelihood", y = "")
p <- a + plot_spacer() + b + c & theme(plot.margin = margin(0,0,0,0, "cm"))
ggsave(filename = "plots/mnp_likelihoods_bycluster_byset.png", plot = p, height = 6, width = 6*1, units = "in", dpi = 400)

# [10] M7 plot ----

num_top = 2
geneset_use = "M7"
head(scaled_marker_scores)
plot_scores <- scaled_marker_scores %>% filter(geneset == geneset_use)
plot_scores$set_cluster <- paste0(plot_scores$set, "_", plot_scores$group)
passing_clusters <- plot_scores %>%
  group_by(set, group) %>%
  summarize(mean = median(score), set_cluster = unique(set_cluster)) %>%
  top_n(num_top, mean) %>%
  select(set_cluster)
plot_scores <- plot_scores %>% filter(set_cluster %in% passing_clusters$set_cluster)
plot_scores <- plot_scores %>% filter(source %in% sources)
score_vlnplot <- ggplot(plot_scores,
  aes(x = source, y = score, fill = disease_classification)) +
  geom_jitter(shape = 21, alpha = 0.25, size = 0.2, position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.1)) +
  #geom_boxplot(alpha = 0.75, width = 0.75, position = position_dodge(width = 0.5), outlier.color = "Gray60", outlier.alpha = 0.25) +
  geom_violin(alpha = 0.6, width = 0.75, position = position_dodge(width = 0.5), draw_quantiles = c(0.25, 0.5, 0.75)) +
  scale_fill_manual(values = c("#d55e00", "#0072b2"), name = "Condition", labels = c("Disease", "Health")) +
  scale_x_discrete(labels = gsub("_", " ", unique(plot_scores$source))) +
  scale_y_continuous(limits = c(min(plot_scores$score),
    ceiling(max(plot_scores$score)))) +
  # annotate(geom = "text", label = glue("P = {format.pval(unique(kpis$plogit))}, \u03b2 = {round(unique(kpis$blogit), 2)}"), x = "health",
  #   y = ceiling(max(scores_to_plot$score)), hjust = 0) +
  labs(x = "Dataset", y = "M7 Score") +
  GeneralTheme(12) #+ theme(axis.title.y = element_text(angle = 0, hjust = 0.5, vjust = 1))
score_vlnplot
#SavePlot(filename = "plots/m5_violin_plot.xpdf", plot = m5_violin_plot, base_height = 1.9*2, base_asp = 3.15/1.9)

# check statistics if needed
unique(plot_scores$source)
source_test <- unique(plot_scores$source)[3]
wilcox.test(plot_scores$score[plot_scores$source == source_test & plot_scores$disease_classification == "disease"],
  plot_scores$score[plot_scores$source == source_test & plot_scores$disease_classification == "health"])

names(marker_results) # check names
stats <- imap(marker_results, ~ {
  df <- .x$stats
  df$source_geneset <- .y
  return(df)
})
stats <- bind_rows(stats)
subset_stats <- stats %>% filter(source_geneset %in% paste0(sources, "_M7"))
percentile_plot <- ggplot(subset_stats,
  aes(x = source_geneset, y = diff)) +
  geom_jitter(alpha = 0.1, width = 0.15) +
  #tidybayes::stat_eye(side = "right") +
  #tidybayes::stat_dots(scale = 1, width = 1, alpha = 0.8, size = 20) +
  #geom_point() +
  tidybayes::stat_pointinterval(position = position_nudge(x = 0.2)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "Gray90", size = 1) +
  scale_color_brewer(name = "Median QI") +
  scale_y_continuous(limits = c(-2, 2)) +
  # annotate(geom = "text", label = glue("P = {format.pval(unique(kpis$plogit))}, \u03b2 = {round(unique(kpis$blogit), 2)}"), x = "health",
  #   y = ceiling(max(scores_to_plot$score)), hjust = 0) +
  labs(x = "", y = "90th quantile\n difference") +
  GeneralTheme(12) +
  guides(color = guide_legend(override.aes = list(size = 10))) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks = element_blank()
    #axis.title.y = element_text(angle = 0, hjust = 0.5, vjust = 1)
  )
percentile_plot

pw_b <- (percentile_plot + theme(plot.margin = margin(0,0,0,0))) / score_vlnplot + plot_layout(heights = c(1, 2))
pw_b
SavePlot(filename = "plots/m7_plot.pdf", pw_b, base_height = 1.9*4, base_asp = 3.15/1.9)
SavePlot(filename = "plots/m7_plot.png", pw_b, base_height = 1.9*4, base_asp = 3.15/1.9)

# [10] All marker scores, all program scores ----

head(scaled_marker_scores)
sum_marker_scores <- scaled_marker_scores %>%
  group_by(set, group, geneset) %>%
  summarize(mean_score = mean(score), mean_percent = mean(percents))
sum_marker_scores$set_group <- paste0(sum_marker_scores$set, "_", sum_marker_scores$group)
sum_marker_scores$max <- FALSE
sum_marker_scores <- sum_marker_scores %>% group_by(set, geneset) %>% mutate(max = ifelse(mean_score %in% tail(sort(mean_score), 2), TRUE, FALSE))
sum_marker_scores$geneset <- factor(sum_marker_scores$geneset,
  levels = unique(sum_marker_scores$geneset[order(as.numeric(gsub("M", "", sum_marker_scores$geneset)))]))

marker_plot <- ggplot(sum_marker_scores, aes(x = geneset, y = set_group, fill = mean_score, size = mean_percent, color = max)) +
  geom_point(shape = 21, alpha = 0.8, stroke = 0.5) +
  scale_y_discrete(labels = gsub("_", " ", unique(sum_marker_scores$set_group))) +
  scale_color_manual(values = c("gray80", "#D55E00"), name = "Top 2") +
  scale_fill_gradient2(low = "white", high = "black", name = "Average Score", limits = c(-0.5, 2), na.value = "black") +
  #colorspace::scale_fill_continuous_sequential("Blues 3", name = "Average Score") +
  scale_size_continuous(range = c(0.5, 2), name = "Average % detected") +
  labs(x = "Geneset", y = "Dataset and cluster") +
  GeneralTheme(12) + theme(axis.text.y = element_text(size = 6), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
marker_plot

head(scaled_program_scores)
sum_program_scores <- scaled_program_scores %>%
  group_by(set, group, geneset) %>%
  summarize(mean_score = mean(score), mean_percent = mean(percents))
sum_program_scores$set_group <- paste0(sum_program_scores$set, "_", sum_program_scores$group)
sum_program_scores$max <- FALSE
sum_program_scores <- sum_program_scores %>% group_by(set, geneset) %>% mutate(max = ifelse(mean_score %in% tail(sort(mean_score), 2), TRUE, FALSE))
sum_program_scores$geneset <- factor(sum_program_scores$geneset,
  levels = unique(sum_program_scores$geneset[order(as.numeric(gsub("P", "", sum_program_scores$geneset)))]))

program_plot <- ggplot(sum_program_scores, aes(x = geneset, y = set_group, fill = mean_score, size = mean_percent, color = max)) +
  geom_point(shape = 21, alpha = 0.8, stroke = 0.5) +
  scale_y_discrete(labels = gsub("_", " ", unique(sum_program_scores$set_group))) +
  scale_color_manual(values = c("gray80", "#D55E00"), name = "Top 2") +
  scale_fill_gradient2(low = "white", high = "black", name = "Average Score", limits = c(-0.5, 2), na.value = "black") +
  #colorspace::scale_fill_continuous_sequential("Blues 3", name = "Average Score") +
  scale_size_continuous(range = c(0.5, 2), name = "Average % detected") +
  labs(x = "Geneset", y = "Dataset and cluster") +
  GeneralTheme(12) + theme(axis.text.y = element_text(size = 6), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
program_plot

sum_scores <- rbind(sum_marker_scores, sum_program_scores)
sum_scores$geneset <- factor(sum_scores$geneset,
  levels = c(as.character(unique(sum_marker_scores$geneset[order(as.numeric(gsub("M", "", sum_marker_scores$geneset)))])),
    as.character(unique(sum_program_scores$geneset[order(as.numeric(gsub("P", "", sum_program_scores$geneset)))]))))

combined_plot <- program_plot <- ggplot(sum_scores, aes(x = geneset, y = set_group, fill = mean_score, size = mean_percent, color = max)) +
  geom_point(shape = 21, alpha = 0.8, stroke = 0.5) +
  scale_y_discrete(labels = gsub("_", " ", unique(sum_scores$set_group))) +
  scale_color_manual(values = c("gray80", "#D55E00"), name = "Top 2") +
  scale_fill_gradient2(low = "white", high = "black", name = "Average Score", limits = c(-0.5, 2), na.value = "black") +
  #colorspace::scale_fill_continuous_sequential("Blues 3", name = "Average Score") +
  scale_size_continuous(range = c(0.5, 2), name = "Average % detected") +
  labs(x = "Geneset", y = "Dataset and cluster") +
  GeneralTheme(12) + theme(axis.text.y = element_text(size = 6), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
combined_plot

pw_figs4 <- ((percentile_plot + theme(plot.margin = margin(0,0,0,0))) /
    (score_vlnplot + theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)))) /
  (combined_plot) + plot_layout(heights = c(1, 2, 12))
SavePlot(filename = "plots/figs4.png", pw_figs4, base_height = 14, base_asp = 0.77)
SavePlot(filename = "plots/figs4.pdf", pw_figs4, base_height = 14, base_asp = 0.77)

# [12] Plot individual gene expression ----

object_names

sets_touse <- as.character(cdf$set[cdf$source == "Reyfman_2019"])
markers_topull <- markers_ls[["M10"]]

norm_data <- map(sets_touse, ~ {

  # pull data and convert to dataframe
  norm_data <- GetAssayData(object_ls[[.x]], slot = "counts")
  norm_data <- norm_data[intersect(markers_topull, rownames(norm_data)), ]
  norm_data <- as.data.frame(as.matrix(norm_data))
  norm_data$gene <- rownames(norm_data)
  norm_data <- norm_data %>% pivot_longer(names_to = "cell", -gene)
  norm_data$set <- .x

  umis <- object_ls[[.x]]$nCount_RNA
  umis <- qdapTools::vect2df(umis)
  colnames(umis) <- c("cell", "nUMI")
  norm_data <- merge(norm_data, umis, by = "cell")

  return(norm_data)

})

norm_data <- bind_rows(norm_data)

markers_topull[c(1:10)]

for (i in markers_topull) {
  gg <- ggplot(norm_data %>% filter(gene == i),
    aes(x = value, y = nUMI, fill = set)) +
    geom_point(alpha = 0.5, shape = 21, color = "black") +
    colorspace::scale_fill_discrete_diverging("Blue-Red 2") +
    scale_y_continuous(trans = scales::log10_trans(),
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    labs(x = "Expr.", y = "nUMI") +
    GeneralTheme(12) +
    NoLegend()
  ggsave(plot = gg, filename = glue("plots/m10_plots_reyfman/{i}.png"))
}


merged <- merge(object_ls[[sets_touse[1]]], object_ls[[sets_touse[2]]])
DotPlot(merged, features = markers_topull, group.by = "disease_classification")

merged <- SetIdent(merged, value = "disease_classification")
markers <- FindAllMarkers(merged, features = markers_topull, slot = "counts", test.use = "LR",
  logfc.threshold = 0)
fm <- markers %>% filter(p_val_adj <= 0.05 & avg_logFC >= 0.05)
table(fm$cluster)

# [13] Calculate M1 and M2 score from Xue et al. ----

m1m2 <- readxl::read_xlsx("data/m1m2.xlsx")
msets <- list(m1m2$Gene[m1m2$Classification == "M1"], m1m2$Gene[m1m2$Classification == "M2"])
names(msets) <- c("M1", "M2")
msets$M1 <- UpdateSymbolList(msets$M1)
msets$M2 <- UpdateSymbolList(msets$M2)
msets$M2 <- c("SEPP1", "LPAR6", "SLC4A7", "CA2", "TGFBI", "LTA4H", "SLCO2B1",
  "TPST2", "CERK", "HS3ST2", "IGF1", "ADK", "HNMT", "GAS7", "P2RY13",
  "MAF", "MS4A4A", "MRC1", "CTSC", "HEXB", "LIPA", "FGL2", "CCL13",
  "CD36", "MS4A6A", "CLEC4F", "CD209", "FN1", "P2RY14", "SLC38A6",
  "CLEC7A", "CCL18", "CXCR4", "MSR1", "HS3ST1", "CD302", "CCL23",
  "CHN2", "TLR5", "ALOX15", "EGR2", "DHX8", "TGFBR2")

cors <- pbapply::pblapply(cdf$set[cdf$compare], function(x) {
  ui_info("Plotting {x}")
  object <- object_ls[[x]]

  # score object
  ui_todo("Scoring object...")
  object <- AddScore(object = object, features = msets, nbin = 30, ctrl = 30, name = "Geneset", seed = 1)
  colnames(object@meta.data)[grep("Geneset[[:digit:]]", colnames(object@meta.data))] <- names(msets)
  ui_done("Object scored")

  df <- object@meta.data[, c("nCount_RNA", "disease_classification", "set", "merged_leiden", "batch", names(msets))]
  df <- df %>% rownames_to_column("cell_barcode")
  df <- df %>% gather(key = "program", value = "score", names(msets))
  head(df)
  colnames(df)[5] <- "group"

  # extract scores
  ps <- scaled_program_scores %>% filter(set == x) %>% select(-percents, -percents_geneset)
  ps$cell_barcode <- gsub("(d|h)_", "", ps$cell_barcode)

  # merge scores
  ps <- merge(ps, df, all.x = TRUE, by = "cell_barcode")

  # calculate correlation
  cors <- ps %>%
    group_by(geneset, program) %>%
    summarize(
      cors = cor(score.x, score.y))
  cors$set <- x
  return(cors)
})
cors <- bind_rows(cors)

cors$geneset<- factor(cors$geneset,
  levels = unique(cors$geneset[order(as.numeric(gsub("P", "", cors$geneset)))]))
corrplot <- ggplot(cors, aes(x = geneset, y = cors, fill = program)) +
  geom_point(shape = 21, size = 1, alpha = 0.5, position = position_dodge(0.5)) +
  geom_boxplot(alpha = 0.5, width = 0.4, position = position_dodge(0.5)) +
  scale_fill_manual(values = c("#D55E00", "#009E73"), name = "Polarization") +
  labs(x = "Program", y = "Pearson correlation") +
  GeneralTheme(18) + theme(
    panel.grid.major.x = element_line(color = "gray90", size = 0.5)
  )
SavePlot(filename = "plots/m1m2corr.pdf", plot = corrplot, base_asp = 2.1, base_height = 1.9*3)
SavePlot(filename = "plots/m1m2corr.png", plot = corrplot, base_asp = 2.2, base_height = 1.9*3)


