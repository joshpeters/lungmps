# ---
# Description: Identify conserved, efficient MNP markers
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

suffix <- "mnp_final.rds"
dir <- "data/objects"
load_mnp_files <- list.files(path = dir, pattern = suffix, full.names = TRUE)
ids <- str_match(string = load_mnp_files, pattern = glue("(data/objects/)(.*_(health|disease))_{suffix}"))[, 3, drop = TRUE]
ui_todo("Loading in {length(load_mnp_files)} files")
names(load_mnp_files) <- ids
load_anno_files <- list.files(path = dir, pattern = "annotated.rds", full.names = TRUE)
names(load_anno_files) <- ids
ids <- SelfName(ids)

# load data
Clean()
mnp_ls <- pblapply(ids, function(x) readRDS(file = load_mnp_files[x]))
names(mnp_ls) <- ids

# takes longer, add progress bar
anno_ls <- pblapply(ids, function(x) readRDS(file = load_anno_files[x]))
names(anno_ls) <- ids
Clean()

# add final mnp metadata
anno_ls <- pblapply(X = ids, FUN = function(x) {

  anno_ls[[x]]$final_mp <- FALSE
  anno_ls[[x]]$final_mp[Cells(mnp_ls[[x]])] <- TRUE

  anno_ls[[x]]$mnp_withclusters <- "NonMP"
  #assertthat::assert_that(all.equal(names(anno_ls[[x]]$mnp_withclusters[Cells(mnp_ls[[x]])]), names(mnp_ls[[x]]$merged_leiden)))
  anno_ls[[x]]$mnp_withclusters[Cells(mnp_ls[[x]])] <- mnp_ls[[x]]$merged_leiden
  anno_ls[[x]] <- SetIdent(anno_ls[[x]], value = "mnp_withclusters")

  return(anno_ls[[x]])
})

# [2] Determine markers ----
# determine all markers
Clean()
StartFutureLapply()
StopFutureLapply()
dfs <- pblapply(X = ids, FUN = function(x) {

  ui_info("Finding markers for {x}")
  clusters <- as.character(unique(mnp_ls[[x]]$merged_leiden))
  de_rna <- map(clusters, CalcAUC, anno_ls[[x]])
  markers <- bind_rows(de_rna)
  anno_ls[[x]]@misc$mp_wa_markers <- markers
  filtered_wa <- markers %>%
    group_by(group) %>%
    filter(logFC > log(2) & auc > 0.6 & padj < 1E-5)
  filtered_wa$set <- x
  Clean()
  return(filtered_wa)

})

# save all markers
saveRDS(dfs, file = "data/interim/conserved_markers_df.rds")

# compile markers
dfs <- bind_rows(dfs)
dfs_sum <- dfs %>%
  group_by(feature) %>%
  summarize(n_gene = n(),
    n_set = n_distinct(set),
    median_logFC = median(logFC, na.rm = TRUE),
    median_AUC = median(auc, na.rm = TRUE))
hist(dfs_sum$n_set)
auc_genes <- dfs_sum %>% filter(n_set >= 18) %>% top_n(20, median_AUC)
saveRDS(auc_genes, file = "data/interim/auc_genes_df.rds")
logfc_genes <- dfs_sum %>% filter(n_set >= 18) %>% top_n(20, median_logFC)
features_to_plot <- unique(c(auc_genes$feature, logfc_genes$feature))

# compile expression
expr <- pblapply(X = ids, FUN = function(x) {

  anno_ls[[x]] <- SetIdent(anno_ls[[x]], value = "final_mp")
  avg_expr <- AverageExpression(object = anno_ls[[x]],
    features = auc_genes$feature, verbose = TRUE)
  avg_expr <- as.matrix(avg_expr$RNA)
  avg_expr <- as.data.frame(t(avg_expr))
  avg_expr$mnp <- rownames(avg_expr)
  mnp_expr <- gather(avg_expr, key = "feature", value = "expr", -mnp)
  mnp_expr$set <- x
  return(mnp_expr)

})

expr_collapsed <- bind_rows(expr)
expr_collapsed <- expr_collapsed %>% mutate(label = glue("{gsub('_', ' ', set)} {mnp}"))
expr_collapsed <- expr_collapsed %>%
  group_by(feature) %>%
  #filter(feature %in% auc_genes$feature) %>%
  mutate(nor_expr = expr/max(expr), scale_expr = (expr - mean(expr))/sd(expr))
typeof(expr_collapsed$nor_expr)

saveRDS(expr_collapsed, file = "data/processed/mp_consensusmarkers_fig1.rds")

LoadColors()
marker_plot <- ggplot(expr_collapsed, aes(x = fct_reorder(label, nor_expr, .desc = TRUE),
  y = fct_rev(feature), fill = scale_expr)) +
  geom_tile(color = "black", size = 0.5) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  #colorspace::scale_fill_continuous_diverging("Blue-Red 3", name = "Gene\nstandardized\nexpression\n") +
  scale_fill_gradient2(low = colorspace::lighten("#0072B2", amount = 0.2),
    mid = "Gray95", high = colorspace::lighten("#D55E00", amount = 0.2), name = "Gene\nstandardized\nexpression\n") +
  labs(x = "", y = "Feature") +
  theme_bw(base_size = 12) +
  ggeasy::easy_all_text_color("black") +
  theme(
    plot.subtitle = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 6),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(face = "bold", margin = margin(t = 2, unit = "pt")),
    axis.title.y = element_text(face = "bold", margin = margin(r = 12, unit = "pt")),
    plot.margin = margin(4, 4, 4, 4, unit = "mm"),
    axis.ticks.x = element_blank(),
    panel.border = element_rect(size = 1, color = "black", fill = "transparent"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    legend.position = "right",
    legend.justification = "top",
    legend.key.size = unit(1, "line"),
  )
marker_plot
cowplot::save_plot(filename = glue("plots/mnp_marker.pdf"),
  plot = marker_plot, base_height = 6, base_asp = 1.333)

# [3] Determine thresholds ----

anno_ls <- pbapply::pblapply(ids, function(x) {
  anno_ls[[x]] <- AddModuleScore(anno_ls[[x]], features = list("reduced" = auc_genes$feature), nbin = 30, ctrl = 30, name = "mp_reduced")
  return(anno_ls[[x]])
})

# [2] Pull scores
scores_ls <- pbapply::pblapply(ids, function(x) {
  scores <- anno_ls[[x]]@meta.data %>% select(mp_diff_score, mp_reduced1, final_mp, batch)
  scores$set <- x
  return(scores)
})
scores <- Reduce(rbind, scores_ls)

red <- GetCBPalette()[7]
blue <- GetCBPalette()[6]

aucs <- map(ids, ~ {
  scores.use <- scores %>% filter(set == .x)
  positive <- scores.use$mp_diff_score[scores.use$ final_mp == TRUE]
  negative <- scores.use$mp_diff_score[scores.use$ final_mp == FALSE]
  predictions <- ROCR::prediction(predictions = c(positive, negative),
    labels = c(rep(x = 1, length(x = positive)), rep(x = 0, length(x = negative))),
    label.ordering = 0:1)
  perf.use <- ROCR::performance(prediction.obj = predictions, measure = "auc")
  auc.use <- round(x = perf.use@y.values[[1]], digits = 3)
  return(auc.use)
})

aucs <- map(ids, ~ {
  scores.use <- scores %>% filter(set == .x)
  positive <- scores.use$mp_reduced1[scores.use$ final_mp == TRUE]
  negative <- scores.use$mp_reduced1[scores.use$ final_mp == FALSE]
  predictions <- ROCR::prediction(predictions = c(positive, negative),
    labels = c(rep(x = 1, length(x = positive)), rep(x = 0, length(x = negative))),
    label.ordering = 0:1)
  perf.use <- ROCR::performance(prediction.obj = predictions, measure = "auc")
  auc.use <- signif(x = perf.use@y.values[[1]], digits = 3)
  return(auc.use)
})

ids
bar_labels <- gsub("_", " ", ids)
bar_labels <- snakecase::to_title_case(bar_labels)
bar_labels <- bar_labels[order(bar_labels)]
bar_labels

c <- ggplot(scores, aes(x = mp_reduced1, y = set, fill = final_mp)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray80", size = 0.5) +
  ggridges::geom_density_ridges(alpha = 0.8) +
  scale_fill_manual(values = c(shades::lightness(blue, shades::delta(20)), shades::lightness(red, shades::delta(20))),
    name = "Filtered\nMononuclear\nPhagocyte", labels = c("False", "True")) +
  scale_y_discrete(labels = paste0(bar_labels, "   (", aucs, ")  ")) +
  labs(x = "MP Score", y = "Dataset") +
  theme_classic(base_size = 14) +
  ggeasy::easy_all_text_color("black") +
  theme(
    axis.title.y = element_text(angle = 90, vjust = 0.5, margin = margin(0, 12, 0, 0), face = "bold"),
    axis.title.x = element_text(margin = margin(12, 0, 0, 0), face = "bold"),
    axis.text.y = element_text(margin = margin(0, 4, 0, 0)),
    axis.line = element_blank(),
    #axis.ticks.x = element_blank(),
    plot.margin = margin(4, 4, 4, 4, unit = "mm"),
    panel.border = element_rect(size = 1, color = "black", fill = "transparent"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.position = "right",
    legend.justification = "top",
    legend.key.size = unit(1, "line")
  )
c

cowplot::save_plot(filename = "plots/type_score.png",
  plot = c, base_asp = 1.4, base_height = 5)
cowplot::save_plot(filename = "plots/type_score.pdf",
  plot = c, base_asp = 1.4, base_height = 5)
saveRDS(list(scores, bar_labels), file = "data/interim/mp_signaturescores.rds")

# [4] Plot proportions of MNPs ----
names(anno_ls)
anno_ls <- anno_ls[-c(8,9)]
ids <- ids[-c(8,9)]

i <- 1
set <- ids[[i]]
ui_info("{set}")
table(anno_ls[[set]]$leukocyte)

anno_ls[[set]]$leukocyte <- FALSE
anno_ls[[set]]$leukocyte[anno_ls[[set]]$jmp_class == "immune"] <- TRUE
#anno_ls[[set]]$leukocyte[anno_ls[[set]]$base_clustering %in% c(6)] <- FALSE
table(anno_ls[[set]]$leukocyte)

# 4, 17
# 5, 18
# 17, 6

props <- lapply(X = ids, FUN = function(x) {
  res_props <- anno_ls[[x]]@meta.data %>%
    dplyr::group_by(batch) %>%
    dplyr::summarise(class = sum(leukocyte == TRUE),
      mnp = sum(final_mp == TRUE)) %>%
    dplyr::mutate(prop = mnp / class)
  res_props$set <- x
  return(res_props)
})
props <- do.call(rbind, props)

sets <- readxl::read_excel("data/study_metadata.xlsx")
lung_sets <- sets %>% filter(Tissue == "Lung")
lung_ids <- lung_sets$ID
disease_sets <- lung_sets %>% filter(Type != "Healthy") %>% pull(ID)
healthy_sets <- lung_sets %>% filter(Type == "Healthy") %>% pull(ID)
names(disease_sets) <- disease_sets
names(healthy_sets) <- healthy_sets
dataset_healthy <- data.frame(set = healthy_sets, disease = rep("Control", length(healthy_sets)))
dataset_disease <- data.frame(set = c("Valenzi_2019", "Morse_2019", "Habermann_2019", "Reyfman_2019",
  "Lambrechts_2018", "Laughney_2020", "Zilionis_2019", "Mayr_2020"),
  disease = c("ILD", "ILD", "ILD", "ILD", "Cancer", "Cancer", "Cancer", "ILD"))
conditions <- stringr::str_match(ids, "([a-zA-Z]*_[0-9]{4})_(.*)$")
conditions <- conditions[, c(1,2,3)]
conditions <- as.data.frame(conditions)
colnames(conditions) <- c("ID", "Study", "Condition")
# conditions$Type <- c("Healthy", "Chronic", "Healthy", "Cancer", "Healthy",
#   "Cancer", "Healthy", "Healthy", "Chronic", "Healthy",
#   "Chronic", "Healthy", "Healthy", "Chronic", "Healthy",
#   "Healthy", "Chronic", "Cancer")
conditions$Disease <- "Health"
matched <- match(conditions$Study[conditions$Condition != "health"], dataset_disease$set)
conditions$Disease[conditions$Condition != "health"] <- as.character(dataset_disease$disease[matched])

props <- merge(props, conditions, by.x = "set", by.y = "ID")
unique(as.factor(props$set))
labels <- gsub("_", " ", sort(unique(props$Study)))
med_props <- props %>% group_by(Disease) %>% summarise(median_props = median(prop))

head(props)
plot <- ggplot(props, aes(x = Study, y = prop, fill = Disease)) +
  geom_hline(yintercept = med_props$median_props[1], size = 1, color = "#E69F00", alpha = 0.5, linetype = "dotted") +
  geom_hline(yintercept = med_props$median_props[2], size = 1, color = "#CC79A7", alpha = 0.5, linetype = "dotted") +
  geom_hline(yintercept = med_props$median_props[3], size = 1, color = "#0072B2", alpha = 0.5, linetype = "dotted") +
  geom_point(alpha = 0.6, position = position_dodge(width = 0.5), shape = 21, color = "black", size = 2.5) +
  stat_summary(fun = median, fun.min = median, fun.max = median,
    geom = "crossbar", width = 0.5, position = position_dodge(width = 0.5)) +
  #geom_violin(position = "dodge", width = 0.5) +
  scale_x_discrete(labels = labels) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = ggthemes::colorblind_pal()(8)[c(6,8,4)]) +
  labs(x = "", y = "% phagocytes of immune cells") +
  theme_classic(base_size = 14) +
  ggeasy::easy_all_text_color("black") +
  theme(
    axis.title.y = element_text(angle = 90, vjust = 0.5, margin = margin(0, 12, 0, 0), face = "bold"),
    axis.title.x = element_text(margin = margin(12, 0, 0, 0), face = "bold"),
    axis.text.y = element_text(margin = margin(0, 4, 0, 0)),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.line = element_blank(),
    #axis.ticks.x = element_blank(),
    plot.margin = margin(2, 2, 2, 2, unit = "mm"),
    panel.border = element_rect(size = 1, color = "black", fill = "transparent"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.position = "right",
    legend.justification = "top",
    legend.key.size = unit(1, "line")
  )
plot

cowplot::save_plot(plot = plot, filename = glue("plots/mnp_props.png"), base_height = 5.5, base_asp = 1.2)
cowplot::save_plot(plot = plot, filename = glue("plots/mnp_props.pdf"), base_height = 5.5, base_asp = 1.2)
saveRDS(list(props, med_props, labels), "data/interim/props_data.rds")

# [5] Resave annotated and MNP objects ----
pblapply(X = ids, function(x) {
  ui_info("\nSaving {x} ({round(as.numeric(object.size(anno_ls[[x]])/1E6), 2)} mb)")
  saveRDS(anno_ls[[x]],
    file = glue::glue("data/objects/{x}_annotated.rds"))
  Clean()
  ui_info("\nSaving {x} ({round(as.numeric(object.size(mnp_ls[[x]])/1E6), 2)} mb)")
  saveRDS(mnp_ls[[x]],
    file = glue::glue("data/objects/{x}_mnp_final.rds"))
  Clean()
  ui_done(("\nCompleted {x}"))
})

