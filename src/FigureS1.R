# ---
# Description: Figure S1 plotting
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
verbose <- TRUE

load_files <- list.files(path = dir, pattern = suffix, full.names = TRUE)
load_files <- load_files[grep(paste0(ids, collapse = "|"), x = load_files)]
ui_todo("Loading in {length(load_files)} files")
names(load_files) <- ids

object_ls <- pblapply(ids, function(x) readRDS(glue("data/objects/{x}_base_object.rds")))
names(object_ls) <- ids
Clean()

# [2] A ----

gene_stats <- lapply(ids, function(x) {
  genes <- object_ls[[x]]@assays$RNA@meta.features
  total_genes <- nrow(genes)
  no_ensembl <- sum(is.na(genes$ensembl_id))
  not_mapped <- sum(is.na(genes$ensembl_id) & is.na(genes$hgnc_id))
  hgnc_mapped <- sum(!is.na(genes$hgnc_id))
  total_changes <- sum(genes$orig_name != genes$name)
  gene_stats <- data.frame(total = total_genes, num_noensembl = no_ensembl, num_nohgnc = not_mapped, changes = total_changes, set = x)
  return(gene_stats)
})
gene_stats <- Reduce(rbind, gene_stats)
gene_stats <- gather(gene_stats, key = "stat", value = "count", -set)
gene_stats$stat <- factor(gene_stats$stat, levels = c("total", "num_noensembl", "num_nohgnc", "changes"))

labels <- gsub("_", " ", ids)
labels <- snakecase::to_title_case(labels)
labels <- labels[order(labels)]

gene_stats$set <- factor(gene_stats$set, levels = as.character(unique(gene_stats$set))[order(as.character(unique(gene_stats$set)))])
unique(gene_stats$stat)
cols <- ggthemes::colorblind_pal()(4)
cols[1] <- "#999999"
names(cols) <- unique(gene_stats$stat)
geneplot <- ggplot(gene_stats, aes(x = set, y = count, fill = stat)) +
  geom_col(position = position_dodge(0.75), width = 0.75, color = "black") +
  geom_text(aes(label = count), position = position_dodge(0.75), hjust = -0.5, size = 3) +
  scale_x_discrete(labels = labels) +
  #scale_y_continuous(expand = c(0, 0), limits = c(0, 5E4)) +
  scale_y_sqrt(expand = c(0, 0), limits = c(0, 7E4), breaks = c(0, 1000, 10000, 40000)) +
  scale_fill_manual(values = cols, name = "Number of genes", labels = c("Total", "No Ensembl mapping", "No HGNC mapping", "Name changes")) +
  labs(x = "Dataset", y = "Number of genes") +
  theme_classic(base_size = 14) +
  ggeasy::easy_all_text_color("black") +
  theme(
    axis.title.y = element_text(angle = 90, vjust = 0.5, margin = margin(0, 12, 0, 0), face = "bold"),
    axis.title.x = element_text(margin = margin(12, 0, 0, 0), face = "bold"),
    axis.text.y = element_text(margin = margin(0, 4, 0, 0)),
    axis.line = element_blank(),
    panel.border = element_rect(size = 1, color = "black", fill = "transparent"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.justification = "top",
    legend.key.size = unit(1, "line")
  ) + coord_flip()
geneplot
cowplot::save_plot(filename = "plots/genechanges.png", plot = geneplot,
  base_height = 6, base_asp = 1.2)

# [3] B ----

features <- c("nCount_RNA", "nFeature_RNA", "percent_mit", "percent_rib", "percent_hgb", "percent_hsp")
feature_df <- lapply(ids, function(x) {
  df <- object_ls[[x]][[]] %>% select(!!c("batch", features))
  df$set <- x
  return(df)
})
feature_df <- Reduce(rbind, feature_df)

a <- ggplot(feature_df, aes(x = set, y = !!sym(features[1]))) +
  geom_hline(yintercept = 200, size = 1, linetype = "dotted", color = "#D55E00") +
  geom_violin(draw_quantiles = c(0.1, 0.5, 0.9), width = 0.5, alpha = 0.8, fill = "#999999") +
  scale_x_discrete(labels = labels) +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  guides(fill = FALSE) +
  labs(x = "Dataset", y = "UMI counts") +
  theme_classic(base_size = 18) +
  ggeasy::easy_all_text_color("black") +
  theme(
    axis.title.y = element_text(angle = 90, vjust = 0.5, margin = margin(0, 12, 0, 0), face = "bold"),
    axis.title.x = element_text(margin = margin(12, 0, 0, 0), face = "bold"),
    axis.text.y = element_text(margin = margin(0, 4, 0, 0)),
    axis.line = element_blank(),
    panel.border = element_rect(size = 1, color = "black", fill = "transparent"),
  ) +
  coord_flip()

b <- ggplot(feature_df, aes(x = set, y = !!sym(features[2]))) +
  geom_hline(yintercept = 100, size = 1, linetype = "dotted", color = "#D55E00") +
  geom_violin(draw_quantiles = c(0.1, 0.5, 0.9), width = 0.5, alpha = 0.8, fill = "#999999") +
  scale_x_discrete(labels = labels) +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  guides(fill = FALSE) +
  labs(x = "", y = "Genes detected") +
  theme_classic(base_size = 18) +
  ggeasy::easy_all_text_color("black") +
  coord_flip()

c <- ggplot(feature_df, aes(x = set, y = !!sym(features[3]))) +
  geom_hline(yintercept = 20, size = 1, linetype = "dotted", color = "#D55E00") +
  geom_violin(draw_quantiles = c(0.1, 0.5, 0.9), width = 0.5, alpha = 0.8, fill = "#999999") +
  scale_x_discrete(labels = labels) +
  scale_y_continuous(limits = c(0, 60), labels = scales::percent_format(accuracy = 1, scale = 1, suffix = "%")) +
  guides(fill = FALSE) +
  labs(x = "", y = "% Mitochondrial") +
  theme_classic(base_size = 18) +
  ggeasy::easy_all_text_color("black") +
  coord_flip()

d <- ggplot(feature_df, aes(x = set, y = !!sym(features[4]))) +
  geom_hline(yintercept = 50, size = 1, linetype = "dotted", color = "#D55E00") +
  geom_violin(draw_quantiles = c(0.1, 0.5, 0.9), width = 0.5, alpha = 0.8, fill = "#999999") +
  scale_x_discrete(labels = labels) +
  scale_y_continuous(limits = c(0, 60), labels = scales::percent_format(accuracy = 1, scale = 1, suffix = "%")) +
  guides(fill = FALSE) +
  labs(x = "", y = "% Ribosomal") +
  theme_classic(base_size = 18) +
  ggeasy::easy_all_text_color("black") +
  coord_flip()

e <- ggplot(feature_df, aes(x = set, y = !!sym(features[5]))) +
  geom_hline(yintercept = 5, size = 1, linetype = "dotted", color = "#D55E00") +
  geom_violin(draw_quantiles = c(0.1, 0.5, 0.9), width = 0.5, alpha = 0.8, fill = "#999999") +
  scale_x_discrete(labels = labels) +
  scale_y_continuous(limits = c(0, 10), labels = scales::percent_format(accuracy = 1, scale = 1, suffix = "%")) +
  guides(fill = FALSE) +
  labs(x = "", y = "% Hemoglobin") +
  theme_classic(base_size = 18) +
  ggeasy::easy_all_text_color("black") +
  coord_flip()

f <- ggplot(feature_df, aes(x = set, y = !!sym(features[6]))) +
  geom_hline(yintercept = 5, size = 1, linetype = "dotted", color = "#D55E00") +
  geom_violin(draw_quantiles = c(0.1, 0.5, 0.9), width = 0.5, alpha = 0.8, fill = "#999999") +
  scale_x_discrete(labels = labels) +
  scale_y_continuous(limits = c(0, 10), labels = scales::percent_format(accuracy = 1, scale = 1, suffix = "%")) +
  guides(fill = FALSE) +
  labs(x = "", y = "% HSP Family") +
  theme_classic(base_size = 18) +
  ggeasy::easy_all_text_color("black") +
  coord_flip()

grid <- cowplot::plot_grid(a, b + RemoveY(), c + RemoveY(), d + RemoveY(), e + RemoveY(), f + RemoveY(),
  rel_widths = c(1.75, 1, 1, 1, 1, 1), align = "h", nrow = 1)
cowplot::save_plot(filename = "plots/qc_thresholds.png", plot = grid, base_height = 8, base_asp = 3)

# [4] C ----

lapply(object_ls, ncol)
object_ls[[ids[11]]] <- SetIdent(object_ls[[ids[11]]], value = "batch")
qc_stats <- lapply(ids, function(x) {
  object_ls[[x]] <- FilterCells(
    object_ls[[x]],
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
  qc_stats <- as.data.frame(table(object_ls[[x]]$passing_qc))
  qc_stats$set <- x
  colnames(qc_stats) <- c("passing", "freq", "set")
  return(qc_stats)
})
qc_stats <- Reduce(rbind, qc_stats)

gene_stats$set <- factor(gene_stats$set, levels = as.character(unique(gene_stats$set))[order(as.character(unique(gene_stats$set)))])
unique(gene_stats$stat)
scales::show_col(ggthemes::colorblind_pal()(8))
cols <- ggthemes::colorblind_pal()(8)[c(2, 6)]
names(cols) <- unique(qc_stats$passing)
qc_plot <- ggplot(qc_stats, aes(x = set, y = freq, fill = passing)) +
  geom_col(position = position_dodge(0.75), width = 0.75, color = "black") +
  geom_text(aes(label = freq), position = position_dodge(0.75), hjust = -0.5, size = 3) +
  scale_x_discrete(labels = labels) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.5E5), breaks = c(0, 75000, 150000)) +
  scale_fill_manual(values = cols, name = "Passing\nfilters", labels = c("False", "True")) +
  labs(x = "Dataset", y = "Number of cells") +
  theme_classic(base_size = 14) +
  ggeasy::easy_all_text_color("black") +
  theme(
    axis.title.y = element_text(angle = 90, vjust = 0.5, margin = margin(0, 12, 0, 0), face = "bold"),
    axis.title.x = element_text(margin = margin(12, 0, 0, 0), face = "bold"),
    axis.text.y = element_text(margin = margin(0, 4, 0, 0)),
    axis.line = element_blank(),
    axis.ticks.x = element_blank(),
    panel.border = element_rect(size = 1, color = "black", fill = "transparent"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.justification = "top",
    legend.key.size = unit(1, "line")
  ) + coord_flip()
qc_plot
cowplot::save_plot(filename = "plots/passing_qc_counts.png", plot = qc_plot,
  base_height = 6, base_asp = 1.2)

# [5] D ----

suffix <- "annotated.rds"
dir <- "data/objects"
verbose <- TRUE

load_files <- list.files(path = dir, pattern = suffix, full.names = TRUE)
ids <- str_match(string = load_files, pattern = "(data/objects/)(.*_(health|disease))_annotated.rds")[, 3, drop = TRUE]
ui_todo("Loading in {length(load_files)} files")
names(load_files) <- ids

anno_ls <- pblapply(ids, function(x) { readRDS(file = load_files[x]) })
names(anno_ls) <- ids
ids <- SelfName(ids)
Clean()

df <- pblapply(X = ids, function(x) {
  harmony_vars <- apply(anno_ls[[x]]@reductions$harmony@cell.embeddings, 2, var)
  percent <- harmony_vars/sum(harmony_vars)
  df <- data.frame(component = seq(1:50), var = harmony_vars, percent = percent, set = x)
  df$mle <- "No"
  df$mle[df$component == unique(anno_ls[[x]]@misc$harmony_heuristics$tp_pc)] <- "Yes"
  Clean()
  return(df)
})

df <- Reduce(rbind, df)

varplot <- ggplot(df %>% filter(component <= 30), aes(x = component, y = percent, color = set)) +
  geom_line(size = 1, alpha = 0.6) +
  geom_line(data = df %>% filter(component <=30) %>% group_by(component) %>% summarize(mean = mean(percent)), aes(x = component, y = mean), color = "black", size = 1, alpha = 0.8) +
  colorspace::scale_color_discrete_qualitative("Cold", name = "Dataset") +
  scale_x_continuous(expand = c(0.01, 0), breaks = c(0, 5, 10, 15, 20, 25, 30)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Component", y = "Fraction variance explained") +
  guides(color = FALSE) +
  theme_classic(base_size = 14) +
  ggeasy::easy_all_text_color("black") +
  theme(
    axis.title.y = element_text(angle = 90, vjust = 0.5, margin = margin(0, 12, 0, 0), face = "bold"),
    axis.title.x = element_text(margin = margin(12, 0, 0, 0), face = "bold"),
    axis.text.y = element_text(margin = margin(0, 4, 0, 0)),
    axis.line = element_blank(),
    panel.border = element_rect(size = 1, color = "black", fill = "transparent"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.justification = "top",
    legend.key.size = unit(2, "line")
  )

cowplot::save_plot(filename = "plots/varplot.png", plot = varplot,
  base_height = 6, base_asp = 1)
cowplot::save_plot(filename = "plots/varplot.pdf", plot = varplot,
  base_height = 6, base_asp = 1)
# [6] G ----

objects_with_celltype <- sapply(object_ls, function(x) "celltype" %in% colnames(x[[]]))
objects_with_celltype <- ids[objects_with_celltype]

stats <- list()

unique(object_ls[[objects_with_celltype[1]]]$celltype)
unique(object_ls[[objects_with_celltype[2]]]$celltype)
unique(object_ls[[objects_with_celltype[3]]]$celltype)
unique(object_ls[[objects_with_celltype[4]]]$celltype)
unique(object_ls[[objects_with_celltype[5]]]$celltype)
unique(object_ls[[objects_with_celltype[6]]]$celltype)
unique(object_ls[[objects_with_celltype[7]]]$celltype)
unique(object_ls[[objects_with_celltype[8]]]$celltype)
unique(object_ls[[objects_with_celltype[9]]]$celltype)
unique(object_ls[[objects_with_celltype[10]]]$celltype)

mp_idx <- list(c(5), #1
  c(1, 7, 10, 11, 29),
  c(3, 5, 9, 25, 28), #3
  c(4, 6, 15, 20, 24),
  c(1, 9, 17, 29, 31), #5
  c(4, 5, 8, 11, 14, 17, 19, 21, 25),
  c(3, 19, 21), #7
  c(2, 21, 22),
  c(6, 7, 8, 10, 11, 12, 13, 18, 19, 21, 46, 47),
  c(2, 12))
for (i in seq_len(length(objects_with_celltype))) {
  set <- objects_with_celltype[i]
  mp_names <- unique(anno_ls[[set]]$celltype)[mp_idx[[i]]]
  specificity <- caret::specificity(table(anno_ls[[set]]$final_mp, anno_ls[[set]]$celltype %in% mp_names))
  sensitivity <- caret::sensitivity(table(anno_ls[[set]]$final_mp, anno_ls[[set]]$celltype %in% mp_names))
  ari <- mclust::adjustedRandIndex(anno_ls[[set]]$final_mp, anno_ls[[set]]$celltype %in% mp_names)
  stats[[set]] <- data.frame(sens = sensitivity, spec = specificity, ari = ari, set = set)
}
stats <- bind_rows(stats)
stats <- stats %>% gather(key = "stat", value = "value", -set)

statplot <- ggplot(stats, aes(x = stat, y = value, fill = set)) +
  geom_point(shape = 21, size = 4, position = position_dodge(width = 0.5)) +
  colorspace::scale_fill_discrete_sequential("Grays") +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_discrete(labels = c("ARI", "Sensitivity", "Specificity")) +
  labs(x = "Metric", y = "Value") +
  guides(fill = FALSE) +
  theme_classic(base_size = 18) +
  ggeasy::easy_all_text_color("black") +
  theme(
    axis.title.y = element_text(angle = 90, vjust = 0.5, margin = margin(0, 12, 0, 0), face = "bold"),
    axis.title.x = element_text(face = "bold", margin = margin(12, 0, 0, 0)),
    axis.text.y = element_text(margin = margin(0, 4, 0, 0)),
    axis.ticks.y = element_blank(),
    axis.line = element_blank(),
    panel.grid.major.y = element_line(color = "gray90", size = 0.5, linetype = "solid"),
    panel.border = element_rect(size =1, color = "black", fill = "transparent"),
    legend.text = element_text(size = 10),
    legend.justification = "top",
    legend.key.size = unit(0.5, "line")
  )
statplot
cowplot::save_plot(plot = statplot, filename = "plots/agreementstats.png", base_height = 6, base_asp = 1)
cowplot::save_plot(plot = statplot, filename = "plots/agreementstats.pdf", base_height = 6, base_asp = 1)

# [7] H ----

suffix <- "annotated.rds"
dir <- "data/objects"
verbose <- TRUE

load_anno_files <- list.files(path = dir, pattern = suffix, full.names = TRUE)
ids <- str_match(string = load_anno_files, pattern = "(data/objects/)(.*_(health|disease))_annotated.rds")[, 3, drop = TRUE]
ui_todo("Loading in {length(load_anno_files)} files")
names(load_anno_files) <- ids
load_mnp_files <- list.files(path = dir, pattern = "mnp.rds", full.names = TRUE)
names(load_mnp_files) <- ids

# cells_df <- map_dfr(.x = ids, .f = ~ {
#   anno <- readRDS(file = load_anno_files[.x])
#   mnp <- readRDS(file = load_mnp_files[.x])
#   cells <- data.frame(set = .x, mnp = ncol(mnp), total = ncol(anno), frac = ncol(mnp)/ncol(anno))
#   return(cells)
# })

# if already loaded
cells_df <- map_dfr(.x = ids, .f = ~ {
  cells <- data.frame(set = .x, mnp = ncol(mnp_ls[[.x]]), total = ncol(anno_ls[[.x]]), frac = ncol(mnp_ls[[.x]])/ncol(anno_ls[[.x]]))
  return(cells)
})

sum(cells_df$mnp)
sum(map_int(anno_ls, ~ { ncol(.x) }))
sum(map_int(mnp_ls, ~ { length(unique(.x$donor)) }))
sum(map_int(mnp_ls, ~ { sum(table(.x$batch) > 100) }))

cells_df_g <- gather(cells_df, key = "type", value = "count", -frac, -set)

labels <- gsub("_", " ", ids)
labels <- snakecase::to_title_case(labels)
labels <- labels[order(labels)]

cols <- c(ggthemes::colorblind_pal()(8)[c(4)], "Gray80")
count_plot <- ggplot(cells_df_g, aes(x = set, y = count, fill = type)) +
  geom_col(position = "dodge", width = 0.75, color = "black", alpha = 1) +
  scale_x_discrete(expand = c(0, 0.5), labels = labels) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, max(cells_df_g$count) + 1000), breaks = c(0, 10000, 25000, 50000)) +
  scale_fill_manual(values = cols, name = "", labels = c("MNP", "Total")) +
  labs(x = "Set", y = "Number of cells") +
  theme_classic(base_size = 14) +
  ggeasy::easy_all_text_color("black") +
  theme(
    axis.title.y = element_text(angle = 90, vjust = 0.5, margin = margin(0, 12, 0, 0), face = "bold"),
    axis.title.x = element_text(margin = margin(12, 0, 0, 0), face = "bold"),
    axis.text.y = element_text(margin = margin(0, 4, 0, 0)),
    axis.line = element_blank(),
    plot.margin = margin(4, 4, 4, 4, unit = "mm"),
    panel.border = element_rect(size = 1, color = "black", fill = "transparent"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.position = "right",
    legend.justification = "left",
    legend.key.size = unit(1, "line")
  ) + coord_flip()
count_plot
cowplot::save_plot(plot = count_plot, filename = "plots/mnp_counts.png",
  base_height = 6, base_asp = 1.2)
cowplot::save_plot(plot = count_plot, filename = "plots/mnp_counts.pdf",
  base_height = 6, base_asp = 1.2)
