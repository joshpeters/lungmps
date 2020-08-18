# ---
# Description: Identify correlated surface receptors
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

# [2] Load files ----

# load surface receptors
cs_genes <- readRDS(file = "/Users/jpeters/Projects/sctbgranulomas/data/reference/cell_surface_matrix.rds")

# load Seurat objects
object_ls <- LoadObjectList(dir.path = "/Users/jpeters/Projects/lungmps/data/objects",
  file_pattern = "_mnp_final.rds$", name_pattern = "_mnp_final.rds$", exclude = "Liao_2020", verbose = 1)
object_names <- names(object_ls)
names(object_names) <- object_names
object_ls <- lapply(object_names, function(x) {
  object_ls[[x]]$set <- x
  return(object_ls[[x]])
})

object_names
object_names <- object_names[cdf$set[cdf$use == TRUE]]
object_names

# [3] Calculate correlations ----
genesets <- programs_ls
geneset_names <- programs
slot_touse = "counts"
percent_cutoff = 0.1
corr_method = "spearman"
receptors <- cs_genes$symbol[cs_genes$cspa_num == 1 | cs_genes$cspa_num == 2]
hugo_cd <- read_csv(file = "/Users/jpeters/Projects/sctbgranulomas/data/raw/group-471.csv", skip = 1)
receptors <- c(receptors, hugo_cd$`Approved symbol`)
receptors <- unique(receptors)
"TREM2" %in% receptors
receptors <- c(receptors, "TREM2")

receptor_cors <- pbapply::pblapply(object_names, function(x) {

  # prepare matrices
  object_ls[[x]] <- AddModuleScore(object_ls[[x]], genesets, nbin = 30, ctrl = 30, name = "Geneset")
  assertthat::assert_that(length(grep("Geneset[[:digit:]]", colnames(object_ls[[x]]@meta.data))) == length(genesets))
  colnames(object_ls[[x]]@meta.data)[grep("Geneset[[:digit:]]", colnames(object_ls[[x]]@meta.data))] <- geneset_names
  scores <- object_ls[[x]][[]][, programs]
  norm_counts <- GetAssayData(object = object_ls[[x]], slot = slot_touse)
  cs_genes_touse <- intersect(receptors, rownames(object_ls[[x]]))

  # filter genes
  norm_counts <- norm_counts[cs_genes_touse, ]
  mean_norm_counts <- Matrix::rowMeans(norm_counts, na.rm = TRUE)
  norm_counts <- norm_counts[Matrix::rowSums(norm_counts) > 0, ]
  norm_counts_percent <- apply(norm_counts, 1, function(x) { sum(x > 0)/length(x)})
  norm_counts_percent <- norm_counts_percent[norm_counts_percent > percent_cutoff]
  norm_counts <- norm_counts[names(norm_counts_percent), ]

  # calculate correlations
  receptor_program_cors <- cor(t(as.matrix(norm_counts)), scores, method = corr_method)
  assertthat::assert_that(all(rownames(receptor_program_cors) %in% cs_genes_touse))

  return(receptor_program_cors)
})

# [4] Collate results ----

receptor_cors_format <- lapply(receptor_cors, function(x) {
  x <- as.data.frame(x)
  x$gene <- rownames(x)
  return(x)
})
receptor_cors_rbind <- Reduce(f = rbind, x = receptor_cors_format)
receptor_cors_gather <- gather(receptor_cors_rbind, key = "program", value = "corr", -gene)
hist(receptor_cors_gather$corr)
sum_rcors <- receptor_cors_gather %>%
  group_by(gene, program) %>%
  summarize(n = n(), mean_corr = mean(corr), median_corr = median(corr))
genes_toplot <- sum_rcors %>%
  group_by(gene, program) %>%
  filter(n > 6 & abs(mean_corr) > 0.1 & gene %in% cs_genes$symbol[cs_genes$cspa_num == 1]) %>%
  group_by(program) %>%
  filter(dense_rank(mean_corr) <= 3 | dense_rank(desc(mean_corr)) <= 3)
rcors_toplot <- sum_rcors %>% filter(gene %in% genes_toplot$gene & abs(mean_corr) >= 0.1)

rcors_toplot$program <- factor(rcors_toplot$program,
  levels = stringr::str_sort(unique(rcors_toplot$program), numeric = TRUE))
#reordered_names <- program_gene_names[as.character(unique(rcors_toplot$program[order(as.numeric(gsub("P_", "", rcors_toplot$program)))]))]
rcors_toplot$program

# [5] Visualize results ----
rec_cor_plot <- ggplot(rcors_toplot, aes(x = gene, y = program, fill = mean_corr)) +
  geom_point(shape = 21, alpha = 1, stroke = 1, size = 6, color = "gray80") +
  colorspace::scale_fill_continuous_diverging("Blue-Red 3", name = "Mean\nSpearman\ncorrelation") +
  #scale_y_discrete(labels = paste(gsub("P_", "", levels(rcors_toplot$program)), reordered_names, sep = " ")) +
  labs(x = "Gene", y = "Program") +
  GeneralTheme(base_size = 18) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
    panel.grid.major = element_line(size = 0.5, color = "gray90"))
rec_cor_plot
SavePlot(plot = rec_cor_plot, filename = "plots/receptor_corr.pdf", base_asp = 2.6, base_height = 1.6*4)
SavePlot(plot = rec_cor_plot, filename = "plots/receptor_corr.png", base_asp = 2.6, base_height = 1.6*4)

