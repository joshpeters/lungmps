# ---
# Description: Prepare supporting data
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

# [2] Prepare and save supporting data ----
sheets <- readxl::excel_sheets(path = "data/travaglini_2020_markers.xlsx")
sheets <- sheets[-grep("SS2|SS", sheets)]
degs <- pblapply(sheets, function(x) {
  degs <- readxl::read_xlsx(path = "data/travaglini_2020_markers.xlsx",
    skip = 1, sheet = x)
  degs <- degs %>% filter(avg_logFC > log(2) & p_val_adj < 1E-10) %>% top_n(20, avg_logFC) %>% pull(Gene)
  return(degs)
})
deg_names <- pblapply(sheets, function(x) {
  deg_names <- readxl::read_xlsx(path = "data/travaglini_2020_markers.xlsx",
    n_max = 1, sheet = x, col_names = FALSE)
  deg_names <- as.character(deg_names[1, 1])
  return(deg_names)
})
deg_names <- snakecase::to_snake_case(unlist(deg_names))
deg_names <- paste0(deg_names, "_score")
names(degs) <- deg_names
type_df <- read.csv("data/travaglini_2020_type.csv")
deg_names_df <- data.frame(index = seq(1:length(deg_names)), name = deg_names)
deg_names_df <- merge(deg_names_df, type_df, by.x = "name", by.y = "score_name", sort = FALSE)
score_names_touse <- as.character(deg_names_df$name[deg_names_df$usage == 1])
deg_names_touse <- deg_names_df[deg_names_df$usage == 1, ]
deg_names_touse$index <- seq(1:nrow(deg_names_touse))
score_labels <- gsub("_", " ", score_names_touse)
score_labels <- snakecase::to_title_case(score_labels)
score_labels <- gsub("Score", " ", score_labels)
score_labels <- str_trim(score_labels)

lm22_expr <- readxl::read_xlsx(path = "data/lm22_expr.xlsx")
is_expr <- readxl::read_xlsx(path = "data/immunoStates_expr.xlsx")
is_expr <- janitor::clean_names(is_expr)
colnames(is_expr)
ref_genes <- is_expr$gene
is_expr <- is_expr[, -c(1)]
rownames(is_expr) <- ref_genes
lm22_expr <- janitor::clean_names(lm22_expr)
colnames(lm22_expr)
ref_genes <- lm22_expr$gene_symbol
lm22_expr <- lm22_expr[, -c(1:3)]
rownames(lm22_expr) <- ref_genes

# [3] Save data ----

saveRDS(object = degs, file = "data/interim/tra20_degs.rds")
saveRDS(object = deg_names_touse, file = "data/interim/tra20_degs_names.rds")
saveRDS(object = is_expr, file = "data/interim/is_expr.rds")
saveRDS(object = lm22_expr, file = "data/interim/lm22_expr.rds")
saveRDS(object = type_df, file = "data/interim/tra20_type_df.rds")
saveRDS(object = score_labels, file = "data/interim/tra20_score_labels.rds")
saveRDS(object = score_names_touse, file = "data/interim/tra20_score_names.rds")
