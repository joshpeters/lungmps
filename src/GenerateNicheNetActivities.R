# ---
# Description: Run NicheNet analysis
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

library(nichenetr)

# [2] Load files ----

# load nn files
ligand_target_matrix <- readRDS("~/Projects/lungmps/data/ligand_target_matrix.rds")
sort(ligand_target_matrix[, "IFNG"], decreasing = TRUE)
lr_network <- readRDS("data/lr_network.rds")
weighted_networks <- readRDS("data/weighted_networks.rds")
weighted_networks_lr <- readRDS("data/weighted_networks_lr.rds")

# [3] Prepare inputs ----

object_ls <- LoadObjectList(dir.path = "/Users/jpeters/Projects/lungmps/data/objects",
  file_pattern = "_annotated.rds$", name_pattern = "_annotated.rds$", exclude = "Liao_2020", verbose = 1)
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
loadings <- loadings[all_programs, ]

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

# [4] Test only top 2 clusters with scores ------------------------------------

activities <- pblapply(cdf$set[cdf$use], function(x) {
  sa <- ScoreLigandActivity(object = object_ls[[x]], object_name = x, genesets = programs_ls, geneset_names = programs)
  return(sa)
})
saveRDS(activities, file = "data/interim/nichenet_activities_programs_200720.rds")
#activities <- readRDS("data/interim/nichenet_activities_programs.rds")

activities <- bind_rows(activities)
sum_activities <- activities %>% group_by(program, test_ligand) %>% summarize(mp = mean(pearson), mauroc = mean(auroc))
top_sa <- sum_activities %>% group_by(program) %>% filter(mp > 0.01) %>% top_n(10, mp)
#top_sa$program <- factor(top_sa$program, levels = as.character(order(as.numeric(top_sa$program))))
top_sa$program <- factor(top_sa$program, levels = programs)
#reordered_names <- programs[as.character(unique(top_sa$program[order(as.numeric(gsub("P_", "", top_sa$program)))]))]
unique(top_sa$program)

a <- ggplot(top_sa, aes(x = test_ligand, y = program, fill = mp, size = mauroc)) +
  geom_point(shape = 21, stroke = 1, alpha = 0.75) +
  #colorspace::scale_fill_discrete_qualitative("Dark 2", name = "Program", labels = paste(gsub("P_", "", levels(top_sa$program)), reordered_names, sep = " ")) +
  colorspace::scale_fill_continuous_sequential("Blues 3", name = "Pearson\ncorrelation\ncoefficient") +
  scale_size_continuous(name = "auROC") +
  #scale_y_discrete(labels = paste(programs, program_gene_names, sep = " ")) +
  labs(x = "Ligands", y = "Program") +
  PointTheme(base_size = 18) +
  guides(fill = guide_legend(override.aes = list(size = 3))) +
  theme(
    legend.key.size = unit(1, "line"),
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1, vjust = 1)
    )
a
SavePlot(plot = a, filename = "plots/nichenet_programs_200720.pdf", base_height = 6.5, base_asp = 2.6)