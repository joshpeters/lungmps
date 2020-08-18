# ---
# Description: iNMF script for cluster
# Author: Josh Peters
# ---

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(glue))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(liger))

# grab arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 4) {
  stop("Missing arguments")
}

# set variables and load additional functions
task_id <- as.numeric(args[[1]])
name <- as.character(args[[2]])
k <- as.numeric(args[[3]])
lambda <- as.numeric(args[[4]])

date_var <- format(Sys.Date(), "%Y%m%d")
time_var <- format(Sys.time(), "%Y%m%d%H%M")
set.seed(task_id)

cat(glue("----\nScript initiated at: {time_var}\nArgs: \n{task_id}\n{k}\n{lambda}"),
  file = glue("/broad/blainey_lab/jpeters/projects/lung_biorxiv/final_run/log_{date_var}.txt"), append = TRUE, sep = "\n")


working_dir <- ""
# load genes
genes <- readRDS(glue("{working_dir}/genestouse.rds"))

# load matrices
files <- list.files(working_dir, full.names = TRUE)
files <- files[grep("matrices.rds", files)]
object_names <- c(Braga_2019_health = "Braga_2019_health", Habermann_2019_disease = "Habermann_2019_disease",
  Habermann_2019_health = "Habermann_2019_health", Lambrechts_2018_disease = "Lambrechts_2018_disease",
  Lambrechts_2018_health = "Lambrechts_2018_health", Laughney_2020_disease = "Laughney_2020_disease",
  Laughney_2020_health = "Laughney_2020_health", Madissoon_2019_health = "Madissoon_2019_health",
  Mayr_2020_disease = "Mayr_2020_disease", Mayr_2020_health = "Mayr_2020_health",
  Morse_2019_disease = "Morse_2019_disease", Morse_2019_health = "Morse_2019_health",
  Raredon_2019_health = "Raredon_2019_health", Reyfman_2019_disease = "Reyfman_2019_disease",
  Reyfman_2019_health = "Reyfman_2019_health", Travaglini_2020_health = "Travaglini_2020_health",
  Valenzi_2019_disease = "Valenzi_2019_disease", Zilionis_2019_disease = "Zilionis_2019_disease"
)

object_list <- lapply(files, function(x) {
  obj <- readRDS(file = x)
  num_cells <- lapply(obj, ncol)
  obj <- obj[num_cells > 100]
  return(obj)
})
names(object_list) <- object_names

object_list <- lapply(object_names, function(x) {
  liger_obj <- liger::createLiger(object_list[[x]], take.gene.union = FALSE)
  liger_obj <- liger::normalize(liger_obj)
  all_genes <- lapply(liger_obj@norm.data, rownames)
  all_genes <- Reduce(intersect, all_genes)
  genes <- genes[genes %in% all_genes]
  liger_obj@var.genes <- genes
  liger_obj <- liger::scaleNotCenter(liger_obj)
  liger_obj <- liger::optimizeALS(liger_obj, k = k, lambda = lambda, thresh = 1E-05, max.iters = 100, nrep = 1, rand.seed = task_id)
  liger_obj <- liger::quantile_norm(liger_obj)
  saveRDS(list(W = liger_obj@W, H = liger_obj@H, H_norm = liger_obj@H.norm, V = liger_obj@V, genes = liger_obj@var.genes),
    file = glue("{working_dir}/{name}/{x}_{task_id}_{k}_{time_var}.rds"))
  gc()
  return(liger_obj)
})

